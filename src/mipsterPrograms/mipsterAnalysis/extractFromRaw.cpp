/*
 * extractFromRaw.cpp
 *
 *  Created on: Aug 9, 2017
 *      Author: nick
 */

// MIPWrangler - A library for analyzing sequence data from molecular inversion probe analysis
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of MIPWrangler.
//
// MIPWrangler is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MIPWrangler is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MIPWrangler.  If not, see <http://www.gnu.org/licenses/>.
//
#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"


namespace njhseq {

int mipsterAnalysisRunner::extractFromRaw(const njh::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	extractFromRawParsMultiple pars;
	setUp.setUpExtractFromRawMultiple(pars);

	//check for required external programs

	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.mips_->setAllMinimumExpectedLen(pars.minLen);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	if("" != pars.sampleMetaFnp){
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, njh::files::make_path( mipMaster.directoryMaster_.logsDir_,
			pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);


	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for(const auto & f : files){
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find("_"));
		if(njh::in(sampName, mipMaster.names_->samples_)){
			readPairs[sampName].emplace_back(filename);
		}else{
			readPairsUnrecognized[sampName].emplace_back(filename);
		}
	}
	auto keys = getVectorOfMapKeys(readPairs);
	njh::sort(keys);
	std::vector<bfs::path> emptyFiles;
	VecStr samplesExtracted;
	VecStr samplesEmpty;
	njh::concurrent::LockableQueue<std::string> filesKeys(keys);
	Json::Value logs;
	logs["mainCommand"] = setUp.commands_.commandLine_;
	logs["qualityFilteringPars"] = pars.qFilPars_.toJson();
	std::mutex logsMut;

	uint64_t maxLen = pars.minLen * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();

	auto extractFiles =
			[&pars,&samplesExtracted,&samplesEmpty,&filesKeys,&emptyFiles,&aligners,&logs,&logsMut,&mipMaster](
					const std::string & outputDirectory,
					const std::unordered_map<std::string, VecStr>& readPairs) {
				std::string key = "";
				VecStr currentEmptySamps;
				VecStr currentSamplesExtracted;
				std::vector<bfs::path> currentEmptyFiles;
				std::unordered_map<std::string, Json::Value> currentLogs;
				std::regex r1_pattern{R"((.*)_R1_(.*).fastq.gz$)"};
				std::regex r2_pattern{R"((.*)_R2_(.*).fastq.gz$)"};
				auto alignerObjForFamilyDet = aligners.popAligner();
				while(filesKeys.getVal(key)){
					Json::Value log;
					njh::stopWatch watch;
					std::stringstream out;
					std::stringstream err;
					bool success = true;
					bool allEmpty = true;
					out << key << std::endl;
					auto sampDir = njh::files::makeDir(outputDirectory, njh::files::MkdirPar(key));
					std::vector<bfs::path> r1_files;
					std::vector<bfs::path> r2_files;
					struct PairedFnps{
						bfs::path r1_fnp = "";
						bfs::path r2_fnp = "";
					};
					std::unordered_map<std::string, PairedFnps> pairedFnpForSample;
					for(const auto & f : readPairs.at(key)){
						auto fqName = bfs::path(f);
						auto inputName = njh::files::make_path(pars.dir,fqName);
						//if the file is empty copy to an empty directory to avoid concatenating errors
						if(0 == bfs::file_size(inputName)){
							emptyFiles.emplace_back(inputName);
							continue;
						}
						allEmpty = false;
						std::smatch r1Match;
						std::smatch r2Match;

						bool matchesR1 = std::regex_match(fqName.string(), r1Match, r1_pattern);
						bool matchesR2 = std::regex_match(fqName.string(), r2Match, r2_pattern);
						if(matchesR1 && matchesR2){
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, " << fqName << ", matched both R1 or R2 pattern" << "\n";
							throw std::runtime_error{ss.str()};
						}else if(matchesR1){
							std::string r1Name = njh::pasteAsStr(r1Match[1], r1Match[2]);
							if(njh::has(pairedFnpForSample, r1Name) && "" != pairedFnpForSample[r1Name].r1_fnp ){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r1Name].r1_fnp << " for " << f << " for " << r1Name << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								pairedFnpForSample[r1Name].r1_fnp = inputName;
							}
						}else if(matchesR2){
							std::string r2Name = njh::pasteAsStr(r2Match[1], r2Match[2]);
							if(njh::has(pairedFnpForSample, r2Name) && "" != pairedFnpForSample[r2Name].r2_fnp ){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r2Name].r2_fnp << " for " << f << " for " << r2Name << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								pairedFnpForSample[r2Name].r2_fnp = inputName;
							}
						}else{
							std::stringstream ss;
							ss << __PRETTY_FUNCTION__ << ", error, " << fqName << ", didn't match either R1 or R2 pattern" << "\n";
							throw std::runtime_error{ss.str()};
						}
					}

					if(allEmpty){
						currentEmptySamps.emplace_back(key);
						success = false;
					}
					std::vector<SeqIOOptions> sampleOpts;
					if(pairedFnpForSample.size() >= 1 && success && !allEmpty){
						for(const auto & pairedFnp : pairedFnpForSample){
							if("" == pairedFnp.second.r1_fnp  ||
									"" == pairedFnp.second.r2_fnp){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, " << key << ", " << pairedFnp.first << ", had an empty r1 or r2 file list" << "\n";
								ss << "R1:" << pairedFnp.second.r1_fnp << "\n";
								ss << "R2:" << pairedFnp.second.r2_fnp << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								sampleOpts.emplace_back(SeqIOOptions::genPairedInGz(pairedFnp.second.r1_fnp, pairedFnp.second.r2_fnp));
								sampleOpts.back().outFormat_ = SeqIOOptions::outFormats::FASTQPAIRED;
								sampleOpts.back().out_.outExtention_ = "_R1.fastq";
								sampleOpts.back().revComplMate_ = true;
							}
						}
					}
					if(success){
						currentSamplesExtracted.emplace_back(key);
						auto samplePars = pars.createForSample(key);
						MipExtractor mExtractor(pars.verbose_);
						mExtractor.extractFilterSampleForMipsPaired(sampleOpts,
								mipMaster,
								*alignerObjForFamilyDet,
								samplePars);
					}
					log["sample"] = key;
					log["time"] = watch.totalTime();
					log["out"] = out.str();
					log["err"] = err.str();
					log["success"] = njh::json::toJson(success);
					currentLogs[key] = log;
				}
				{
					std::lock_guard<std::mutex> lock(logsMut);
					addOtherVec(samplesEmpty, currentEmptySamps);
					addOtherVec(samplesExtracted, currentSamplesExtracted);
					addOtherVec(emptyFiles, currentEmptyFiles);
					for(const auto & keyLog : currentLogs){
						logs[keyLog.first] = keyLog.second;
					}
				}
			};
	std::vector<std::thread> threads;
	for (uint32_t tNum = 0; tNum < pars.numThreads; ++tNum) {
		threads.emplace_back(
				std::thread(extractFiles,
						std::cref(mipMaster.directoryMaster_.masterDir_.string()),
						std::cref(readPairs)));
	}

	njh::concurrent::joinAllThreads(threads);

	logFile << logs << std::endl;


	//copy over resources;

	bfs::copy_file(pars.mipArmsFileName,njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
	njh::files::makeDirP(njh::files::MkdirPar(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo")));
	std::ofstream outSamplesFoundFile;
	openTextFile(outSamplesFoundFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesFound.tab.txt")));
	std::ofstream outSamplesEmptyfile;
	openTextFile(outSamplesEmptyfile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesEmpty.tab.txt")));
	std::ofstream gzCatSamplesFile;
	openTextFile(gzCatSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/gzCatSamples.tab.txt")));
	std::ofstream pearExtractSamplesFile;
	openTextFile(pearExtractSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/emptyFiles.tab.txt")));
	std::ofstream samplesMissingFile;
	openTextFile(samplesMissingFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/samplesMissing.tab.txt")));
	auto allFoundSamps = concatVecs(samplesExtracted, samplesEmpty);
	njh::sort(allFoundSamps);
	njh::sort(samplesEmpty);
	njh::sort(samplesExtracted);
	njh::sort(emptyFiles);
	if(!allFoundSamps.empty()) outSamplesFoundFile << njh::conToStr(allFoundSamps, "\n") << std::endl;
	if(!samplesEmpty.empty()) outSamplesEmptyfile << njh::conToStr(samplesEmpty, "\n") << std::endl;
	if(!samplesExtracted.empty()) gzCatSamplesFile << njh::conToStr(samplesExtracted, "\n") << std::endl;
	if(!emptyFiles.empty()) pearExtractSamplesFile << njh::conToStr(emptyFiles, "\n") << std::endl;

	std::ofstream allMipsSamplesFile;
	openTextFile(allMipsSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "allMipsSamplesNames.tab.txt")));
	MipsSamplesNames goodSamples = *mipMaster.names_;
	goodSamples.setSamples(samplesExtracted);
	goodSamples.write(allMipsSamplesFile);
	bfs::copy_file(pars.mipsSamplesFile,njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "original_allMipsSamplesNames.tab.txt"));
	mipMaster.createPopClusMipDirs(pars.numThreads);

	if (nullptr != mipMaster.meta_) {
		bfs::copy_file(mipMaster.meta_->groupingsFile_,
				njh::files::make_path(mipMaster.directoryMaster_.masterDir_,
						"resources", "samplesMeta.tab.txt"));
	}

	return 0;
}

}  // namespace njhseq
