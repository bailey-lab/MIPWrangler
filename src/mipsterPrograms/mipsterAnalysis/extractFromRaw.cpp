/*
 * extractFromRaw.cpp
 *
 *  Created on: Aug 9, 2017
 *      Author: nick
 */
#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"


namespace bibseq {

int mipsterAnalysisRunner::extractFromRaw(const bib::progutils::CmdArgs & inputCommands) {
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
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	pars.dir = bib::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << bib::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, bib::files::make_path( mipMaster.directoryMaster_.logsDir_,
			pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);


	auto files = bib::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for(const auto & f : files){
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find("_"));
		if(bib::in(sampName, mipMaster.names_->samples_)){
			readPairs[sampName].emplace_back(filename);
		}else{
			readPairsUnrecognized[sampName].emplace_back(filename);
		}
	}
	auto keys = getVectorOfMapKeys(readPairs);
	bib::sort(keys);
	std::vector<bfs::path> emptyFiles;
	VecStr samplesExtracted;
	VecStr samplesEmpty;
	bib::concurrent::LockableQueue<std::string> filesKeys(keys);
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
					bib::stopWatch watch;
					std::stringstream out;
					std::stringstream err;
					bool success = true;
					bool allEmpty = true;
					out << key << std::endl;
					auto sampDir = bib::files::makeDir(outputDirectory, bib::files::MkdirPar(key));
					std::vector<bfs::path> r1_files;
					std::vector<bfs::path> r2_files;
					struct PairedFnps{
						bfs::path r1_fnp = "";
						bfs::path r2_fnp = "";
					};
					std::unordered_map<std::string, PairedFnps> pairedFnpForSample;
					for(const auto & f : readPairs.at(key)){
						auto fqName = bfs::path(f);
						auto inputName = bib::files::make_path(pars.dir,fqName);
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
							std::string r1Name = bib::pasteAsStr(r1Match[1], r1Match[2]);
							if(bib::has(pairedFnpForSample, r1Name) && "" != pairedFnpForSample[r1Name].r1_fnp ){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r1Name].r1_fnp << " for " << f << " for " << r1Name << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								pairedFnpForSample[r1Name].r1_fnp = inputName;
							}
						}else if(matchesR2){
							std::string r2Name = bib::pasteAsStr(r2Match[1], r2Match[2]);
							if(bib::has(pairedFnpForSample, r2Name) && "" != pairedFnpForSample[r2Name].r2_fnp ){
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
						mExtractor.extractFilterSampleForMipsPaired(sampleOpts, mipMaster, *alignerObjForFamilyDet,
																samplePars);
					}
					log["sample"] = key;
					log["time"] = watch.totalTime();
					log["out"] = out.str();
					log["err"] = err.str();
					log["success"] = bib::json::toJson(success);
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

	bib::concurrent::joinAllThreads(threads);

	logFile << logs << std::endl;


	//copy over resources;

	bfs::copy(pars.mipArmsFileName,bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
	bib::files::makeDirP(bib::files::MkdirPar(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo")));
	std::ofstream outSamplesFoundFile;
	openTextFile(outSamplesFoundFile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesFound.tab.txt")));
	std::ofstream outSamplesEmptyfile;
	openTextFile(outSamplesEmptyfile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesEmpty.tab.txt")));
	std::ofstream gzCatSamplesFile;
	openTextFile(gzCatSamplesFile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/gzCatSamples.tab.txt")));
	std::ofstream pearExtractSamplesFile;
	openTextFile(pearExtractSamplesFile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/emptyFiles.tab.txt")));
	std::ofstream samplesMissingFile;
	openTextFile(samplesMissingFile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/samplesMissing.tab.txt")));
	auto allFoundSamps = concatVecs(samplesExtracted, samplesEmpty);
	bib::sort(allFoundSamps);
	bib::sort(samplesEmpty);
	bib::sort(samplesExtracted);
	bib::sort(emptyFiles);
	if(!allFoundSamps.empty()) outSamplesFoundFile << bib::conToStr(allFoundSamps, "\n") << std::endl;
	if(!samplesEmpty.empty()) outSamplesEmptyfile << bib::conToStr(samplesEmpty, "\n") << std::endl;
	if(!samplesExtracted.empty()) gzCatSamplesFile << bib::conToStr(samplesExtracted, "\n") << std::endl;
	if(!emptyFiles.empty()) pearExtractSamplesFile << bib::conToStr(emptyFiles, "\n") << std::endl;

	std::ofstream allMipsSamplesFile;
	openTextFile(allMipsSamplesFile,OutOptions(bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "allMipsSamplesNames.tab.txt")));
	MipsSamplesNames goodSamples = *mipMaster.names_;
	goodSamples.setSamples(samplesExtracted);
	goodSamples.write(allMipsSamplesFile);
	bfs::copy(pars.mipsSamplesFile,bib::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "original_allMipsSamplesNames.tab.txt"));
	mipMaster.createPopClusMipDirs(pars.numThreads);

	if (nullptr != mipMaster.meta_) {
		bfs::copy_file(mipMaster.meta_->groupingsFile_,
				bib::files::make_path(mipMaster.directoryMaster_.masterDir_,
						"resources", "samplesMeta.tab.txt"));
	}

	return 0;
}

}  // namespace bibseq
