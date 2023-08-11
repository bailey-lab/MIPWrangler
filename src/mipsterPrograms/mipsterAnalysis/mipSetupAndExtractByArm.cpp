/*
 * mipSetupAndExtractByArm.cpp
 *
 *  Created on: Sep 5, 2018
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

int mipsterAnalysisRunner::mipSetupAndExtractByArm(const njh::progutils::CmdArgs & inputCommands) {
	extractFromRawParsMultiple pars;
	std::string mipServerName;
	bool runRest = false;
	int32_t stitchMatchScore = 2;
	int32_t stitchMismatchScore = -2;
	uint32_t stitchGapOpen = 10;
	uint32_t stitchGapExtend = 1;
	mipsterAnalysisSetUp setUp(inputCommands);

	setUp.setOption(stitchMatchScore, "--stitchMatchScore", "Match Score for stitching alignments");
	setUp.setOption(stitchMismatchScore, "--stitchMismatchScore", "Mismatch penalty for stitching alignments");
	setUp.setOption(stitchGapOpen, "--stitchGapOpen", "Gap opening penalty for stitching alignments, the penalty will be the negative value given here");
	setUp.setOption(stitchGapExtend, "--stitchGapExtend", "Gap extending penalty for stitching alignments, the penalty will be the negative value given here");


	setUp.setOption(runRest, "--runRest", "Run the rest of the analysis as well with defaults");
	setUp.setOption(mipServerName, "--mipServerNumber", "Name of the mip server, e.g. 1", true);
	setUp.setUpExtractFromRawMultiple(pars);

	//check for required external programs

	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setServerName(njh::pasteAsStr("mip", mipServerName));
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	if (!pars.sampleMetaFnp.empty()) {
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, njh::files::make_path(mipMaster.directoryMaster_.logsDir_,
			pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);
	setUp.startARunLog(njh::files::make_path(mipMaster.directoryMaster_.logsDir_).string());

	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for (const auto &f: files) {
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find('_'));
		if (njh::in(sampName, mipMaster.names_->samples_)) {
			readPairs[sampName].emplace_back(filename);
		} else {
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

	uint64_t maxLen = pars.minCaptureLength * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();
	/*
	 * 	pars_.generalMatch_,
			pars_.generalMismatch_,
			pars_.degenScoring_,
			pars_.lessNScoring_,
			pars_.caseInsensitiveScoring_
	 */
	auto scoringForStitching = substituteMatrix::createScoreMatrix(
			stitchMatchScore,
			stitchMismatchScore,
			false,
			true,
			true);
	aligner alignerObjForStitching(maxLen, gapScoringParameters(stitchGapOpen,stitchGapExtend,0,0,0,0), scoringForStitching, false);
	alignerObjForStitching.qScorePars_.qualThresWindow_ = 0;
	concurrent::AlignerPool alignersForStitching(alignerObjForStitching, pars.numThreads);

	alignersForStitching.initAligners();
//	{
//		auto alignerObjForStitching = alignersForStitching.popAligner();
//		std::cout << alignerObjForStitching->parts_.gapScores_.toJson() << std::endl;
//
//	}

	auto extractFiles =
			[&pars,&samplesExtracted,&samplesEmpty,&filesKeys,&emptyFiles,&aligners,&alignersForStitching,&logs,&logsMut,&mipMaster](
					const std::string & outputDirectory,
					const std::unordered_map<std::string, VecStr>& readPairs) {
				std::string key;
				VecStr currentEmptySamps;
				VecStr currentSamplesExtracted;
				std::vector<bfs::path> currentEmptyFiles;
				std::unordered_map<std::string, Json::Value> currentLogs;
				std::regex r1_pattern{R"((.*)_R1_(.*).fastq.gz$)"};
				std::regex r2_pattern{R"((.*)_R2_(.*).fastq.gz$)"};
				auto alignerObjForFamilyDet = aligners.popAligner();
				auto alignerObjForStitching = alignersForStitching.popAligner();
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
							if(njh::has(pairedFnpForSample, r1Name) && !pairedFnpForSample[r1Name].r1_fnp.empty() ){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r1Name].r1_fnp << " for " << f << " for " << r1Name << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								pairedFnpForSample[r1Name].r1_fnp = inputName;
							}
						}else if(matchesR2){
							std::string r2Name = njh::pasteAsStr(r2Match[1], r2Match[2]);
							if(njh::has(pairedFnpForSample, r2Name) && !pairedFnpForSample[r2Name].r2_fnp.empty() ){
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
					if(!pairedFnpForSample.empty() && success && !allEmpty){
						for(const auto & pairedFnp : pairedFnpForSample){
							if(pairedFnp.second.r1_fnp.empty()  ||
									pairedFnp.second.r2_fnp.empty()){
								std::stringstream ss;
								ss << __PRETTY_FUNCTION__ << ", error, " << key << ", " << pairedFnp.first << ", had an empty r1 or r2 file list" << "\n";
								ss << "R1:" << pairedFnp.second.r1_fnp << "\n";
								ss << "R2:" << pairedFnp.second.r2_fnp << "\n";
								throw std::runtime_error{ss.str()};
							}else{
								sampleOpts.emplace_back(SeqIOOptions::genPairedInGz(pairedFnp.second.r1_fnp, pairedFnp.second.r2_fnp));
								sampleOpts.back().outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
								sampleOpts.back().out_.outExtention_ = "_R1.fastq.gz";
								sampleOpts.back().revComplMate_ = true;
							}
						}
					}
					if(success){
						currentSamplesExtracted.emplace_back(key);
						auto samplePars = pars.createForSample(key);
						MipExtractor mExtractor(pars.verbose_);
						mExtractor.extractFilterSampleForMipsPairedStitch(sampleOpts,
								mipMaster,
								*alignerObjForFamilyDet,
								*alignerObjForStitching,
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
				extractFiles,
						std::cref(mipMaster.directoryMaster_.masterDir_.string()),
						std::cref(readPairs));
	}
	njh::concurrent::joinAllThreads(threads);
	logFile << logs << std::endl;

	//copy over resources;
	bfs::copy_file(pars.mipArmsFileName, njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
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


	// set up scripts

	// run barcode correction, mipBarcodeCorrectionMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipBarcodeCorrectionMultiple --masterDir "
				<< njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
		<< " --numThreads " << pars.numThreads
		<< " --logFile mipBarcodeCorrecting_run1";
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run clustering step, mipClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipClusteringMultiple --masterDir "
				<< njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
		<< " --numThreads " << pars.numThreads
		<< " --logFile mipClustering_run1" ;
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run population clustering, mipPopulationClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipPopulationClusteringMultiple --masterDir "
				<< njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
		<< " --numThreads " << pars.numThreads
		<< " --logFile mipPopClustering_run1";
		if(!pars.refDir.empty()){
			mipScriptOut << " --refDir " << njh::files::normalize(pars.refDir);
		}
		if (pars.keepIntermediateFiles) {
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run setup for viewer, mipAnalysisServerSetUp
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then 
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipAnalysisServerSetUp --masterDir "
				<< njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
		<< " --numThreads " << pars.numThreads
		<< " --name mip$1 --verbose";
		mipScriptOut << std::endl;
	}
	// run viewer, mav
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_viewer.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then 
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << "nohup " << setUp.commands_.masterProgramRaw_ << " mav --masterDir "
				<< njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
		<< " --numThreads " << pars.numThreads
		<< " --port $((10000+$1)) --name mip$1 --verbose &";
		mipScriptOut << std::endl;
	}

	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh")) << " " << mipServerName << std::endl;
		mipScriptOut << std::endl;
	}

	if(runRest){
		std::stringstream cmdSs;
		cmdSs << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh")) << " " << mipServerName << std::endl;
		auto runLog = njh::sys::run(VecStr{cmdSs.str()});
		OutOptions restOfAnalysisRunLogOpts(njh::files::make_path(mipMaster.directoryMaster_.logsDir_, "run_restOfAnalysis_log.json"));
		OutputStream restOfAnalysisRunLogOut(restOfAnalysisRunLogOpts);
		restOfAnalysisRunLogOut << runLog.toJson() << std::endl;
	}



	return 0;
}

int mipsterAnalysisRunner::mipSetup(const njh::progutils::CmdArgs & inputCommands) {
	extractFromRawParsMultiple pars;
	std::string mipServerName;
	mipsterAnalysisSetUp setUp(inputCommands);
	setUp.setOption(mipServerName, "--mipServerNumber", "Name of the mip server, e.g. 1", true);
	setUp.setUpExtractFromRawMultiple(pars);

	//check for required external programs

	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setServerName(njh::pasteAsStr("mip", mipServerName));
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	if (!pars.sampleMetaFnp.empty()) {
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
						<< "Error in directory structure, make sure you are in the correct analysis directory"
						<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	setUp.startARunLog(njh::files::make_path(mipMaster.directoryMaster_.logsDir_).string());

	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	std::vector<bfs::path> emptyFiles;
	for (const auto &f: files) {
		auto filename = f.first.filename().string();
		if (njh::files::isFileEmpty(f.first)) {
			emptyFiles.emplace_back(f.first);
		} else {
			std::string sampName = filename.substr(0, filename.find('_'));
			if (njh::in(sampName, mipMaster.names_->samples_)) {
				readPairs[sampName].emplace_back(filename);
			} else {
				readPairsUnrecognized[sampName].emplace_back(filename);
			}
		}
	}
	auto keys = getVectorOfMapKeys(readPairs);
	njh::sort(keys);

	{
		auto sampleInputFilesFilesOpts = OutOptions(
						njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "resources/sampleInputFiles.tab.txt"));
		OutputStream sampleInputFilesFilesOut(sampleInputFilesFilesOpts);
		sampleInputFilesFilesOut << "sample\tfiles" << std::endl;
		for (const auto &inputFnps: readPairs) {
			sampleInputFilesFilesOut << inputFnps.first << "\t"
															 << njh::conToStr(inputFnps.second, ",") << std::endl;
		}

		auto unrecognizedSampleInputFilesOpts = OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "resources/unrecognizedSampleInputFiles.tab.txt"));
		OutputStream unrecognizedSampleInputFilesOut(unrecognizedSampleInputFilesOpts);

		unrecognizedSampleInputFilesOut << "sample\tfiles" << std::endl;
		for (const auto &inputFnps: readPairsUnrecognized) {
			unrecognizedSampleInputFilesOut << inputFnps.first << "\t"
															 << njh::conToStr(inputFnps.second, ",") << std::endl;
		}
	}

	Json::Value logs;
	logs["mainCommand"] = setUp.commands_.commandLine_;
	std::mutex logsMut;



	//copy over resources;
	bfs::copy_file(pars.mipArmsFileName, njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
	OutputStream allMipsSamplesFile(OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "allMipsSamplesNames.tab.txt")));
	MipsSamplesNames goodSamples = *mipMaster.names_;
	goodSamples.setSamples(keys);
	goodSamples.write(allMipsSamplesFile);
	bfs::copy_file(pars.mipsSamplesFile,njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "original_allMipsSamplesNames.tab.txt"));
	mipMaster.createPopClusMipDirs(pars.numThreads);

	if (nullptr != mipMaster.meta_) {
		bfs::copy_file(mipMaster.meta_->groupingsFile_,
									 njh::files::make_path(mipMaster.directoryMaster_.masterDir_,
																				 "resources", "samplesMeta.tab.txt"));
	}
	// set up scripts

	// run barcode correction, mipBarcodeCorrectionMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipBarcodeCorrectionMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipBarcodeCorrecting_run1";
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run clustering step, mipClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipClustering_run1" ;
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run population clustering, mipPopulationClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipPopulationClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipPopClustering_run1";
		if(!pars.refDir.empty()){
			mipScriptOut << " --refDir " << njh::files::normalize(pars.refDir);
		}
		if (pars.keepIntermediateFiles) {
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run setup for viewer, mipAnalysisServerSetUp
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipAnalysisServerSetUp --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --name mip$1 --verbose";
		mipScriptOut << std::endl;
	}
	// run viewer, mav
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_viewer.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << "nohup " << setUp.commands_.masterProgramRaw_ << " mav --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --port $((10000+$1)) --name mip$1 --verbose &";
		mipScriptOut << std::endl;
	}

	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh")) << " " << mipServerName << std::endl;
		mipScriptOut << std::endl;
	}
	return 0;
}

int mipsterAnalysisRunner::mipExtractByArm(const njh::progutils::CmdArgs & inputCommands) {
	extractFromRawParsMultiple pars;
	std::string mipServerName;
	bool runRest = false;
	int32_t stitchMatchScore = 2;
	int32_t stitchMismatchScore = -2;
	uint32_t stitchGapOpen = 10;
	uint32_t stitchGapExtend = 1;
	mipsterAnalysisSetUp setUp(inputCommands);

	setUp.setOption(stitchMatchScore, "--stitchMatchScore", "Match Score for stitching alignments");
	setUp.setOption(stitchMismatchScore, "--stitchMismatchScore", "Mismatch penalty for stitching alignments");
	setUp.setOption(stitchGapOpen, "--stitchGapOpen", "Gap opening penalty for stitching alignments, the penalty will be the negative value given here");
	setUp.setOption(stitchGapExtend, "--stitchGapExtend", "Gap extending penalty for stitching alignments, the penalty will be the negative value given here");


	setUp.setOption(runRest, "--runRest", "Run the rest of the analysis as well with defaults");
	setUp.setOption(mipServerName, "--mipServerNumber", "Name of the mip server, e.g. 1", true);
	setUp.setUpExtractFromRawMultiple(pars);

	//check for required external programs

	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setServerName(njh::pasteAsStr("mip", mipServerName));
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	if (!pars.sampleMetaFnp.empty()) {
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
						<< "Error in directory structure, make sure you are in the correct analysis directory"
						<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, njh::files::make_path(mipMaster.directoryMaster_.logsDir_,
																							pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);
	setUp.startARunLog(njh::files::make_path(mipMaster.directoryMaster_.logsDir_).string());

	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for(const auto & f : files){
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find('_'));
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

	uint64_t maxLen = pars.minCaptureLength * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
																	 setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();
	/*
	 * 	pars_.generalMatch_,
			pars_.generalMismatch_,
			pars_.degenScoring_,
			pars_.lessNScoring_,
			pars_.caseInsensitiveScoring_
	 */
	auto scoringForStitching = substituteMatrix::createScoreMatrix(
					stitchMatchScore,
					stitchMismatchScore,
					false,
					true,
					true);
	aligner alignerObjForStitching(maxLen, gapScoringParameters(stitchGapOpen,stitchGapExtend,0,0,0,0), scoringForStitching, false);
	alignerObjForStitching.qScorePars_.qualThresWindow_ = 0;
	concurrent::AlignerPool alignersForStitching(alignerObjForStitching, pars.numThreads);

	alignersForStitching.initAligners();
//	{
//		auto alignerObjForStitching = alignersForStitching.popAligner();
//		std::cout << alignerObjForStitching->parts_.gapScores_.toJson() << std::endl;
//
//	}

	auto extractFiles =
					[&pars,&samplesExtracted,&samplesEmpty,&filesKeys,&emptyFiles,&aligners,&alignersForStitching,&logs,&logsMut,&mipMaster](
									const std::string & outputDirectory,
									const std::unordered_map<std::string, VecStr>& readPairs) {
						std::string key;
						VecStr currentEmptySamps;
						VecStr currentSamplesExtracted;
						std::vector<bfs::path> currentEmptyFiles;
						std::unordered_map<std::string, Json::Value> currentLogs;
						std::regex r1_pattern{R"((.*)_R1_(.*).fastq.gz$)"};
						std::regex r2_pattern{R"((.*)_R2_(.*).fastq.gz$)"};
						auto alignerObjForFamilyDet = aligners.popAligner();
						auto alignerObjForStitching = alignersForStitching.popAligner();
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
									if(njh::has(pairedFnpForSample, r1Name) && !pairedFnpForSample[r1Name].r1_fnp.empty() ){
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r1Name].r1_fnp << " for " << f << " for " << r1Name << "\n";
										throw std::runtime_error{ss.str()};
									}else{
										pairedFnpForSample[r1Name].r1_fnp = inputName;
									}
								}else if(matchesR2){
									std::string r2Name = njh::pasteAsStr(r2Match[1], r2Match[2]);
									if(njh::has(pairedFnpForSample, r2Name) && !pairedFnpForSample[r2Name].r2_fnp.empty() ){
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
							if(!pairedFnpForSample.empty() && success && !allEmpty){
								for(const auto & pairedFnp : pairedFnpForSample){
									if(pairedFnp.second.r1_fnp.empty()  ||
										 pairedFnp.second.r2_fnp.empty()){
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error, " << key << ", " << pairedFnp.first << ", had an empty r1 or r2 file list" << "\n";
										ss << "R1:" << pairedFnp.second.r1_fnp << "\n";
										ss << "R2:" << pairedFnp.second.r2_fnp << "\n";
										throw std::runtime_error{ss.str()};
									}else{
										sampleOpts.emplace_back(SeqIOOptions::genPairedInGz(pairedFnp.second.r1_fnp, pairedFnp.second.r2_fnp));
										sampleOpts.back().outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
										sampleOpts.back().out_.outExtention_ = "_R1.fastq.gz";
										sampleOpts.back().revComplMate_ = true;
									}
								}
							}
							if(success){
								currentSamplesExtracted.emplace_back(key);
								auto samplePars = pars.createForSample(key);
								MipExtractor mExtractor(pars.verbose_);
								mExtractor.extractFilterSampleForMipsPairedStitch(sampleOpts,
																																	mipMaster,
																																	*alignerObjForFamilyDet,
																																	*alignerObjForStitching,
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
						extractFiles,
						std::cref(mipMaster.directoryMaster_.masterDir_.string()),
						std::cref(readPairs));
	}
	njh::concurrent::joinAllThreads(threads);
	logFile << logs << std::endl;

	//copy over resources;
	bfs::copy_file(pars.mipArmsFileName, njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
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


	// set up scripts

	// run barcode correction, mipBarcodeCorrectionMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipBarcodeCorrectionMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipBarcodeCorrecting_run1";
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run clustering step, mipClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipClustering_run1" ;
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run population clustering, mipPopulationClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipPopulationClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipPopClustering_run1";
		if(!pars.refDir.empty()){
			mipScriptOut << " --refDir " << njh::files::normalize(pars.refDir);
		}
		if (pars.keepIntermediateFiles) {
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run setup for viewer, mipAnalysisServerSetUp
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipAnalysisServerSetUp --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --name mip$1 --verbose";
		mipScriptOut << std::endl;
	}
	// run viewer, mav
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_viewer.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << "nohup " << setUp.commands_.masterProgramRaw_ << " mav --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --port $((10000+$1)) --name mip$1 --verbose &";
		mipScriptOut << std::endl;
	}

	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh")) << " " << mipServerName << std::endl;
		mipScriptOut << std::endl;
	}

	if(runRest){
		std::stringstream cmdSs;
		cmdSs << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh")) << " " << mipServerName << std::endl;
		auto runLog = njh::sys::run(VecStr{cmdSs.str()});
		OutOptions restOfAnalysisRunLogOpts(njh::files::make_path(mipMaster.directoryMaster_.logsDir_, "run_restOfAnalysis_log.json"));
		OutputStream restOfAnalysisRunLogOut(restOfAnalysisRunLogOpts);
		restOfAnalysisRunLogOut << runLog.toJson() << std::endl;
	}



	return 0;
}

int mipsterAnalysisRunner::mipExtractByArmMultiple(const njh::progutils::CmdArgs & inputCommands) {
	extractFromRawParsMultiple pars;
	std::string mipServerName;
	bool runRest = false;
	int32_t stitchMatchScore = 2;
	int32_t stitchMismatchScore = -2;
	uint32_t stitchGapOpen = 10;
	uint32_t stitchGapExtend = 1;
	mipsterAnalysisSetUp setUp(inputCommands);

	setUp.setOption(stitchMatchScore, "--stitchMatchScore", "Match Score for stitching alignments");
	setUp.setOption(stitchMismatchScore, "--stitchMismatchScore", "Mismatch penalty for stitching alignments");
	setUp.setOption(stitchGapOpen, "--stitchGapOpen", "Gap opening penalty for stitching alignments, the penalty will be the negative value given here");
	setUp.setOption(stitchGapExtend, "--stitchGapExtend", "Gap extending penalty for stitching alignments, the penalty will be the negative value given here");


	setUp.setOption(runRest, "--runRest", "Run the rest of the analysis as well with defaults");
	setUp.setOption(mipServerName, "--mipServerNumber", "Name of the mip server, e.g. 1", true);
	setUp.setUpExtractFromRawMultiple(pars);

	//check for required external programs

	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setServerName(njh::pasteAsStr("mip", mipServerName));
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	if (!pars.sampleMetaFnp.empty()) {
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
						<< "Error in directory structure, make sure you are in the correct analysis directory"
						<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, njh::files::make_path(mipMaster.directoryMaster_.logsDir_,
																							pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);
	setUp.startARunLog(njh::files::make_path(mipMaster.directoryMaster_.logsDir_).string());

	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for(const auto & f : files){
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find('_'));
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

	uint64_t maxLen = pars.minCaptureLength * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
																	 setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();
	/*
	 * 	pars_.generalMatch_,
			pars_.generalMismatch_,
			pars_.degenScoring_,
			pars_.lessNScoring_,
			pars_.caseInsensitiveScoring_
	 */
	auto scoringForStitching = substituteMatrix::createScoreMatrix(
					stitchMatchScore,
					stitchMismatchScore,
					false,
					true,
					true);
	aligner alignerObjForStitching(maxLen, gapScoringParameters(stitchGapOpen,stitchGapExtend,0,0,0,0), scoringForStitching, false);
	alignerObjForStitching.qScorePars_.qualThresWindow_ = 0;
	concurrent::AlignerPool alignersForStitching(alignerObjForStitching, pars.numThreads);

	alignersForStitching.initAligners();
//	{
//		auto alignerObjForStitching = alignersForStitching.popAligner();
//		std::cout << alignerObjForStitching->parts_.gapScores_.toJson() << std::endl;
//
//	}

	auto extractFiles =
					[&pars,&samplesExtracted,&samplesEmpty,&filesKeys,&emptyFiles,&aligners,&alignersForStitching,&logs,&logsMut,&mipMaster](
									const std::string & outputDirectory,
									const std::unordered_map<std::string, VecStr>& readPairs) {
						std::string key;
						VecStr currentEmptySamps;
						VecStr currentSamplesExtracted;
						std::vector<bfs::path> currentEmptyFiles;
						std::unordered_map<std::string, Json::Value> currentLogs;
						std::regex r1_pattern{R"((.*)_R1_(.*).fastq.gz$)"};
						std::regex r2_pattern{R"((.*)_R2_(.*).fastq.gz$)"};
						auto alignerObjForFamilyDet = aligners.popAligner();
						auto alignerObjForStitching = alignersForStitching.popAligner();
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
									if(njh::has(pairedFnpForSample, r1Name) && !pairedFnpForSample[r1Name].r1_fnp.empty() ){
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error, already have " <<  pairedFnpForSample[r1Name].r1_fnp << " for " << f << " for " << r1Name << "\n";
										throw std::runtime_error{ss.str()};
									}else{
										pairedFnpForSample[r1Name].r1_fnp = inputName;
									}
								}else if(matchesR2){
									std::string r2Name = njh::pasteAsStr(r2Match[1], r2Match[2]);
									if(njh::has(pairedFnpForSample, r2Name) && !pairedFnpForSample[r2Name].r2_fnp.empty() ){
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
							if(!pairedFnpForSample.empty() && success && !allEmpty){
								for(const auto & pairedFnp : pairedFnpForSample){
									if(pairedFnp.second.r1_fnp.empty()  ||
										 pairedFnp.second.r2_fnp.empty()){
										std::stringstream ss;
										ss << __PRETTY_FUNCTION__ << ", error, " << key << ", " << pairedFnp.first << ", had an empty r1 or r2 file list" << "\n";
										ss << "R1:" << pairedFnp.second.r1_fnp << "\n";
										ss << "R2:" << pairedFnp.second.r2_fnp << "\n";
										throw std::runtime_error{ss.str()};
									}else{
										sampleOpts.emplace_back(SeqIOOptions::genPairedInGz(pairedFnp.second.r1_fnp, pairedFnp.second.r2_fnp));
										sampleOpts.back().outFormat_ = SeqIOOptions::outFormats::FASTQPAIREDGZ;
										sampleOpts.back().out_.outExtention_ = "_R1.fastq.gz";
										sampleOpts.back().revComplMate_ = true;
									}
								}
							}
							if(success){
								currentSamplesExtracted.emplace_back(key);
								auto samplePars = pars.createForSample(key);
								MipExtractor mExtractor(pars.verbose_);
								mExtractor.extractFilterSampleForMipsPairedStitch(sampleOpts,
																																	mipMaster,
																																	*alignerObjForFamilyDet,
																																	*alignerObjForStitching,
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
						extractFiles,
						std::cref(mipMaster.directoryMaster_.masterDir_.string()),
						std::cref(readPairs));
	}
	njh::concurrent::joinAllThreads(threads);
	logFile << logs << std::endl;

	//copy over resources;
	bfs::copy_file(pars.mipArmsFileName, njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
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


	// set up scripts

	// run barcode correction, mipBarcodeCorrectionMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipBarcodeCorrectionMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipBarcodeCorrecting_run1";
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run clustering step, mipClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipClustering_run1" ;
		if(pars.keepIntermediateFiles){
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run population clustering, mipPopulationClusteringMultiple
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipPopulationClusteringMultiple --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --logFile mipPopClustering_run1";
		if(!pars.refDir.empty()){
			mipScriptOut << " --refDir " << njh::files::normalize(pars.refDir);
		}
		if (pars.keepIntermediateFiles) {
			mipScriptOut << " --keepIntermediateFiles";
		}
		mipScriptOut << std::endl;
	}
	// run setup for viewer, mipAnalysisServerSetUp
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << setUp.commands_.masterProgramRaw_ << " mipAnalysisServerSetUp --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --name mip$1 --verbose";
		mipScriptOut << std::endl;
	}
	// run viewer, mav
	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_viewer.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << R"(if [[ $# -ne 1 ]]; then
    echo "Illegal number of parameters, needs 1 argument, 1) name of mip server number"
    exit
fi
)";

		mipScriptOut << "nohup " << setUp.commands_.masterProgramRaw_ << " mav --masterDir "
								 << njh::files::normalize(mipMaster.directoryMaster_.masterDir_)
								 << " --numThreads " << pars.numThreads
								 << " --port $((10000+$1)) --name mip$1 --verbose &";
		mipScriptOut << std::endl;
	}

	{
		OutOptions mipScriptOpts(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh"));
		std::ofstream mipScriptOut;
		mipScriptOpts.openExecutableFile(mipScriptOut);
		mipScriptOut << "#!/usr/bin/env bash" << "\n";
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipBarcodeCorrectionMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipPopulationClusteringMultiple.sh")) << std::endl;
		mipScriptOut << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_mipAnalysisServerSetUp.sh")) << " " << mipServerName << std::endl;
		mipScriptOut << std::endl;
	}

	if(runRest){
		std::stringstream cmdSs;
		cmdSs << njh::files::normalize(njh::files::make_path(mipMaster.directoryMaster_.scriptsDir_, "run_restOfAnalysis.sh")) << " " << mipServerName << std::endl;
		auto runLog = njh::sys::run(VecStr{cmdSs.str()});
		OutOptions restOfAnalysisRunLogOpts(njh::files::make_path(mipMaster.directoryMaster_.logsDir_, "run_restOfAnalysis_log.json"));
		OutputStream restOfAnalysisRunLogOut(restOfAnalysisRunLogOpts);
		restOfAnalysisRunLogOut << runLog.toJson() << std::endl;
	}



	return 0;
}

//


}  // namespace njhseq
