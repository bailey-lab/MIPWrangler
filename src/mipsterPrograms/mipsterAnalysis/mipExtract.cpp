/*
 * mipExtract.cpp
 *
 *  Created on: Feb 5, 2016
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

void extractFilterSampleForMips(const SeqIOOptions & sampleIOOpts,
		const SetUpMaster & mipMaster, aligner & alignerObjForFamilyDet,
		const mipIllumArmExtractionPars & pars, const QualFilteringPars & qFilPars,
		bool verbose) {
	njh::stopWatch watch;
	//std::cout << "On Thread: " << std::this_thread::get_id() << std::endl;
	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
	sampDirMaster.createExtractDirectory(pars.overWriteDirs);
	alignerObjForFamilyDet.resetAlnCache();
	alignerObjForFamilyDet.processAlnInfoInputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose);
	//set up sub directories
	bfs::path filteredOffDir = njh::files::makeDir(sampDirMaster.extractDir_.string(), njh::files::MkdirPar("filteredOff"));
	//create out files
	MultiSeqOutCache<seqInfo> mipOuts;
	mipOuts.setOpenLimit(pars.fileOpenLimit_);

	mipOuts.addReader("indeterminate",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "indeterminate").string(), sampleIOOpts.outFormat_,sampleIOOpts.out_));
	mipOuts.addReader("unmatchedReads",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "unmatchedReads").string(), sampleIOOpts.outFormat_,sampleIOOpts.out_));
	mipOuts.addReader("smallFragment",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "smallFragment").string(), sampleIOOpts.outFormat_,sampleIOOpts.out_));
	VecStr filterOutNames = { "_failedQuality", "_failedLigation",
			"_failedMinLen", "_containsNs" };
	VecStr allMipTargets = mipMaster.getAllMipTargets();
	for (const auto & mip : allMipTargets) {
		mipOuts.addReader(mip, SeqIOOptions(njh::files::join(VecStr {
			sampDirMaster.extractDir_.string(), mip, mip }), sampleIOOpts.outFormat_,sampleIOOpts.out_));
		for (const auto & outName : filterOutNames) {
			mipOuts.addReader(mip + outName,
					SeqIOOptions(
							njh::files::join(
									VecStr { sampDirMaster.extractDir_.string(), mip, mip
											+ outName }), sampleIOOpts.outFormat_,sampleIOOpts.out_));
		}
	}

	MipExtractionStats allExtractStats(pars.sampleName);
	//read in reads
	SeqIO readerOpt(sampleIOOpts);
	readerOpt.openIn();
	seqInfo read;
	uint32_t readCount = 1;
	while (readerOpt.readNextRead(read)) {
		if (readCount % 100 == 0 && verbose) {
			std::cout << "\r" << "currently on " << readCount;
			std::cout.flush();
		}
		++readCount;
		if (len(read) < pars.smallFragmentLength) {
			allExtractStats.increaseSmallFragment();
			mipOuts.add("smallFragment", read);
			continue;
		}
		//get length and resize aligner vector if needed
		uint64_t maxLen = alignerObjForFamilyDet.parts_.maxSize_;
		readVec::getMaxLength(read, maxLen);
		if (maxLen > alignerObjForFamilyDet.parts_.maxSize_) {
			alignerObjForFamilyDet.parts_.setMaxSize(maxLen);
		}
		bool found = false;
		std::unordered_map<std::string, std::pair<std::vector<Mip::ArmPosScore>, std::vector<Mip::ArmPosScore>>> possibleArms;
		std::unordered_map<std::string, std::vector<Mip::ArmPosScore>> possibleExtArms;

		for (const auto & mKey : allMipTargets) {
			const auto & mip = mipMaster.mips_->mips_.at(mKey);
			//check arm
			auto armPosMotif = mip.getPossibleExtArmPos(read);
			//if found arm and at the right position continue on;
			if (!armPosMotif.empty()) {
				found = true;
				possibleExtArms.emplace(mip.name_, armPosMotif);
				auto ligArmPosMotif = mip.getPossibleLigArmPos(read);
				if(!ligArmPosMotif.empty()){
					possibleArms.emplace(mip.name_, std::make_pair(armPosMotif, ligArmPosMotif));
				}
			}
		}
		if(possibleArms.empty()){
			if(possibleExtArms.empty()){
				//no matches found
				allExtractStats.increaseUnmatched();
				mipOuts.add("unmatchedReads", read);
			}else if(possibleExtArms.size() == 1){
				const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
				//quality control
				SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(read, qFilPars);

				std::string failedQaulifierName = MipExtractionStats::getNameForCase(
						eCase);
				//log and write read
				if(!allExtractStats.haveStatFor(mip.name_)){
					njh::files::makeDirP(sampDirMaster.extractDir_.string(), njh::files::MkdirPar(mip.name_));
				}
				mipOuts.add(mip.name_ + failedQaulifierName, read);
				allExtractStats.increaseCount(mip.name_, eCase);
			}else{
				const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
				auto currentMip = mipMaster.mips_->determineBestMipInFamily(read, mip, alignerObjForFamilyDet);
				//quality control
				SinlgeMipExtractInfo::extractCase eCase = currentMip.checkRead(read, qFilPars);

				std::string failedQaulifierName = MipExtractionStats::getNameForCase(
						eCase);
				//log and write read
				if(!allExtractStats.haveStatFor(currentMip.name_)){
					njh::files::makeDirP(sampDirMaster.extractDir_.string(), njh::files::MkdirPar(currentMip.name_));
				}
				mipOuts.add(currentMip.name_ + failedQaulifierName, read);
				allExtractStats.increaseCount(currentMip.name_, eCase);
			}

		}else if(1 == possibleArms.size()){
			//only one possible match
			const auto & mip = mipMaster.mips_->mips_.at(possibleArms.begin()->first);
			//quality control
			SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(read, qFilPars);

			std::string failedQaulifierName = MipExtractionStats::getNameForCase(
					eCase);
			//log and write read
			if(!allExtractStats.haveStatFor(mip.name_)){
				njh::files::makeDirP(sampDirMaster.extractDir_.string(), njh::files::MkdirPar(mip.name_));
			}
			mipOuts.add(mip.name_ + failedQaulifierName, read);
			allExtractStats.increaseCount(mip.name_, eCase);
		}else{
			//multiple possible matches
			uint32_t bestScore = 0;
			std::vector<std::string> bestMips;
			for(const auto & possible : possibleArms){
				uint32_t bestExtScore = 0;
				uint32_t bestLigScore = 0;
				for(const auto & ext : possible.second.first){
					if(ext.score_ > bestExtScore){
						bestExtScore = ext.score_;
					}
				}
				for(const auto & lig : possible.second.second){
					if(lig.score_ > bestLigScore){
						bestLigScore = lig.score_;
					}
				}
				uint32_t currentScore = bestLigScore + bestExtScore;
				if (currentScore > bestScore) {
					bestScore = currentScore;
					bestMips.clear();
					bestMips.emplace_back(possible.first);
				} else if (0 != currentScore && currentScore == bestScore) {
					bestMips.emplace_back(possible.first);
				}
			}
			if(bestMips.size() == 1){
				const auto & mip = mipMaster.mips_->mips_.at(bestMips.front());

				//quality control
				SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(read, qFilPars);

				std::string failedQaulifierName = MipExtractionStats::getNameForCase(
						eCase);
				//log and write read
				if(!allExtractStats.haveStatFor(mip.name_)){
					njh::files::makeDirP(sampDirMaster.extractDir_.string(), njh::files::MkdirPar(mip.name_));
				}
				mipOuts.add(mip.name_ + failedQaulifierName, read);
				allExtractStats.increaseCount(mip.name_, eCase);

			}else if(bestMips.size() > 1){
				//no matches found

				allExtractStats.increaseIndeterminate();
				std::stringstream appName;
				appName << "[";
				for(const auto & best : bestMips){
					const auto & possible = possibleArms.at(best);
					uint32_t bestExtScore = 0;
					uint32_t bestLigScore = 0;

					for(const auto & ext : possible.first){
						if(ext.score_ > bestExtScore){
							bestExtScore = ext.score_;
						}
					}
					for(const auto & lig : possible.second){
						if(lig.score_ > bestLigScore){
							bestLigScore = lig.score_;
						}
					}
					appName << "name=" << best << ";"
							<< "extScore=" << bestExtScore << ";"
							<< "ligScore=" << bestLigScore << ";";
				}
				appName << "]";
				read.name_.append(appName.str());
				mipOuts.add("indeterminate", read);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << std::endl;
				ss << "Best Mips vector is empty, which should be able to happen" << std::endl;
				ss << "For : " << sampleIOOpts.firstName_ << std::endl;
				throw std::runtime_error{ss.str()};
			}
		}
	}
	if (verbose) {
		std::cout << std::endl;
	}
	mipOuts.closeOutAll();
	VecStr extracOnlyColNames {"sampleName", "mipTarget", "mipFamily", "readNumber",
			"goodReads", "failedLigationArm", "failedMinLen" };
	if (qFilPars.checkingQFrac_) {
		extracOnlyColNames.emplace_back(
				"failed_q" + estd::to_string(qFilPars.qualCheck_) + "<"
						+ estd::to_string(qFilPars.qualCheckCutOff_));
	} else if (qFilPars.checkingQWindow) {
		extracOnlyColNames.emplace_back("failed_qW" + qFilPars.qualWindow_);
	} else {
		extracOnlyColNames.emplace_back("failed_quality(noneDone)");
	}
	extracOnlyColNames.emplace_back("containsNs");
	table infoTabByTarget(extracOnlyColNames);
	infoTabByTarget.content_ = allExtractStats.outputContents(*mipMaster.mips_, "\t");
	infoTabByTarget.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoByTarget.txt", ".txt"), "\t", infoTabByTarget.hasHeader_));
	std::ofstream sampLog;
	openTextFile(sampLog, sampDirMaster.extractDir_.string() + "log.txt", ".txt", false, true);
	sampLog << "Ran on: " << njh::getCurrentDate() << std::endl;
	sampLog << "Number of Alignments Done: "
			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	sampLog << "Run Time: " << watch.timeLapFormatted(6) << std::endl;
	alignerObjForFamilyDet.processAlnInfoOutputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose);
}

int mipsterAnalysisRunner::mipIllumExtractByArmAndFilter(
		const njh::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	mipIllumArmExtractionPars pars;
	setUp.setUpMipIllumArmExtraction(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (!mipMaster.names_->hasSample(pars.sampleName)) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", unrecognized sample name: "
				<< pars.sampleName << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	uint64_t maxLen = pars.minCaptureLength * 1.5;
	aligner alignerObjForFamilyDet = aligner(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_);
	//read in mip info and convert to mip class
	MipCollection mips(pars.mipArmsFileName, pars.allowableErrors);
	extractFilterSampleForMips(setUp.pars_.ioOptions_, mipMaster,
			alignerObjForFamilyDet, pars, setUp.pars_.qFilPars_,
			setUp.pars_.verbose_);

	setUp.rLog_ << "Number of Alignments Done: "
			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	return 0;
}

void extractMultiSamples(const SetUpMaster & mipMaster,
		const mipIllumArmExtractionParsMultiple& pars,
		const SeqSetUpPars & setUpPars, concurrent::AlignerPool & aligners,
		njh::concurrent::LockableQueue<std::string>& sampsQueue) {
	std::string sampleName = "";
	auto alignerObjForFamilyDet = aligners.popAligner();
	while (sampsQueue.getVal(sampleName)) {
		SeqIOOptions sampOpts = SeqIOOptions::genFastqIn(
				njh::files::make_path(mipMaster.directoryMaster_.masterDir_.string(), sampleName,
								sampleName + pars.seqFileSuffix).string());
		mipIllumArmExtractionPars samplePars = pars.createForSample(sampleName);
		extractFilterSampleForMips(sampOpts, mipMaster, *alignerObjForFamilyDet,
				samplePars, setUpPars.qFilPars_, setUpPars.verbose_);
	}
}

int mipsterAnalysisRunner::mipIllumExtractByArmAndFilterMultiple(
		const njh::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	mipIllumArmExtractionParsMultiple pars;
	setUp.setUpMipIllumArmExtractionMultiple(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	mipMaster.mips_->setAllMinCaptureLength(pars.minCaptureLength);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	std::ofstream logFile;
	openTextFile(logFile,
			njh::files::make_path(mipMaster.directoryMaster_.logsDir_, pars.logFilename),
			".txt", pars.overWriteLog, true);
	logFile << "Ran on: " << getCurrentDate() << std::endl;
	logFile << "Ran from: " << inputCommands.workingDir_ << std::endl;
	logFile << "Command: " << inputCommands.commandLine_ << std::endl;
	njh::concurrent::LockableQueue<std::string> sampsQueue(mipMaster.names_->samples_);
	uint64_t maxLen = pars.minCaptureLength * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();
	//read in mip info and convert to mip class
	std::vector<std::thread> threads;
	for (uint32_t threadNum = 0; threadNum < pars.numThreads; ++threadNum) {
		threads.emplace_back(extractMultiSamples, std::cref(mipMaster), std::cref(pars),
				std::cref(setUp.pars_), std::ref(aligners), std::ref(sampsQueue));
	}
	for (auto & t : threads) {
		t.join();
	}
	logFile << "Total Run Time: " << setUp.timer_.totalTimeFormatted(6) << std::endl;
	return 0;
}

}

