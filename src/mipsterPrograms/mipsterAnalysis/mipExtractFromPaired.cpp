/*
 * mipExtractFromPaired.cpp
 *
 *  Created on: May 10, 2017
 *      Author: nick
 */


#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"

namespace bibseq {



//void extractFilterSampleForMipsPaired(const std::vector<SeqIOOptions> & sampleIOOpts,
//		const SetUpMaster & mipMaster,
//		aligner & alignerObjForFamilyDet,
//		const mipIllumArmExtractionPars & pars,
//		const QualFilteringPars & qFilPars,
//		bool verbose) {
//
//	bib::stopWatch watch;
//	//std::cout << "On Thread: " << std::this_thread::get_id() << std::endl;
//	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
//	sampDirMaster.createExtractDirectory(pars.overWriteDirs);
//	alignerObjForFamilyDet.resetAlnCache();
//	alignerObjForFamilyDet.processAlnInfoInputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose);
//	//set up sub directories
//	bfs::path filteredOffDir = bib::files::makeDir(sampDirMaster.extractDir_.string(), bib::files::MkdirPar("filteredOff"));
//	//create out files
//	MultiSeqOutCache<PairedRead> mipOuts;
//	mipOuts.setOpenLimit(pars.fileOpenLimit_);
//
//	mipOuts.addReader("indeterminate",
//			SeqIOOptions(bib::files::make_path(filteredOffDir, "indeterminate").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
//	mipOuts.addReader("unmatchedReads",
//			SeqIOOptions(bib::files::make_path(filteredOffDir, "unmatchedReads").string(), sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
//	mipOuts.addReader("smallFragment",
//			SeqIOOptions(bib::files::make_path(filteredOffDir, "smallFragment").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
//	VecStr filterOutNames = { "_failedQuality", "_failedLigation",
//			"_failedMinLen", "_containsNs" };
//	VecStr allMipTargets = mipMaster.getAllMipTargets();
//	for (const auto & mip : allMipTargets) {
//		mipOuts.addReader(mip, SeqIOOptions(bib::files::join(VecStr {
//			sampDirMaster.extractDir_.string(), mip, mip }), sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
//		for (const auto & outName : filterOutNames) {
//			mipOuts.addReader(mip + outName,
//					SeqIOOptions(
//							bib::files::join(
//									VecStr { sampDirMaster.extractDir_.string(), mip, mip
//											+ outName }), sampleIOOpts.front().outFormat_,
//											sampleIOOpts.front().out_) );
//		}
//	}
//
//	MipExtractionStats allExtractStats(pars.sampleName);
//	for(const auto & sampleIOOpt : sampleIOOpts){
//		//read in reads
//		SeqIO readerOpt(sampleIOOpt);
//		readerOpt.openIn();
//		PairedRead seq;
//		uint32_t readCount = 1;
//		while (readerOpt.readNextRead(seq)) {
//			if (readCount % 100 == 0 && verbose) {
//				std::cout << "\r" << "currently on " << readCount;
//				std::cout.flush();
//			}
//			++readCount;
//			//get length and resize aligner vector if needed
//			uint64_t maxLen = alignerObjForFamilyDet.parts_.maxSize_;
//			readVec::getMaxLength(seq, maxLen);
//			if (maxLen > alignerObjForFamilyDet.parts_.maxSize_) {
//				alignerObjForFamilyDet.parts_.setMaxSize(maxLen);
//			}
//			bool found = false;
//			std::unordered_map<std::string, std::pair<std::vector<Mip::ArmPosScore>, std::vector<Mip::ArmPosScore>>> possibleArms;
//			std::unordered_map<std::string, std::vector<Mip::ArmPosScore>> possibleExtArms;
//
//	//		allExtractStats.increaseUnmatched();
//	//		mipOuts.add("unmatchedReads", seq);
//	//
//			for (const auto & mKey : allMipTargets) {
//				const auto & mip = mipMaster.mips_->mips_.at(mKey);
//				//check arm
//				auto armPosMotif = mip.getPossibleExtArmPos(seq.seqBase_);
//				//if found arm and at the right position continue on;
//				if (!armPosMotif.empty()) {
//					found = true;
//					possibleExtArms.emplace(mip.name_, armPosMotif);
//					auto ligArmPosMotif = mip.getPossibleLigArmPos(seq.mateSeqBase_);
//					if(!ligArmPosMotif.empty()){
//						possibleArms.emplace(mip.name_, std::make_pair(armPosMotif, ligArmPosMotif));
//					}
//				}
//			}
//			if(possibleArms.empty()){
//				if(possibleExtArms.empty()){
//					//no matches found
//					allExtractStats.increaseUnmatched();
//					mipOuts.add("unmatchedReads", seq);
//
//				}else if(possibleExtArms.size() == 1){
//					const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
//					//quality control
//					SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, qFilPars);
//
//					std::string failedQaulifierName = MipExtractionStats::getNameForCase(
//							eCase);
//					//log and write read
//					if(!allExtractStats.haveStatFor(mip.name_)){
//						bib::files::makeDirP(sampDirMaster.extractDir_.string(), bib::files::MkdirPar(mip.name_));
//					}
//					mipOuts.add(mip.name_ + failedQaulifierName, seq);
//					//mipOuts.add("unmatchedReads", seq);
//					allExtractStats.increaseCount(mip.name_, eCase);
//				}else{
//					const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
//					auto currentMip = mipMaster.mips_->determineBestMipInFamily(seq, mip, alignerObjForFamilyDet);
//					//quality control
//					SinlgeMipExtractInfo::extractCase eCase = currentMip.checkRead(seq, qFilPars);
//
//					std::string failedQaulifierName = MipExtractionStats::getNameForCase(
//							eCase);
//					//log and write read
//					if(!allExtractStats.haveStatFor(currentMip.name_)){
//						bib::files::makeDirP(sampDirMaster.extractDir_.string(), bib::files::MkdirPar(currentMip.name_));
//					}
//					mipOuts.add(currentMip.name_ + failedQaulifierName, seq);
//					//mipOuts.add("unmatchedReads", seq);
//					allExtractStats.increaseCount(currentMip.name_, eCase);
//				}
//			}else if(1 == possibleArms.size()){
//				//only one possible match
//				const auto & mip = mipMaster.mips_->mips_.at(possibleArms.begin()->first);
//				//quality control
//				SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, qFilPars);
//
//				std::string failedQaulifierName = MipExtractionStats::getNameForCase(
//						eCase);
//				//log and write read
//				if(!allExtractStats.haveStatFor(mip.name_)){
//					bib::files::makeDirP(sampDirMaster.extractDir_.string(), bib::files::MkdirPar(mip.name_));
//				}
//				mipOuts.add(mip.name_ + failedQaulifierName, seq);
//				//mipOuts.add("unmatchedReads", seq);
//				allExtractStats.increaseCount(mip.name_, eCase);
//			}else{
//				//multiple possible matches
//				uint32_t bestScore = 0;
//				std::vector<std::string> bestMips;
//				for(const auto & possible : possibleArms){
//					uint32_t bestExtScore = 0;
//					uint32_t bestLigScore = 0;
//					for(const auto & ext : possible.second.first){
//						if(ext.score_ > bestExtScore){
//							bestExtScore = ext.score_;
//						}
//					}
//					for(const auto & lig : possible.second.second){
//						if(lig.score_ > bestLigScore){
//							bestLigScore = lig.score_;
//						}
//					}
//					uint32_t currentScore = bestLigScore + bestExtScore;
//					if (currentScore > bestScore) {
//						bestScore = currentScore;
//						bestMips.clear();
//						bestMips.emplace_back(possible.first);
//					} else if (0 != currentScore && currentScore == bestScore) {
//						bestMips.emplace_back(possible.first);
//					}
//				}
//				if(bestMips.size() == 1){
//					const auto & mip = mipMaster.mips_->mips_.at(bestMips.front());
//
//					//quality control
//					SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, qFilPars);
//
//					std::string failedQaulifierName = MipExtractionStats::getNameForCase(
//							eCase);
//					//log and write read
//					if(!allExtractStats.haveStatFor(mip.name_)){
//						bib::files::makeDirP(sampDirMaster.extractDir_.string(), bib::files::MkdirPar(mip.name_));
//					}
//					mipOuts.add(mip.name_ + failedQaulifierName, seq);
//					//mipOuts.add("unmatchedReads", seq);
//					allExtractStats.increaseCount(mip.name_, eCase);
//
//				}else if(bestMips.size() > 1){
//					//no matches found
//
//					allExtractStats.increaseIndeterminate();
//
//					std::stringstream appName;
//					appName << "[";
//					for(const auto & best : bestMips){
//						const auto & possible = possibleArms.at(best);
//						uint32_t bestExtScore = 0;
//						uint32_t bestLigScore = 0;
//
//						for(const auto & ext : possible.first){
//							if(ext.score_ > bestExtScore){
//								bestExtScore = ext.score_;
//							}
//						}
//						for(const auto & lig : possible.second){
//							if(lig.score_ > bestLigScore){
//								bestLigScore = lig.score_;
//							}
//						}
//						appName << "name=" << best << ";"
//								<< "extScore=" << bestExtScore << ";"
//								<< "ligScore=" << bestLigScore << ";";
//					}
//					appName << "]";
//					seq.seqBase_.name_.append(appName.str());
//					seq.mateSeqBase_.name_.append(appName.str());
//					mipOuts.add("indeterminate", seq);
//					//mipOuts.add("unmatchedReads", seq);
//				}else{
//					std::stringstream ss;
//					ss << __PRETTY_FUNCTION__ << std::endl;
//					ss << "Best Mips vector is empty, which shouldn't be able to happen" << std::endl;
//					ss << "For : " << sampleIOOpt.firstName_ << std::endl;
//					throw std::runtime_error{ss.str()};
//				}
//			}
//		}
//	}
//
//	if (verbose) {
//		std::cout << std::endl;
//	}
//	mipOuts.closeOutAll();
//	VecStr extracOnlyColNames {"sampleName", "mipTarget", "mipFamily", "readNumber",
//			"goodReads", "failedLigationArm", "failedMinLen(<"
//					+ estd::to_string(pars.minLen) + ")" };
//	if (qFilPars.checkingQFrac_) {
//		extracOnlyColNames.emplace_back(
//				"failed_q" + estd::to_string(qFilPars.qualCheck_) + "<"
//						+ estd::to_string(qFilPars.qualCheckCutOff_));
//	} else if (qFilPars.checkingQWindow) {
//		extracOnlyColNames.emplace_back("failed_qW" + qFilPars.qualWindow_);
//	} else {
//		extracOnlyColNames.emplace_back("failed_quality(noneDone)");
//	}
//	extracOnlyColNames.emplace_back("containsNs");
//	table infoTabByTarget(extracOnlyColNames);
//	infoTabByTarget.content_ = allExtractStats.outputContents(*mipMaster.mips_, "\t");
//	infoTabByTarget.outPutContents(
//			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoByTarget.txt", ".txt"), "\t", infoTabByTarget.hasHeader_));
//	std::ofstream sampLog;
//	openTextFile(sampLog, sampDirMaster.extractDir_.string() + "log.txt", ".txt", false, true);
//	sampLog << "Ran on: " << bib::getCurrentDate() << std::endl;
//	sampLog << "Number of Alignments Done: "
//			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
//	sampLog << "Run Time: " << watch.timeLapFormatted(6) << std::endl;
//	alignerObjForFamilyDet.processAlnInfoOutputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose);
//}




int mipsterAnalysisRunner::mipIllumExtractByArmAndFilterPaired(
		const bib::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	mipIllumArmExtractionPars pars;
	setUp.setUpMipIllumArmExtractionPaired(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (!mipMaster.names_->hasSample(pars.sampleName)) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", unrecognized sample name: "
				<< pars.sampleName << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllMinimumExpectedLen(pars.minLen);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	uint64_t maxLen = pars.minLen * 2;
	aligner alignerObjForFamilyDet = aligner(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_);
	//read in mip info and convert to mip class
	MipCollection mips(pars.mipArmsFileName, pars.allowableErrors);
	MipExtractor mExtractor(setUp.pars_.verbose_);
	mExtractor.extractFilterSampleForMipsPaired({setUp.pars_.ioOptions_}, mipMaster,
			alignerObjForFamilyDet, pars);

//	setUp.rLog_ << "Number of Alignments Done: "
//			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	//
	return 0;
}



void extractMultiSamplesPaired(const SetUpMaster & mipMaster,
		const mipIllumArmExtractionParsMultiple& pars,
		const SeqSetUpPars & setUpPars, concurrent::AlignerPool & aligners,
		bib::concurrent::LockableQueue<std::string>& sampsQueue) {
	std::string sampleName = "";
	auto alignerObjForFamilyDet = aligners.popAligner();
	while (sampsQueue.getVal(sampleName)) {
		SeqIOOptions sampOpts = SeqIOOptions::genPairedIn(
				mipMaster.pathSampleRawDataFirstRead(MipFamSamp("", sampleName)),
				mipMaster.pathSampleRawDataSecondRead(MipFamSamp("", sampleName)));
		sampOpts.revComplMate_ = true;
		mipIllumArmExtractionPars samplePars = pars.createForSample(sampleName);
		MipExtractor mExtractor(setUpPars.verbose_);
		mExtractor.extractFilterSampleForMipsPaired({sampOpts}, mipMaster, *alignerObjForFamilyDet,
				samplePars);
	}
}

int mipsterAnalysisRunner::mipIllumExtractByArmAndFilterMultiplePaired(
		const bib::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	mipIllumArmExtractionParsMultiple pars;
	setUp.setUpMipIllumArmExtractionMultiplePaired(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	mipMaster.mips_->setAllMinimumExpectedLen(pars.minLen);
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	std::ofstream logFile;
	openTextFile(logFile,
			bib::files::make_path(mipMaster.directoryMaster_.logsDir_, pars.logFilename),
			".txt", pars.overWriteLog, true);
	logFile << "Ran on: " << getCurrentDate() << std::endl;
	logFile << "Ran from: " << inputCommands.workingDir_ << std::endl;
	logFile << "Command: " << inputCommands.commandLine_ << std::endl;
	bib::concurrent::LockableQueue<std::string> sampsQueue(mipMaster.names_->samples_);
	uint64_t maxLen = pars.minLen * 2;
	concurrent::AlignerPool aligners(maxLen, setUp.pars_.gapInfo_,
			setUp.pars_.scoring_, pars.numThreads);
	aligners.initAligners();
	//read in mip info and convert to mip class
	std::vector<std::thread> threads;
	for (uint32_t threadNum = 0; threadNum < pars.numThreads; ++threadNum) {
		threads.emplace_back(extractMultiSamplesPaired, std::cref(mipMaster), std::cref(pars),
				std::cref(setUp.pars_), std::ref(aligners), std::ref(sampsQueue));
	}
	for (auto & t : threads) {
		t.join();
	}
	logFile << "Total Run Time: " << setUp.timer_.totalTimeFormatted(6) << std::endl;
	return 0;
}

}


