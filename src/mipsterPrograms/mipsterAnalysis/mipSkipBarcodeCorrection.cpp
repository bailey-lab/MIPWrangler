/*
 * mipSkipBarcodeCorrection.cpp
 *
 *  Created on: May 8, 2017
 *      Author: nick
 */

/*
 * mipBarcodeCorrect.cpp
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







void runBarCorSkipForMipFamForSamp(const MipFamSamp &mipSampName,
		const SetUpMaster & mipMaster, aligner & alignerObj,
		const SeqSetUpPars & setUpPars,
		const mipBarcodeCorrectionPars & pars,
		const SampleDirectoryMaster & sampDirMaster) {
	if(setUpPars.debug_){
		std::cout << "Starting analysis for mip: " << mipSampName.mipFam_
				<< " in samp: " << mipSampName.samp_ << std::endl;
	}
	if(setUpPars.debug_){
		std::cout << "Found reads for mip: " << mipSampName.mipFam_
				<< " in samp: " << mipSampName.samp_ << std::endl;
	}
	auto mipsForFam = mipMaster.mips_->getMipsForFamily(mipSampName.mipFam_);
	BarcodeFilterStats filterStats(mipSampName.samp_);
	njh::stopWatch watch;
	if(setUpPars.debug_){
		std::cout << "Making directory for mip: " << mipSampName.mipFam_
				<< " in samp: " << mipSampName.samp_ << std::endl;
	}
	bfs::path mipFamilyDir = njh::files::makeDir(
			sampDirMaster.barCorDir_.string(),njh::files::MkdirPar( mipSampName.mipFam_,
			pars.overWriteDirs));
	if (setUpPars.debug_) {
		std::cout << "Done making directory for mip: " << mipSampName.mipFam_
				<< " in samp: " << mipSampName.samp_ << std::endl;
	}
	uint32_t mipCorrectedNumber = 0;
	SeqIOOptions outOpts = SeqIOOptions::genFastqOut( mipFamilyDir.string() + mipSampName.mipFam_
					+ "_all.fastq");
	SeqOutput writer(outOpts);
	writer.openOut();
	bool foundNone = true;
	alignerObj.resetAlnCache();
	alignerObj.processAlnInfoInputNoCheck(
			njh::files::make_path(sampDirMaster.barCorAlnCacheDir_, mipSampName.mipFam_).string(),
			setUpPars.debug_);
	for (const auto & mipName : mipsForFam) {
		if(setUpPars.debug_){
			std::cout << "Creating seqIoOptions for mip: " << mipSampName.mipFam_
				<< " in samp: " << mipSampName.samp_ << std::endl;
		}
		SeqIOOptions options = SeqIOOptions::genFastqIn(njh::files::join(
				VecStr { sampDirMaster.extractDir_.string(), mipName, mipName + ".fastq" }));
		if(setUpPars.debug_){
			std::cout << options.firstName_ << std::endl;
		}
		if (options.inExists()) {
			foundNone = false;
			charCounter ligBarCounter(std::vector<char> { 'A', 'C', 'G', 'T' });
			charCounter extBarCounter(std::vector<char> { 'A', 'C', 'G', 'T' });
			std::ofstream barcodeFile;
			openTextFile(barcodeFile, njh::files::make_path(mipFamilyDir, mipName ).string() + "_barcodes",
					".tab.txt", false, true);
			//barcodeFile << "Barcode\treadCnt\tfilePos" << std::endl;
			barcodeFile << "Barcode\treadCnt" << "\n";
			BarcodeFilterStats::BarcodeFilterStat tarStat(mipName, mipSampName.mipFam_);
			SeqInput reader(options);
			auto reads = reader.readAllReadsPtrs<MippedRead>();
			tarStat.initial_ = reads.size();
			//key1 ext bar, key2 lig bar
			std::unordered_map<std::string,
					std::unordered_map<std::string,
							std::vector<std::shared_ptr<MippedRead>>> >sameBarcodes;
			uint32_t readCount = 0;
			//determine barcodes
			for (const auto & read : reads) {
				++readCount;
				if (setUpPars.debug_) {
					std::cout << "\rCurrently on " << readCount << " of " << reads.size()
							<< " for " << mipName;
					std::cout.flush();
				}
				//this will both determine the barcodes and trim off the barcodes and the arms
				read->barInfo_ = std::make_shared<BarcodeInfo>(
						mipMaster.mips_->mips_.at(mipName).determineBarcodesTrim(read->seqBase_));
				sameBarcodes[read->barInfo_->extBar_][read->barInfo_->ligBar_].push_back(
						read);
			}
			if (setUpPars.debug_) {
				std::cout << std::endl;
			}
			std::vector<std::vector<std::shared_ptr<MippedRead>>>readsPerBarcode;
			for (const auto & extBar : sameBarcodes) {
				if (setUpPars.debug_) {
					std::cout << "Currently on extension barcode: " << extBar.first
							<< " which has " << extBar.second.size() << " ligation barcode(s)"
							<< std::endl;
				}
				// if all second barcodes are the same just emplace the reads
//				filterOnMultipleLigBar(extBar.second, readsPerBarcode, tarStat);
				//skipping here
				for(const auto & ext : extBar.second){
					readsPerBarcode.push_back(ext.second);
				}
			}
			for (const auto & barReads : readsPerBarcode) {
				for(const auto & barRead : barReads){
					std::shared_ptr<MippedRead> correctedRead = barRead;
					BarcodeInfo currentBarInfo = *(correctedRead->barInfo_);
					//skipping here
//					auto correctedRead = filterWithBarcodeCoverage(barReads, alignerObj,
//							setUpPars, pars, tarStat);
					tarStat.final_ += correctedRead->seqBase_.cnt_;
					tarStat.barCoverage_.push_back(std::round(correctedRead->seqBase_.cnt_));
					std::shared_ptr<MippedRead> topRead = std::make_shared<MippedRead>(
							seqInfo(correctedRead->seqBase_.name_, correctedRead->seqBase_.seq_,
									correctedRead->seqBase_.qual_));
					topRead->barInfo_ = std::make_shared<BarcodeInfo>(
							*(correctedRead->barInfo_));
					topRead->seqBase_.name_ = njh::pasteAsStr(mipCorrectedNumber,
							"[samp=",pars.sampleName, ";",
							"mipTar=", mipName, ";",
							"mipFam=", mipSampName.mipFam_,";",
							"bar=", topRead->barInfo_->fullBar_, ";",
							"readCnt=",correctedRead->seqBase_.cnt_, "]");
					topRead->updateName();
					topRead->appendName(
							"_R" + estd::to_string(correctedRead->seqBase_.cnt_));
					++mipCorrectedNumber;
					if (!topRead->barInfo_->ligBar_.empty()) {
						ligBarCounter.increaseCountByString(topRead->barInfo_->ligBar_,
								correctedRead->seqBase_.cnt_);
					}
					extBarCounter.increaseCountByString(topRead->barInfo_->extBar_,
							correctedRead->seqBase_.cnt_);
					barcodeFile << topRead->barInfo_->fullBar_ << "\t"
							<< correctedRead->seqBase_.cnt_ << "\n";
					/*barcodeFile << topRead->barInfo_->fullBar_ << "\t"
							<< correctedRead->seqBase_.cnt_ << "\t" << writer.tellpPri()
							<< std::endl;*/
					writer.writeNoCheck(topRead);
					if (setUpPars.debug_) {
						std::cout << "Done Adding barcode: "
								<< correctedRead->barInfo_->fullBar_ << std::endl;
					}
				}
			}
			if (mipMaster.mips_->mips_.at(mipName).ligBarcodeLen_ > 0) {
				std::ofstream ligBarCompOutFile;
				openTextFile(ligBarCompOutFile,
						njh::files::make_path(mipFamilyDir,  mipName + "_ligBarNucComp").string(), ".tab.txt", false, true);
				ligBarCounter.resetAlphabet(true);
				ligBarCounter.setFractions();
				ligBarCounter.outPutInfo(ligBarCompOutFile, false);
			}
			std::ofstream extBarCompOutFile;
			openTextFile(extBarCompOutFile, njh::files::make_path(mipFamilyDir, mipName).string() + "_extBarNucComp",
					".tab.txt", false, true);
			extBarCounter.resetAlphabet(true);
			extBarCounter.setFractions();
			extBarCounter.outPutInfo(extBarCompOutFile, false);
			filterStats.addFilterStat(tarStat);
		}
	}
	if (foundNone) {
		//didn't find any extracted reads remove the directory created so it doesn't look like we did
		writer.closeOut();
		njh::files::rmDirForce(mipFamilyDir);
	} else {
		std::ofstream filterInfoFile;
		openTextFile(filterInfoFile,
				OutOptions(njh::files::make_path(mipFamilyDir , "barcodeFilterStats.tab.txt")));
		filterStats.printInfo(filterInfoFile, "\t");

		alignerObj.processAlnInfoOutputNoCheck(
				njh::files::join(sampDirMaster.barCorAlnCacheDir_.string(),
						mipSampName.mipFam_).string(), setUpPars.debug_);
		std::ofstream logfile;
		openTextFile(logfile, OutOptions(njh::files::make_path(mipFamilyDir, "log.txt")));
		logfile << "Ran on: " << njh::getCurrentDate() << std::endl;
		logfile << "Number of Alignments Done: "
				<< alignerObj.numberOfAlingmentsDone_ << "\n";
		logfile << "Run Length: " << watch.totalTimeFormatted(6) << std::endl;
	}
}


int mipsterAnalysisRunner::mipSkipBarcodeCorrection(const njh::progutils::CmdArgs & inputCommands){
	mipsterAnalysisSetUp setUp(inputCommands);
	mipBarcodeCorrectionPars pars;
	setUp.setUpMipBarcodeCorrection(pars);
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
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	if (!mipMaster.names_->hasSample(pars.sampleName)) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", unrecognized sample name: "
				<< pars.sampleName << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	//check for extraction directory and the rawDir and the extractInfo file
	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
	bfs::path extractInfoFilename = njh::files::make_path(sampDirMaster.extractDir_, "extractInfoByTarget.txt");

	sampDirMaster.checkForExtractDirectoryThrow();
	checkExistenceThrow(extractInfoFilename, __PRETTY_FUNCTION__);

	sampDirMaster.ensureBarCorDirectoryExist(pars.cacheAlignments);
	aligner alignerObj = aligner(800, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_),
			setUp.pars_.qScorePars_, setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);

	for (const auto & mipFam : mipMaster.mips_->mipFamilies_) {
		runBarCorSkipForMipFamForSamp(MipFamSamp(mipFam, pars.sampleName), mipMaster,
				alignerObj, setUp.pars_, pars, sampDirMaster);
	}

	TableIOOpts barFilOpts = TableIOOpts(
			OutOptions(sampDirMaster.barCorDir_.string() + "barcodeFilterStats.tab.txt", ".tab.txt", "tab"
			,false,true, false), "\t", true);
	barFilOpts.hasHeader_ = true;
	MasterTableCache masterBarFilterTab(barFilOpts, sampDirMaster.barCorDir_.string(), true,
			2, std::regex(".*barcodeFilterStats.tab.txt"));
	masterBarFilterTab.writeTab();
	return 0;
}


int mipsterAnalysisRunner::mipSkipBarcodeCorrectionMultiple(const njh::progutils::CmdArgs & inputCommands){
	mipsterAnalysisSetUp setUp(inputCommands);
	mipBarcodeCorrectionParsMultiple pars;
	setUp.setUpMipBarcodeCorrectionMultiple(pars);

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
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	Json::Value logInfo;
	std::ofstream logFile;
	openTextFile(logFile,njh::files::make_path(
			mipMaster.directoryMaster_.logsDir_, pars.logFilename), ".json",
			pars.overWriteLog, true);
	logInfo["date"] = getCurrentDate();
	logInfo["workingDir"] = inputCommands.workingDir_;
	logInfo["command"] = inputCommands.commandLine_;
	if(setUp.pars_.debug_){
		std::cout << "Making Barcode correction directories" << std::endl;
	}
	mipMaster.makeBarcodeCorDirs(pars.cacheAlignments);
	if(setUp.pars_.debug_){
		std::cout << "Done making Barcode correction directories" << std::endl;
		std::cout << "Creating mip and sample parings" << std::endl;
	}
	std::vector<MipFamSamp> allPairings = mipMaster.getPairsWithExtracted(pars.numThreads);
	njh::concurrent::LockableQueue<MipFamSamp> pairQueue(allPairings);
	if(setUp.pars_.debug_){
		std::cout << "Creating aligner pool" << std::endl;
	}
	aligner alignerObj = aligner(800, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_), setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	concurrent::AlignerPool alnPool(alignerObj, pars.numThreads);
	alnPool.initAligners();
	auto & pairingLog = logInfo["PerMipSample"];
	std::mutex logMut;
	auto runMipSamp = [&alnPool,&pairQueue,&pars,&setUp,&logMut,&pairingLog](const SetUpMaster & mipMaster){
		MipFamSamp mipSamp("", "");
		auto curAlignerObj = alnPool.popAligner();
		Json::Value currrentPairingLog;
		while(pairQueue.getVal(mipSamp)){
			njh::stopWatch watch;
			auto sampPars = pars.createForSample(mipSamp.samp_);
			if(setUp.pars_.debug_){
				std::cout << "Checking for extraction directory and if mip: " << mipSamp.mipFam_
						<< " has reads in samp: " << mipSamp.samp_  << std::endl;
			}
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, mipSamp);
			if(setUp.pars_.debug_){
				std::cout << "Mip " <<  mipSamp.mipFam_ << " does have reads, correcting with barcodes" << std::endl;
			}
			runBarCorSkipForMipFamForSamp(mipSamp, mipMaster, *curAlignerObj, setUp.pars_, sampPars,sampDirMaster);
			if(setUp.pars_.debug_){
				std::cout << "Done correcting with barcodes for mip: " << mipSamp.mipFam_
						<< " in samp: " << mipSamp.samp_ <<  std::endl;
			}
			auto & currentPair = currrentPairingLog[mipSamp.samp_ + "_" + mipSamp.mipFam_];
			currentPair["sample"] = mipSamp.samp_;
			currentPair["mip"] = mipSamp.mipFam_;
			currentPair["runTime"] = watch.totalTimeFormatted(6);
		}
		{
			std::lock_guard<std::mutex> logLock(logMut);
			pairingLog.append(currrentPairingLog);
		}
	};
	std::vector<std::thread> threads;
	if(setUp.pars_.debug_){
		std::cout << "Starting Analysis" << std::endl;
	}
	for(uint32_t t = 0; t < pars.numThreads; ++t){
		threads.emplace_back(std::thread(runMipSamp, std::cref(mipMaster)));
	}
	for(auto & t : threads){
		t.join();
	}

	//std::unordered_map<std::string, MasterTableCache> sampleExtractTables;
	for(const auto & samp : mipMaster.names_->samples_){
		bfs::path sampBarCorDir = njh::files::join(VecStr{samp, samp + "_mipBarcodeCorrection/"});
		TableIOOpts barFilOpts = TableIOOpts(OutOptions(njh::files::make_path(
				sampBarCorDir, "barcodeFilterStats.tab.txt"), ".tab.txt", "tab",
				false,true, false), "\t", true);
		barFilOpts.hasHeader_ = true;
		MasterTableCache masterBarFilterTab(barFilOpts, sampBarCorDir.string(), true,
				2, std::regex(".*barcodeFilterStats.tab.txt"));
		masterBarFilterTab.writeTab();
	}
	logInfo["totalRunTime"] = setUp.timer_.totalTimeFormatted(6);
	logFile << logInfo;
	return 0;
}

}




