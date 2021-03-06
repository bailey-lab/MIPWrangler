/*
 * mipClustering.cpp
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



std::vector<std::shared_ptr<readObject>> getBeginningClusters(
		const cluster & clus, std::vector<identicalCluster> & identicalClusters,
		const std::unordered_map<std::string, uint32_t> & idenToIndex) {
	std::vector<std::shared_ptr<readObject>> ret;
	for (const auto & read : clus.reads_) {
		auto identPosSearch = idenToIndex.find(read->getStubName(true));
		if (idenToIndex.end() == identPosSearch) {
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< ", could not find index for: " << read->seqBase_.name_
					<< " in " << clus.seqBase_.name_ << std::endl;
			throw std::runtime_error{ss.str()};
		} else {
			addOtherVec(ret,identicalClusters[identPosSearch->second].reads_);
		}
	}
	return ret;
}



std::vector<identicalCluster> collapseIdenticalReadsForBarCorrectedReads(
     const std::vector<std::shared_ptr<readObject>> &reads,
		const std::string &repQual){
  std::vector<identicalCluster> ret;
  uint32_t count = 0;
  for (const auto &read : reads) {
  	if(!getSeqBase(read).on_){
  		continue;
  	}
    ++count;
    if (count == 1) {
    	ret.emplace_back(getSeqBase(read));
      continue;
    }
    bool foundMatch = false;
    for (auto &clusterIter : ret) {
      if (getSeqBase(read).seq_ == clusterIter.seqBase_.seq_) {
        clusterIter.addRead(getSeqBase(read));
        foundMatch = true;
        break;
      }
    }
    if (!foundMatch) {
    	ret.emplace_back(getSeqBase(read));
    }
  }


	if ("median" == repQual || "mean" == repQual || "average" == repQual) {
		std::function<uint32_t(std::vector<uint32_t>&)> qualCalcFunc =
				[](std::vector<uint32_t> & quals) {
					return round(vectorMean(quals));
				};
		if ("median" == repQual) {
			qualCalcFunc = [](std::vector<uint32_t> & quals) {
				return round(vectorMedianRef(quals));
			};
		}
		for(auto & seq : ret){
			std::vector<uint32_t> readCountsPerReads;
			uint32_t readTotal = 0;
			for(const auto & subSeq : seq.reads_){
				MetaDataInName subSeqMeta(subSeq->seqBase_.name_);
				readCountsPerReads.emplace_back(subSeqMeta.getMeta<uint32_t>("readCnt"));
				readTotal+= readCountsPerReads.back();
			}
			for(const auto & qualPos : iter::range(seq.seqBase_.qual_.size())){
				std::vector<uint32_t> quals;
				quals.reserve(readTotal);
				for(const auto & subSeqPos : iter::range(seq.reads_.size())){
					addOtherVec(quals, std::vector<uint32_t>(readCountsPerReads[subSeqPos], seq.reads_[subSeqPos]->seqBase_.qual_[qualPos]));
				}
				seq.seqBase_.qual_[qualPos] = qualCalcFunc(quals);
			}
		}
	} else if (repQual == "bestSeq") {
		for(auto & seq : ret){
			std::vector<uint32_t> readCountsPerReads;
			uint32_t readTotal = 0;
			std::vector<uint32_t> bestSeqPositions;
			uint32_t maxReadCount = 0;
			for(const auto & subSeqPos : iter::range(seq.reads_.size())){
				const auto & subSeq = seq.reads_[subSeqPos];
				MetaDataInName subSeqMeta(subSeq->seqBase_.name_);
				readCountsPerReads.emplace_back(subSeqMeta.getMeta<uint32_t>("readCnt"));
				readTotal+= readCountsPerReads.back();
				if(maxReadCount > readCountsPerReads.back()){
					maxReadCount = readCountsPerReads.back();
					bestSeqPositions.clear();
					bestSeqPositions.emplace_back(subSeqPos);
				}else if(maxReadCount == readCountsPerReads.back()){
					bestSeqPositions.emplace_back(subSeqPos);
				}
			}
			if(bestSeqPositions.size() == 1){
				seq.seqBase_.qual_ = seq.reads_[bestSeqPositions.front()]->seqBase_.qual_;
			}else{
				std::vector<std::shared_ptr<readObject>> bestSeqs;
				for(const auto & best : bestSeqPositions){
					bestSeqs.push_back(seq.reads_[best]);
				}
				readVec::allSetQualCheck(bestSeqs, 30);
				readVecSorter::sortByQualCheck(bestSeqs, true);
				seq.seqBase_.qual_ = bestSeqs.front()->seqBase_.qual_;
			}
		}
	} else if (repQual == "bestQual") {
		for(auto & seq : ret){
			std::vector<uint32_t> readCountsPerReads;
			uint32_t readTotal = 0;
			for(const auto & subSeq : seq.reads_){
				MetaDataInName subSeqMeta(subSeq->seqBase_.name_);
				readCountsPerReads.emplace_back(subSeqMeta.getMeta<uint32_t>("readCnt"));
				readTotal+= readCountsPerReads.back();
			}
			for(const auto & qualPos : iter::range(seq.seqBase_.qual_.size())){
				std::vector<uint32_t> quals;
				quals.reserve(readTotal);
				for(const auto & subSeqPos : iter::range(seq.reads_.size())){
					quals.emplace_back(seq.reads_[subSeqPos]->seqBase_.qual_[qualPos]);
				}
				seq.seqBase_.qual_[qualPos] = vectorMaximum(quals);
			}
		}
	} else {
		std::stringstream ss;
		ss << "Unrecognized qualRep: " << repQual << std::endl;
		ss << "Needs to be median, average, bestSeq, bestQual, or worst"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
  //identicalCluster::setIdneticalClusterQual(ret, repQual);
  readVec::allSetLetterCount(ret);
  readVec::allUpdateName(ret);
  return ret;
}



void runClusteringForMipFamForSamp(const MipFamSamp &mipSampName,
		const SampleDirectoryMaster & sampDirMaster,
		const mipClusteringPars & pars,
		const SeqSetUpPars & seqPars,
		aligner & alignerObj,
		collapser collapserObj){
	SeqIOOptions options = SeqIOOptions::genFastqInGz(
			njh::files::join(
					VecStr { sampDirMaster.barCorDir_.string(), mipSampName.mipFam_,
							mipSampName.mipFam_ + "_all.fastq.gz" }));
	if (options.inExists()) {
		njh::stopWatch watch;
		SeqInput reader(options);
		reader.openIn();
		std::vector<std::shared_ptr<readObject>> allReads = reader.readAllReadsPtrs<
				readObject>();
		//process reads and set on_ to false if they have remove set to true
		for(const auto & read : allReads){
			read->processNameForMeta();
			if(read->containsMeta("remove") && "true" == read->getMeta("remove")){
				read->seqBase_.on_ = false;
			}
		}
		//now to cluster these all the barcode corrected reads;
		auto identicalClusters = collapseIdenticalReadsForBarCorrectedReads(allReads, pars.qualRep);

		//now check to see if there are any reads to cluster
		//if all of them were contamination this should be zero
		if(!identicalClusters.empty()){
			bfs::path mipFamilyDir = njh::files::makeDir(
					sampDirMaster.clusDir_.string(), njh::files::MkdirPar(mipSampName.mipFam_,
					pars.overWriteDirs));
			alignerObj.resetAlnCache();
			if(pars.cacheAlignments){
				alignerObj.processAlnInfoInputNoCheck(
						njh::files::join(sampDirMaster.clusAlnCacheDir_.string(), mipSampName.mipFam_).string(),seqPars.debug_);
			}
			alignerObj.numberOfAlingmentsDone_ = 0;
			std::unordered_map<std::string, uint32_t> idenToIndex;
			for (const auto & iPos : iter::range(identicalClusters.size())) {
				idenToIndex[identicalClusters[iPos].getStubName(true)] = iPos;
			}
			std::vector<cluster> clusters = baseCluster::convertVectorToClusterVector<
					cluster>(identicalClusters);
			if (clusters.size() > 1) {
				if (clusters.size() > 10) {
					collapserObj.opts_.verboseOpts_.verbose_ = seqPars.verbose_;
				} else {
					collapserObj.opts_.verboseOpts_.verbose_ = false;
				}
				readVecSorter::sort(clusters);
				auto cutOff = processRunCutoff(collapserObj.opts_.kmerOpts_.runCutOffString_,
						readVec::getTotalReadCount(clusters));
				auto currentKMaps = indexKmers(clusters, collapserObj.opts_.kmerOpts_.kLength_, cutOff,
						collapserObj.opts_.kmerOpts_.kmersByPosition_, seqPars.expandKmerPos_,
						seqPars.expandKmerSize_);
				alignerObj.kMaps_ = currentKMaps;
				std::pair<std::vector<cluster>, std::vector<cluster>> clusSplit =
						readVecSplitter::splitVectorOnReadCount(clusters, 1);
				clusters = clusSplit.first;
				/**@todo add functionality offered by runFullClustering of bin and such */
				collapserObj.runFullClustering(clusters, pars.iterMap,
						pars.iterMap, alignerObj, "", seqPars.ioOptions_,
						seqPars.refIoOptions_, SnapShotsOpts(false, ""));
				{
					addOtherVec(clusters, clusSplit.second);
					collapserObj.runFullClustering(clusters, pars.iterMap,
											pars.iterMap, alignerObj, "", seqPars.ioOptions_,
											seqPars.refIoOptions_, SnapShotsOpts(false, ""));
				}
			}
			if(pars.cacheAlignments){
				alignerObj.processAlnInfoOutputNoCheck(
						njh::files::join(sampDirMaster.clusAlnCacheDir_.string(),
								mipSampName.mipFam_).string(), seqPars.debug_);
			}
			SeqIOOptions opts = SeqIOOptions::genFastqOutGz(
					njh::files::make_path(mipFamilyDir,  mipSampName.mipFam_).string() + "_clustered");

			renameReadNames(clusters, "[samp=" + mipSampName.samp_ + ";" + "mipFam=" + mipSampName.mipFam_ +"]",
					true, true, true);
			bfs::path mipFamilyAllClustersDir = "";
			if(pars.writeOutClusters){
				mipFamilyAllClustersDir = njh::files::makeDir(mipFamilyDir, njh::files::MkdirPar("allInputClusters"));
			}

			OutputStream outInfoFile(OutOptions(njh::files::make_path(mipFamilyDir,"info.tab.txt")));
			outInfoFile
					<< "clusterName\tbarcodes\tbarcodeFraction\treads\treadsFraction"
					<< std::endl;
			std::unordered_map<std::string, uint32_t> readAmounts;
			uint32_t totalReads = 0;
			uint32_t barcodeTotal = 0;
			for (auto & clus : clusters) {
				auto begReads = getBeginningClusters(clus, identicalClusters, idenToIndex);
				uint32_t readAmount = 0;
				for(auto & read : begReads){
					read->processNameForMeta();
					readAmount += read->getMeta<uint32_t>("readCnt");
				}
				barcodeTotal+= clus.seqBase_.cnt_;
				totalReads += readAmount;
				clus.appendName("_R" + estd::to_string(readAmount));
				MetaDataInName meta(clus.seqBase_.name_);
				meta.addMeta("readCnt", readAmount, true);
				meta.addMeta("barcodeCnt", clus.seqBase_.cnt_, true);
				meta.resetMetaInName(clus.seqBase_.name_);
				readAmounts[clus.seqBase_.name_] = readAmount;
				if(pars.writeOutClusters){
					SeqIOOptions outOpts = opts;
					outOpts.out_.outFilename_ = njh::files::make_path(mipFamilyAllClustersDir, clus.seqBase_.name_);
					SeqOutput clusWriter(outOpts);
					clusWriter.openWrite(begReads);
				}
			}
			for (auto & clus : clusters) {
				outInfoFile << clus.seqBase_.name_
						<< "\t" << clus.seqBase_.cnt_
						<< "\t" << clus.seqBase_.cnt_ / static_cast<double>(barcodeTotal)
						<< "\t" << readAmounts[clus.seqBase_.name_]
						<< "\t" << readAmounts[clus.seqBase_.name_] / static_cast<double>(totalReads)
						<< std::endl;
			}
			SeqOutput::write(clusters, opts);
			if(pars.writeOutClusters){
				auto mipFamilyClustersDir = njh::files::makeDir(mipFamilyDir, njh::files::MkdirPar("clusters"));
				clusterVec::allWriteClustersInDir(clusters, mipFamilyClustersDir.string(), opts);
			}
			OutputStream logfile(OutOptions(njh::files::make_path(mipFamilyDir, "log.txt")));
			logfile << "Ran on: " << njh::getCurrentDate() << std::endl;
			logfile << "Number of Alignments Done: " << alignerObj.numberOfAlingmentsDone_ << "\n";
			logfile << "Run Length: " << watch.totalTimeFormatted(6) << std::endl;
		}
	}
	if (bfs::exists(options.firstName_.parent_path()) && !pars.keepIntermediateFiles) {
		njh::files::rmDirForce(options.firstName_.parent_path());
	}
}


int mipsterAnalysisRunner::mipClustering(
		const njh::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	mipClusteringPars pars;
	setUp.setUpMipClustering(pars);
	//read in mip files
	MipCollection mips(pars.mipArmsFileName, pars.allowableErrors);
	MipsSamplesNames mipSamps(pars.mipsSamplesFile);
	if (!mipSamps.hasSample(pars.sampleName)) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", unrecognized sample name: "
				<< pars.sampleName << std::endl;
		throw std::runtime_error { ss.str() };
	}
	MipAnalysisDirectoryMaster directoryMaster(pars.masterDir);
	SampleDirectoryMaster sampDirMaster(directoryMaster,
			MipFamSamp("", pars.sampleName));
	//check for folder created by previous steps to see if we can run
	sampDirMaster.checkForExtractDirectoryThrow();
	sampDirMaster.checkForBarCorDirectoryThrow();
	//create folders neccessary for this sample
	sampDirMaster.ensureClusDirectoryExist();
	//set up collapser
	collapser collapserObj = collapser(setUp.pars_.colOpts_);
	aligner alignerObj = aligner(800, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_), setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_, setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	for (const auto & mipFam : mips.mipFamilies_) {
		runClusteringForMipFamForSamp(MipFamSamp(mipFam, pars.sampleName),
				sampDirMaster, pars, setUp.pars_, alignerObj, collapserObj);
	}
	return 0;
}


int mipsterAnalysisRunner::mipClusteringMultiple(const njh::progutils::CmdArgs & inputCommands){
	mipsterAnalysisSetUp setUp(inputCommands);
	mipClusteringParsMultiple pars;
	setUp.setUpMipClusteringMultiple(pars);

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

	Json::Value logInfo;
	std::ofstream logFile;
	openTextFile(logFile,
			njh::files::make_path(mipMaster.directoryMaster_.logsDir_, pars.logFilename), ".json",
			pars.overWriteLog, true);
	logInfo["date"] = getCurrentDate();
	logInfo["workingDir"] = inputCommands.workingDir_;
	logInfo["command"] = inputCommands.commandLine_;

	mipMaster.makeClusteringDirs();


	//set up collapser
	collapser collapserObj = collapser(setUp.pars_.colOpts_);
	aligner alignerObj = aligner(800, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			KmerMaps(setUp.pars_.colOpts_.kmerOpts_.kLength_), setUp.pars_.qScorePars_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_,
			setUp.pars_.colOpts_.iTOpts_.weighHomopolyer_);
	concurrent::AlignerPool alnPool(alignerObj, pars.numThreads);
	alnPool.initAligners();

	std::vector<MipFamSamp> allPairings = mipMaster.getPairsWithBarCor(pars.numThreads);
	njh::concurrent::LockableQueue<MipFamSamp> pairQueue(allPairings);

	auto & pairingLog = logInfo["PerMipSample"];
	std::mutex logMut;
	auto runMipSamp = [&alnPool,&pairQueue,&pars,&setUp,&logMut,&pairingLog,&collapserObj](const SetUpMaster & mipMaster){
		MipFamSamp mipSamp("", "");
		//std::cout << std::this_thread::get_id() << " " << __LINE__ << std::endl;

		auto curAlignerObj = alnPool.popAligner();
		Json::Value currrentPairingLog;
//		std::cout << std::this_thread::get_id() << " " << __LINE__ << std::endl;

		while(pairQueue.getVal(mipSamp)){
//			std::cout << std::this_thread::get_id() << " samp: " << mipSamp.samp_ << " mipFam: " << mipSamp.mipFam_  << __LINE__ << std::endl;
			njh::stopWatch watch;
			auto sampPars = pars.createForSample(mipSamp.samp_);
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, mipSamp);
			sampDirMaster.checkForClusDirectoryThrow();
			runClusteringForMipFamForSamp(mipSamp, sampDirMaster, pars, setUp.pars_, *curAlignerObj, collapserObj);

			auto & currentPair = currrentPairingLog[mipSamp.samp_ + "_" + mipSamp.mipFam_];
			currentPair["sample"] = mipSamp.samp_;
			currentPair["mip"] = mipSamp.mipFam_;
			currentPair["runTime"] = watch.totalTimeFormatted(6);
//			std::cout << std::this_thread::get_id() << " samp: " << mipSamp.samp_ << " mipFam: " << mipSamp.mipFam_  << __LINE__ << std::endl;
		}
//		std::cout << std::this_thread::get_id() << " " << __LINE__ << std::endl;
		{
			std::lock_guard<std::mutex> logLock(logMut);
			pairingLog.append(currrentPairingLog);
		}
//		std::cout << std::this_thread::get_id() << " " << __LINE__ << std::endl;

	};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars.numThreads; ++t){
		threads.emplace_back(std::thread(runMipSamp, std::cref(mipMaster)));
	}
	for(auto & t : threads){
		t.join();
	}

	logInfo["totalRunTime"] = setUp.timer_.totalTimeFormatted(6);
	logFile << logInfo;
	return 0;
}


}
