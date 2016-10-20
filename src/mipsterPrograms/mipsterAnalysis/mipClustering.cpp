/*
 * mipClustering.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: nick
 */



#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"


namespace bibseq {



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



void runClusteringForMipFamForSamp(const MipFamSamp &mipSampName,
		const SampleDirectoryMaster & sampDirMaster,
		const mipClusteringPars & pars,
		const SeqSetUpPars & seqPars,
		aligner & alignerObj,
		collapser collapserObj){
	SeqIOOptions options = SeqIOOptions::genFastqIn(
			bib::files::join(
					VecStr { sampDirMaster.barCorDir_.string(), mipSampName.mipFam_,
							mipSampName.mipFam_ + "_all.fastq" }));
	if (options.inExists()) {
		bib::stopWatch watch;
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
		auto identicalClusters = clusterCollapser::collapseIdenticalReads(allReads,
				"median");
		//now check to see if there are any reads to cluster
		//if all of them were contamination this should be zero
		if(!identicalClusters.empty()){
			std::string mipFamilyDir = bib::files::makeDir(
					sampDirMaster.clusDir_.string(), bib::files::MkdirPar(mipSampName.mipFam_,
					pars.overWriteDirs));
			alignerObj.resetAlnCache();
			alignerObj.processAlnInfoInputNoCheck(
					bib::files::join(sampDirMaster.clusAlnCacheDir_.string(),
							mipSampName.mipFam_),seqPars.debug_);
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
			alignerObj.processAlnInfoOutputNoCheck(
					bib::files::join(sampDirMaster.clusAlnCacheDir_.string(),
							mipSampName.mipFam_), seqPars.debug_);
			SeqIOOptions opts = SeqIOOptions::genFastqOut(
					mipFamilyDir + mipSampName.mipFam_ + "_clustered");

			renameReadNames(clusters, "[samp=" + mipSampName.samp_ + ";" + "mipFam=" + mipSampName.mipFam_ +"]",
					true, true, true);
			auto mipFamilyAllClustersDir = bib::files::makeDir(mipFamilyDir,
					bib::files::MkdirPar("allInputClusters"));
			std::ofstream outInfoFile;
			openTextFile(outInfoFile, OutOptions(mipFamilyDir + "info.tab.txt"));
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
				readAmounts[clus.seqBase_.name_] = readAmount;
				SeqIOOptions outOpts = opts;
				outOpts.out_.outFilename_ = mipFamilyAllClustersDir + clus.seqBase_.name_;
				SeqOutput clusWriter(outOpts);
				clusWriter.openWrite(begReads);
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
			auto mipFamilyClustersDir = bib::files::makeDir(mipFamilyDir, bib::files::MkdirPar("clusters"));
			clusterVec::allWriteClustersInDir(clusters, mipFamilyClustersDir, opts);
			std::ofstream logfile;
			openTextFile(logfile, OutOptions(mipFamilyDir + "log.txt"));
			logfile << "Ran on: " << bib::getCurrentDate() << std::endl;
			logfile << "Number of Alignments Done: "
					<< alignerObj.numberOfAlingmentsDone_ << "\n";
			logfile << "Run Length: " << watch.totalTimeFormatted(6) << std::endl;
		}
	}
}


int mipsterAnalysisRunner::mipClustering(
		const bib::progutils::CmdArgs & inputCommands) {
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


int mipsterAnalysisRunner::mipClusteringMultiple(const bib::progutils::CmdArgs & inputCommands){
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
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	Json::Value logInfo;
	std::ofstream logFile;
	openTextFile(logFile,
			mipMaster.directoryMaster_.logsDir_.string() + pars.logFilename, ".json",
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
	bib::concurrent::LockableQueue<MipFamSamp> pairQueue(allPairings);

	auto & pairingLog = logInfo["PerMipSample"];
	std::mutex logMut;
	auto runMipSamp = [&alnPool,&pairQueue,&pars,&setUp,&logMut,&pairingLog,&collapserObj](const SetUpMaster & mipMaster){
		MipFamSamp mipSamp("", "");
		auto curAlignerObj = alnPool.popAligner();
		Json::Value currrentPairingLog;
		while(pairQueue.getVal(mipSamp)){
			bib::stopWatch watch;
			auto sampPars = pars.createForSample(mipSamp.samp_);
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, mipSamp);
			sampDirMaster.checkForClusDirectoryThrow();
			runClusteringForMipFamForSamp(mipSamp, sampDirMaster, pars, setUp.pars_, *curAlignerObj, collapserObj);
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
