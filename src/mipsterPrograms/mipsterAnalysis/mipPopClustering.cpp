/*
 * mipPopClustering.cpp
 *
 *  Created on: Feb 9, 2016
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




void runPopClusForMip(const MipFamSamp & mipSamp,
		const MipAnalysisDirectoryMaster & directoryMaster,
		const MipsSamplesNames & mipSamps,
		const mipPopulationClusteringPars & pars,
		const SeqSetUpPars & seqPars) {
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	bfs::path mipFamilyDir = njh::files::makeDir(
			directoryMaster.populationClusteringDir_.string() + mipSamp.mipFam_,
			njh::files::MkdirPar("analysis/", pars.overWriteDirs));
	njh::stopWatch watch;
	uint64_t maxSize = 0;
	VecStr missingSamples;
	VecStr foundSamples;
	SeqIOOptions genOpts = SeqIOOptions::genFastqInGz("", true);

	for (const auto& samp : mipSamps.samples_) {
		SampleDirectoryMaster sampDirMaster(directoryMaster, MipFamSamp("", samp));
		SeqIOOptions opts = SeqIOOptions::genFastqInGz(sampDirMaster.getClusteredHapFnp(pars.mipName).string(), true);
		if(bfs::exists(opts.firstName_)){
			SeqInput reader(opts);
			reader.openIn();
			seqInfo seq;
			while(reader.readNextRead(seq)){
				readVec::getMaxLength(seq, maxSize);
			}
			foundSamples.emplace_back(samp);
		}else{
			missingSamples.emplace_back(samp);
		}
	}

	if(foundSamples.empty()){
		if(seqPars.verbose_){
			std::cout << "No reads for mip: " << mipSamp.mipFam_ << std::endl;
		}
		return;
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	// reading expected sequences to compare to
	bool checkingExpected = !pars.refIoOptions.firstName_.empty();
	std::vector<readObject> expectedSeqs;
	if (checkingExpected) {
		expectedSeqs = SeqInput::getReferenceSeq(pars.refIoOptions, maxSize);
	}
	// create aligner class object
	// Hard-coded parameters for alignment parameters follow follow:
	KmerMaps emptyMaps;
	gapScoringParameters gapPars(seqPars.gapInfo_);
	aligner alignerObj = aligner(maxSize, gapPars, seqPars.scoring_, emptyMaps,
			seqPars.qScorePars_,
			seqPars.colOpts_.alignOpts_.countEndGaps_,
			seqPars.colOpts_.iTOpts_.weighHomopolyer_);
	bfs::path alnCacheDir = njh::files::join(
			VecStr { directoryMaster.populationClusteringDir_.string(), pars.mipName,
					"alnCache" });
	alignerObj.processAlnInfoInputNoCheck(alnCacheDir.string(), seqPars.debug_);

	collapser collapserObj = collapser(seqPars.colOpts_);
	collapse::SampleCollapseCollection sampColl(genOpts, pars.masterDir,
			mipFamilyDir,
			PopNamesInfo(pars.mipName, foundSamples, VecStr{}),
			pars.clusteringCutOffs);

	if(!pars.sampleMetaFnp.empty()){
		sampColl.addGroupMetaData(pars.sampleMetaFnp);
	}
	std::vector<sampleCluster> allSamples;
	//
	collapserObj.opts_.kmerOpts_.checkKmers_ = false;
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	for (const auto & samp : foundSamples) {
		if(seqPars.verbose_){
			std::cout << "Starting: " << samp << std::endl;
		}
		SampleDirectoryMaster sampDirMaster(directoryMaster, MipFamSamp("", samp));
		auto inputFiles = std::vector<collapse::SampleCollapseCollection::RepFile>{collapse::SampleCollapseCollection::RepFile{samp,sampDirMaster.getClusteredHapFnp(pars.mipName).string() }};
		for(auto & input : inputFiles){
			input.reNameInput_ = false;
		}
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		sampColl.setUpSample(samp,
				inputFiles,
				alignerObj, collapserObj, seqPars.chiOpts_);
		//sampColl.clusterSample(samp, alignerObj, collapserObj, pars.iteratorMap);
		//exclude clusters that don't have the necessary replicate number
		//defaults to the number of input replicates if none supplied
		if (0 != pars.runsRequired) {
			sampColl.sampleCollapses_[samp]->excludeBySampNum(pars.runsRequired, false);
		} else {
			sampColl.sampleCollapses_[samp]->excludeBySampNum(
					sampColl.sampleCollapses_[samp]->input_.info_.infos_.size(), false);
		}
		if (!expectedSeqs.empty()) {
			bool oldWeighHomopolymers = alignerObj.weighHomopolymers_;
			alignerObj.weighHomopolymers_ = false;
			sampColl.sampleCollapses_[samp]->collapsed_.checkAgainstExpected(
					expectedSeqs, alignerObj, false);
			alignerObj.weighHomopolymers_ = oldWeighHomopolymers;
		}
		for(auto & clus : sampColl.sampleCollapses_[samp]->collapsed_.clusters_){
			updateNameWithBarinfo(clus);
		}

		for(auto & clus : sampColl.sampleCollapses_[samp]->excluded_.clusters_){
			updateNameWithBarinfo(clus);
		}

		if (!pars.keepChimeras) {
			sampColl.sampleCollapses_[samp]->excludeChimeras(false);
		}
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		sampColl.sampleCollapses_[samp]->excludeFraction(pars.fracCutoff, true);

		std::string sortBy = "fraction";
		addOtherVec(allSamples, sampColl.sampleCollapses_[samp]->createOutput(false, sortBy));

		sampColl.dumpSample(samp);
		if(seqPars.verbose_){
			std::cout << "Ending: " << samp << std::endl;
		}
		if (bfs::exists(sampDirMaster.getClusteredHapFnp(pars.mipName).parent_path()) && !pars.keepIntermediateFiles) {
			njh::files::rmDirForce(sampDirMaster.getClusteredHapFnp(pars.mipName).parent_path() );
		}
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(!allSamples.empty()){
		if(seqPars.verbose_){
			std::cout << njh::bashCT::boldGreen("Pop Clustering") << std::endl;
		}
		sampColl.doPopulationClustering(allSamples,
						alignerObj, collapserObj,
						pars.popIteratorMap);
		if (seqPars.verbose_) {
			std::cout << "Ref file for " << pars.mipName << ": "
					<< pars.previousPopFilename << std::endl;
		}
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if (!pars.previousPopFilename.empty()) {
			sampColl.addRefMetaToName(getSeqs<readObject>(pars.previousPopFilename.string()), pars.previousPopErrors);
		}
		//auto sampTab = sampColl.genSampleCollapseInfo(std::set<std::string>{foundSamples.begin(), foundSamples.end()});
		//auto popTab = sampColl.genPopulationCollapseInfo();


		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		sampColl.dumpPopulation(false);

		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(seqPars.verbose_){
			std::cout << njh::bashCT::boldBlack("Printing info...") << std::endl;
		}
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto popTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(mipFamilyDir, "population", "populationCluster.tab.txt"), true);
		auto sampTabOpts = TableIOOpts::genTabFileOut(njh::files::make_path(mipFamilyDir, "selectedClustersInfo.tab.txt.gz"), true);
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		auto tabs = printMipSampleCollapseInfo(sampColl, true, pars.mipName);
		//auto tabs = printMipSampleCollapseInfo(sampColl, !expectedSeqs.empty(), pars.mipName);
		//popTab.outPutContents(popTabOpts);
		//sampTab.outPutContents(sampTabOpts);
		tabs.popTab_.outPutContents(popTabOpts);
		tabs.sampTab_.outPutContents(sampTabOpts);
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		if(!sampColl.popCollapse_){
			sampColl.loadInPreviousPop();
		}
		std::cout << __FILE__ << " " << __LINE__ << std::endl;
		for (const auto & samp : foundSamples) {
			sampColl.setUpSampleFromPrevious(samp);
			/*
			std::cout << sampColl.popCollapse_->collapsed_.subClustersPositions_.size() << std::endl;
			std::cout << samp << std::endl;
			std::cout << mipSamp.mipFam_ << std::endl;*/
			std::cout << __FILE__ << " " << __LINE__ << std::endl;
			auto sampResults = sampColl.sampleCollapses_.at(samp);
			for (auto & clus : sampResults->collapsed_.clusters_) {
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				std::cout << clus.seqBase_.name_ << std::endl;

				clus.processNameForMeta();
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				auto popUID =
						sampColl.popCollapse_->collapsed_.clusters_[sampColl.popCollapse_->collapsed_.subClustersPositions_.at(
								clus.seqBase_.getStubName(true))].seqBase_.getStubName(true);
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				if(MetaDataInName::nameHasMetaData(popUID)){
					MetaDataInName::removeMetaDataInName(popUID);
				}
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
				clus.addMeta("h_popUID", popUID, true);
				//clus.addMeta("region", "", true);
				clus.resetMetaInName();
				std::cout << __FILE__ << " " << __LINE__ << std::endl;
			}
			sampColl.dumpSample(samp);
		}
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	if(pars.cacheAlignments){
		alignerObj.processAlnInfoOutputNoCheck(alnCacheDir.string(), seqPars.debug_);
	}
	if(!pars.keepIntermediateFiles){
		auto clustersDir = njh::files::make_path(mipFamilyDir, "population", "clusters");
		if (bfs::exists(clustersDir) ) {
			njh::files::rmDirForce(clustersDir);
		}
		for (const auto & samp : foundSamples) {
			auto excludedClustersDir = njh::files::make_path(mipFamilyDir, "samplesOutput", samp, "excluded");
			if (bfs::exists(excludedClustersDir) ) {
				njh::files::rmDirForce(excludedClustersDir);
			}
			auto finalInputClustersDir = njh::files::make_path(mipFamilyDir, "samplesOutput", samp, "final", "clusters");
			if (bfs::exists(finalInputClustersDir) ) {
				njh::files::rmDirForce(finalInputClustersDir);
			}
		}
	}
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::ofstream logfile;
	openTextFile(logfile, OutOptions(njh::files::make_path(mipFamilyDir,"log.txt")));
	logfile << "Ran on: " << njh::getCurrentDate() << std::endl;
	logfile << "Number of Alignments Done: "
			<< alignerObj.numberOfAlingmentsDone_ << "\n";
	logfile << "Run Length: " << watch.totalTimeFormatted(6) << std::endl;
	/*
	for (const auto& readsIter : outputReads) {
		if (seqPars.verbose_) {
			std::cout << readsIter.first << " clustersNum: "
					<< readsIter.second.size() << " readsNum: "
					<< readVec::getTotalReadCount(readsIter.second) << std::endl;
		}
		logfile << readsIter.first << " clustersNum: "
				<< readsIter.second.size() << " readsNum: "
				<< readVec::getTotalReadCount(readsIter.second) << std::endl;
	}*/
	std::cout << __FILE__ << " " << __LINE__ << std::endl;
}

int mipsterAnalysisRunner::mipPopulationClustering(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterAnalysisSetUp setUp(inputCommands);
	mipPopulationClusteringPars pars;
	setUp.setUpMipPopulationClustering(pars);
	MipCollection mips(pars.mipArmsFileName, pars.allowableErrors);
	MipAnalysisDirectoryMaster directoryMaster(pars.masterDir);
	MipsSamplesNames mipSamps(pars.mipsSamplesFile);
	if (!mipSamps.hasMip(pars.mipName)) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", unrecognized mip name: "
				<< pars.mipName << std::endl;
		throw std::runtime_error { ss.str() };
	}
	runPopClusForMip(MipFamSamp(pars.mipName, ""), directoryMaster, mipSamps,
			pars, setUp.pars_);
	if (setUp.pars_.verbose_) {
		setUp.logRunTime(std::cout);
	}
	return 0;
}


int mipsterAnalysisRunner::mipPopulationClusteringMultiple(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterAnalysisSetUp setUp(inputCommands);
	mipPopulationClusteringParsMultiple pars;
	setUp.setUpMipPopulationClusteringMultiple(pars);
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
			njh::files::make_path(mipMaster.directoryMaster_.logsDir_ , pars.logFilename), ".json",
			pars.overWriteLog, true);
	logInfo["date"] = getCurrentDate();
	logInfo["workingDir"] = inputCommands.workingDir_;
	logInfo["command"] = inputCommands.commandLine_;
	std::vector<std::thread> threads;
	njh::concurrent::LockableQueue<std::string> mipsQueue(
			mipMaster.names_->mips_);
	std::mutex logMut;
	auto & pairingLog = logInfo["PerMip"];
	auto runPopClusMultiMip = [&logMut,&pairingLog](njh::concurrent::LockableQueue<std::string>& mipsQueue,
			const mipPopulationClusteringParsMultiple& pars,
			const SeqSetUpPars & setUpPars,
			const SetUpMaster & mipMaster){
		std::string mipName;
		while (mipsQueue.getVal(mipName)) {
			njh::stopWatch watch;
			auto currentSampPars = pars.createForMip(mipName);
			runPopClusForMip(MipFamSamp(mipName, ""), mipMaster.directoryMaster_, *mipMaster.names_,
					currentSampPars, setUpPars);
			{
				std::lock_guard<std::mutex> logLock(logMut);
				auto & currentPair = pairingLog[mipName];
				currentPair["mip"] = mipName;
				currentPair["runTime"] = watch.totalTimeFormatted(6);
			}
		}
	};

	for (uint32_t t = 0; t < pars.numThreads; ++t) {
		threads.emplace_back(
				runPopClusMultiMip, std::ref(mipsQueue), std::cref(pars),
						std::cref(setUp.pars_), std::cref(mipMaster));
	}

	for (auto & t : threads) {
		t.join();
	}

	TableIOOpts tabOpts(InOptions(), "\t",
			OutOptions(
					njh::files::make_path(
							mipMaster.directoryMaster_.populationClusteringDir_,
							"allInfo.tab.txt.gz")), "\t", true);
	tabOpts.out_.overWriteFile_ = true;
	std::vector<bfs::path> infofilePaths;
	for(const auto & mip : mipMaster.names_->mips_){
		auto filePath = mipMaster.pathMipPopClusSampInfo(MipFamSamp(mip,""));
		if(bfs::exists(filePath)){
			infofilePaths.emplace_back(filePath);
		}
	}
	MasterTableStaticCache allInfoTab(tabOpts, infofilePaths, true);

	allInfoTab.writeTab();
	logInfo["totalRunTime"] = setUp.timer_.totalTimeFormatted(6);
	logFile << logInfo;
	return 0;
}

}  // namespace njhseq


