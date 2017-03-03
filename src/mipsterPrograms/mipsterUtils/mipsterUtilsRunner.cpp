
//  mipsterUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterUtilsRunner.hpp"
    
namespace bibseq {

mipsterUtilsRunner::mipsterUtilsRunner()
    : bib::progutils::programRunner({addFunc("alignTarget", alignTargets, false),
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 addFunc("processProcessedMips", processProcessedMips, false),
																		 addFunc("scanForContam", scanForContam, false),
																		 addFunc("processMipOverlapGraph", processMipOverlapGraph, false),
																		 addFunc("processMipOverlapGraphSingle", processMipOverlapGraphSingle, false),
																		 addFunc("rearmTargetsAndCombine", rearmTargetsAndCombine, false)
},//
                    "mipsterUtils") {}


int mipsterUtilsRunner::processMipOverlapGraph(
		const bib::progutils::CmdArgs & inputCommands) {

	mipCorePars pars;
	comparison noErrors;
	uint32_t minOverlap = 5;
	mipsterUtilsSetUp setUp(inputCommands);
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.setOption(minOverlap, "--minOverlap",
			"Minimum Overlap Allowed Between targets");
	setUp.processComparison(noErrors);
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gap_ = "5,1";
	setUp.processAlignerDefualts();
	setUp.processDirectoryOutputName("processMipOverlapGraph_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);




	auto regions = mipMaster.mips_->getMipRegionsForFams(mipMaster.names_->mips_);
	for(const auto & samp : mipMaster.names_->samples_){
		for(const auto & region : regions){
			MipOverlapGraph mog(region, noErrors, minOverlap);
			auto mipFams = mipMaster.mips_->getMipFamsForRegion(region);
			uint64_t maxSize = 0;
			for(const auto & mipFam : mipFams){
				if(bib::in(mipFam,mipMaster.names_->mips_ )){
					//std::cout << mipMaster.pathPopClusFinalHaplo(MipFamSamp(mipFam, samp)) << std::endl;

					auto seqOpts = SeqIOOptions::genFastqIn(mipMaster.pathPopClusFinalHaplo(MipFamSamp(mipFam, samp)).string(), true);
					seqInfo seq;
					SeqInput reader(seqOpts);
					reader.openIn();
					while(reader.readNextRead(seq)){
						readVec::getMaxLength(seq, maxSize);
						mog.addMipSeqToMipSeqMap(seq);
					}
				}
			}
			aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
					setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
			alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
			mog.genMipOverlapGraph(alignerObj);
			std::stringstream ss;
			mog.overlapGraph_->writePaths(ss);
			auto sampDirPath = bib::files::make_path(setUp.pars_.directoryName_, samp, region);
			bib::files::makeDirP(bib::files::MkdirPar(sampDirPath.string()));
			std::ofstream possibleHapsFile;
			std::ofstream lociInfoFile;
			openTextFile(possibleHapsFile,bib::files::join(sampDirPath.string(), "possibleHapsFile"),".txt",false, true );
			openTextFile(lociInfoFile,bib::files::join(sampDirPath.string(), "lociInfoFile"),".txt",false, true );
			OutOptions lociAlleNameKeyOpts(bib::files::make_path(sampDirPath, "nameKey.tab.txt"));
			std::ofstream lociAlleNameKeyFile;
			lociAlleNameKeyOpts.openFile(lociAlleNameKeyFile);
			lociAlleNameKeyFile << "seqName\tLociName" << std::endl;
			auto overlapPaths = streamToVecStr(ss);

			for(const auto & oPath : overlapPaths){
				auto toks = tokenizeString(oPath, " -> ");
				for(const auto & tok : toks){
					seqInfo tempSeq(tok);
					auto alleleNum = bib::lexical_cast<uint32_t>(tempSeq.getReadId());
					std::unordered_map<std::string, std::string> meta;
					tempSeq.processNameForMeta(meta);
					uint32_t lociNum = 0;
					if (bib::has(meta, "mipFam")) {
						auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
						lociNum = bib::lexical_cast<uint32_t>(
								bib::replaceString(mipFamToks.back(), "mip", ""));
					} else {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam"
								<< std::endl;
						throw std::runtime_error { ss.str() };
					}
					if(tok != toks.front()){
						possibleHapsFile << " ";
					}

					possibleHapsFile << "L" << lociNum + 1 << ".A" << alleleNum + 1;
				}
				possibleHapsFile << std::endl;
			}
			lociInfoFile << "LOCI " << mog.seqsByMipNum_.size() << std::endl;
			for(const auto & mipSubRegion : mog.seqsByMipNum_){
				lociInfoFile << "L" << mipSubRegion.first + 1 << " " << mipSubRegion.second.size() << std::endl;;

			}
			for (const auto & mipSubRegion : mog.seqsByMipNum_) {
				for (const auto & seq : mipSubRegion.second) {
					auto alleleNum = bib::lexical_cast<uint32_t>(seq->getReadId());
					lociInfoFile << "L" << mipSubRegion.first + 1 << " A" << alleleNum + 1
							<< " " << roundDecPlaces(seq->frac_, 3) << std::endl;
					lociAlleNameKeyFile << seq->name_ << "\t" << "L"
							<< mipSubRegion.first + 1 << ".A" << alleleNum + 1 << std::endl;
				}
			}
			alignerObj.processAlnInfoOutput(setUp.pars_.alnInfoDirName_,
					setUp.pars_.verbose_);
		}
	}



	return 0;
}

int mipsterUtilsRunner::processMipOverlapGraphSingle(const bib::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gap_ = "5,1";
	comparison noErrors;
	uint32_t minOverlap = 5;
	std::string regionName = "mip";
	setUp.setOption(minOverlap, "--minOverlap",
			"Minimum Overlap Allowed Between targets");
	setUp.setOption(regionName, "--regionName", "Region Name");

	setUp.processComparison(noErrors);
	setUp.processAlignerDefualts();
	setUp.processDefaultReader(true);
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "json";
	setUp.pars_.ioOptions_.out_.outExtention_ = ".json";
	setUp.processDebug();
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);


	uint64_t maxSize = 0;
	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto reads = reader.readAllReadsPtrs<seqInfo>();
	MipOverlapGraph mog(regionName, noErrors, minOverlap);
	for(const auto & seq : reads){
		readVec::getMaxLength(seq, maxSize);
		mog.addMipSeqToMipSeqMap(seq);
	}

	if(!setUp.pars_.ioOptions_.processed_){
		mog.resetFractions();
	}

	aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
			setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
	mog.genMipOverlapGraph(alignerObj);

	if(setUp.pars_.debug_){
		mog.overlapGraph_->writePaths(std::cout);
	}
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_);
	Json::FastWriter jWriter;
	auto graphJson = mog.genOverlapJson();
	outFile << jWriter.write(graphJson) << std::endl;
	alignerObj.processAlnInfoOutput(setUp.pars_.alnInfoDirName_,
			setUp.pars_.verbose_);
	return 0;
}


int mipsterUtilsRunner::scanForContam(
		const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterUtilsSetUp setUp(inputCommands);
	std::string mipFamName = "";
	mipCorePars pars;
	pars.processDefaults(setUp);
	std::string outFile = "out";
	setUp.setOption(outFile, "--outFilename", "Name of a file to print to");
	setUp.processWritingOptions();
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(mipFamName, "--mipName", "Name of mip to investigate", true);
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	table sharedTab(VecStr {"mipTarName", "SharedBar", "sample", "readCnt" });
	auto mipFamTars = mipMaster.mips_->getMipsForFamily(mipFamName);
	//key1 = samp, key2 = mipFam, val = vec of barcodes to remove
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> barcodesToRemove;
	for(const auto & mipTarName : mipFamTars){
		std::vector<MipFamSamp> barcodesExists;
		for (const auto & samp : mipMaster.names_->samples_) {
			MipFamSamp mipSamp(mipTarName, samp);
			if (bfs::exists(mipMaster.pathMipSampBarCorBars(MipFamSamp(mipFamName, samp), mipTarName))) {
				barcodesExists.emplace_back(mipSamp);
			}
		}

		std::unordered_map<std::string, strCounter> barCounts;
		for (const auto & bars : barcodesExists) {
			std::ifstream barFile(mipMaster.pathMipSampBarCorBars(MipFamSamp(mipFamName, bars.samp_), bars.mipFam_).string());
			if (!barFile) {
				std::stringstream ss;
				ss << "Error in opening " << mipMaster.pathMipSampBarCorBars(MipFamSamp(mipFamName, bars.samp_), bars.mipFam_)
						<< std::endl;
				throw std::runtime_error { ss.str() };
			}
			std::string line = "";
			while (bib::files::crossPlatGetline(barFile, line)) {
				if (!bib::beginsWith(line, "Barcode")) {
					auto toks = tokenizeString(line, "\t");
					if (2 != toks.size()) {
						std::stringstream ss;
						ss << "Error in reading barcode line: " << std::endl;
						ss << line << std::endl;
						ss << "lines should have two columns, this line had " << toks.size()
								<< std::endl;
						throw std::runtime_error { ss.str() };
					}
					barCounts[bars.samp_].increaseCountByString(toks[0],
							bib::lexical_cast<uint32_t>(toks[1]));
				}
			}
		}
		std::unordered_map<std::string, VecStr> sharedBars;
		for (const auto & sampCount : barCounts) {
			for (const auto & bars : sampCount.second.counts_) {
				sharedBars[bars.first].emplace_back(sampCount.first);
			}
		}
		for (const auto & shared : sharedBars) {
			if (shared.second.size() > 1) {
				//std::cout << "Shared Bar: " << shared.first << std::endl;
				for (const auto & samp : shared.second) {
					//std::cout << "\t" << samp << " " << barCounts[samp].counts_[shared.first] << std::endl;
					sharedTab.content_.emplace_back(
							toVecStr(mipTarName,shared.first, samp,
									barCounts[samp].counts_[shared.first]));
				}
			}
		}
	}


	sharedTab.outPutContents(
			TableIOOpts(
					OutOptions(outFile, ".tab.txt", "tab", false,
							setUp.pars_.ioOptions_.out_.overWriteFile_, false), "\t",
					sharedTab.hasHeader_));

	return 0;
}

int mipsterUtilsRunner::processProcessedMips(
		const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterUtilsSetUp setUp(inputCommands);
	std::string dirName = "";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapInfo_ = gapScoringParameters(setUp.pars_.gap_,
			setUp.pars_.gapLeft_, setUp.pars_.gapRight_);
	setUp.processAlignerDefualts();
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	setUp.processSeq(true);
	setUp.processDebug();
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	uint64_t maxLen = 0;
	readVec::getMaxLength(inReads, maxLen);
	readVec::getMaxLength(setUp.pars_.seqObj_, maxLen);

	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
	alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_);

	std::unordered_map<std::string, mipTargetReads> mReads = processMipReads(
			inReads, setUp.pars_.seqObj_, alignerObj, true, false);

	auto targetNames = getVectorOfMapKeys(mReads);
	bib::sort(targetNames);
	for (auto & targetReadName : targetNames) {
		auto & targetReads = mReads.at(targetReadName);
		if (setUp.pars_.debug_) {
			std::cout << targetReads.targetName_ << std::endl;
			std::cout << "\t" << "start: " << targetReads.start_ << std::endl;
			std::cout << "\t" << "stop: " << targetReads.stop_ << std::endl;
			std::cout << "\t" << "ReverseStrand: "
					<< bib::colorBool(targetReads.reverseStrand_) << std::endl;
			VecStr subReadsNames;
			for (const auto & read : targetReads.reads_) {
				subReadsNames.emplace_back(read->seqBase_.name_);
			}
			printVector(subReadsNames, ", ");
			std::cout << std::endl;
		}
	}
	comparison comp;
	mippedGene currentGene(setUp.pars_.seqObj_.seqBase_, mReads);
	currentGene.setUpGraph(alignerObj, comp);
	currentGene.setGroups();
	//currentGene.printEdges(std::cout);
	currentGene.printAllPaths(std::cout);
	currentGene.printLociInfo(std::cout);
	currentGene.printPossibleHaps(std::cout);
	//currentGene.printGroupInfo();

	if (setUp.pars_.writingOutAlnInfo_) {
		alignerObj.alnHolder_.write(setUp.pars_.alnInfoDirName_);
	}

	return 0;
}

int mipsterUtilsRunner::alignTargets(
		const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterUtilsSetUp setUp(inputCommands);
	std::string dirName = "";
	setUp.pars_.gapRight_ = "0,0";
	setUp.pars_.gapLeft_ = "0,0";
	setUp.pars_.gapInfo_ = gapScoringParameters(setUp.pars_.gap_,
			setUp.pars_.gapLeft_, setUp.pars_.gapRight_);
	setUp.processAlignerDefualts();
	setUp.pars_.ioOptions_.out_.outFilename_ = "out";
	setUp.processDefaultReader(true);
	setUp.processSeq(true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	reader.openIn();
	auto inReads = reader.readAllReads<readObject>();
	uint64_t maxLen = 0;
	readVec::getMaxLength(inReads, maxLen);
	readVec::getMaxLength(setUp.pars_.seqObj_, maxLen);
	aligner alignerObj(maxLen, setUp.pars_.gapInfo_, setUp.pars_.scoring_);
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_.outFilename_, ".fastq",
			setUp.pars_.ioOptions_.out_);

	setUp.pars_.seqObj_.seqBase_.outPutFastq(outFile);
	for (const auto & read : inReads) {
		alignerObj.alignCache(setUp.pars_.seqObj_, read, false);
		auto firstAlign = alignerObj.alignObjectB_;
		int32_t bestScore = alignerObj.parts_.score_;
		auto readComp = read;
		readComp.seqBase_.reverseComplementRead(true);
		alignerObj.alignCache(setUp.pars_.seqObj_, readComp, false);
		if (alignerObj.parts_.score_ > bestScore) {
			alignerObj.alignObjectB_.seqBase_.outPutFastq(outFile);
		} else {
			firstAlign.seqBase_.outPutFastq(outFile);
		}
	}
	return 0;
}

int mipsterUtilsRunner::rearmTargetsAndCombine(
		const bib::progutils::CmdArgs & inputCommands) {
	bfs::path inputDir = "";
	bfs::path outputDir = "";
	mipCorePars corePars;
	bool overWriteDir = false;
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.setOption(corePars.mipArmsFileName, "--mipArmsFilename",
				"Name of the mip arms file", true);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(inputDir, "--inputDir", "Input Directory", true);
	setUp.setOption(outputDir, "--outputDir", "Output Directory", true);
	setUp.setOption(overWriteDir, "--overWriteDir",
			"Over write output directory if it already exists");
	setUp.finishSetUp(std::cout);
	if (!bfs::exists(inputDir)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, input directory "
				<< bib::bashCT::boldRed(inputDir.string()) << " doesn't exist" << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (bfs::exists(outputDir) && !overWriteDir) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, output directory "
				<< bib::bashCT::boldRed(outputDir.string())
				<< " already exists, use --overWriteDir to over write " << "\n";
		throw std::runtime_error { ss.str() };
	}

	MipCollection mips(corePars.mipArmsFileName, corePars.allowableErrors);

	bib::files::MkdirPar outputPars(outputDir);
	outputPars.overWriteDir_ = true;
	bib::files::makeDir(outputPars);
	setUp.startARunLog(bib::appendAsNeededRet(outputDir.string(), "/"));

	auto files = bib::files::gatherFiles(inputDir, ".fasta", false);
	if(setUp.pars_.debug_){
		printVector(files, "\n");
	}

	std::unordered_map<std::string, std::vector<bfs::path>> inputByRegion;

	for(const auto & f : files){
		inputByRegion[f.filename().string().substr(0, f.filename().string().find("_"))].emplace_back(f);
	}

	uint64_t maxLen = 0;
	//region, genome, subRegion
	std::unordered_map<std::string,std::unordered_map<std::string, std::map<std::string, std::string>>> seqsByGenomeAndRegion;
	for(const auto & reg : inputByRegion){
		for(const auto & f : reg.second){
			SeqInput reader(SeqIOOptions::genFastaIn(f));
			auto seqs = reader.readAllReads<seqInfo>();
			for(auto & seq : seqs){
				seq.prepend(mips.mips_[bfs::basename(f)].extentionArmObj_.seqBase_.seq_);
				seq.append(mips.mips_[bfs::basename(f)].ligationArmObj_.seqBase_.seq_);
				auto genomes = tokenizeString(seq.name_, "-");
				readVec::getMaxLength(seq, maxLen);
				for(const auto & genome : genomes){
					seqsByGenomeAndRegion[reg.first][genome][bfs::basename(f)] = seq.seq_;
				}
			}
		}
	}
	std::unordered_map<std::string, std::vector<seqInfo>> seqsByGenome;
	for(const auto & region : seqsByGenomeAndRegion){
		for(const auto & genome : region.second){
			for(const auto & subReg : genome.second){
				seqsByGenome[genome.first].emplace_back(seqInfo(subReg.first, subReg.second));
			}
		}
	}

	for(auto & genome : seqsByGenome){
		MipNameSorter::sort(genome.second);
		SeqOutput::write(genome.second, SeqIOOptions::genFastaOut(bib::files::make_path(outputDir, genome.first)));
	}


	return 0;
}




                    
} // namespace bibseq
