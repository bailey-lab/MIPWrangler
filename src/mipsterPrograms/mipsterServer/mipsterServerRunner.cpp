
//  mipsterServerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterServerRunner.hpp"


namespace bibseq {

mipsterServerRunner::mipsterServerRunner()
    : bib::progutils::programRunner({addFunc("mipServerSetUp", mipServerSetUp, false),
																		 addFunc("mipAnalysisServerSetUp", mipAnalysisServerSetUp, false),
																		 addFunc("mav", mavRunner, false)},
                    "mipsterServer") {}//



int mipsterServerRunner::mipAnalysisServerSetUp(const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterServerSetUp setUp(inputCommands);
	mipCorePars pars;
	std::string name = "";
	setUp.setOption(pars.mipArmsFileName, "--mipArmsFilename", "Name of the mip arms file", true);
	setUp.setOption(pars.masterDir, "--masterDir","Name of main analysis directory", true);
	setUp.setOption(pars.mipsSamplesFile, "--mipSampleFile",
			"Mip sample filename, two columns, one column is mips, other is samples",
			true);
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(pars.seqFileSuffix, "--seqFileSuffix", "The raw data suffix");
	setUp.setOption(name, "--name", "Name of root of the server", true);
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setServerName(name);
	mipMaster.setRawDataSuffix(pars.seqFileSuffix);

	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}
	mipMaster.prepareMipAnalysisServer(pars.numThreads);
	return 0;
}

int mipsterServerRunner::mipServerSetUp(const bib::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterServerSetUp setUp(inputCommands);
	std::string dirName = "";
	std::string outDirName = "";
	std::string refSeqGenesFilename = "";
	std::string twoBitFileName = "";
	std::string genomicDir = "";
	std::string sampNamefile = "";
	setUp.setOption(sampNamefile, "--file", "Filename of mips and samples", true);
	setUp.setOption(refSeqGenesFilename, "-refSeqGeneFile", "Name of the ref seq gene file (zero based, non-inclusive end positions)");
	setUp.setOption(genomicDir, "-genomicDir", "Name of the genomic directory");
	setUp.setOption(twoBitFileName, "-twoBitFile", "Name of the two bit file to extract dna from using the seqGeneFile ");
	setUp.setOption(outDirName, "--outDirName", "Directory Name Of An Existing Output directory", true);
	setUp.setOption(dirName, "--dirName", "Directory Name Containing Analysis", true);
	setUp.processVerbose();
	bool human = false;
	setUp.setOption(human, "--human", "Human Mips");
	setUp.finishSetUp(std::cout);
	if(!bib::files::bfs::exists(outDirName)){
		std::cerr << "Error, --outDirName should be an exisitng directory, " << outDirName << "doesn't exist" << std::endl;
		exit(1);
	}
  table mipSampInfo(sampNamefile, "\t", true);
  auto mips = mipSampInfo.getColumn("mips");
  auto samps = mipSampInfo.getColumn("samples");
  removeElement(samps, std::string(""));
  removeElement(mips, std::string(""));
	std::unordered_map<std::string, bfs::path> mipAnalysisFolders;
	std::unordered_map<std::string, bfs::path> sampAnalysisFolders;
	std::unordered_map<std::string, std::unordered_map<std::string, bfs::path>> sampleMipAnalysisFolders;
	std::unordered_map<std::string, bibseq::table> allInfoBySample;
	table stats;
	table mipStats;
	bib::stopWatch watch;
	bib::appendAsNeeded(dirName, "/");
	std::regex analysisPat { ".*Analysis\\b" };
	auto files = bib::files::listAllFiles(dirName, false, { analysisPat });
	printOutMapContents(files, "\t", std::cout);
	bfs::path analysisFolder = "";

	if (files.empty()) {
		std::cout
				<< "Error, need a folder that ends with Analysis that contains the analysis folders"
				<< std::endl;
		exit(1);
	} else {
		analysisFolder = files.begin()->first;
	}
	std::regex mipFolderPat;
	if(human){
		mipFolderPat = std::regex{ ".*_C\\d+\\b" };
	}else{
		mipFolderPat = std::regex{ ".*_(mip)*\\d+\\b" };
	}

	if(setUp.pars_.verbose_){
		std::cout << "Analysis folder: " << analysisFolder << std::endl;
	}
	watch.startNewLap("Find mip analysis folders");
	auto mipFinalAnalysisFolders = bib::files::listAllFiles(
			analysisFolder.string(), false, { mipFolderPat });
	for (const auto & mip : mipFinalAnalysisFolders) {
		if (mip.second) {
			if(setUp.pars_.verbose_){
				std::cout << "Currently on mip folder: " << mip.first << std::endl;
			}
			//load analysis
			if (bfs::is_directory(bib::appendAsNeededRet(mip.first.string(), "/") + "analysis")) {
				mipAnalysisFolders[mip.first.filename().string()] = bfs::path(
						bib::appendAsNeededRet(mip.first.string(), "/") + "analysis");
			} else {
				std::cerr << "Error, did not find a folder called analysis under "
						<< bib::appendAsNeededRet(mip.first.string(), "/") << " for mip target "
						<< mip.first.filename().string() << std::endl;
				exit(1);
			}
		}
	}
	watch.startNewLap("Gather extraction stats");
	if (!bib::files::bfs::exists(outDirName + "perMipStats.tab.txt")) {
		mipStats = getSampleMipStats(dirName, setUp.pars_.verbose_, samps);
		mipStats.outPutContents(
				TableIOOpts(OutOptions(outDirName + "perMipStats.tab.txt", ".tab.txt", "tab",
						false,true, false), "\t", mipStats.hasHeader_));
	} else {
		mipStats = table(outDirName + "perMipStats.tab.txt", "\t", true);
	}

	if (!bib::files::bfs::exists(outDirName + "stats.tab.txt")) {
		stats = getSampleStats(dirName, setUp.pars_.verbose_, samps);
		stats.outPutContents(
				TableIOOpts(OutOptions(outDirName + "stats.tab.txt", ".tab.txt", "tab", false, true, false), "\t", stats.hasHeader_));
	} else {
		stats = table(outDirName + "stats.tab.txt", "\t", true);
	}
	auto readsUsed = stats.getColumn("readsUsed");
	auto sampNames = stats.getColumn("sampleName");
	watch.startNewLap("Read in extraction info files");
	for (const auto & e : iter::enumerate(readsUsed)) {
		if (estd::stou(e.element) > 0) {
			std::string sampName = sampNames[e.index];
			if(setUp.pars_.verbose_){
				std::cout << "Currently on sample: " << sampName << std::endl;
			}
			auto extractionDirs = bib::files::listAllFiles(dirName + sampName,
					false, VecStr { sampName + ".assembled_mip" });
			std::map<bfs::path, bool, std::greater<bfs::path>> eDirs;
			for (const auto & ed : extractionDirs) {
				if (ed.second) {
					eDirs.emplace(ed);
				}
			}
			if (!eDirs.empty()) {
				auto resultsDir = bib::appendAsNeededRet(eDirs.begin()->first.string(), "/");
				sampAnalysisFolders[sampName] = (eDirs.begin()->first);
				auto mipFolders = bib::files::listAllFiles(
						eDirs.begin()->first.string(), false, { mipFolderPat });
				for (const auto & m : mipFolders) {
					if (m.second) {
						sampleMipAnalysisFolders[sampName][m.first.filename().string()] =
								m.first;
					}
				}
			} else {
				std::cout << "Couldn't find latest analysis file for sample "
						<< sampName << std::endl;
			}
		}
	}

	std::unordered_map<std::string, std::set<std::string>> sampNamesForGeneSet;
	watch.startNewLap("Read in analysis results");
	bfs::path allInfoBySampleDir = bib::files::makeDir(outDirName, bib::files::MkdirPar("allInfoBySample"));
	bfs::path namesDir = bib::files::makeDir(outDirName, bib::files::MkdirPar("names"));
	std::ofstream mipAnalysisFoldersFile;
	openTextFile(mipAnalysisFoldersFile,bib::files::make_path(namesDir, "mipAnalysisFolders").string(),".tab.txt", false, true);

	for (const auto & mipAnalysis : mipAnalysisFolders) {
		mipAnalysisFoldersFile << mipAnalysis.first << "\t" << mipAnalysis.second.string() << "\n";
		if (!bfs::exists(
				bfs::path(
						mipAnalysis.second.string() + "/selectedClustersInfo.tab.txt"))) {
			continue;
		}
		if(setUp.pars_.verbose_){
			std::cout << "Currently on mip : " << mipAnalysis.first << std::endl;
		}
		auto tab = bibseq::table(
				mipAnalysis.second.string() + "/selectedClustersInfo.tab.txt", "\t",
				true);
		auto expNames = tab.getColumn("h_popUID");
		if (expNames.empty()) {
			continue;
		}
		auto expName = expNames.front().substr(0, expNames.front().find("_"));
		auto targetName = expNames.front().substr(0, expNames.front().find("."));
		tab.columnNames_.insert(tab.columnNames_.begin(), "geneName");
		tab.columnNames_.insert(tab.columnNames_.begin(), "mipName");
		for (auto & row : tab.content_) {
			row.insert(row.begin(), expName);
			row.insert(row.begin(), targetName);
		}
		auto split = tab.splitTableOnColumn("s_sName");
		for (const auto & s : split) {
			auto search = allInfoBySample.find(s.first);
			if (search == allInfoBySample.end()) {
				allInfoBySample[s.first] = s.second;
			} else {
				allInfoBySample[s.first].rbind(s.second, true);
			}
		}
	}
	std::unordered_map<std::string, VecStr> sampNamesForMip;
	std::unordered_map<std::string, VecStr> sampNamesForGene;
	std::unordered_map<std::string, VecStr> mipNamesForSamp;
	std::unordered_map<std::string, VecStr> geneNamesForSamp;
	for(auto & samp: allInfoBySample){
		samp.second.sortTable("mipName", false);
		samp.second.outPutContents(TableIOOpts(OutOptions(bib::files::make_path(allInfoBySampleDir , samp.first).string(), ".tab.txt"), "\t", samp.second.hasHeader_));

		auto mipNames = samp.second.getColumnLevels("mipName");
		auto geneNames = samp.second.getColumnLevels("geneName");
		for (const auto & mip : mipNames) {
			sampNamesForMip[mip].emplace_back(samp.first);
			mipNamesForSamp[samp.first].emplace_back(mip);
		}
		for (const auto & gene : geneNames) {
			sampNamesForGene[gene].emplace_back(samp.first);
			geneNamesForSamp[samp.first].emplace_back(gene);
		}
	}

	watch.startNewLap("Output");


	std::ofstream sampAnalysisFoldersFile;
	openTextFile(sampAnalysisFoldersFile,bib::files::make_path(namesDir, "sampAnalysisFolders").string(),".tab.txt", false, true);
	for(const auto & sampAnalysis : sampAnalysisFolders){
		sampAnalysisFoldersFile << sampAnalysis.first << "\t" << sampAnalysis.second.string() << "\n";
	}
	std::ofstream sampleMipAnalysisFoldersFile;
	openTextFile(sampleMipAnalysisFoldersFile,bib::files::make_path(namesDir, "sampleMipAnalysisFolders").string(),".tab.txt", false, true);
	for(const auto & sampAnalysis : sampleMipAnalysisFolders){
		for(const auto & mipAnalysis : sampAnalysis.second){
			sampleMipAnalysisFoldersFile << sampAnalysis.first << "\t" << mipAnalysis.first << "\t" << mipAnalysis.second.string() << "\n";
		}
	}
	std::ofstream sampNamesForMipFile;
	openTextFile(sampNamesForMipFile, bib::files::make_path(namesDir, "sampNamesForMip").string(),".tab.txt", false, true);
	for (auto & names : sampNamesForMip) {
		bib::sort(names.second);
		sampNamesForMipFile << names.first << "\t" << vectorToString(names.second, ",") << "\n";
	}
	std::ofstream mipNamesForSampFile;
	openTextFile(mipNamesForSampFile,bib::files::make_path(namesDir , "mipNamesForSamp").string(),".tab.txt", false, true);
	for (auto & names : mipNamesForSamp) {
		bib::sort(names.second);
		mipNamesForSampFile << names.first << "\t" << vectorToString(names.second, ",") << "\n";
	}
	std::ofstream sampNamesForGeneFile;
	openTextFile(sampNamesForGeneFile,bib::files::make_path(namesDir , "sampNamesForGene").string(),".tab.txt", false, true);
	for (auto & names : sampNamesForGene) {
		bib::sort(names.second);
		sampNamesForGeneFile << names.first << "\t" << vectorToString(names.second, ",") << "\n";
	}
	std::ofstream geneNamesForSampFile;
	openTextFile(geneNamesForSampFile,bib::files::make_path(namesDir , "geneNamesForSamp").string(),".tab.txt", false, true);
	for (auto & names : geneNamesForSamp) {
		bib::sort(names.second);
		geneNamesForSampFile << names.first << "\t" << vectorToString(names.second, ",") << "\n";
	}

	watch.startNewLap("snp_info");
	if(!bfs::exists(outDirName + "allSeqSnps.txt") && refSeqGenesFilename != ""){
		std::ofstream outFile(bib::files::findNonexitantFile(outDirName + "allSeqSnps.txt").string());
		std::ofstream outFileProt(bib::files::findNonexitantFile(outDirName + "allProteinAAChange.tab.txt").string());

		auto mipNames = bibseq::getVectorOfMapKeys(mipAnalysisFolders);
		bib::sort(mipNames);
		std::set<std::string> geneNames;
		for (const auto & mPos : iter::range(mipNames.size())) {
			geneNames.emplace(
					mipNames[mPos].substr(0, mipNames[mPos].find_last_of("_")));
		}
		VecStr geneNamesVec(geneNames.begin(), geneNames.end());
		uint32_t snpTabCount = 0;
		uint32_t snpProtTabCount = 0;
		for (const auto & geneName : geneNamesVec) {
			std::cout << geneName << std::endl;
			for (const auto & sampName : sampNamesForGene[geneName]) {
				std::cout << sampName << std::endl;
				auto gSearch = sampNamesForGene.find(geneName);
				if (gSearch == sampNamesForGene.end()) {
					std::cout << "oneGeneOneSampSnpData: couldn't find gene: " << geneName
							<< " redirecting" << std::endl;
				} else {
					if (!bib::in(sampName, gSearch->second)) {
						std::cout << "oneGeneOneSampSnpData: couldn't find samp: " << sampName
								<< " for gene name " << geneName << " redirecting" << std::endl;
					} else {
						auto geneRecords = getRefSeqRecs(refSeqGenesFilename, VecStr { });
						if (!bib::contains(geneRecords, geneName)) {
							//if(false){
							std::cout
									<< "oneGeneOneSampSnpData: couldn't find genomic seq info for "
									<< geneName << " in " << refSeqGenesFilename << std::endl;
						} else {
							auto record = geneRecords.at(geneName);
							uint64_t maxLen = 0;
							//find genome data
							//std::cout << "Reading in two bit file: " << twoBitFileGenome_ << std::endl;
							//TwoBit::TwoBitFile twoBitFile(twoBitFileGenome_);
							//auto seqNames = twoBitFile.sequenceNames();
							//printVector(seqNames);
							std::string refGenome = genomicDir + geneName + "_genomic.fasta";
							//std::cout << record->toStrLine() << std::endl;
							if (!bfs::exists(refGenome)) {
								//if(!bib::in(record->chrom_, seqNames)){
								//if(false){
								//std::cout << "showOneGeneOneSamp: couldn't find chrom: " << record->chrom_ << " in " << vectorToString(seqNames, ",") << std::endl;
								std::cout << "showOneGeneOneSamp: couldn't refGenome file: "
										<< refGenome << " in " << genomicDir << std::endl;
							} else {
								std::vector<readObject> allReads;
								//find all mip target data
								for (const auto & mipAnalysis : mipAnalysisFolders) {
									if (bib::beginsWith(mipAnalysis.first, geneName)) {
										std::string seqFilename = mipAnalysis.second.string()
												+ "/final/" + sampName + ".fastq";
										if (bfs::exists(seqFilename)) {
											/**@todo fix this */
											SeqIOOptions opts =SeqIOOptions::genFastqIn(seqFilename);
											addOtherVec(allReads,
													SeqInput::getReferenceSeq(opts, maxLen));
										}
									}
								}
								//readVec::allPrintSeqs(allReads, std::cout);
								if (allReads.empty()) {
									std::cout << "showOneGeneOneSamp: "
											<< "couldn't find any reads for samp: " << sampName
											<< " for geneName: " << geneName << std::endl;
								} else {
									SeqIOOptions opts =SeqIOOptions::genFastaIn(refGenome);
									SeqInput genomeReader(opts);
									genomeReader.openIn();
									readObject refSeqObj;
									genomeReader.readNextRead(refSeqObj);
									//align and return
									readVecSorter::sortReadVector(allReads, "name", false);
									readVec::getMaxLength(refSeqObj, maxLen);
									aligner alignerObj(maxLen,
											gapScoringParameters(5, 1, 0, 0, 0, 0),
											substituteMatrix(2, -2));
									//alignerObj.parts_.scoring_.printScores(std::cout);
									auto mReads = processMipReads(allReads, refSeqObj, alignerObj,
											false, true);
									mippedGene currentGene(refSeqObj.seqBase_, mReads);
									currentGene.refGeneRecord_ = record;


									auto alns = currentGene.getAlignedTargets(alignerObj);
									auto snps = currentGene.findSnps(alignerObj);
									auto sKeys = getVectorOfMapKeys(snps);
									bib::sort(sKeys);

									for (const auto & k : sKeys) {
										//std::cout << k << std::endl;
										if (snps[k].content_.empty()) {
											continue;
										}
										if(k != "all"){
											continue;
										}
										auto outTab = snps[k];
										outTab.addColumn(VecStr { sampName }, "Sample");
										outTab.addColumn(VecStr { geneName }, "Gene");
										outTab.setColNamePositions();
										if(snpTabCount != 0){
											outTab.hasHeader_ = false;

										}
										++snpTabCount;
										outTab.outPutContents(outFile, "\t");
									}
									auto psnps = currentGene.findProteinSnps(alignerObj);
									auto pKeys = getVectorOfMapKeys(psnps);
									bib::sort(pKeys);
									for (const auto & k : pKeys) {
										if (psnps[k].content_.empty()) {
											continue;
										}
										if(k != "all"){
											continue;
										}
										//std::cout << k << std::endl;
										auto outTab = psnps[k];
										outTab.addColumn(VecStr { sampName }, "Sample");
										outTab.addColumn(VecStr { geneName }, "Gene");
										outTab.setColNamePositions();

										if(snpProtTabCount != 0){
											outTab.hasHeader_ = false;

										}
										++snpProtTabCount;
										outTab.outPutContents(outFileProt, "\t");
									}
								}
							}
						}
					}
				}
			}
		}
		outFile.close();
		outFileProt.close();
	}


	watch.logLapTimes(std::cout, true, 6, true);
	return 0;
}
                    
} // namespace bibseq
