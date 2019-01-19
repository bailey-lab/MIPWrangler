
//  mipsterUtilsRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

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
#include "mipsterUtilsRunner.hpp"
//#include <elucidator/seqToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>


#include <unordered_map>


namespace njhseq {




mipsterUtilsRunner::mipsterUtilsRunner()
    : njh::progutils::ProgramRunner({addFunc("alignTarget", alignTargets, true),
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 addFunc("processProcessedMips", processProcessedMips, true),
																		 addFunc("scanForContam", scanForContam, true),
																		 addFunc("processMipOverlapGraph", processMipOverlapGraph, false),
																		 addFunc("processMipOverlapGraphSingle", processMipOverlapGraphSingle, false),
																		 addFunc("rearmTargetsAndCombine", rearmTargetsAndCombine, false),
																		 addFunc("createExtArmFastas", createExtArmFastas, false),
																		 addFunc("createLigArmFastas", createLigArmFastas, false),
																		 addFunc("mipFastasToSeqTable", mipFastasToSeqTable, false),
																		 addFunc("extractPossibleMipCapturesFromGenome", extractPossibleMipCapturesFromGenome, true),
																		 addFunc("creatingMipArmsFromSeqs", creatingMipArmsFromSeqs, true),
																		 addFunc("fixingMipBedFiles", fixingMipBedFiles, true),
																		 addFunc("writeOutPossibleHaplotypes", writeOutPossibleHaplotypes, false),
																		 addFunc("creatingSeqTableFromDirectory", creatingSeqTableFromDirectory, true),
																		 addFunc("createMipArmFromSelectedMips", createMipArmFromSelectedMips, true),
																		 addFunc("typeFinalHaplotypes", typeFinalHaplotypes, false),
																		 addFunc("createPrimerFileFromArmFile", createPrimerFileFromArmFile, false),
																		 addFunc("ExtractTargetsFromGenomes", ExtractTargetsFromGenomes, true),
																		 addFunc("benchmarkingForControlMixtures", benchmarkingForControlMixtures, false),
},
                    "MIPWranglerUtils") {}

//


//std::vector<std::vector<SeqOverlapGraph::node>> getPossibleHaplotypes(const MipOverlapGraph & mog){
//	/**@todo add amount of overlap in edge
//	 *
//	 */
//	std::vector<std::shared_ptr<SeqOverlapGraph::node>> nodesToProcess;
//	for(const auto & n : mog.overlapGraph_->nodes_){
//		if(n.second->headless() || n.second->headEdges_.size() > 1){
//			nodesToProcess.emplace_back(n);
//		}
//	}
//	std::vector<std::vector<SeqOverlapGraph::node>> ret;
//	for(const auto & n : nodesToProcess){
//		std::vector<SeqOverlapGraph::node> currentPath;
//		auto currentNode = n;
//		while(n->tailEdges_.size() == 1){
//			currentPath.emplace_back(currentNode);
//			currentNode = n->tailEdges_.front()->tail_.lock();
//		}
//		ret.emplace_back(currentPath);
//	}
//	return ret;
//}




int mipsterUtilsRunner::createPrimerFileFromArmFile(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path mipArmsFnp;
	OutOptions outOpts(bfs::path(""));
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.setOption(mipArmsFnp, "--mipArmsFilename",
				"Name of the mip arms file", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	MipCollection mips(mipArmsFnp, 6);
	OutputStream out(outOpts);
	out << "target\tforward\treverse" << std::endl;
	auto tarNames = mips.getMipTars();
	for(const auto & mipName : tarNames){

		out << mipName
				<< "\t" << mips.mips_[mipName].extentionArm_
				<< "\t" << mips.mips_[mipName].ligationArm5_to_3prime_ << std::endl;
	}

	return 0;
}
int mipsterUtilsRunner::writeOutPossibleHaplotypes(
		const njh::progutils::CmdArgs & inputCommands) {
	double freqCutOff = 0;
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
	setUp.setOption(freqCutOff, "--freqCutOff", "freq Cut Off to exclude sequences");
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
		ss << njh::conToStr(warnings, "\n") << std::endl;
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
				if(njh::in(mipFam,mipMaster.names_->mips_ )){
					//std::cout << mipMaster.pathPopClusFinalHaplo(MipFamSamp(mipFam, samp)) << std::endl;

					auto seqOpts = SeqIOOptions::genFastqIn(mipMaster.pathPopClusFinalHaplo(MipFamSamp(mipFam, samp)).string(), true);
					seqInfo seq;
					SeqInput reader(seqOpts);
					reader.openIn();
					while(reader.readNextRead(seq)){
						if(seq.frac_ >freqCutOff){
							readVec::getMaxLength(seq, maxSize);
							mog.addMipSeqToMipSeqMap(seq);
						}
					}
				}
			}
			aligner alignerObj(maxSize, setUp.pars_.gapInfo_, setUp.pars_.scoring_,
					setUp.pars_.colOpts_.alignOpts_.countEndGaps_);
			alignerObj.processAlnInfoInput(setUp.pars_.alnInfoDirName_, setUp.pars_.verbose_);
			mog.genMipOverlapGraph(alignerObj);

			std::stringstream ss;
			mog.overlapGraph_->writePaths(ss);
			auto sampDirPath = njh::files::make_path(setUp.pars_.directoryName_, samp, region);
			njh::files::makeDirP(njh::files::MkdirPar(sampDirPath.string()));
			std::ofstream possibleHapsFile;
			std::ofstream lociInfoFile;
			openTextFile(possibleHapsFile,njh::files::join(sampDirPath.string(), "possibleHapsFile"),".txt",false, true );
			openTextFile(lociInfoFile,njh::files::join(sampDirPath.string(), "lociInfoFile"),".txt",false, true );
			OutOptions lociAlleNameKeyOpts(njh::files::make_path(sampDirPath, "nameKey.tab.txt"));
			std::ofstream lociAlleNameKeyFile;
			lociAlleNameKeyOpts.openFile(lociAlleNameKeyFile);
			lociAlleNameKeyFile << "seqName\tLociName" << std::endl;
			auto overlapPaths = streamToVecStr(ss);
			auto seqOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(sampDirPath, "allOverlapsStitched.fasta"));
			SeqOutput writer(seqOutOpts);
			writer.openOut();
			for(const auto & oPath : overlapPaths){
				auto toks = tokenizeString(oPath, " -> ");
				if(toks.size() > 1){
					seqInfo buildingSeq = *(mog.overlapGraph_->nodes_.at(toks[0])->val_);
					MetaDataInName meta(buildingSeq.name_);
					buildingSeq.name_ = meta.getMeta("h_popUID");
					for(const auto & tokPos : iter::range<size_t>(0,toks.size()  - 1)){
						alignerObj.alignCacheGlobal(mog.overlapGraph_->nodes_.at(toks[tokPos])->val_,
								mog.overlapGraph_->nodes_.at(toks[tokPos + 1])->val_);

						auto lastGap = alignerObj.alignObjectA_.seqBase_.seq_.find_last_not_of('-');
						if(std::string::npos != lastGap ){
							auto nextSeqPos = alignerObj.getSeqPosForAlnBPos(lastGap);
							buildingSeq.append(mog.overlapGraph_->nodes_.at(toks[tokPos + 1])->val_->getSubRead(nextSeqPos));

							MetaDataInName meta(mog.overlapGraph_->nodes_.at(toks[tokPos + 1])->val_->name_);
							buildingSeq.name_.append("-" + meta.getMeta("h_popUID"));
						}
					}
					writer.write(buildingSeq);
				}
			}

			for(const auto & oPath : overlapPaths){
				auto toks = tokenizeString(oPath, " -> ");
				for(const auto & tok : toks){

					seqInfo tempSeq(tok);
					auto alleleNum = estd::stou(tempSeq.getReadId());
					std::unordered_map<std::string, std::string> meta;
					tempSeq.processNameForMeta(meta);
					uint32_t lociNum = 0;
					if (njh::has(meta, "mipFam")) {
						auto mipFamToks = njh::tokenizeString(meta.at("mipFam"), "_");
						lociNum = estd::stou(
								njh::replaceString(mipFamToks.back(), "mip", ""));
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
					auto alleleNum = estd::stou(seq->getReadId());
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


int mipsterUtilsRunner::processMipOverlapGraph(
		const njh::progutils::CmdArgs & inputCommands) {

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
		ss << njh::conToStr(warnings, "\n") << std::endl;
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
				if(njh::in(mipFam,mipMaster.names_->mips_ )){
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
			auto sampDirPath = njh::files::make_path(setUp.pars_.directoryName_, samp, region);
			njh::files::makeDirP(njh::files::MkdirPar(sampDirPath.string()));
			std::ofstream possibleHapsFile;
			std::ofstream lociInfoFile;
			openTextFile(possibleHapsFile,njh::files::join(sampDirPath.string(), "possibleHapsFile"),".txt",false, true );
			openTextFile(lociInfoFile,njh::files::join(sampDirPath.string(), "lociInfoFile"),".txt",false, true );
			OutOptions lociAlleNameKeyOpts(njh::files::make_path(sampDirPath, "nameKey.tab.txt"));
			std::ofstream lociAlleNameKeyFile;
			lociAlleNameKeyOpts.openFile(lociAlleNameKeyFile);
			lociAlleNameKeyFile << "seqName\tLociName" << std::endl;
			auto overlapPaths = streamToVecStr(ss);

			for(const auto & oPath : overlapPaths){
				auto toks = tokenizeString(oPath, " -> ");
				for(const auto & tok : toks){
					seqInfo tempSeq(tok);
					auto alleleNum = estd::stou(tempSeq.getReadId());
					std::unordered_map<std::string, std::string> meta;
					tempSeq.processNameForMeta(meta);
					uint32_t lociNum = 0;
					if (njh::has(meta, "mipFam")) {
						auto mipFamToks = njh::tokenizeString(meta.at("mipFam"), "_");
						lociNum = estd::stou(
								njh::replaceString(mipFamToks.back(), "mip", ""));
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
					auto alleleNum = estd::stou(seq->getReadId());
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

int mipsterUtilsRunner::processMipOverlapGraphSingle(const njh::progutils::CmdArgs & inputCommands) {
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
	auto graphJson = mog.genOverlapJson();
	outFile << njh::json::writeAsOneLine(graphJson) << std::endl;
	alignerObj.processAlnInfoOutput(setUp.pars_.alnInfoDirName_,
			setUp.pars_.verbose_);
	return 0;
}


int mipsterUtilsRunner::scanForContam(
		const njh::progutils::CmdArgs & inputCommands) {
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
		ss << njh::conToStr(warnings, "\n") << std::endl;
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
			while (njh::files::crossPlatGetline(barFile, line)) {
				if (!njh::beginsWith(line, "Barcode")) {
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
							estd::stou(toks[1]));
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
		const njh::progutils::CmdArgs & inputCommands) {
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
	njh::sort(targetNames);
	for (auto & targetReadName : targetNames) {
		auto & targetReads = mReads.at(targetReadName);
		if (setUp.pars_.debug_) {
			std::cout << targetReads.targetName_ << std::endl;
			std::cout << "\t" << "start: " << targetReads.start_ << std::endl;
			std::cout << "\t" << "stop: " << targetReads.stop_ << std::endl;
			std::cout << "\t" << "ReverseStrand: "
					<< njh::colorBool(targetReads.reverseStrand_) << std::endl;
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
		const njh::progutils::CmdArgs & inputCommands) {
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
		const njh::progutils::CmdArgs & inputCommands) {
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
				<< njh::bashCT::boldRed(inputDir.string()) << " doesn't exist" << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (bfs::exists(outputDir) && !overWriteDir) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error, output directory "
				<< njh::bashCT::boldRed(outputDir.string())
				<< " already exists, use --overWriteDir to over write " << "\n";
		throw std::runtime_error { ss.str() };
	}

	MipCollection mips(corePars.mipArmsFileName, corePars.allowableErrors);

	njh::files::MkdirPar outputPars(outputDir);
	outputPars.overWriteDir_ = true;
	njh::files::makeDir(outputPars);
	setUp.startARunLog(njh::appendAsNeededRet(outputDir.string(), "/"));

	auto files = njh::files::gatherFiles(inputDir, ".fasta", false);
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
				seq.prepend(mips.mips_[bfs::basename(f)].extentionArmObj_.seq_);
				seq.append(mips.mips_[bfs::basename(f)].ligationArmObj_.seq_);
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
		SeqOutput::write(genome.second, SeqIOOptions::genFastaOut(njh::files::make_path(outputDir, genome.first)));
	}


	return 0;
}


int mipsterUtilsRunner::createLigArmFastas(
		const njh::progutils::CmdArgs & inputCommands) {
	mipCorePars corePars;
	bool overWrite = false;
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(overWrite, "--overWrite", "Over write files");
	setUp.setOption(corePars.mipArmsFileName, "--mipArmsFilename",
				"Name of the mip arms file", true);

	setUp.finishSetUp(std::cout);


	MipCollection mips(corePars.mipArmsFileName, corePars.allowableErrors);

	for(const auto & m : mips.mips_){
		auto ligOpts = SeqIOOptions::genFastaOut(m.first + "-lig");
		ligOpts.out_.overWriteFile_ = overWrite;
		auto ligArmObj = m.second.ligationArmObj_;
		ligArmObj.reverseComplementRead(false, true);
		SeqOutput::write(std::vector<seqInfo>{ligArmObj}, ligOpts);
	}
	return 0;
}

int mipsterUtilsRunner::createExtArmFastas(
		const njh::progutils::CmdArgs & inputCommands) {
	mipCorePars corePars;
	bool overWrite = false;
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(overWrite, "--overWrite", "Over write files");
	setUp.setOption(corePars.mipArmsFileName, "--mipArmsFilename",
				"Name of the mip arms file", true);

	setUp.finishSetUp(std::cout);


	MipCollection mips(corePars.mipArmsFileName, corePars.allowableErrors);

	for(const auto & m : mips.mips_){
		auto extOpts = SeqIOOptions::genFastaOut(m.first + "-ext");
		extOpts.out_.overWriteFile_ = overWrite;
		auto extArmObj = m.second.extentionArmObj_;
		SeqOutput::write(std::vector<seqInfo>{extArmObj}, extOpts);
	}
	return 0;
}


int mipsterUtilsRunner::mipFastasToSeqTable(const njh::progutils::CmdArgs & inputCommands){
	std::string filePat = "^trimmed_.*.fasta$";
	bfs::path inputDir = "./";
	auto outOpts = TableIOOpts::genTabFileOut("out.tab.txt", true);
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.processWritingOptions(outOpts.out_);
	setUp.setOption(filePat, "--filePat", "File pattern to search for");
	setUp.setOption(inputDir, "--inputDir", "Input directory to search");
	setUp.finishSetUp(std::cout);

	auto files = njh::files::listAllFiles(inputDir, false, {std::regex{filePat}});
	table outTab(VecStr{"target", "name","length", "seq"});
	for(const auto & f : files){
		if(setUp.pars_.verbose_){
			std::cout << f.first << std::endl;
		}
		auto seqOpts = SeqIOOptsWithTime(SeqIOOptions::genFastaIn(f.first, false));
		auto seqs = seqOpts.get<seqInfo>();
		for(const auto & seq : seqs){
			if(setUp.pars_.verbose_ && setUp.pars_.debug_){
				std::cout << seq.name_ << std::endl;
			}
			auto bName = bfs::basename(f.first);
			bName = njh::replaceString(bName, "trimmed_", "");
			outTab.addRow(bName, seq.name_,len(seq), seq.seq_);
		}
	}
	outTab.outPutContents(outOpts);
	return 0;
}

int mipsterUtilsRunner::extractPossibleMipCapturesFromGenome(const njh::progutils::CmdArgs & inputCommands){
	comparison comp;
	bfs::path ligBamFnp = "";
	bfs::path extBamFnp = "";
	size_t insertSize = 1000;
	bfs::path twoBitFnp = "";
	OutOptions outOpts(bfs::path(""));

	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(twoBitFnp, "--twoBitFnp", "Two Bit file path", true);
	setUp.setOption(ligBamFnp, "--ligBamFnp", "lig Bam file path", true);
	setUp.setOption(extBamFnp, "--extBamFnp", "ext Bam file path", true);
	setUp.processComparison(comp);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	std::vector<std::shared_ptr<AlignmentResults>> alnResultsExt = gatherMapResults(
			extBamFnp, twoBitFnp, comp);
	std::vector<std::shared_ptr<AlignmentResults>> alnResultsLig = gatherMapResults(
			ligBamFnp, twoBitFnp, comp);

	auto extractResults = getPossibleGenomeExtracts(alnResultsExt, alnResultsLig, insertSize);
	if(setUp.pars_.debug_){
		std::cout << "Found " << extractResults.size() << " possible extractions" << std::endl;
	}

	std::ofstream outfile;
	std::ostream out(determineOutBuf(outfile, outOpts));

	for(const auto & extract : extractResults){
		out << extract.gRegion_->genBedRecordCore().toDelimStr() << std::endl;
	}
	return 0;
}

int mipsterUtilsRunner::creatingMipArmsFromSeqs(const njh::progutils::CmdArgs & inputCommands){
	uint32_t extensionArmSize = 15;
	uint32_t ligationArmSize = 15;
	uint32_t extensionBarSize = 10;
	uint32_t ligationBarSize = 4;
	std::string mipSet = "";
	std::string nameSuffix = "";
	seqSetUp setUp(inputCommands);
	setUp.setOption(extensionArmSize, "--extensionArmSize","Extension Arm Size");
	setUp.setOption(ligationArmSize,  "--ligationArmSize", "Ligation Arm Size");
	setUp.setOption(extensionBarSize, "--extensionBarSize","Extension Bar Size");
	setUp.setOption(ligationBarSize,  "--ligationBarSize", "Ligation Bar size");
	setUp.setOption(nameSuffix, "--nameSuffix","Name Suffix to the name in the seq file");
	setUp.setOption(mipSet, "--mipSet","Mip Set Name", true);
	setUp.processDefaultReader(true);
	setUp.finishSetUp(std::cout);

	SeqInput reader(setUp.pars_.ioOptions_);
	std::string mipNameSuffix = "_S0_Sub0_mip0";
	table ret(VecStr{"mip_family","mip_id","extension_arm","ligation_arm","extension_barcode_length","ligation_barcode_length","gene_name","mipset"});
	seqInfo seq;
	reader.openIn();
	while(reader.readNextRead(seq)){
		ret.addRow(seq.name_ + nameSuffix + mipNameSuffix,
				seq.name_ + nameSuffix + mipNameSuffix,
				seq.seq_.substr(0, extensionArmSize),
				seqUtil::reverseComplement(seq.seq_.substr(len(seq) - ligationArmSize), "DNA"),
				extensionBarSize,
				ligationBarSize,
				seq.name_ + nameSuffix,
				mipSet);
	}
	auto outOpts = TableIOOpts::genTabFileOut(setUp.pars_.ioOptions_.out_.outFilename_, true);
	outOpts.out_.overWriteFile_ = setUp.pars_.ioOptions_.out_.overWriteFile_;
	ret.outPutContents(outOpts);
	return 0;
}

int mipsterUtilsRunner::fixingMipBedFiles(const njh::progutils::CmdArgs & inputCommands){
	bfs::path bedFile = "";
	OutOptions bedOut(bfs::path("out.bed"));
	seqSetUp setUp(inputCommands);
	setUp.setOption(bedFile, "--bed", "Bed file to fix", true);
	setUp.processWritingOptions(bedOut);
	setUp.finishSetUp(std::cout);

	BioDataFileIO<Bed6RecordCore> reader(IoOptions(InOptions(bedFile), bedOut));

	reader.openIn();

	Bed6RecordCore bRecord;
	std::vector<Bed6RecordCore> bRecords;
	while(reader.readNextRecord(bRecord)){
		MetaDataInName meta(bRecord.name_);
		auto mipName = meta.getMeta("mipTar");
		bRecord.name_ = mipName;
		bRecords.emplace_back(bRecord);
	}
	reader.closeIn();
	reader.openOut();

	for(const auto & record : bRecords){
		reader.write(record, [](const Bed6RecordCore & record, std::ostream & out){
			out << record.toDelimStr() << std::endl;
		});
	}



	return 0;
}


int mipsterUtilsRunner::creatingSeqTableFromDirectory(const njh::progutils::CmdArgs & inputCommands){
	bfs::path directory = "";
	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	setUp.setOption(directory, "--directory", "Input directory", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	std::ofstream outFile;
	std::ostream out(outOpts.determineOutBuf(outFile));
	out << "target\tgenome\tseq" << std::endl;
	auto inputFiles = njh::files::gatherFiles(directory, ".fasta", false);
	for(const auto & inputFile : inputFiles){
		auto seqIn = SeqIOOptions::genFastaIn(inputFile);
		SeqInput reader(seqIn);
		auto seqs = reader.readAllReads<seqInfo>();
		auto targetName = bfs::basename(inputFile);
		for(const auto & seq : seqs){
			auto toks = tokenizeString(seq.name_, "-");
			for(const auto & tok : toks){
				out << targetName
						<<  "\t" << tok
						<< "\t" << seq.seq_ << std::endl;
			}
		}
	}



	return 0;
}
int mipsterUtilsRunner::createMipArmFromSelectedMips(const njh::progutils::CmdArgs & inputCommands){
	mipCorePars pars;

	OutOptions outOpts(bfs::path(""));
	seqSetUp setUp(inputCommands);
	pars.processDefaults(setUp);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);


	MipCollection mips(pars.mipArmsFileName, pars.allowableErrors);

	MipsSamplesNames names(pars.mipsSamplesFile);
	OutputStream out(outOpts);
	out << njh::conToStr(Mip::writeInfoLineHeader(),"\t") << "\n";
	for(const auto & mipFam : names.mips_){
		for(const auto & mip : mips.getMipsForFamily(mipFam)){
			mips.mips_.at(mip).writeInfoLine(out);
		}
	}

	return 0;
}




                    
} // namespace njhseq
