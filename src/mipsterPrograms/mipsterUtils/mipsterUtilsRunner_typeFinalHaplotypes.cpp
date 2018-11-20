/*
 * mipsterUtilsRunner_typeFinalHaplotypes.cpp
 *
 *  Created on: Oct 4, 2017
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

#include "mipsterUtilsRunner.hpp"
#include <elucidator/seqToolsUtils.h>
#include <elucidator/objects/BioDataObject.h>
#include <njhseq/objects/Gene/GeneFromGffs.hpp>
#include <elucidator/BamToolsUtils.h>

#include <njhseq/objects/Gene/GenomicAminoAcidPositionTyper.hpp>

#include <TwoBit.h>

#include <unordered_map>


namespace njhseq {




int mipsterUtilsRunner::typeFinalHaplotypes(
		const njh::progutils::CmdArgs & inputCommands) {
	mipCorePars pars;
	bfs::path genomeFnp = "";
	bfs::path gffFnp = "";
	bfs::path proteinMutantTypingFnp = "";
	bool zeroBased = false;
	//bool onlyInput = false;
	mipsterUtilsSetUp setUp(inputCommands);
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();
	//setUp.setOption(onlyInput, "--onlyInput", "Output information only on the input amino acid");
	setUp.setOption(proteinMutantTypingFnp, "--proteinMutantTypingFnp", "Protein Mutant Typing Fnp, columns should be ID=gene id in gff, AAPosition=amino acid position", true);
	setUp.setOption(zeroBased, "--zeroBased", "If the positions in the proteinMutantTypingFnp are zero based");
	setUp.setOption(genomeFnp, "--genomeFnp", "Genome to align to", true);
	setUp.setOption(gffFnp, "--gffFnp", "GFF path for gene annotation", true);
	setUp.processDirectoryOutputName("haplotypeTyping_TODAY", true);
	setUp.setOption(pars.numThreads, "--numThreads", "Number of jobs to run in parallel");
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

	njh::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);
	BioCmdsUtils bRunner(setUp.pars_.verbose_);
	bRunner.RunFaToTwoBit(genomeFnp);
	bRunner.RunBowtie2Index(genomeFnp);

	auto gprefix = bfs::path(genomeFnp).replace_extension("");
	auto twoBitFnp = gprefix.string() + ".2bit";

	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	if("" != pars.sampleMetaFnp){
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}

	TwoBit::TwoBitFile tReader(twoBitFnp);


	//initiate typer
	GenomicAminoAcidPositionTyper aaTyper(proteinMutantTypingFnp, zeroBased);


	//get all the population haplotypes
	auto popMips = mipMaster.getMipFamsWithPopClustered(pars.numThreads);
	std::vector<bfs::path> popMipsFnps;
	for(const auto & mipName : popMips){
		popMipsFnps.emplace_back(mipMaster.pathMipPopClusHaplo(mipName));
	}
	OutOptions allPopHapsOpts(njh::files::make_path(setUp.pars_.directoryName_, "allPopSeqs.fastq"));
	concatenateFiles(popMipsFnps, allPopHapsOpts);

	//align to genome
	auto alignOpts = SeqIOOptions::genFastqIn(allPopHapsOpts.outName());
	alignOpts.out_.outFilename_ = njh::files::make_path(setUp.pars_.directoryName_, "allPopSeqs.sorted.bam");
	alignOpts.out_.outExtention_= ".sorted.bam";
	auto bowtieRunOut = bRunner.bowtie2Align(alignOpts, genomeFnp, njh::pasteAsStr("-p ", pars.numThreads));


	//get overlapped gene ids
	auto gCounter = GenomicRegionCounter::countRegionsInBam(alignOpts.out_.outName());

	auto idsFromData = gCounter.getIntersectingGffIds(gffFnp);


	// get gene information
	auto geneInfoDir = njh::files::make_path(setUp.pars_.directoryName_, "geneInfos");
	njh::files::makeDir(njh::files::MkdirPar{geneInfoDir});

	OutOptions outOpts(njh::files::make_path(geneInfoDir, "gene"));
	std::set<std::string> allIds = aaTyper.getGeneIds();
	allIds.insert(idsFromData.begin(), idsFromData.end());

	std::unordered_map<std::string, VecStr> idToTranscriptName;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, allIds);
	std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> geneInfos;
	OutOptions geneIdsOpts(njh::files::make_path(geneInfoDir, "geneIds.txt"));
	OutputStream geneIdsOut(geneIdsOpts);
	uint64_t proteinMaxLen = 0;
	std::unordered_map<std::string, std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>>> genesByChrom;

	for(const auto & gene : genes){
		genesByChrom[gene.second->gene_->seqid_].emplace(gene.first, gene.second);
		geneIdsOut << gene.first << std::endl;
		for(const auto & transcript : gene.second->mRNAs_){
			idToTranscriptName[gene.second->gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
		}
		geneInfos[gene.first] = gene.second->generateGeneSeqInfo(tReader, false).begin()->second;
		readVec::getMaxLength(geneInfos[gene.first]->protein_, proteinMaxLen);
		gene.second->writeOutGeneInfo(tReader, outOpts);
	}

	//check for multiple transcripts
	/**@todo add inputing what transcript to avoid this check, cannot be handled now because the input positions can only refere to one transcript*/
	bool failed = false;
	VecStr idsWithMoreThanOneTranscript;
	for(const auto & idTrans : idToTranscriptName){
		if(idTrans.second.size() > 1){
			failed = true;
			idsWithMoreThanOneTranscript.emplace_back(idTrans.first);
		}
	}
	if (failed) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error the following ids were found to have more than 1 transcript which can't be handled in the current version"
				<< "\n";
		ss << "ids: " << njh::conToStr(idsWithMoreThanOneTranscript, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}
	//check for missing ids from the input query
	VecStr idsMissing;
	for(const auto & id : aaTyper.getGeneIds()){
		if(!njh::in(id, idToTranscriptName)){
			idsMissing.emplace_back(id);
			failed = true;
		}
	}

	if (failed) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", the following ids were not found in the gff file, " << gffFnp
				<< "\n";
		ss << "ids: " << njh::conToStr(idsMissing, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}
	//check to make sure the asked for positons are not out of range
	bool failedPosition = false;
	std::stringstream ssAminoAcidPosCheck;
	ssAminoAcidPosCheck << __PRETTY_FUNCTION__ <<  "\n" ;
	for (const auto & geneId : aaTyper.aminoPositionsForTyping_) {
		std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> gsInfos = genes[geneId.first]->generateGeneSeqInfo(tReader, false);
		auto gsInfo = gsInfos.begin()->second;
		std::map<uint32_t, char> refAminoAcids;
		for (const auto & aaPos : geneId.second) {
			if (aaPos >= gsInfo->protein_.seq_.size()) {
				ssAminoAcidPosCheck << "amino acid position "
						<< (zeroBased ? aaPos : aaPos + 1)
						<< " is out of range of the gene id " << geneId.first
						<< ", max position: "
						<< (zeroBased ?
								gsInfo->protein_.seq_.size() - 1 : gsInfo->protein_.seq_.size())
						<< '\n';
			} else {
				refAminoAcids[aaPos] = gsInfo->protein_.seq_[aaPos];
			}
		}
		aaTyper.aminoPositionsForTypingWithInfo_.emplace(geneId.first, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo(geneId.first, refAminoAcids));
		if(aaTyper.altNamesForIds_[geneId.first].size() == 1){
			aaTyper.aminoPositionsForTypingWithInfo_.at(geneId.first).altId_ = *aaTyper.altNamesForIds_[geneId.first].begin();
		}
	}
	if(failedPosition){
		throw std::runtime_error{ssAminoAcidPosCheck.str()};
	}

	aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));


	std::unordered_map<std::string, std::unordered_map<std::string, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo>> popHapsTyped;
	std::unordered_map<std::string, std::set<std::string>> regionsToGeneIds;
	//targetName, GeneID, AA Position
	std::map<std::string, std::map<std::string, std::set<uint32_t>>> targetNameToAminoAcidPositions;
	std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
	for(const auto & gCount : gCounter.counts_){
		for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
			if (gCount.second.region_.overlaps(*g.second->gene_)) {
				alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
			}
		}
	}
	{
		//typing by protein changes
		BamTools::BamReader bReader;
		bReader.Open(alignOpts.out_.outName().string());
		checkBamOpenThrow(bReader, alignOpts.out_.outName());
		auto refData = bReader.GetReferenceData();

		uint32_t numOfThreadsToUse = std::min<uint32_t>(refData.size(), pars.numThreads);
		//create bam readers
		concurrent::BamReaderPool bamPool(alignOpts.out_.outName(), numOfThreadsToUse);
		bamPool.openBamFile();
		//create alingers
		concurrent::AlignerPool alnPool(alignObj, numOfThreadsToUse);
		alnPool.initAligners();

		MultiSeqOutCache<seqInfo> proteinSeqOuts;
		for(const auto & g :genes){
			proteinSeqOuts.addReader(g.first, SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, g.first)));
		}
		std::mutex proteinSeqOutsMut;
		std::mutex transferInfoMut;
		auto chromRegions = genGenRegionsFromRefData(refData);
		njh::concurrent::LockableVec<GenomicRegion> chromRegionsVec(chromRegions);

		auto typeOnChrom = [&bamPool,&alnPool,&proteinSeqOuts,&proteinSeqOutsMut,
												&chromRegionsVec,&refData,&mipMaster,&genes,
												&geneInfos,&alnRegionToGeneIds,&transferInfoMut,
												&setUp,&twoBitFnp,&zeroBased,&aaTyper,
												&targetNameToAminoAcidPositions,&popHapsTyped,&regionsToGeneIds](){
			GenomicRegion currentChrom;
			auto curAligner = alnPool.popAligner();
			auto curBReader = bamPool.popReader();
			TwoBit::TwoBitFile tReader(twoBitFnp);
			BamTools::BamAlignment bAln;

			std::unordered_map<std::string, std::unordered_map<std::string, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo>> currentPopHapsTyped;
			//targetName, GeneID, AA Position
			std::map<std::string, std::map<std::string, std::set<uint32_t>>> currentTargetNameToAminoAcidPositions;
			std::unordered_map<std::string, std::set<std::string>> currentRegionsToGeneIds;

			while(chromRegionsVec.getVal(currentChrom)){
				setBamFileRegionThrow(*curBReader, currentChrom);
				while (curBReader->GetNextAlignment(bAln)) {
					if (bAln.IsMapped() && bAln.IsPrimaryAlignment()) {
						auto balnGenomicRegion = GenomicRegion(bAln, refData);

						if (setUp.pars_.verbose_) {
							std::cout << balnGenomicRegion.genBedRecordCore().toDelimStr() << std::endl;
						}
						if(!njh::in(balnGenomicRegion.createUidFromCoords(), alnRegionToGeneIds)){
							continue;
						}
						for (const auto & g : alnRegionToGeneIds.at(balnGenomicRegion.createUidFromCoords())) {
							bAln.Name = bAln.Name.substr(0, bAln.Name.rfind("_f"));
							auto targetName = bAln.Name.substr(0, bAln.Name.find("."));
							auto regionName = mipMaster.getGroupForMipFam(targetName);
							currentRegionsToGeneIds[regionName].emplace(g);

							const auto & currentGene = genes.at(g);
							const auto & currentGeneInfo = geneInfos.at(g);

							auto aminoTyping = aaTyper.typeAlignment(bAln, *currentGene, *currentGeneInfo, tReader, *curAligner, refData);

							for(const auto & aType : aminoTyping){
								if(' ' != aType.second){
									currentTargetNameToAminoAcidPositions[targetName][g].emplace(aType.first);
								}
							}
							if (!aminoTyping.empty() && njh::in(g, aaTyper.aminoPositionsForTyping_)) {
								currentPopHapsTyped[bAln.Name].emplace(g, GenomicAminoAcidPositionTyper::GeneAminoTyperInfo(g, aminoTyping));
								auto typeMeta = MetaDataInName::mapToMeta(aminoTyping);
								curAligner->alignObjectB_.seqBase_.name_.append(typeMeta.createMetaName());
							}
							if(setUp.pars_.debug_){
								std::lock_guard<std::mutex> lock(proteinSeqOutsMut);
								proteinSeqOuts.add(g, curAligner->alignObjectA_.seqBase_);
								proteinSeqOuts.add(g, curAligner->alignObjectB_.seqBase_);
							}
						}
					}
				}
				{
					std::lock_guard<std::mutex> lock(transferInfoMut);
					for(const auto & tar : currentTargetNameToAminoAcidPositions){
						targetNameToAminoAcidPositions.emplace(tar.first, tar.second);
					}
					for(const auto & popHapTyped : currentPopHapsTyped){
						popHapsTyped.emplace(popHapTyped.first, popHapTyped.second);
					}
					for(const auto & regionToGeneId : currentRegionsToGeneIds){
						regionsToGeneIds.emplace(regionToGeneId.first, regionToGeneId.second);
					}
				}
			}
		};

		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < numOfThreadsToUse; ++t){
			threads.emplace_back(typeOnChrom);
		}
		njh::concurrent::joinAllJoinableThreads(threads);
	}

	OutOptions targetToAminoAcidsCoveredOpt(njh::files::make_path(setUp.pars_.directoryName_, "targetToAminoAcidsCovered.tab.txt"));
	OutputStream targetToAminoAcidsCoveredOut(targetToAminoAcidsCoveredOpt);
	targetToAminoAcidsCoveredOut << "targetName\tGeneID\taltName\taminoAcidPosition\trefAminoAcid" << std::endl;

	for(const auto & tar : targetNameToAminoAcidPositions){
		for(const auto & gId : tar.second){
			for(const auto & pos : gId.second){
				targetToAminoAcidsCoveredOut << tar.first
						<< "\t" << gId.first
						<< "\t" << aaTyper.aminoPositionsForTypingWithInfo_.at(gId.first).altId_
						<< "\t" << pos
						<< "\t" << aaTyper.aminoPositionsForTypingWithInfo_.at(gId.first).aminos_.at(zeroBased ? pos : pos - 1)
						<< std::endl;
			}
		}
	}

	TableReader allInfoReader(TableIOOpts::genTabFileIn(mipMaster.pathToAllPopInfo()));
	VecStr row;
	std::unordered_map<std::string, std::unique_ptr<OutputStream>> outputs;

	VecStr cols       {            "s_Sample",  "p_geneName",   "p_targetName","h_popUID", "s_usedTotalClusterCnt","s_usedTotalBarcodeCnt","c_barcodeCnt","c_barcodeFrac"};
	VecStr renamedCols{"GeneID",   "Sample",    "GenomicRegion","TargetName",  "h_popUID", "TargetCOI",            "TotalBarcodes",        "Barcodes",    "BarcodesFraction"};
	std::vector<uint32_t> colPos;
	for(const auto & col : cols){
		colPos.emplace_back(allInfoReader.header_.getColPos(col));
	}
	if(nullptr != mipMaster.meta_){
		auto metaFields = getVectorOfMapKeys(mipMaster.meta_->groupData_);
		addOtherVec(renamedCols, metaFields);
	}



	while(allInfoReader.getNextRow(row)){
		auto outRow = getTargetsAtPositions(row, colPos);
		auto sample =     outRow[0];
		auto regionName = outRow[1];
		auto targetName = outRow[2];
		auto popName =    outRow[3];
		std::string outPutName = "others";
		bool aminoPosTyping = njh::in(regionName, regionsToGeneIds) &&
				std::any_of(regionsToGeneIds[regionName].begin(),
										regionsToGeneIds[regionName].end(), [&aaTyper](const std::string & gId){
									return njh::in(gId, aaTyper.aminoPositionsForTypingWithInfo_);
								});
		if (aminoPosTyping) {
			outPutName = regionName + "-" + targetName;
		}
		if(!njh::in(outPutName, outputs)){
			outputs.emplace(outPutName, std::make_unique<OutputStream>(OutOptions(njh::files::make_path(setUp.pars_.directoryName_, outPutName + ".tab.txt"))));
			(*outputs.at(outPutName)) << njh::conToStr(renamedCols, "\t");
			if(aminoPosTyping){
				VecStr additionalColumns;
				for(const auto & g : regionsToGeneIds.at(regionName)){
					if(njh::in(g, aaTyper.aminoPositionsForTyping_)){
						for(const auto & aaPos : aaTyper.aminoPositionsForTyping_.at(g)){
							uint32_t currentPos = zeroBased ? aaPos : aaPos + 1;
							if(njh::in(currentPos, targetNameToAminoAcidPositions[targetName][g])){
								additionalColumns.emplace_back(njh::pasteAsStr(g, "-", currentPos));
							}
						}
					}
				}
				if(!additionalColumns.empty()){
					(*outputs.at(outPutName)) << "\t"<< njh::conToStr(additionalColumns, "\t");
				}
			}
			(*outputs.at(outPutName)) << std::endl;
		}
		if(aminoPosTyping){
			(*outputs.at(outPutName)) << njh::conToStr(regionsToGeneIds.at(regionName), ",") << "\t"<< njh::conToStr(outRow, "\t");
			VecStr additionalColumns;
			for(const auto & g : regionsToGeneIds.at(regionName)){
				if(njh::in(g, aaTyper.aminoPositionsForTyping_)){
					for(const auto & aaPos : aaTyper.aminoPositionsForTyping_.at(g)){
						uint32_t currentPos = zeroBased ? aaPos : aaPos + 1;
						if(njh::in(currentPos, targetNameToAminoAcidPositions[targetName][g])){
							if(njh::in(popName, popHapsTyped)){
								if(njh::in(g, popHapsTyped[popName])){
									additionalColumns.emplace_back(std::string(1, popHapsTyped[popName].at(g).aminos_[currentPos]));
								}else{
									additionalColumns.emplace_back("");
								}
							}else{
								additionalColumns.emplace_back("");
							}
						}
					}
				}
			}
			if(nullptr != mipMaster.meta_){
				for(const auto & group : mipMaster.meta_->groupData_){
					(*outputs.at(outPutName)) << "\t" << group.second->getGroupForSample(sample);
				}
			}
			if(!additionalColumns.empty()){
				(*outputs.at(outPutName)) << "\t"<< njh::conToStr(additionalColumns, "\t");
			}
			(*outputs.at(outPutName)) << std::endl;
		} else {
			std::string geneId = "";
			if(njh::in(regionName, regionsToGeneIds)){
				geneId = *regionsToGeneIds.at(regionName).begin();
			}
			(*outputs.at(outPutName)) << geneId
					<< "\t" << njh::conToStr(outRow, "\t");
			if(nullptr != mipMaster.meta_){
				for(const auto & group : mipMaster.meta_->groupData_){
					(*outputs.at(outPutName))
							<< "\t" << group.second->getGroupForSample(sample);
				}
			}
			(*outputs.at(outPutName)) << std::endl;
		}
	}
	return 0;
}


} //namespace njhseq

