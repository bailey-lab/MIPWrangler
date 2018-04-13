/*
 * mipsterUtilsRunner_typeFinalHaplotypes.cpp
 *
 *  Created on: Oct 4, 2017
 *      Author: nick
 */


#include "mipsterUtilsRunner.hpp"
#include <elucidator/seqToolsUtils.h>
#include <elucidator/objects/BioDataObject.h>
#include <elucidator/objects/Gene/GeneFromGffs.hpp>
#include <elucidator/BamToolsUtils.h>

#include <TwoBit.h>

#include <unordered_map>


namespace bibseq {


int mipsterUtilsRunner::typeFinalHaplotypes(
		const bib::progutils::CmdArgs & inputCommands) {
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
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}

	bib::files::checkExistenceThrow(genomeFnp, __PRETTY_FUNCTION__);
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

	table proteinMutantTypingTab(proteinMutantTypingFnp, "\t", true);
	proteinMutantTypingTab.checkForColumnsThrow(VecStr { "ID", "AAPosition" }, __PRETTY_FUNCTION__);

	struct GeneAminoTyperInfo {
		GeneAminoTyperInfo(const std::string & geneId) :
				geneId_(geneId), altId_(geneId) {
		}
		GeneAminoTyperInfo(const std::string & geneId,
				const std::map<uint32_t, char> & aminos) :
				geneId_(geneId), altId_(geneId), aminos_(aminos) {
		}
		std::string geneId_;
		std::string altId_;
		std::map<uint32_t, char> aminos_;
	};

	std::unordered_map<std::string, std::vector<uint32_t>> aminoPositionsForTyping;
	std::unordered_map<std::string, std::set<std::string>> altNamesForIds;
	std::unordered_map<std::string, GeneAminoTyperInfo> aminoPositionsForTypingWithInfo;

	for (const auto & row : proteinMutantTypingTab.content_) {
		uint32_t aaPos = zeroBased ?
				bib::StrToNumConverter::stoToNum<uint32_t>(
						row[proteinMutantTypingTab.getColPos("AAPosition")]) :
				bib::StrToNumConverter::stoToNum<uint32_t>(
						row[proteinMutantTypingTab.getColPos("AAPosition")]) - 1;
		if(bib::in(aaPos, aminoPositionsForTyping[row[proteinMutantTypingTab.getColPos("ID")]])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in reading in "
					<< proteinMutantTypingFnp << " gene id " << row[proteinMutantTypingTab.getColPos("ID")]
					<< " already has position: " << aaPos  << "\n";
			throw std::runtime_error { ss.str() };
		}
		aminoPositionsForTyping[row[proteinMutantTypingTab.getColPos("ID")]].emplace_back(
				aaPos);
		if(proteinMutantTypingTab.containsColumn("Gene")){
			altNamesForIds[row[proteinMutantTypingTab.getColPos("ID")]].emplace(row[proteinMutantTypingTab.getColPos("Gene")]);
		}
	}
	for (const auto & geneId : aminoPositionsForTyping) {
		if (altNamesForIds[geneId.first].size() > 1) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error in reading in "
					<< proteinMutantTypingFnp << " gene id " << geneId.first
					<< " had more than one alternative ID in Gene column, found: "
					<< bib::conToStr(altNamesForIds[geneId.first], ", ") << "\n";
			throw std::runtime_error { ss.str() };
		}
	}

	for(auto & positions : aminoPositionsForTyping){
		bib::sort(positions.second);
	}

	auto idsVec = proteinMutantTypingTab.getColumnLevels("ID");
	std::set<std::string> ids(idsVec.begin(), idsVec.end());



	auto popMips = mipMaster.getMipFamsWithPopClustered(pars.numThreads);
	std::vector<bfs::path> popMipsFnps;
	for(const auto & mipName : popMips){
		popMipsFnps.emplace_back(mipMaster.pathMipPopClusHaplo(mipName));
	}
	OutOptions allPopHapsOpts(bib::files::make_path(setUp.pars_.directoryName_, "allPopSeqs.fastq"));
	concatenateFiles(popMipsFnps, allPopHapsOpts);
	auto alignOpts = SeqIOOptions::genFastqIn(allPopHapsOpts.outName());
	alignOpts.out_.outFilename_ = bib::files::make_path(setUp.pars_.directoryName_, "allPopSeqs.sorted.bam");
	alignOpts.out_.outExtention_= ".sorted.bam";
	auto bowtieRunOut = bRunner.bowtie2Align(alignOpts, genomeFnp, bib::pasteAsStr("-p ", pars.numThreads));
	GenomicRegionCounter gCounter;


	{
		BamTools::BamReader bReader;
		bReader.Open(alignOpts.out_.outName().string());
		checkBamOpenThrow(bReader, alignOpts.out_.outName());

		BamTools::BamAlignment bAln;
		auto refData = bReader.GetReferenceData();
		while(bReader.GetNextAlignment(bAln)){
			if(bAln.IsMapped()){
				gCounter.increaseCount(GenomicRegion(bAln, refData), 1);
			}
		}
		BioDataFileIO<GFFCore> reader{IoOptions(InOptions(gffFnp))};
		reader.openIn();
		uint32_t count = 0;
		std::string line = "";
		std::shared_ptr<GFFCore> gRecord = reader.readNextRecord();
		while (nullptr != gRecord) {
			if("gene" == gRecord->type_){
				for(const auto & gCount : gCounter.counts_){
					if(gCount.second.region_.overlaps(*gRecord)){
						ids.emplace(gRecord->getIDAttr());
						break;
					}
				}
			}
			bool end = false;
			while ('#' == reader.inFile_->peek()) {
				if (bib::files::nextLineBeginsWith(*reader.inFile_, "##FASTA")) {
					end = true;
					break;
				}
				bib::files::crossPlatGetline(*reader.inFile_, line);
			}
			if (end) {
				break;
			}
			gRecord = reader.readNextRecord();
			++count;
		}
	}
	auto geneInfoDir = bib::files::make_path(setUp.pars_.directoryName_, "geneInfos");
	bib::files::makeDir(bib::files::MkdirPar{geneInfoDir});
	OutOptions outOpts(bib::files::make_path(geneInfoDir, "gene"));

	std::unordered_map<std::string, VecStr> idToTranscriptName;
	std::unordered_map<std::string, std::shared_ptr<GeneFromGffs>> genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, ids);
	std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> geneInfos;
	OutOptions geneIdsOpts(bib::files::make_path(geneInfoDir, "geneIds.txt"));
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
		ss << "ids: " << bib::conToStr(idsWithMoreThanOneTranscript, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}

	VecStr idsMissing;
	for(const auto & id : ids){
		if(!bib::in(id, idToTranscriptName)){
			idsMissing.emplace_back(id);
			failed = true;
		}
	}

	if (failed) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", the following ids were not found in the gff file, " << gffFnp
				<< "\n";
		ss << "ids: " << bib::conToStr(idsMissing, ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}
	bool failedPosition = false;
	std::stringstream ssAminoAcidPosCheck;
	ssAminoAcidPosCheck << __PRETTY_FUNCTION__ <<  "\n" ;
	for (const auto & geneId : aminoPositionsForTyping) {
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
		aminoPositionsForTypingWithInfo.emplace(geneId.first, GeneAminoTyperInfo(geneId.first, refAminoAcids));
		if(altNamesForIds[geneId.first].size() == 1){
			aminoPositionsForTypingWithInfo.at(geneId.first).altId_ = *altNamesForIds[geneId.first].begin();
		}
	}

	if(failedPosition){
		throw std::runtime_error{ssAminoAcidPosCheck.str()};
	}

	aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));


	std::unordered_map<std::string, std::unordered_map<std::string, GeneAminoTyperInfo>> popHapsTyped;
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
			proteinSeqOuts.addReader(g.first, SeqIOOptions::genFastaOut(bib::files::make_path(setUp.pars_.directoryName_, g.first)));
		}
		std::mutex proteinSeqOutsMut;
		std::mutex transferInfoMut;
		auto chromRegions = genGenRegionsFromRefData(refData);
		bib::concurrent::LockableVec<GenomicRegion> chromRegionsVec(chromRegions);

		auto typeOnChrom = [&bamPool,&alnPool,&proteinSeqOuts,&proteinSeqOutsMut,
												&chromRegionsVec,&refData,&mipMaster,&genes,
												&geneInfos,&alnRegionToGeneIds,&transferInfoMut,
												&setUp,&twoBitFnp,&zeroBased,&aminoPositionsForTyping,
												&targetNameToAminoAcidPositions,&popHapsTyped,&regionsToGeneIds](){
			GenomicRegion currentChrom;
			auto curAligner = alnPool.popAligner();
			auto curBReader = bamPool.popReader();
			TwoBit::TwoBitFile tReader(twoBitFnp);
			BamTools::BamAlignment bAln;

			std::unordered_map<std::string, std::unordered_map<std::string, GeneAminoTyperInfo>> currentPopHapsTyped;
			//targetName, GeneID, AA Position
			std::map<std::string, std::map<std::string, std::set<uint32_t>>> currentTargetNameToAminoAcidPositions;
			std::unordered_map<std::string, std::set<std::string>> currentRegionsToGeneIds;

			while(chromRegionsVec.getVal(currentChrom)){
				setBamFileRegionThrow(*curBReader, currentChrom);
				while (curBReader->GetNextAlignment(bAln)) {
					if (bAln.IsMapped()) {
						auto results = std::make_shared<AlignmentResults>(bAln, refData, true);
						if (setUp.pars_.verbose_) {
							std::cout << results->gRegion_.genBedRecordCore().toDelimStr() << std::endl;
						}
						if(!bib::in(results->gRegion_.createUidFromCoords(), alnRegionToGeneIds)){
							continue;
						}
						results->setRefSeq(tReader);
						results->setComparison(true);

						for (const auto & g : alnRegionToGeneIds.at(results->gRegion_.createUidFromCoords())) {
							const auto & currentGene = genes.at(g);
							bAln.Name = bAln.Name.substr(0, bAln.Name.rfind("_f"));
							auto targetName = bAln.Name.substr(0, bAln.Name.find("."));
							auto regionName = mipMaster.getGroupForMipFam(targetName);
							currentRegionsToGeneIds[regionName].emplace(g);

							bool endsAtStopCodon = false;
							uint32_t transStart = 0;
							std::unordered_map<size_t, alnInfoLocal> balnAlnInfo = bamAlnToAlnInfoLocal(bAln);
							auto genePosInfoByGDna = geneInfos.at(g)->getInfosByGDNAPos();
							const auto & transcript = currentGene->mRNAs_.front();
							seqInfo balnSeq(bAln.Name);
							std::vector<GFFCore> cDNAIntersectedWith;
							for (const auto & cDna : currentGene->CDS_.at(transcript->getIDAttr())) {
								if (results->gRegion_.overlaps(*cDna)) {
									cDNAIntersectedWith.emplace_back(*cDna);
								}
							}
							if (cDNAIntersectedWith.size() == 1
									&& results->gRegion_.start_
											>= cDNAIntersectedWith.front().start_ - 1
									&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
								balnSeq = *(results->alnSeq_);
								if (currentGene->gene_->isReverseStrand()) {
									if (genePosInfoByGDna.at(results->gRegion_.start_).cDNAPos_
											== geneInfos.at(g)->cDna_.seq_.size() - 1) {
										endsAtStopCodon = true;
									}
									uint32_t gPos = results->gRegion_.end_ - 1;
									auto codon = genePosInfoByGDna.at(gPos).codonPos_;
									while (0 != codon) {
										--gPos;
										codon = genePosInfoByGDna.at(gPos).codonPos_;
										++transStart;
									}
								} else {
									if (genePosInfoByGDna.at(results->gRegion_.end_ - 1).cDNAPos_
											== geneInfos.at(g)->cDna_.seq_.size() - 1) {
										endsAtStopCodon = true;
									}
									uint32_t gPos = results->gRegion_.start_;
									uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
									while (0 != codon) {
										++gPos;
										codon = genePosInfoByGDna.at(gPos).codonPos_;
										++transStart;
									}
								}
							} else {
								bib::sort(cDNAIntersectedWith,
										[](const GenomicRegion & reg1, const GenomicRegion & reg2) {
											if(reg1.start_ < reg2.start_) {
												return true;
											}
											return false;
										});

								if (currentGene->gene_->isReverseStrand()) {
									auto cDnaStop = cDNAIntersectedWith.back().end_;
									uint32_t gPos = std::min(cDnaStop, results->gRegion_.end_) - 1;
									auto codon = genePosInfoByGDna.at(gPos).codonPos_;
									while (0 != codon) {
										--gPos;
										codon = genePosInfoByGDna.at(gPos).codonPos_;
										++transStart;
									}
								} else {
									auto cDnaStart = cDNAIntersectedWith.front().start_ - 1;
									uint32_t gPos = std::max(cDnaStart, results->gRegion_.start_);
									uint32_t codon = genePosInfoByGDna.at(gPos).codonPos_;
									while (0 != codon) {
										++gPos;
										codon = genePosInfoByGDna.at(gPos).codonPos_;
										++transStart;
									}
								}
								std::vector<uint32_t> starts;
								std::vector<uint32_t> ends;
								for (const auto & cDna : cDNAIntersectedWith) {
									auto cDnaStart = cDna.start_ - 1;
									auto detStart = std::max(cDnaStart, results->gRegion_.start_);
									auto detStop = std::min(cDna.end_, results->gRegion_.end_);
									ends.emplace_back(detStop);
									starts.emplace_back(detStart);
									detStart -= results->gRegion_.start_;
									detStop -= results->gRegion_.start_;
									auto alnStart = getAlnPosForRealPos(results->refSeqAligned_->seq_,
											detStart);
									auto alnStop = getAlnPosForRealPos(results->refSeqAligned_->seq_,
											detStop - 1);
									balnSeq.append(
											results->alnSeqAligned_->getSubRead(alnStart,
													alnStop - alnStart + 1));
								}
								uint32_t cDnaStart = *std::min_element(starts.begin(),
										starts.end());
								uint32_t cDnaStop = *std::max_element(ends.begin(), ends.end());
								if (currentGene->gene_->isReverseStrand()) {
									if (genePosInfoByGDna.at(cDnaStart).cDNAPos_
											== geneInfos.at(g)->cDna_.seq_.size() - 1) {
										endsAtStopCodon = true;
									}
								} else {
									if (genePosInfoByGDna.at(cDnaStop - 1).cDNAPos_
											== geneInfos.at(g)->cDna_.seq_.size() - 1) {
										endsAtStopCodon = true;
									}
								}
								balnSeq.removeGaps();
							}
							if (currentGene->gene_->isReverseStrand()) {
								balnSeq.reverseComplementRead(false, true);
							}
							auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
							curAligner->alignCacheGlobal(geneInfos.at(g)->protein_, balnSeqTrans);
							curAligner->profilePrimerAlignment(geneInfos.at(g)->protein_,
									balnSeqTrans);
							std::map<uint32_t, char> aminoTyping;
							if (bib::in(g, aminoPositionsForTyping)) {
								uint32_t proteinAlnStart =
										curAligner->alignObjectB_.seqBase_.seq_.find_first_not_of('-');
								uint32_t proteinAlnStop =
										curAligner->alignObjectB_.seqBase_.seq_.find_last_not_of('-');
								uint32_t proteinStart = getRealPosForAlnPos(
										curAligner->alignObjectA_.seqBase_.seq_, proteinAlnStart);
								uint32_t proteinStop = getRealPosForAlnPos(
										curAligner->alignObjectA_.seqBase_.seq_, proteinAlnStop);

								for (const auto & pos : aminoPositionsForTyping.at(g)) {
									if (pos < proteinStart || pos > proteinStop) {
										if (!zeroBased) {
											aminoTyping[pos + 1] = ' ';
										} else {
											aminoTyping[pos] = ' ';
										}
									} else {
										auto posAln = getAlnPosForRealPos(
												curAligner->alignObjectA_.seqBase_.seq_, pos);
										if (!zeroBased) {
											currentTargetNameToAminoAcidPositions[targetName][g].emplace(pos + 1);
											aminoTyping[pos + 1] =
													curAligner->alignObjectB_.seqBase_.seq_[posAln];
										} else {
											currentTargetNameToAminoAcidPositions[targetName][g].emplace(pos);
											aminoTyping[pos] =
													curAligner->alignObjectB_.seqBase_.seq_[posAln];
										}
									}
								}
							}
							if (!aminoTyping.empty() && bib::in(g, aminoPositionsForTyping)) {
								currentPopHapsTyped[bAln.Name].emplace(g, GeneAminoTyperInfo(g, aminoTyping));
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
		bib::concurrent::joinAllJoinableThreads(threads);
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	OutOptions targetToAminoAcidsCoveredOpt(bib::files::make_path(setUp.pars_.directoryName_, "targetToAminoAcidsCovered.tab.txt"));
	OutputStream targetToAminoAcidsCoveredOut(targetToAminoAcidsCoveredOpt);
	targetToAminoAcidsCoveredOut << "targetName\tGeneID\taltName\taminoAcidPosition\trefAminoAcid" << std::endl;

	for(const auto & tar : targetNameToAminoAcidPositions){
		for(const auto & gId : tar.second){
			for(const auto & pos : gId.second){
				targetToAminoAcidsCoveredOut << tar.first
						<< "\t" << gId.first
						<< "\t" << aminoPositionsForTypingWithInfo.at(gId.first).altId_
						<< "\t" << pos
						<< "\t" << aminoPositionsForTypingWithInfo.at(gId.first).aminos_.at(zeroBased ? pos : pos - 1)
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
		bool aminoPosTyping = bib::in(regionName, regionsToGeneIds) &&
				std::any_of(regionsToGeneIds[regionName].begin(),
										regionsToGeneIds[regionName].end(), [&aminoPositionsForTypingWithInfo](const std::string & gId){
									return bib::in(gId, aminoPositionsForTypingWithInfo);
								});
		if (aminoPosTyping) {
			outPutName = regionName + "-" + targetName;
		}
		if(!bib::in(outPutName, outputs)){
			outputs.emplace(outPutName, std::make_unique<OutputStream>(OutOptions(bib::files::make_path(setUp.pars_.directoryName_, outPutName + ".tab.txt"))));
			(*outputs.at(outPutName)) << bib::conToStr(renamedCols, "\t");
			if(aminoPosTyping){
				VecStr additionalColumns;
				for(const auto & g : regionsToGeneIds.at(regionName)){
					if(bib::in(g, aminoPositionsForTyping)){
						for(const auto & aaPos : aminoPositionsForTyping.at(g)){
							uint32_t currentPos = zeroBased ? aaPos : aaPos + 1;
							if(bib::in(currentPos, targetNameToAminoAcidPositions[targetName][g])){
								additionalColumns.emplace_back(bib::pasteAsStr(g, "-", currentPos));
							}
						}
					}
				}
				if(!additionalColumns.empty()){
					(*outputs.at(outPutName)) << "\t"<< bib::conToStr(additionalColumns, "\t");
				}
			}
			(*outputs.at(outPutName)) << std::endl;
		}
		if(aminoPosTyping){
			(*outputs.at(outPutName)) << bib::conToStr(regionsToGeneIds.at(regionName), ",") << "\t"<< bib::conToStr(outRow, "\t");
			VecStr additionalColumns;
			for(const auto & g : regionsToGeneIds.at(regionName)){
				if(bib::in(g, aminoPositionsForTyping)){
					for(const auto & aaPos : aminoPositionsForTyping.at(g)){
						uint32_t currentPos = zeroBased ? aaPos : aaPos + 1;
						if(bib::in(currentPos, targetNameToAminoAcidPositions[targetName][g])){
							if(bib::in(popName, popHapsTyped)){
								if(bib::in(g, popHapsTyped[popName])){
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
				(*outputs.at(outPutName)) << "\t"<< bib::conToStr(additionalColumns, "\t");
			}
			(*outputs.at(outPutName)) << std::endl;
		} else {
			std::string geneId = "";
			if(bib::in(regionName, regionsToGeneIds)){
				geneId = *regionsToGeneIds.at(regionName).begin();
			}
			(*outputs.at(outPutName)) << geneId
					<< "\t" << bib::conToStr(outRow, "\t");
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


} //namespace bibseq

