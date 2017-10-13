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



template<typename MAP>
MetaDataInName mapToMeta(const MAP & m){
	MetaDataInName ret;
	for(const auto & p : m){
		ret.addMeta(estd::to_string(p.first), p.second);
	}
	return ret;
}


class TableReader {
public:

	TableReader(const TableIOOpts & tabOpts): tabOpts_(tabOpts){
		//inital header reader
		bib::files::checkExistenceThrow(tabOpts_.in_.inFilename_);
		std::string currentLine = bib::files::getFirstLine(
				tabOpts_.in_.inFilename_);
		auto toks = tokenizeString(currentLine, tabOpts_.inDelim_, true);
		VecStr columnNames;
		if (!tabOpts_.hasHeader_) {
			for (const auto i : iter::range(toks.size())) {
				columnNames.emplace_back("col." + leftPadNumStr(i, toks.size()));
			}
		} else {
			columnNames = toks;
		}
		header_ = table(columnNames);
		in_ = std::make_unique<InputStream>(tabOpts_.in_);
		if(tabOpts_.hasHeader_){
			bib::files::crossPlatGetline(*in_, currentLine);
		}
	}

	const TableIOOpts tabOpts_;
	table header_;

	std::unique_ptr<InputStream> in_;

	bool getNextRow(VecStr & row){
		std::string currentLine = "";
		row.clear();
		if(bib::files::crossPlatGetline(*in_, currentLine)){
			row = tokenizeString(currentLine, tabOpts_.inDelim_, true);
			if(row.size() != header_.nCol()){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error the row has a different number of columns than the first line" << "\n";
				ss << "rowSize: " << row.size() << ", firstLineSize: " << header_.nCol() << "\n";
				ss << "row: " << currentLine << "\n";
				throw std::runtime_error{ss.str()};
			}
			return true;
		}
		return false;
	}

	VecStr extractCols(const VecStr & row, const VecStr & cols) const{
		VecStr ret;
		ret.reserve(cols.size());
		for(const auto & col : cols){
			ret.emplace_back(row[header_.getColPos(col)]);
		}
		return ret;
	}

};


int mipsterUtilsRunner::typeFinalHaplotypes(
		const bib::progutils::CmdArgs & inputCommands) {
	mipCorePars pars;
	bfs::path genomeFnp = "";
	bfs::path gffFnp = "";
	bfs::path proteinMutantTypingFnp = "";
	bool zeroBased = false;
	mipsterUtilsSetUp setUp(inputCommands);
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();

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
	proteinMutantTypingTab.checkForColumnsThrow(VecStr { "ID", "AAPosition" },
			__PRETTY_FUNCTION__);

	std::unordered_map<std::string, std::vector<uint32_t>> aminoPositionsForTyping;
	for (const auto & row : proteinMutantTypingTab.content_) {
		aminoPositionsForTyping[row[proteinMutantTypingTab.getColPos("ID")]].emplace_back(
				zeroBased ?
						bib::StrToNumConverter::stoToNum<uint32_t>(
								row[proteinMutantTypingTab.getColPos("AAPosition")]) :
						bib::StrToNumConverter::stoToNum<uint32_t>(
								row[proteinMutantTypingTab.getColPos("AAPosition")]) - 1);
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
	std::unordered_map<std::string, GeneFromGffs> genes = GeneFromGffs::getGenesFromGffForIds(gffFnp, ids);
	std::unordered_map<std::string, std::shared_ptr<GeneSeqInfo>> geneInfos;
	OutOptions geneIdsOpts(bib::files::make_path(geneInfoDir, "geneIds.txt"));
	OutputStream geneIdsOut(geneIdsOpts);
	uint64_t proteinMaxLen = 0;
	std::unordered_map<std::string, std::unordered_map<std::string, GeneFromGffs>> genesByChrom;

	for(const auto & gene : genes){
		genesByChrom[gene.second.gene_->seqid_].emplace(gene.first, gene.second);
		geneIdsOut << gene.first << std::endl;
		for(const auto & transcript : gene.second.mRNAs_){
			idToTranscriptName[gene.second.gene_->getIDAttr()].emplace_back(transcript->getIDAttr());
		}
		geneInfos[gene.first] = gene.second.generateGeneSeqInfo(tReader, false);
		readVec::getMaxLength(geneInfos[gene.first]->protein_, proteinMaxLen);
		gene.second.writeOutGeneInfo(tReader, outOpts);
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

	aligner alignObj(proteinMaxLen, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(1,-1));
	struct GeneAminoTyperInfo{
		GeneAminoTyperInfo(const std::string & geneId):geneId_(geneId){

		}
		GeneAminoTyperInfo(const std::string & geneId,
				const std::map<uint32_t, char> & aminos):geneId_(geneId),
						aminos_(aminos){

		}
		std::string geneId_;
		std::map<uint32_t, char> aminos_;
	};
	std::unordered_map<std::string, std::unordered_map<std::string,GeneAminoTyperInfo>> popHapsTyped;
	std::unordered_map<std::string, std::set<std::string>> regionsToGeneIds;
	std::unordered_map<std::string, std::vector<std::string>> alnRegionToGeneIds;
	for(const auto & gCount : gCounter.counts_){
		for (const auto & g : genesByChrom[gCount.second.region_.chrom_]) {
			if (gCount.second.region_.overlaps(*g.second.gene_)) {
				alnRegionToGeneIds[gCount.second.region_.createUidFromCoords()].emplace_back(g.first);
			}
		}
	}
	{
		//typing by protein changes
		BamTools::BamReader bReader;
		bReader.Open(alignOpts.out_.outName().string());
		checkBamOpenThrow(bReader, alignOpts.out_.outName());

		BamTools::BamAlignment bAln;
		auto refData = bReader.GetReferenceData();
		OutOptions proteinOpts(bfs::path("protein_temp.fasta"));
		proteinOpts.overWriteFile_ = true;
		OutputStream proteinOut(proteinOpts);
		OutOptions cdnaOpts(bfs::path("cDna_temp.fasta"));
		cdnaOpts.overWriteFile_ = true;
		OutputStream cdnaOut(cdnaOpts);
		MultiSeqOutCache<seqInfo> proteinSeqOuts;
		for(const auto & g :genes){
			proteinSeqOuts.addReader(g.first, SeqIOOptions::genFastaOut(bib::files::make_path(setUp.pars_.directoryName_, g.first)));
		}
		while (bReader.GetNextAlignment(bAln)) {
			if (bAln.IsMapped()) {
				auto results = std::make_shared<AlignmentResults>(bAln, refData, true);

				if (setUp.pars_.verbose_) {
					std::cout << results->gRegion_.genBedRecordCore().toDelimStr()
							<< std::endl;
				}
				results->setRefSeq(tReader);
				results->setComparison(true);

				for (const auto & g : alnRegionToGeneIds[results->gRegion_.createUidFromCoords()]) {
					const auto & currentGene = genes.at(g);
					bAln.Name = bAln.Name.substr(0, bAln.Name.rfind("_f"));
					auto targetName = bAln.Name.substr(0, bAln.Name.find("."));
					auto regionName = mipMaster.getGroupForMipFam(targetName);
					regionsToGeneIds[regionName].emplace(g);
					bool endsAtStopCodon = false;
					uint32_t transStart = 0;
					std::unordered_map<size_t, alnInfoLocal> balnAlnInfo =
							bamAlnToAlnInfoLocal(bAln);
					auto genePosInfoByGDna = geneInfos[g]->getInfosByGDNAPos();
					const auto & transcript = currentGene.mRNAs_.front();
					seqInfo balnSeq(bAln.Name);
					std::vector<GFFCore> cDNAIntersectedWith;
					for (const auto & cDna : currentGene.CDS_.at(transcript->getIDAttr())) {
						if (results->gRegion_.overlaps(*cDna)) {
							cDNAIntersectedWith.emplace_back(*cDna);
						}
					}
					if (cDNAIntersectedWith.size() == 1
							&& results->gRegion_.start_
									>= cDNAIntersectedWith.front().start_ - 1
							&& results->gRegion_.end_ <= cDNAIntersectedWith.front().end_) {
						balnSeq = *(results->alnSeq_);
						if (currentGene.gene_->isReverseStrand()) {
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

						if (currentGene.gene_->isReverseStrand()) {
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
						if (currentGene.gene_->isReverseStrand()) {
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
					if (currentGene.gene_->isReverseStrand()) {
						balnSeq.reverseComplementRead(false, true);
					}
					auto balnSeqTrans = balnSeq.translateRet(false, false, transStart);
					alignObj.alignCacheGlobal(geneInfos.at(g)->protein_, balnSeqTrans);
					alignObj.profilePrimerAlignment(geneInfos.at(g)->protein_,
							balnSeqTrans);
					std::map<uint32_t, char> aminoTyping;
					if (bib::in(g, aminoPositionsForTyping)) {
						uint32_t proteinAlnStart =
								alignObj.alignObjectB_.seqBase_.seq_.find_first_not_of('-');
						uint32_t proteinAlnStop =
								alignObj.alignObjectB_.seqBase_.seq_.find_last_not_of('-');
						uint32_t proteinStart = getRealPosForAlnPos(
								alignObj.alignObjectA_.seqBase_.seq_, proteinAlnStart);
						uint32_t proteinStop = getRealPosForAlnPos(
								alignObj.alignObjectA_.seqBase_.seq_, proteinAlnStop);
						for (const auto & pos : aminoPositionsForTyping[g]) {
							if (pos < proteinStart | pos > proteinStop) {
								if (!zeroBased) {
									aminoTyping[pos + 1] = ' ';
								} else {
									aminoTyping[pos] = ' ';
								}
							} else {
								auto posAln = getAlnPosForRealPos(
										alignObj.alignObjectA_.seqBase_.seq_, pos);
								if (!zeroBased) {
									aminoTyping[pos + 1] =
											alignObj.alignObjectB_.seqBase_.seq_[posAln];
								} else {
									aminoTyping[pos] =
											alignObj.alignObjectB_.seqBase_.seq_[posAln];
								}
							}
						}
					}

					if (!aminoTyping.empty() && bib::in(g, aminoPositionsForTyping)) {
						popHapsTyped[bAln.Name].emplace(g,
								GeneAminoTyperInfo(g, aminoTyping));
						auto typeMeta = mapToMeta(aminoTyping);
						alignObj.alignObjectB_.seqBase_.name_.append(
								typeMeta.createMetaName());
					}
					proteinSeqOuts.add(g, alignObj.alignObjectA_.seqBase_);
					proteinSeqOuts.add(g, alignObj.alignObjectB_.seqBase_);
				}
			}
		}
	}

	TableReader allInfoReader(TableIOOpts::genTabFileIn(mipMaster.pathToAllPopInfo()));
	VecStr row;
	std::unordered_map<std::string, std::unique_ptr<OutputStream>> outputs;
//	;
	VecStr cols{"s_Sample", "p_geneName","p_targetName","h_popUID", "s_usedTotalClusterCnt","s_usedTotalBarcodeCnt","c_barcodeCnt","c_barcodeFrac"};
	VecStr renamedCols{"GeneID", "Sample", "Region","TargetName","h_popUID", "COI","TotalBarcodes","Barcodes","BarcodesFraction"};
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
		auto regionName = outRow[1];
		auto popName = outRow[3];
		auto sample = outRow[0];
		if(!bib::in(regionName, outputs)){
			outputs.emplace(regionName, std::make_unique<OutputStream>(OutOptions(bib::files::make_path(setUp.pars_.directoryName_, regionName + ".tab.txt"))));
			(*outputs.at(regionName)) << bib::conToStr(renamedCols, "\t");
			if(bib::in(regionName, regionsToGeneIds)){
				VecStr additionalColumns;
				for(const auto & g : regionsToGeneIds.at(regionName)){
					if(bib::in(g, aminoPositionsForTyping)){
						for(const auto & aaPos : aminoPositionsForTyping.at(g)){
							additionalColumns.emplace_back(bib::pasteAsStr(g, "-", aaPos));
						}
					}
				}
				if(!additionalColumns.empty()){
					(*outputs.at(regionName)) << "\t"<< bib::conToStr(additionalColumns, "\t");
				}
				(*outputs.at(regionName)) << std::endl;
			}else{
				(*outputs.at(regionName)) << "\t" << bib::conToStr(renamedCols, "\t") << std::endl;
			}
		}
		if(bib::in(regionName, regionsToGeneIds)){
			(*outputs.at(regionName)) << bib::conToStr(regionsToGeneIds.at(regionName), ",") << "\t"<< bib::conToStr(outRow, "\t");
			VecStr additionalColumns;

			for(const auto & g : regionsToGeneIds.at(regionName)){
				if(bib::in(g, aminoPositionsForTyping)){
					for(const auto & aaPos : aminoPositionsForTyping.at(g)){
						if(bib::in(popName, popHapsTyped)){
							if(bib::in(g, popHapsTyped[popName])){
								additionalColumns.emplace_back(std::string(1, popHapsTyped[popName].at(g).aminos_[zeroBased ? aaPos : aaPos + 1]));
							}else{
								additionalColumns.emplace_back("");
							}
						}else{
							additionalColumns.emplace_back("");
						}
					}
				}
			}
			if(nullptr != mipMaster.meta_){
				for(const auto & group : mipMaster.meta_->groupData_){
					(*outputs.at(regionName)) << "\t" << group.second->getGroupForSample(sample);
				}
			}
			if(!additionalColumns.empty()){
				(*outputs.at(regionName)) << "\t"<< bib::conToStr(additionalColumns, "\t");
			}

			(*outputs.at(regionName)) << std::endl;
		}else{
			(*outputs.at(regionName)) << "\t" << bib::conToStr(outRow, "\t");
			if(nullptr != mipMaster.meta_){
				for(const auto & group : mipMaster.meta_->groupData_){
					(*outputs.at(regionName)) << "\t" << group.second->getGroupForSample(sample);
				}
			}
			(*outputs.at(regionName)) << std::endl;
		}
	}

	return 0;
}


} //namespace bibseq

