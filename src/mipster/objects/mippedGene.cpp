/*
 * mippedGene.cpp
 *
 *  Created on: Apr 15, 2015
 *      Author: nickhathaway
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
#include "mippedGene.hpp"

namespace njhseq {



mippedGene::mippedGene(const seqInfo & seqBase,
		const std::unordered_map<std::string, mipTargetReads> & mReads) :
		genomicDna_(seqBase), processedReads_(mReads) {

}
std::unordered_map<std::string, mipTargetReads> processMipReads(
		const std::vector<readObject> & reads, const readObject & genomicDna,
		aligner & alignerObj, bool processed, bool finalAnalysisNameFormat) {
	std::unordered_map<std::string, mipTargetReads> mReads;
	for (const auto & read : reads) {
		//read.seqBase_.outPutSeq(std::cout);
		auto mRead = std::make_shared<mipRead>(read.seqBase_, processed, finalAnalysisNameFormat);
		auto search = mReads.find(mRead->mipTargetName_);
		alignerObj.alignCacheGlobal(genomicDna, read);
		auto firstAlign = alignerObj.alignObjectB_;
		int32_t bestScore = alignerObj.parts_.score_;
		auto readComp = read;
		readComp.seqBase_.reverseComplementRead(true);
		alignerObj.alignCacheGlobal(genomicDna, readComp);
		bool reverseStrand = false;
		uint32_t start = 0;
		uint32_t stop = 0;
		if (alignerObj.parts_.score_ > bestScore) {
			reverseStrand = true;
			if (alignerObj.alignObjectB_.seqBase_.seq_.front() == '-') {
				start = countBeginChar(alignerObj.alignObjectB_.seqBase_.seq_);
			}
			if (alignerObj.alignObjectB_.seqBase_.seq_.back() == '-') {
				stop = genomicDna.seqBase_.seq_.length()
						- countEndChar(alignerObj.alignObjectB_.seqBase_.seq_);
			}
		} else {
			reverseStrand = false;
			if (firstAlign.seqBase_.seq_.front() == '-') {
				start = countBeginChar(firstAlign.seqBase_.seq_);
			}
			if (firstAlign.seqBase_.seq_.back() == '-') {
				stop = genomicDna.seqBase_.seq_.length()
						- countEndChar(firstAlign.seqBase_.seq_);
			}
		}
		if (reverseStrand) {
			mRead->seqBase_.reverseComplementRead(true, true);
		}
		if (search == mReads.end()) {
			mReads.emplace(mRead->mipTargetName_,
					mipTargetReads(mRead->mipTargetName_, reverseStrand));
			mReads.at(mRead->mipTargetName_).addRead(mRead, start, stop);
		} else {
			search->second.addRead(mRead, start, stop);
		}
	}
	auto targetNames = getVectorOfMapKeys(mReads);
	njh::sort(targetNames);
	for(auto & targetReadName : targetNames){
		auto & targetReads = mReads.at(targetReadName);
		targetReads.setGeneralStartStop();
		targetReads.setFractionsOfReads();
	}
	return mReads;
}


void mippedGene::setUpGraph(aligner & alignerObj, comparison comp){
	auto targetNames = getVectorOfMapKeys(processedReads_);
	njh::sort(targetNames);

	for(const auto & tName : targetNames){
		const auto & reads = processedReads_.at(tName);
		for(const auto & read : reads.reads_){
			overlapGraph_.addNode(read);
		}
	}

	for(const auto & tarPos : iter::range(targetNames.size() - 1)){
		auto & targetReadFirst = processedReads_.at(targetNames[tarPos]);
		auto & targetReadSecond = processedReads_.at(targetNames[tarPos + 1]);

		if(targetReadFirst.overLap(targetReadSecond)){
			//could just try a simple suffix - prefix overlap to avoid cost of alignment
			for(const auto & firstRead : targetReadFirst.reads_){
				for(const auto & secondRead : targetReadSecond.reads_){
					alignerObj.alignCacheGlobal(firstRead->seqBase_, secondRead->seqBase_);
					auto currentComp = alignerObj.compareAlignment(firstRead->seqBase_, secondRead->seqBase_,false);
					if(comp.passErrorProfile(currentComp)){
						overlapGraph_.connectNodes(firstRead->seqBase_.name_, secondRead->seqBase_.name_);
					}
				}
			}
		}
	}
}

void mippedGene::printEdges(std::ostream & out)const {
	for(const auto & n : overlapGraph_.nodes_){
		std::cout << n->value_->seqBase_.name_ << std::endl;
		std::cout <<"\tParents" << std::endl;
		for(const auto & e : n->parentEdges_){
			auto sptr = e->par_.lock();
			std::cout << "\t" << sptr->value_->seqBase_.name_ << std::endl;
		}
		std::cout <<"\tChildren" << std::endl;
		for(const auto & e : n->childrenEdges_){
			auto sptr = e->child_.lock();
			std::cout << "\t" << sptr->value_->seqBase_.name_ << std::endl;
		}
	}
}

void mippedGene::printAllPaths(std::ostream & out)const {
	for(const auto & n : overlapGraph_.nodes_){
		if(n->parentEdges_.empty()){
			n->printDownwardPaths(out, "");
		}
	}
}

void mippedGene::printLociInfo(std::ostream & out)const {
	out << "LOCI " << processedReads_.size() << "\n";
	auto keys = getVectorOfMapKeys(processedReads_);
	njh::sort(keys);
	for(const auto & kPos : iter::range(keys.size())){
		out << "L" << kPos + 1 << " " << processedReads_.at(keys[kPos]).reads_.size()  << "\n";
	}
	for(const auto & kPos : iter::range(keys.size())){
		for(const auto & read : processedReads_.at(keys[kPos]).reads_){
			out << "L" << kPos + 1 << " A" << read->clusterId_ + 1 << " "  << roundDecPlaces(read->seqBase_.frac_, 3)<< "\n";
		}
	}
}

void mippedGene::printPossibleHaps(std::ostream & out)const {
	for(const auto & n : overlapGraph_.nodes_){
		if(n->parentEdges_.empty()){
			n->printPossibleHap(out, "", " ");
		}
	}
}

void mippedGene::setGroups(){
	uint32_t parGroup = 0;
	uint32_t childGroup = 0;
	auto targetNames = getVectorOfMapKeys(processedReads_);
	njh::sort(targetNames);
	for(const auto & targetName : targetNames){
		auto & target = processedReads_.at(targetName);
		std::unordered_map<std::string, VecStr> parentInfo;
		std::unordered_map<std::string, VecStr> childInfo;
		for(const auto & read : target.reads_){
			auto pNames = overlapGraph_.nameToNode_[read->seqBase_.name_]->getChildrenNames();
			njh::sort(pNames);
			parentInfo[vectorToString(overlapGraph_.nameToNode_[read->seqBase_.name_]->getParentNames(), ",")].emplace_back(read->seqBase_.name_);
			auto cNames = overlapGraph_.nameToNode_[read->seqBase_.name_]->getChildrenNames();
			njh::sort(cNames);
			childInfo[vectorToString(cNames, ",")].emplace_back(read->seqBase_.name_);
		}
		for(const auto & pInfo : parentInfo){
			for(const auto & name : pInfo.second){
				overlapGraph_.nameToNode_[name]->groupByParents_ = parGroup;
			}
			++parGroup;
		}
		for(const auto & cInfo : childInfo){
			for(const auto & name : cInfo.second){
				overlapGraph_.nameToNode_[name]->groupByChildren_ = childGroup;
			}
			++childGroup;
		}
	}
}

void mippedGene::printGroupInfo()const{
	std::unordered_map<uint32_t, VecStr> groups;
	for(const auto & n : overlapGraph_.nodes_){
		groups[n->groupByChildren_].emplace_back(n->value_->seqBase_.name_);
	}
	auto groupNumbers = getVectorOfMapKeys(groups);
	njh::sort(groupNumbers);
	for(const auto & g : groupNumbers){
		std::cout << "Group: " << g << "\n";
		std::cout << "\t" << vectorToString(groups[g], ", ") << "\n";
		double parTotal = 0;
		double childTotal = 0;
		for(const auto & name : groups[g]){
			parTotal += overlapGraph_.nameToNode_.at(name)->value_->seqBase_.frac_;
		}
		std::cout << "\t" << parTotal << "\n";
		auto childNames = overlapGraph_.nameToNode_.at(groups[g].front())->getChildrenNames();
		std::cout << "\t" << vectorToString(childNames, ", ") << "\n";
		for(const auto & name : childNames){
			childTotal += overlapGraph_.nameToNode_.at(name)->value_->seqBase_.frac_;
		}
		std::cout << '\t' << childTotal << std::endl;
		std::cout << std::endl;
	}
}

std::vector<readObject> mippedGene::getAlignedTargets(aligner & alignerObj) {
	std::vector<readObject> ret;
	ret.emplace_back(genomicDna_);
	auto keys = getVectorOfMapKeys(processedReads_);
	njh::sort(keys);
	for(const auto & key : keys){
		const auto & reads = processedReads_.at(key);
		for(const auto & read :reads.reads_){
			alignerObj.alignCacheGlobal(genomicDna_, read->seqBase_);
			auto firstAlign = alignerObj.alignObjectB_;
			int32_t bestScore = alignerObj.parts_.score_;
			auto readComp = read->seqBase_;
			readComp.reverseComplementRead(true,true);
			alignerObj.alignCacheGlobal(genomicDna_, readComp);
			if(alignerObj.parts_.score_ > bestScore){
				ret.emplace_back(alignerObj.alignObjectB_.seqBase_);
			}else{
				ret.emplace_back(firstAlign.seqBase_);
			}
		}
	}
	return ret;
}


std::unordered_map<std::string, table> mippedGene::findSnps(aligner & alignerObj){
	if(refGeneRecord_ == nullptr){
		std::cout << "no ref gene recrod" << std::endl;
		return std::unordered_map<std::string, table>{{"",table(VecStr{})}};

	}
	uint32_t refStart = refGeneRecord_->txStart_;;
	bool reverseStrand = false;
	if(refGeneRecord_->strand_ == '-'){
		reverseStrand = true;
		refStart = refGeneRecord_->txEnd_ - 1;
	}
	std::unordered_map<std::string, std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, mismatch>>>> allSnps;
	for(const auto & reads : processedReads_){
		std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, mismatch>>> snps;
		for(const auto & read : reads.second.reads_){
			alignerObj.alignCacheGlobal(genomicDna_, read->seqBase_);
			alignerObj.weighHomopolymers_ = false;
			alignerObj.profilePrimerAlignment(genomicDna_, read->seqBase_);
			for(const auto & m : alignerObj.comp_.distances_.mismatches_){
				uint32_t pos = m.second.refBasePos + refStart;
				if(reverseStrand){
					pos = refStart - m.second.refBasePos;
				}

				auto search = snps[pos][m.second.refBase].find(m.second.seqBase);
				if(search == snps[pos][m.second.refBase].end()){
					auto mCopy = m.second;
					mCopy.frac_ = read->seqBase_.frac_;
					snps[pos][m.second.refBase].emplace(mCopy.seqBase,mCopy);
				}else{
					search->second.frac_ += read->seqBase_.frac_;
				}
			}
		}
		allSnps[reads.first] = snps;
	}
	std::unordered_map<std::string, table> ret;
	for (const auto & targetSnps : allSnps) {
		table currentTab(VecStr {"GenomicPos", "RefChar", "SeqChar", "Chrom","inCoding" ,
				"frac"});
		for (const auto & refPos : targetSnps.second) {
			for (const auto & refBase : refPos.second) {
				for (const auto & seqBase : refBase.second) {
					currentTab.content_.emplace_back(
							toVecStr(refPos.first,refBase.first, seqBase.first, refGeneRecord_->chrom_
									, "not-checked", seqBase.second.frac_ ));
				}
			}
		}
		currentTab.sortTable("GenomicPos", false);
		ret[targetSnps.first] = currentTab;
	}

	table allAveraged(VecStr {"GenomicPos", "RefChar", "SeqChar", "Chrom","inCoding" ,
				"frac"});
	std::unordered_map<std::string, std::vector<double>> allInfos;
	for(const auto & target : ret){
		for(const auto & row : target.second.content_){
			allInfos[vectorToString(VecStr{row[target.second.getColPos("GenomicPos")],
				row[target.second.getColPos("RefChar")],
				row[target.second.getColPos("SeqChar")],
				row[target.second.getColPos("Chrom")],
				row[target.second.getColPos("inCoding")]}, "DELIM")].emplace_back(njh::lexical_cast<double>(row[target.second.getColPos("frac")]));
		}
	}
	for(const auto & info : allInfos){
		VecStr row = tokenizeString(info.first, "DELIM");
		row.emplace_back(estd::to_string(vectorMean(info.second)));
		allAveraged.content_.emplace_back(row);
	}
	allAveraged.sortTable("GenomicPos", false);
	ret["all"] = allAveraged;
	return ret;
}

std::unordered_map<std::string, table> mippedGene::findProteinSnps(aligner & alignerObj){
	if(refGeneRecord_ == nullptr){
		std::cout << "no ref gene recrod" << std::endl;
		return std::unordered_map<std::string, table>{{"",table(VecStr{})}};

	}

	std::unordered_map<std::string, std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, mismatch>>>> allSnps;
	std::unordered_map<std::string, double> barcodeCounts;
	for(const auto & reads : processedReads_){
		std::unordered_map<uint32_t, std::unordered_map<char, std::unordered_map<char, mismatch>>> snps;
		double cnt = 0;
		for(const auto & read : reads.second.reads_){
			cnt += read->seqBase_.cnt_;
			alignerObj.alignCacheGlobal(genomicDna_, read->seqBase_);
			alignerObj.weighHomopolymers_ = false;
			alignerObj.profilePrimerAlignment(genomicDna_, read->seqBase_);
			auto proteinMismatches = getProteinSnps(alignerObj, alignerObj.comp_.distances_.mismatches_);
			for(const auto & m : alignerObj.comp_.distances_.mismatches_){
				//std::cout << m.second.outputJson() << std::endl;
				uint32_t pos = m.second.refBasePos;
				auto search = snps[pos][m.second.refBase].find(m.second.seqBase);
				if(search == snps[pos][m.second.refBase].end()){
					auto mCopy = m.second;
					mCopy.frac_ = read->seqBase_.frac_;
					mCopy.freq = read->seqBase_.cnt_;
					snps[pos][m.second.refBase].emplace(mCopy.seqBase,mCopy);
				}else{
					search->second.frac_ += read->seqBase_.frac_;
					search->second.freq += read->seqBase_.cnt_;
				}
			}
		}
		barcodeCounts[reads.first] = cnt;
		allSnps[reads.first] = snps;
	}
	std::unordered_map<std::string, table> ret;
	for (const auto & targetSnps : allSnps) {
		table currentTab(VecStr {"ProteinPos", "RefProtein", "VariantProtein",
				"frac"});
		for (const auto & refPos : targetSnps.second) {
			for (const auto & refBase : refPos.second) {
				for (const auto & seqBase : refBase.second) {
					currentTab.content_.emplace_back(
							toVecStr(refPos.first,refBase.first, seqBase.first,
									 seqBase.second.frac_));
				}
			}
		}
		currentTab.sortTable("ProteinPos", false);
		ret[targetSnps.first] = currentTab;
	}
	table allAveraged(VecStr {"ProteinPos", "RefProtein", "VariantProtein" ,
				"frac"});
	std::unordered_map<std::string, std::vector<double>> allInfos;
	for(const auto & target : ret){
		for(const auto & row : target.second.content_){
			allInfos[vectorToString(VecStr{row[target.second.getColPos("ProteinPos")],
				row[target.second.getColPos("RefProtein")],
				row[target.second.getColPos("VariantProtein")]}, "DELIM")].emplace_back(njh::lexical_cast<double>(row[target.second.getColPos("frac")]));

		}
	}
	/*
	for (const auto & targetSnps : allSnps) {
		table currentTab(VecStr {"ProteinPos", "RefProtein", "VariantProtein",
				"frac", "freq", "totalFreq"});
		for (const auto & refPos : targetSnps.second) {
			for (const auto & refBase : refPos.second) {
				for (const auto & seqBase : refBase.second) {
					currentTab.content_.emplace_back(
							toVecStr(refPos.first,refBase.first, seqBase.first,
									 seqBase.second.frac_, seqBase.second.freq, barcodeCounts[targetSnps.first]));
				}
			}
		}
		currentTab.sortTable("ProteinPos", false);
		ret[targetSnps.first] = currentTab;
	}
	table allAveraged(VecStr {"ProteinPos", "RefProtein", "VariantProtein" ,
				"frac", "freq"});
	std::unordered_map<std::string, std::vector<double>> allInfos;
	std::unordered_map<std::string, std::vector<uint32_t>> allInfosFreq;
	for(const auto & target : ret){
		for(const auto & row : target.second.content_){
			allInfos[vectorToString(VecStr{row[target.second.getColPos("ProteinPos")],
				row[target.second.getColPos("RefProtein")],
				row[target.second.getColPos("VariantProtein")]}, "DELIM")].emplace_back(njh::lexical_cast<double>(row[target.second.getColPos("frac")]));

			allInfosFreq[vectorToString(VecStr{row[target.second.getColPos("ProteinPos")],
							row[target.second.getColPos("RefProtein")],
							row[target.second.getColPos("VariantProtein")]}, "DELIM")].emplace_back(njh::lexical_cast<double>(row[target.second.getColPos("freq")]));
		}
	}*/
	for(const auto & info : allInfos){
		VecStr row = tokenizeString(info.first, "DELIM");
		row.emplace_back(estd::to_string(vectorMean(info.second)));
		//row.emplace_back(estd::to_string(vectorMean(allInfosFreq[info.first])));
		allAveraged.content_.emplace_back(row);
	}
	allAveraged.sortTable("ProteinPos", false);
	ret["all"] = allAveraged;

	return ret;
}

std::map<uint32_t, mismatch> mippedGene::getProteinSnps(aligner & alignerObj, const std::map<uint32_t, mismatch> & mismatches){
	setProtein();
	auto mutProtein = getMutatedProtein(mismatches);
	//std::cout << "mutProtein: " << mutProtein.seq_ << std::endl;
	substituteMatrix proteinScore;
	proteinScore.setWithPam250();
	//std::cout << "len: " << len(protein_) << std::endl;
	//std::cout << "lenMut: " << len(mutProtein) << std::endl;
	//uint32_t maxSize = std::max(len(protein_), len(mutProtein)) + 20;
	//std::cout << "MaxSize: " << maxSize << std::endl;
	auto gapScoring = gapScoringParameters(5,1,0,0,0,0);

	//aligner proteinAligner(maxSize,gapScoring, proteinScore);
	//aligner proteinAligner(1000,gapScoringParameters(5,1,0,0,0,0), substituteMatrix::createDegenScoreMatrixCaseInsensitive(2,-2));
	//std::cout << proteinAligner.parts_.maxSize_ << std::endl;
	alignerObj.alignCacheGlobal(protein_, mutProtein);
	//std::cout << "here1" << std::endl;
	alignerObj.weighHomopolymers_ = false;
	alignerObj.profilePrimerAlignment(protein_, mutProtein);
	//std::cout << "here2" << std::endl;
	return alignerObj.comp_.distances_.mismatches_;
}

void mippedGene::setProtein(){
	if(!proteinCalculated){
		std::string genomicDna = genomicDna_.seq_;
		std::string cDna = "";
		if(refGeneRecord_->strand_ == '-'){
			genomicDna = seqUtil::reverseComplement(genomicDna, "DNA");
		}
		for(const auto & pos : iter::range(refGeneRecord_->exonStarts_.size())){
			uint32_t eStart = refGeneRecord_->exonStarts_[pos] - refGeneRecord_->txStart_;
			uint32_t eEnd = refGeneRecord_->exonEnds_[pos] - refGeneRecord_->txStart_;
			cDna.append(genomicDna.substr(eStart, eEnd - eStart));

		}

		if(refGeneRecord_->strand_ == '-'){
			cDna = seqUtil::reverseComplement(cDna, "DNA");
		}
		protein_ = seqInfo(genomicDna_.name_, seqUtil::convertToProtein(cDna));
	}

}

seqInfo mippedGene::getMutatedProtein(const std::map<uint32_t, mismatch> & mismatches)const{
	std::string genomicDna = genomicDna_.seq_;
	std::string cDna = "";
	for(const auto & m : mismatches){
		genomicDna[m.second.refBasePos] = m.second.seqBase;
	}
	/**@todo rather than doing a reverse to get the corect positions,
	 *  should just do a positino transform, just eaiser the way it is now but less efficient*/
	if(refGeneRecord_->strand_ == '-'){
		genomicDna = seqUtil::reverseComplement(genomicDna, "DNA");
	}
	for(const auto & pos : iter::range(refGeneRecord_->exonStarts_.size())){
		uint32_t eStart = refGeneRecord_->exonStarts_[pos] - refGeneRecord_->txStart_;
		uint32_t eEnd = refGeneRecord_->exonEnds_[pos] - refGeneRecord_->txStart_;
		cDna.append(genomicDna.substr(eStart, eEnd - eStart));
	}
	if(refGeneRecord_->strand_ == '-'){
		cDna = seqUtil::reverseComplement(cDna, "DNA");
	}

	return seqInfo(genomicDna_.name_, seqUtil::convertToProtein(cDna));
}


} /* namespace njhseq */
