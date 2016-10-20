/*
 * MipOverlapGraph.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */




#include "MipOverlapGraph.hpp"

namespace bibseq {

MipOverlapGraph::MipOverlapGraph(const std::string & regionName,
		const comparison & allowableErrors, uint32_t minOverlap) :
		region_(regionName), allowableErrors_(allowableErrors), minOverlap_(
				minOverlap) {
}

void MipOverlapGraph::resetFractions(){
	for(const auto & mipSubReg :seqsByMipNum_){
		double total = 0;
		for(const auto & seq : mipSubReg.second){
			total += seq->cnt_;
		}
		for(auto & seq : mipSubReg.second){
			seq->frac_ = seq->cnt_/total;
		}
	}
}

void MipOverlapGraph::genMipOverlapGraph(
		aligner & alignerObj) {
	overlapGraph_ = std::make_unique<SeqOverlapGraph>();
	//check to see if target should be reversed
	for (const auto & mips : seqsByMipNum_) {
		auto search = seqsByMipNum_.find(mips.first + 1);
		if (seqsByMipNum_.end() != search) {
			for(const auto & seq1 : mips.second){
				uint32_t forCount = 0;
				uint32_t revCount = 0;
				for(const auto & seq2 : search->second){
					alignerObj.alignCacheGlobal(*seq1, *seq2);
					alignerObj.profilePrimerAlignment(*seq1, *seq2);
					auto forComp = alignerObj.comp_;
					auto seq2Rev = *seq2;
					seq2Rev.reverseComplementRead(true, true);
					alignerObj.alignCacheGlobal(*seq1, seq2Rev);
					alignerObj.profilePrimerAlignment(*seq1, seq2Rev);
					bool passFor = allowableErrors_.passErrorProfile(forComp) && forComp.distances_.query_.covered_ >=minOverlap_;
					bool passRev = allowableErrors_.passErrorProfile(alignerObj.comp_) && alignerObj.comp_.distances_.query_.covered_ >=minOverlap_;
					if(passFor && passRev){
						if(forComp.distances_.query_.covered_ >= alignerObj.comp_.distances_.query_.covered_){
							++forCount;
						}else{
							++revCount;
						}
					}else if(passFor){
						++forCount;
					}else if(passRev){
						++revCount;
					}
				}
				if(revCount > forCount){
					for(auto & seq : search->second){
						seq->reverseComplementRead(true, true);
					}
				}
			}
		}

		//add to graph
		for(const auto & seq : mips.second){
			overlapGraph_->addNode(seq);
		}
	}

	for (const auto & mips : seqsByMipNum_) {
		auto search = seqsByMipNum_.find(mips.first + 1);
		if (seqsByMipNum_.end() != search) {
			for(const auto & seq1 : mips.second){
				for(const auto & seq2 : search->second){
					alignerObj.alignCacheGlobal(*seq1, *seq2);
					alignerObj.profilePrimerAlignment(*seq1, *seq2);
					if(allowableErrors_.passErrorProfile(alignerObj.comp_)
							&& alignerObj.comp_.distances_.query_.covered_ >= minOverlap_){
						bool seq1GapFront = '-' == alignerObj.alignObjectA_.seqBase_.seq_.front();
						bool seq1GapBack  = '-' == alignerObj.alignObjectA_.seqBase_.seq_.back();
						bool seq2GapFront = '-' == alignerObj.alignObjectB_.seqBase_.seq_.front();
						bool seq2GapBack  = '-' == alignerObj.alignObjectB_.seqBase_.seq_.back();
						if ((seq1GapFront && seq1GapBack)
								|| (seq2GapFront && seq2GapBack)) {
							//no overlap
						} else {
							if (seq1GapBack || seq2GapFront) {
								overlapGraph_->addEdge(seq1->name_, seq2->name_);
							} else {
								overlapGraph_->addEdge(seq2->name_, seq1->name_);
							}
						}
					}
				}
			}
		}
	}
}

void MipOverlapGraph::addMipNumToOverlapJson(Json::Value & graphJson) {
	for (auto & node : graphJson["nodes"]) {
		std::unordered_map<std::string, std::string> meta;
		seqInfo seq(node["name"].asString(), "A");
		seq.processNameForMeta(meta);
		if (bib::has(meta, "mipFam")) {
			auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
			uint32_t mipNum = bib::lexical_cast<uint32_t>(
					bib::replaceString(mipFamToks.back(), "mip", ""));
			node["mipNum"] = mipNum;
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam"
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
	}
}

Json::Value MipOverlapGraph::genOverlapJson() const {
	auto sankeyOutput = overlapGraph_->createSankeyOutput();
	addMipNumToOverlapJson(sankeyOutput);
	return sankeyOutput;
}

}  // namespace bibseq
