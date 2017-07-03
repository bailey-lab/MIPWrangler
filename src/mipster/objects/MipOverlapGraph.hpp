#pragma once
/*
 * MipOverlapGraph.hpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */




#include "mipster/common.h"
#include <elucidator/objects/dataContainers/graphs.h>

namespace bibseq {

class MipOverlapGraph {
public:

	MipOverlapGraph(const std::string & regionName,
			const comparison & allowableErrors, uint32_t minOverlap);

	std::string region_;
private:
	comparison allowableErrors_;
	uint32_t minOverlap_;
public:
	std::map<uint32_t, std::vector<std::shared_ptr<seqInfo>>> seqsByMipNum_;
	static void addMipNumToOverlapJson(Json::Value & graphJson);
	std::unique_ptr<SeqOverlapGraph> overlapGraph_;

	void genMipOverlapGraph(aligner & alignerObj);

	void resetFractions();

	template<typename T>
	void addMipSeqToMipSeqMap(const T & seq) {
		std::unordered_map<std::string, std::string> meta;
		getSeqBase(seq).processNameForMeta(meta);
		if(bib::has(meta, "mipFam")) {
			auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
			uint32_t mipNum = estd::stou(bib::replaceString(mipFamToks.back(), "mip", ""));
			seqsByMipNum_[mipNum].emplace_back(std::make_shared<seqInfo>(getSeqBase(seq)));
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam" << std::endl;
			throw std::runtime_error {ss.str()};
		}
	}

	Json::Value genOverlapJson()const;

};







}  // namespace bibseq


