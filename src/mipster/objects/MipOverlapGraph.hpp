#pragma once
/*
 * MipOverlapGraph.hpp
 *
 *  Created on: Jul 26, 2016
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



#include "mipster/common.h"
#include <elucidator/objects/dataContainers/graphs.h>

namespace njhseq {

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
		if(njh::has(meta, "mipFam")) {
			auto mipFamToks = njh::tokenizeString(meta.at("mipFam"), "_");
			uint32_t mipNum = estd::stou(njh::replaceString(mipFamToks.back(), "mip", ""));
			seqsByMipNum_[mipNum].emplace_back(std::make_shared<seqInfo>(getSeqBase(seq)));
		} else {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam" << std::endl;
			throw std::runtime_error {ss.str()};
		}
	}

	Json::Value genOverlapJson()const;

};







}  // namespace njhseq


