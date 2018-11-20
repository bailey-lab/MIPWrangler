#pragma once
/*
 * MipExtractionStats.hpp
 *
 *  Created on: Feb 8, 2016
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
#include "mipster/objects/MipCollection.hpp"
#include "mipster/info/filterStats/SinlgeMipExtractInfo.hpp"

namespace njhseq {


class MipExtractionStats {

	std::unordered_map<std::string, SinlgeMipExtractInfo> stats_;
	uint32_t totalUnmatched_ = 0;
	uint32_t indeterminate_ = 0;

	uint32_t smallFragmentCount_ = 0;

public:


	MipExtractionStats(const std::string & sampName);

	std::string sampName_;

	uint32_t minlen_{0};
	std::string qualCheckStr;

	void increaseCount(const std::string & name,
			SinlgeMipExtractInfo::extractCase eCase);
	void increaseIndeterminate();
	void increaseUnmatched();
	void increaseSmallFragment();
	uint32_t getGoodAmount(const std::string & name);
	VecStr getNames() const;
	bool haveStatFor(const std::string & name) const;
	std::vector<VecStr> outputContents(const MipCollection & mips, const std::string & delim);

	table outputContentsJustTargets(const MipCollection & mips);
	table outputContentsSummary(const MipCollection & mips);
	static std::string getNameForCase(SinlgeMipExtractInfo::extractCase eCase);
};


}
