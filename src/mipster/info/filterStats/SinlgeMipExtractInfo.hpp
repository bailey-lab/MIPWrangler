#pragma once
/*
 * SinlgeMipExtractInfo.hpp
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


namespace njhseq {

class SinlgeMipExtractInfo {
	uint32_t good_ = 0;
	uint32_t failedMinLen_ = 0;
	uint32_t failedQual_ = 0;
	uint32_t failedLig_ = 0;
	uint32_t containsNs_ = 0;
	uint32_t badStitch_ = 0;
public:
	friend class MipExtractionStats;
	enum class extractCase {
		GOOD,
		BADREVERSE,
		CONTAINSNS,
		MINLENBAD,
		MAXLENBAD,
		QUALITYFAILED,
		CONTAMINATION,
		BADSTITCH,
		NONE
	};

	void increaseCount(extractCase eCase);
	static std::string getNameForCase(extractCase eCase);
	std::string toStr(const std::string & delim) const;
	VecStr toVecStr() const;
	static VecStr toVecStrHeader(uint32_t minLen, const std::string & qualCheckStr);
	uint32_t getTotal() const;
};

}



