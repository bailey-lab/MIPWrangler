#pragma once
/*
 * BarcodeFilterStats.hpp
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

class BarcodeFilterStats {
public:
	class BarcodeFilterStat {
	public:
		BarcodeFilterStat(const std::string & mipTarget,
				const std::string & mipFamily);
		std::string mipTarget_;
		std::string mipFamily_;
		uint32_t initial_ = 0;
		uint32_t final_ = 0;
		std::vector<uint32_t> barCoverage_;
		uint32_t ligBarFilter_ = 0;
		uint32_t barFilter_ = 0;

		uint32_t totalFilter() const;
		static std::string printInfoHeader(const std::string & delim);
		void printInfo(std::ostream & out,const std::string & sampName, const std::string & delim) const;
	};
	BarcodeFilterStats(const std::string & sampName);
	std::string sampName_;
	std::unordered_map<std::string, BarcodeFilterStat> allStats_;
	void addFilterStat(const BarcodeFilterStat & stat);
	void printInfo(std::ostream & out, const std::string & delim) const;
};

}  // namespace njhseq

