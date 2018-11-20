/*
 * mipUtils.cpp
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

#include "mipUtils.hpp"

namespace njhseq {

void processNameForBarReadCounts(const std::string & name,
		uint32_t & barNum, uint32_t & readNum){
	auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	auto underPos = name.rfind("_");
	//std::cout << vectorToString(toVecStr(rPos, bPos, underPos, name.substr(rPos + 2, bPos - rPos - 2),name.substr(bPos + 2, underPos - bPos - 2) ), "\n") << std::endl;
	readNum = estd::stou(name.substr(rPos + 2, bPos - rPos - 2));
	barNum = estd::stou(name.substr(bPos + 2, underPos - bPos - 2));
}

uint32_t getReadCntFromMipName(const std::string & name) {
	auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	//auto underPos = name.rfind("_");
	return estd::stou(name.substr(rPos + 2, bPos - rPos - 2));
}

uint32_t getBarcodeCntFromMipName(const std::string & name) {
	//auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	auto underPos = name.rfind("_");
	return estd::stou(name.substr(bPos + 2, underPos - bPos - 2));
}

}  // namespace njhseq
