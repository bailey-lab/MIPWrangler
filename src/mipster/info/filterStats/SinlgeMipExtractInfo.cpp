/*
 * SinlgeMipExtractInfo.cpp
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

#include "SinlgeMipExtractInfo.hpp"

namespace njhseq {



void SinlgeMipExtractInfo::increaseCount(extractCase eCase) {
	std::stringstream ss;
	switch (eCase) {
	case extractCase::GOOD:
		++good_;
		break;
	case extractCase::BADREVERSE:
		++failedLig_;
		break;
	case extractCase::MINLENBAD:
		++failedMinLen_;
		break;
	case extractCase::QUALITYFAILED:
		++failedQual_;
		break;
	case extractCase::CONTAINSNS:
		++containsNs_;
		break;
	case extractCase::BADSTITCH:
		++badStitch_;
		break;
	case extractCase::NONE:
		ss << __PRETTY_FUNCTION__ << ": Error, no case set" << std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << __PRETTY_FUNCTION__ << ": Error, bad case" << std::endl;
		throw std::runtime_error{ss.str()};
		break;
	}
}

std::string SinlgeMipExtractInfo::getNameForCase(extractCase eCase) {
	std::string ret = "";
	std::stringstream ss;
	switch (eCase) {
	case extractCase::GOOD:
		ret = "";
		break;
	case extractCase::BADREVERSE:
		ret = "_failedLigation";
		break;
	case extractCase::MINLENBAD:
		ret = "_failedMinLen";
		break;
	case extractCase::QUALITYFAILED:
		ret = "_failedQuality";
		break;
	case extractCase::CONTAINSNS:
		ret = "_containsNs";
		break;
	case extractCase::BADSTITCH:
		ret = "_badStitch";
		break;
	case extractCase::NONE:
		ss << __PRETTY_FUNCTION__ << ": Error, no case set" << std::endl;
		throw std::runtime_error { ss.str() };
		break;
	default:
		ss << __PRETTY_FUNCTION__ << ": Error, bad case" << std::endl;
		throw std::runtime_error{ss.str()};
		break;
	}
	return ret;
}

std::string SinlgeMipExtractInfo::toStr(const std::string & delim) const {
	return vectorToString(toVecStr(), delim);
}

VecStr SinlgeMipExtractInfo::toVecStrHeader(uint32_t minLen, const std::string & qualCheckStr){
	return VecStr{
		"totalMatched",
		"goodReads",
		"failedLigationArm",
		"failedMinLen(<" + estd::to_string(minLen) + ")",
		 qualCheckStr,
		"containsNs",
		"badStitch"
	};
}

VecStr SinlgeMipExtractInfo::toVecStr() const {
	uint32_t total = getTotal();
	return njhseq::toVecStr(total,
			getPercentageString(good_, total),
			getPercentageString(failedLig_, total),
			getPercentageString(failedMinLen_, total),
			getPercentageString(failedQual_, total),
			getPercentageString(containsNs_, total),
			getPercentageString(badStitch_, total));
}

uint32_t SinlgeMipExtractInfo::getTotal() const {
	return good_ + failedLig_ + failedMinLen_ + failedQual_ + containsNs_ + badStitch_;
}

}

