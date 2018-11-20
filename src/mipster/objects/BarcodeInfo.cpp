/*
 * BarcodeInfo.cpp
 *
 *  Created on: Dec 10, 2015
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
#include "BarcodeInfo.hpp"

namespace njhseq {

BarcodeInfo::BarcodeInfo() {
}

BarcodeInfo::BarcodeInfo(const std::string & mipId, const std::string & extBar,
		const std::string & ligBar) :
		mipId_(mipId), extBar_(extBar), ligBar_(ligBar) {
	setFullBar();
}

BarcodeInfo::BarcodeInfo(const std::string & mipId, const std::string & extBar) :
		mipId_(mipId), extBar_(extBar), ligBar_("") {
	setFullBar();
}
void BarcodeInfo::setFullBar() {
	fullBar_ = extBar_ + ligBar_;
}

Json::Value BarcodeInfo::toJson() const{
	Json::Value ret;
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));

	ret["fullBar_"] = njh::json::toJson(fullBar_);
	ret["extBar_"] = njh::json::toJson(extBar_);
	ret["ligBar_"] = njh::json::toJson(ligBar_);
	ret["mipId_"] = njh::json::toJson(mipId_);

	return ret;
}



}  // namespace njhseq

