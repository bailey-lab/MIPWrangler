#pragma once
/*
 * BarcodeInfo.hpp
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


#include "mipster/common.h"


namespace njhseq {

class BarcodeInfo {
public:
	BarcodeInfo();
	BarcodeInfo(const std::string & mipId, const std::string & extBar,
			const std::string & ligBar);
	BarcodeInfo(const std::string & mipId, const std::string & extBar);
	void setFullBar();

	std::string mipId_;   /**< the mip id the barcode belongs to */
	std::string extBar_;  /**< the extension barcode */
	std::string ligBar_;  /**< the ligation barcode */
	std::string fullBar_; /** the ext and lig barcode pasted together to create just one object to reference as a barcode*/

	Json::Value toJson() const;

};





}  // namespace njhseq

