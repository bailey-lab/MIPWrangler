#pragma once
/*
 * mipUtils.hpp
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
#include "mipster/mipUtils/MipNameSorter.hpp"

namespace njhseq {


/**@brief get barcode number and read number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @param barNum a reference to set to barcode number
 * @param readNum a reference to set to read number
 */
void processNameForBarReadCounts(const std::string & name,
		uint32_t & barNum, uint32_t & readNum);

/**@brief get read number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @return the read count
 */
uint32_t getReadCntFromMipName(const std::string & name);

/**@brief get barcode number, assumes name ends with _RX_BX_
 *
 * @param name name to process
 * @return the barcode count
 */
uint32_t getBarcodeCntFromMipName(const std::string & name);


}  // namespace njhseq


