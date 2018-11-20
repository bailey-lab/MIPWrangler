#pragma once
/*
 * filtering.hpp
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

#include "mipster/info/filterStats/BarcodeFilterStats.hpp"
#include "mipster/info/MipsSamplesNames.hpp"
#include "mipster/objects/MippedRead.hpp"
#include "mipster/objects/Mip.hpp"
#include "mipster/parameters.h"

namespace njhseq {

//filtering on mip arms




//filtering on barcodes
void filterOnMultipleLigBar(
		const std::unordered_map<std::string,
				std::vector<std::shared_ptr<MippedRead>>>& ligBarReads,
				std::vector<std::vector<std::shared_ptr<MippedRead>>>& readsPerBarcode,
				BarcodeFilterStats::BarcodeFilterStat& tarStat);

std::shared_ptr<MippedRead> filterWithBarcodeCoverage(
		const std::vector<std::shared_ptr<MippedRead>> & barReads,aligner & alignerObj,
		const SeqSetUpPars & setUpPars, const mipBarcodeCorrectionPars & pars,
		BarcodeFilterStats::BarcodeFilterStat& tarStat);

}  // namespace njhseq




