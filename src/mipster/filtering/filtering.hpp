#pragma once
/*
 * filtering.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */

#include "mipster/info/filterStats/BarcodeFilterStats.hpp"
#include "mipster/info/MipsSamplesNames.hpp"
#include "mipster/objects/MippedRead.hpp"
#include "mipster/objects/Mip.hpp"
#include "mipster/parameters.h"

namespace bibseq {

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

}  // namespace bibseq




