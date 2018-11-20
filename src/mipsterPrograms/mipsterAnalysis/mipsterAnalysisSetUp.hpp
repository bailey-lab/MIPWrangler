#pragma once
/*
 * mipsterAnalysisSetUp.hpp
 *
 *  Created on: Jan 2, 2015
 *      Author: nickhathaway
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
#include <njhseq/programUtils.h>
#include "mipster.h"
namespace njhseq {




class mipsterAnalysisSetUp : public seqSetUp {
 public:
  // constructors
	using seqSetUp::seqSetUp;

	void setUpExtractFromRawMultiple(extractFromRawParsMultiple & pars);

	void setUpRunGzExtractStitch(runGzExtractStitchPars & pars);

	void setUpMipIllumArmExtraction(mipIllumArmExtractionPars & pars);
	void setUpMipIllumArmExtractionMultiple(
			mipIllumArmExtractionParsMultiple & pars);


	void setUpMipIllumArmExtractionPaired(mipIllumArmExtractionPars & pars);
	void setUpMipIllumArmExtractionMultiplePaired(
			mipIllumArmExtractionParsMultiple & pars);

	void setUpMipBarcodeCorrection(mipBarcodeCorrectionPars & pars);
	void setUpMipBarcodeCorrectionMultiple(
			mipBarcodeCorrectionParsMultiple & pars);

	void setUpMipCorrectForContamWithSameBarcodes(mipCorrectForContamWithSameBarcodesPars & pars);
	void setUpMipCorrectForContamWithSameBarcodesMultiple(mipCorrectForContamWithSameBarcodesParsMultiple & pars);

	void setUpMipClustering(mipClusteringPars & pars);
	void setUpMipClusteringMultiple(mipClusteringParsMultiple & pars);


	void setUpMipPopulationClustering(mipPopulationClusteringPars & pars);
	void setUpMipPopulationClusteringMultiple(mipPopulationClusteringParsMultiple & pars);




	void processDefaults(mipCorePars & pars);

	void processMultipleDefaults(mipCorePars & pars);
};


} /* namespace njhseq */


