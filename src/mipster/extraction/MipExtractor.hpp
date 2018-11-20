#pragma once

/*
 * MipExtractor.hpp
 *
 *  Created on: Sep 9, 2017
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
#include "mipster/parameters.h"
#include "mipster/setUp.h"

namespace njhseq {

class MipExtractor {
public:

	MipExtractor();

	MipExtractor(bool verbose);

	bool verbose_ = false;

	/**@brief
	 *
	 * @param sampleIOOpts a vector of paired end sequence IO options for all samples
	 * @param mipMaster the master of the MIPs with all the info
	 * @param alignerObjForFamilyDet an aligner used to help determine which family a mip goes with
	 * @param pars the parameters for this type of extractions
	 * @param qFilPars the quality parameters for filtering
	 *
	 *
	 * @todo for the paired end reads it assumes the mate has reverse complement for now, this is a carry over when this was ported for single end sequencing, should move to not having to do this
	 *
	 */
	void extractFilterSampleForMipsPaired(const std::vector<SeqIOOptions> & sampleIOOpts,
			const SetUpMaster & mipMaster,
			aligner & alignerObjForFamilyDet,
			const mipIllumArmExtractionPars & pars);

	void extractFilterSampleForMipsPairedStitch(const std::vector<SeqIOOptions> & sampleIOOpts,
			const SetUpMaster & mipMaster,
			aligner & alignerObjForFamilyDet,
			aligner & alignerObjForStitching,
			const mipIllumArmExtractionPars & pars);

};

} /* namespace njhseq */

