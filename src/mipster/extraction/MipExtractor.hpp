#pragma once

/*
 * MipExtractor.hpp
 *
 *  Created on: Sep 9, 2017
 *      Author: nick
 */

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

