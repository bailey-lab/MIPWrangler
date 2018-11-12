#pragma once
/*
 * mipsterAnalysisSetUp.hpp
 *
 *  Created on: Jan 2, 2015
 *      Author: nickhathaway
 */

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


