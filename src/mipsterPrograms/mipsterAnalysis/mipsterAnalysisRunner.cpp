/*
 * mipsterAnalysisRunner.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"

namespace njhseq {

mipsterAnalysisRunner::mipsterAnalysisRunner()
    : njh::progutils::ProgramRunner(
    		{
	 	 	 	   addFunc("runGzExtractStitch", runGzExtractStitch, true),
					 addFunc("mipIllumExtractByArmAndFilter", mipIllumExtractByArmAndFilter, true),
					 addFunc("mipIllumExtractByArmAndFilterMultiple", mipIllumExtractByArmAndFilterMultiple, true),
					 addFunc("mipIllumExtractByArmAndFilterPaired", mipIllumExtractByArmAndFilterPaired, true),
					 addFunc("mipIllumExtractByArmAndFilterMultiplePaired", mipIllumExtractByArmAndFilterMultiplePaired, true),
					 addFunc("mipBarcodeCorrection", mipBarcodeCorrection, true),
					 addFunc("mipBarcodeCorrectionMultiple", mipBarcodeCorrectionMultiple, false),
					 addFunc("mipClustering", mipClustering, true),
					 addFunc("mipClusteringMultiple", mipClusteringMultiple, false),
					 addFunc("mipPopulationClustering", mipPopulationClustering, true),
					 addFunc("mipPopulationClusteringMultiple", mipPopulationClusteringMultiple, false),
					 addFunc("mipCorrectForContamWithSameBarcodes", mipCorrectForContamWithSameBarcodes, true),
					 addFunc("mipCorrectForContamWithSameBarcodesMultiple", mipCorrectForContamWithSameBarcodesMultiple, false),
					 addFunc("mipSkipBarcodeCorrection", mipSkipBarcodeCorrection, true),
					 addFunc("mipSkipBarcodeCorrectionMultiple", mipSkipBarcodeCorrectionMultiple, false),
					 addFunc("extractFromRaw", extractFromRaw, true),
					 addFunc("mipSetupAndExtractByArm", mipSetupAndExtractByArm, false),
           },//
          "MIPWranglerAnalysis") {}
/*
 * 	static int (const njh::progutils::CmdArgs & inputCommands);
	static int (const njh::progutils::CmdArgs & inputCommands);
 */

} /* namespace njhseq */
