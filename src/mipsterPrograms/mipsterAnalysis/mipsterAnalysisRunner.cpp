/*
 * mipsterAnalysisRunner.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"

namespace bibseq {

mipsterAnalysisRunner::mipsterAnalysisRunner()
    : bib::progutils::ProgramRunner(
    		{
	 	 	 	   addFunc("runGzExtractStitch", runGzExtractStitch, false),
					 addFunc("mipIllumExtractByArmAndFilter", mipIllumExtractByArmAndFilter, false),
					 addFunc("mipIllumExtractByArmAndFilterMultiple", mipIllumExtractByArmAndFilterMultiple, false),
					 addFunc("mipIllumExtractByArmAndFilterPaired", mipIllumExtractByArmAndFilterPaired, false),
					 addFunc("mipIllumExtractByArmAndFilterMultiplePaired", mipIllumExtractByArmAndFilterMultiplePaired, false),
					 addFunc("mipBarcodeCorrection", mipBarcodeCorrection, false),
					 addFunc("mipBarcodeCorrectionMultiple", mipBarcodeCorrectionMultiple, false),
					 addFunc("mipClustering", mipClustering, false),
					 addFunc("mipClusteringMultiple", mipClusteringMultiple, false),
					 addFunc("mipPopulationClustering", mipPopulationClustering, false),
					 addFunc("mipPopulationClusteringMultiple", mipPopulationClusteringMultiple, false),
					 addFunc("mipCorrectForContamWithSameBarcodes", mipCorrectForContamWithSameBarcodes, false),
					 addFunc("mipCorrectForContamWithSameBarcodesMultiple", mipCorrectForContamWithSameBarcodesMultiple, false),
					 addFunc("mipSkipBarcodeCorrection", mipSkipBarcodeCorrection, false),
					 addFunc("mipSkipBarcodeCorrectionMultiple", mipSkipBarcodeCorrectionMultiple, false),
					 addFunc("extractFromRaw", extractFromRaw, false)
           },//
          "mipsterAnalysis") {}
/*
 * 	static int (const bib::progutils::CmdArgs & inputCommands);
	static int (const bib::progutils::CmdArgs & inputCommands);
 */

} /* namespace bibseq */
