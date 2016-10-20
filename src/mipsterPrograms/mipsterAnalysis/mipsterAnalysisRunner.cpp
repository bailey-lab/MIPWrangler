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
    : bib::progutils::programRunner(
    		{
	 	 	 	   addFunc("runGzExtractStitch", runGzExtractStitch, false),
					 addFunc("mipIllumExtractByArmAndFilter", mipIllumExtractByArmAndFilter, false),
					 addFunc("mipIllumExtractByArmAndFilterMultiple", mipIllumExtractByArmAndFilterMultiple, false),
					 addFunc("mipBarcodeCorrection", mipBarcodeCorrection, false),
					 addFunc("mipBarcodeCorrectionMultiple", mipBarcodeCorrectionMultiple, false),
					 addFunc("mipClustering", mipClustering, false),
					 addFunc("mipClusteringMultiple", mipClusteringMultiple, false),
					 addFunc("mipPopulationClustering", mipPopulationClustering, false),
					 addFunc("mipPopulationClusteringMultiple", mipPopulationClusteringMultiple, false),
					 addFunc("mipCorrectForContamWithSameBarcodes", mipCorrectForContamWithSameBarcodes, false),
					 addFunc("mipCorrectForContamWithSameBarcodesMultiple", mipCorrectForContamWithSameBarcodesMultiple, false)
           },//
          "mipsterAnalysis") {}


} /* namespace bibseq */
