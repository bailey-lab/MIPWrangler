#pragma once
/*
 * mipsterAnalysisRunner.hpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include <bibseq.h>
#include <bibcpp/progutils/programRunner.hpp>
#include "mipster.h"



namespace bibseq {


class mipsterAnalysisRunner : public bib::progutils::programRunner {
 public:
	mipsterAnalysisRunner();

	static int runGzExtractStitch(const bib::progutils::CmdArgs & inputCommands);

	static int mipIllumExtractByArmAndFilter(const bib::progutils::CmdArgs & inputCommands);
	static int mipIllumExtractByArmAndFilterMultiple(const bib::progutils::CmdArgs & inputCommands);

	static int mipBarcodeCorrection(const bib::progutils::CmdArgs & inputCommands);
	static int mipBarcodeCorrectionMultiple(const bib::progutils::CmdArgs & inputCommands);

	static int mipClustering(const bib::progutils::CmdArgs & inputCommands);
	static int mipClusteringMultiple(const bib::progutils::CmdArgs & inputCommands);

	static int mipPopulationClustering(const bib::progutils::CmdArgs & inputCommands);
	static int mipPopulationClusteringMultiple(const bib::progutils::CmdArgs & inputCommands);

	static int mipCorrectForContamWithSameBarcodes(const bib::progutils::CmdArgs & inputCommands);
	static int mipCorrectForContamWithSameBarcodesMultiple(const bib::progutils::CmdArgs & inputCommands);

};



} /* namespace bibseq */


