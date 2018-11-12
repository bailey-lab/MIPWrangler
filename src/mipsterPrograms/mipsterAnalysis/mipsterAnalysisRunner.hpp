#pragma once
/*
 * mipsterAnalysisRunner.hpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include <njhseq.h>
#include <njhcpp/progutils/programRunner.hpp>
#include "mipster.h"



namespace njhseq {


class mipsterAnalysisRunner : public njh::progutils::ProgramRunner {
 public:
	mipsterAnalysisRunner();


	static int extractFromRaw(const njh::progutils::CmdArgs & inputCommands);
	static int runGzExtractStitch(const njh::progutils::CmdArgs & inputCommands);

	static int mipSetupAndExtractByArm(const njh::progutils::CmdArgs & inputCommands);

	static int mipIllumExtractByArmAndFilter(const njh::progutils::CmdArgs & inputCommands);
	static int mipIllumExtractByArmAndFilterMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int mipIllumExtractByArmAndFilterPaired(const njh::progutils::CmdArgs & inputCommands);
	static int mipIllumExtractByArmAndFilterMultiplePaired(const njh::progutils::CmdArgs & inputCommands);

	static int mipBarcodeCorrection(const njh::progutils::CmdArgs & inputCommands);
	static int mipBarcodeCorrectionMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int mipSkipBarcodeCorrection(const njh::progutils::CmdArgs & inputCommands);
	static int mipSkipBarcodeCorrectionMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int mipClustering(const njh::progutils::CmdArgs & inputCommands);
	static int mipClusteringMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int mipPopulationClustering(const njh::progutils::CmdArgs & inputCommands);
	static int mipPopulationClusteringMultiple(const njh::progutils::CmdArgs & inputCommands);

	static int mipCorrectForContamWithSameBarcodes(const njh::progutils::CmdArgs & inputCommands);
	static int mipCorrectForContamWithSameBarcodesMultiple(const njh::progutils::CmdArgs & inputCommands);

};



} /* namespace njhseq */


