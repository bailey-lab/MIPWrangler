#pragma once
/*
 * mipsterAnalysisRunner.hpp
 *
 *  Created on: Dec 30, 2014
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
	static int mipSetup(const njh::progutils::CmdArgs & inputCommands);
	static int mipExtractByArm(const njh::progutils::CmdArgs & inputCommands);
	static int mipExtractByArmMultiple(const njh::progutils::CmdArgs & inputCommands);

	//mipSetup, mipExtractByArm, mipExtractByArmMultiple
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


