/*
 * mipsterAnalysisRunner.cpp
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
