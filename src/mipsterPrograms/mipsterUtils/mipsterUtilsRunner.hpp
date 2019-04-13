#pragma once
//

//  mipsterUtilsRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
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

#include "mipsterUtilsSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterUtilsRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterUtilsRunner();
  
	static int processProcessedMips(const njh::progutils::CmdArgs & inputCommands);
	static int alignTargets(const njh::progutils::CmdArgs & inputCommands);


	static int scanForContam(const njh::progutils::CmdArgs & inputCommands);

	static int writeOutPossibleHaplotypes(const njh::progutils::CmdArgs & inputCommands);
	static int processMipOverlapGraph(const njh::progutils::CmdArgs & inputCommands);
	static int processMipOverlapGraphSingle(const njh::progutils::CmdArgs & inputCommands);


	static int rearmTargetsAndCombine(const njh::progutils::CmdArgs & inputCommands);

	static int createExtArmFastas(const njh::progutils::CmdArgs & inputCommands);
	static int createLigArmFastas(const njh::progutils::CmdArgs & inputCommands);


	static int mipFastasToSeqTable(const njh::progutils::CmdArgs & inputCommands);
	static int extractPossibleMipCapturesFromGenome(const njh::progutils::CmdArgs & inputCommands);
	static int creatingMipArmsFromSeqs(const njh::progutils::CmdArgs & inputCommands);
	static int fixingMipBedFiles(const njh::progutils::CmdArgs & inputCommands);

	static int createPrimerFileFromArmFile(const njh::progutils::CmdArgs & inputCommands);

	static int creatingSeqTableFromDirectory(const njh::progutils::CmdArgs & inputCommands);

	static int createMipArmFromSelectedMips(const njh::progutils::CmdArgs & inputCommands);


	static int typeFinalHaplotypes(const njh::progutils::CmdArgs & inputCommands);

	static int ExtractTargetsFromGenomes(const njh::progutils::CmdArgs & inputCommands);

	static int benchmarkingForControlMixtures(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
