#pragma once
//

//  mipsterUtilsRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
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

};
} // namespace njhseq
