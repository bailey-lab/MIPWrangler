#pragma once
//

//  mipsterUtilsRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterUtilsSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterUtilsRunner : public bib::progutils::programRunner {
 public:
  mipsterUtilsRunner();
  
	static int processProcessedMips(const bib::progutils::CmdArgs & inputCommands);
	static int alignTargets(const bib::progutils::CmdArgs & inputCommands);


	static int scanForContam(const bib::progutils::CmdArgs & inputCommands);

	static int processMipOverlapGraph(const bib::progutils::CmdArgs & inputCommands);
	static int processMipOverlapGraphSingle(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
