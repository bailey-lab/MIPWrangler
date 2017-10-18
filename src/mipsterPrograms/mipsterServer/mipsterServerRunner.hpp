#pragma once
//

//  mipsterServerRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterServerSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterServerRunner : public bib::progutils::ProgramRunner {
 public:
  mipsterServerRunner();
  
	static int mipServerSetUp(const bib::progutils::CmdArgs & inputCommands);
	static int mipAnalysisServerSetUp(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
