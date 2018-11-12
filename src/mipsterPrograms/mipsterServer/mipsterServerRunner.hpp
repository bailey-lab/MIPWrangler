#pragma once
//

//  mipsterServerRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterServerSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterServerRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterServerRunner();
  
	static int mipAnalysisServerSetUp(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
