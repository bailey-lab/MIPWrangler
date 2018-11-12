#pragma once
//

//  mipsterMipExplorerRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterMipExplorerSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterMipExplorerRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterMipExplorerRunner();
  
  static int setUpViewMipsOnGenome(const njh::progutils::CmdArgs & inputCommands);
  static int viewMipsOnGenome(const njh::progutils::CmdArgs & inputCommands);


  static int mipsAgainstHaplotypes(const njh::progutils::CmdArgs & inputCommands);
};
} // namespace njhseq
