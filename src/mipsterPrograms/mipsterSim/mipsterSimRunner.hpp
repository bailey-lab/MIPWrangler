#pragma once
//

//  mipsterSimRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterSimSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterSimRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterSimRunner();
  static int createArmFastaFiles(const njh::progutils::CmdArgs & inputCommands);
  static int mipSimSetup(const njh::progutils::CmdArgs & inputCommands);
  static int searchForArms(const njh::progutils::CmdArgs & inputCommands);
  static int extractMipCaptureSequence(const njh::progutils::CmdArgs & inputCommands);
  static int simMips(const njh::progutils::CmdArgs & inputCommands);

  static int testMipExtract(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
