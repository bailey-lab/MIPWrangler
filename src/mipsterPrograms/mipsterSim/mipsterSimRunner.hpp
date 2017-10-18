#pragma once
//

//  mipsterSimRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterSimSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterSimRunner : public bib::progutils::ProgramRunner {
 public:
  mipsterSimRunner();
  static int createArmFastaFiles(const bib::progutils::CmdArgs & inputCommands);
  static int mipSimSetup(const bib::progutils::CmdArgs & inputCommands);
  static int searchForArms(const bib::progutils::CmdArgs & inputCommands);
  static int extractMipCaptureSequence(const bib::progutils::CmdArgs & inputCommands);
  static int simMips(const bib::progutils::CmdArgs & inputCommands);

  static int testMipExtract(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
