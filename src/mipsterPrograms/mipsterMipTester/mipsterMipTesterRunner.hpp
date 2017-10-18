#pragma once
//

//  mipsterMipTesterRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterMipTesterSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterMipTesterRunner : public bib::progutils::ProgramRunner {
 public:
  mipsterMipTesterRunner();
  
  static int callMircosateliteSizes(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
