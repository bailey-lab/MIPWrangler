#pragma once
//

//  mipsterMipTesterRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterMipTesterSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterMipTesterRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterMipTesterRunner();
  
  static int callMircosateliteSizes(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
