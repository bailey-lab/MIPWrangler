#pragma once
//

//  mipsterMipExplorerRunner.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include "mipsterMipExplorerSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterMipExplorerRunner : public bib::progutils::programRunner {
 public:
  mipsterMipExplorerRunner();
  
  static int viewMipsOnGenome(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
