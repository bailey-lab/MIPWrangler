#pragma once
//

//  mipsterServerSetUp.hpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

#include <bibseq.h>
#include <bibcpp.h>
#include "mipster.h"
namespace bibseq {

class mipsterServerSetUp : public seqSetUp {

 public:
    using seqSetUp::seqSetUp;
};
} // namespace bibseq