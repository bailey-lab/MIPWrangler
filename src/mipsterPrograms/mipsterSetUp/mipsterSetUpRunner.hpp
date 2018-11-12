#pragma once
//

//  mipsterSetUpRunner.hpp
//
//  Created by Nick Hathaway on 2016/02/13.
//  Copyright (c) 2016 Nick Hathaway. All rights reserved.
//

#include "mipsterSetUpSetUp.hpp"
#include "mipster.h"
namespace njhseq {

class mipsterSetUpRunner : public njh::progutils::ProgramRunner {
 public:
  mipsterSetUpRunner();
  
	static int makeMipPopClusDirectories(const njh::progutils::CmdArgs & inputCommands);

	static int checkDirectoryStructure(const njh::progutils::CmdArgs & inputCommands);
	static int createSkeletonDirectoryStructure(const njh::progutils::CmdArgs & inputCommands);

	static int checkForRawData(const njh::progutils::CmdArgs & inputCommands);
	static int checkForExtracted(const njh::progutils::CmdArgs & inputCommands);
	static int checkForBarCor(const njh::progutils::CmdArgs & inputCommands);
	static int checkForClustered(const njh::progutils::CmdArgs & inputCommands);

	static int cloneAnalysisDirectory(const njh::progutils::CmdArgs & inputCommands);

};
} // namespace njhseq
