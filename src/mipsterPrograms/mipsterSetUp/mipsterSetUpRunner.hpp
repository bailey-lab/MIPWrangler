#pragma once
//

//  mipsterSetUpRunner.hpp
//
//  Created by Nick Hathaway on 2016/02/13.
//  Copyright (c) 2016 Nick Hathaway. All rights reserved.
//

#include "mipsterSetUpSetUp.hpp"
#include "mipster.h"
namespace bibseq {

class mipsterSetUpRunner : public bib::progutils::ProgramRunner {
 public:
  mipsterSetUpRunner();
  
	static int makeMipPopClusDirectories(const bib::progutils::CmdArgs & inputCommands);

	static int checkDirectoryStructure(const bib::progutils::CmdArgs & inputCommands);
	static int createSkeletonDirectoryStructure(const bib::progutils::CmdArgs & inputCommands);

	static int checkForRawData(const bib::progutils::CmdArgs & inputCommands);
	static int checkForExtracted(const bib::progutils::CmdArgs & inputCommands);
	static int checkForBarCor(const bib::progutils::CmdArgs & inputCommands);
	static int checkForClustered(const bib::progutils::CmdArgs & inputCommands);

	static int cloneAnalysisDirectory(const bib::progutils::CmdArgs & inputCommands);

};
} // namespace bibseq
