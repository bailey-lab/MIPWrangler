#pragma once
//

//  mipsterSetUpRunner.hpp
//
//  Created by Nick Hathaway on 2016/02/13.
// MIPWrangler - A library for analyzing sequence data from molecular inversion probe analysis
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of MIPWrangler.
//
// MIPWrangler is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MIPWrangler is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MIPWrangler.  If not, see <http://www.gnu.org/licenses/>.
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
