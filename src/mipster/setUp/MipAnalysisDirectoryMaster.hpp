#pragma once
/*
 * MipAnalysisDirectoryMaster.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */
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
#include "mipster/common.h"

namespace njhseq {

class MipAnalysisDirectoryMaster {
public:
	MipAnalysisDirectoryMaster(const bfs::path & masterDir) ;
	void checkForDirectoriesThrow() const;
	bfs::path masterDir_;
	bfs::path resourceDir_;
	bfs::path serResourceDir_;
	bfs::path populationClusteringDir_;
	bfs::path logsDir_;
	bfs::path scriptsDir_;

	static std::string initMipAnalysisDirectoryStructure(const bfs::path & masterDirPath);

	bfs::path getMipSerDir(const std::string & mipServerName) const;
};

}  // namespace njhseq


