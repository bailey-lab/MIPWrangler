/*
 * MipAnalysisDirectoryMaster.cpp
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


#include "MipAnalysisDirectoryMaster.hpp"

namespace njhseq {
MipAnalysisDirectoryMaster::MipAnalysisDirectoryMaster(
		const bfs::path & masterDir) :
		masterDir_(njh::appendAsNeededRet(masterDir.string(), "/")),
				resourceDir_(
				njh::files::join(masterDir.string(), "resources/")),
				serResourceDir_(
				njh::files::join(masterDir.string(), "serverResources/")),
				populationClusteringDir_(
				njh::files::join(masterDir.string(), "populationClustering/")),
				logsDir_(
				njh::files::join(masterDir.string(), "logs/")),
				scriptsDir_(
				njh::files::join(masterDir.string(), "scripts/")){

	//checkForDirectoriesThrow();
}

void MipAnalysisDirectoryMaster::checkForDirectoriesThrow() const {
	try {
		checkExistenceThrow(masterDir_.string());
		checkExistenceThrow(resourceDir_.string());
		checkExistenceThrow(serResourceDir_.string());
		checkExistenceThrow(populationClusteringDir_.string());
		checkExistenceThrow(logsDir_.string());
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis, not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}

std::string MipAnalysisDirectoryMaster::initMipAnalysisDirectoryStructure(const bfs::path & masterDirPath){
	njh::files::makeDir(njh::files::MkdirPar(masterDirPath.string()));
	njh::files::makeDir(masterDirPath.string(), njh::files::MkdirPar("resources/"));
	njh::files::makeDir(masterDirPath.string(), njh::files::MkdirPar("serverResources/"));
	njh::files::makeDir(masterDirPath.string(), njh::files::MkdirPar("populationClustering/"));
	njh::files::makeDir(masterDirPath.string(), njh::files::MkdirPar("logs/"));
	njh::files::makeDir(masterDirPath.string(), njh::files::MkdirPar("scripts/"));


	return njh::appendAsNeededRet(masterDirPath.string(), "/");
}

bfs::path MipAnalysisDirectoryMaster::getMipSerDir(
		const std::string & mipServerName) const {
	return njh::files::make_path(serResourceDir_,
			njh::appendAsNeededRet(mipServerName, "/"));
}


}  // namespace njhseq

