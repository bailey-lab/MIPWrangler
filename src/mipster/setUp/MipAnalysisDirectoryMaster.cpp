/*
 * MipAnalysisDirectoryMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */




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

