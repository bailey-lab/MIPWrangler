/*
 * MipAnalysisDirectoryMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */




#include "MipAnalysisDirectoryMaster.hpp"

namespace bibseq {
MipAnalysisDirectoryMaster::MipAnalysisDirectoryMaster(
		const bfs::path & masterDir) :
		masterDir_(bib::appendAsNeededRet(masterDir.string(), "/")),
				resourceDir_(
				bib::files::join(masterDir.string(), "resources/")),
				serResourceDir_(
				bib::files::join(masterDir.string(), "serverResources/")),
				populationClusteringDir_(
				bib::files::join(masterDir.string(), "populationClustering/")),
				logsDir_(
				bib::files::join(masterDir.string(), "logs/")),
				scriptsDir_(
				bib::files::join(masterDir.string(), "scripts/")){

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
	bib::files::makeDir(bib::files::MkdirPar(masterDirPath.string()));
	bib::files::makeDir(masterDirPath.string(), bib::files::MkdirPar("resources/"));
	bib::files::makeDir(masterDirPath.string(), bib::files::MkdirPar("serverResources/"));
	bib::files::makeDir(masterDirPath.string(), bib::files::MkdirPar("populationClustering/"));
	bib::files::makeDir(masterDirPath.string(), bib::files::MkdirPar("logs/"));
	bib::files::makeDir(masterDirPath.string(), bib::files::MkdirPar("scripts/"));


	return bib::appendAsNeededRet(masterDirPath.string(), "/");
}

bfs::path MipAnalysisDirectoryMaster::getMipSerDir(
		const std::string & mipServerName) const {
	return bib::files::make_path(serResourceDir_,
			bib::appendAsNeededRet(mipServerName, "/"));
}


}  // namespace bibseq

