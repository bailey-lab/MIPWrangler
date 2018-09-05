#pragma once
/*
 * MipAnalysisDirectoryMaster.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */

#include "mipster/common.h"

namespace bibseq {

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

}  // namespace bibseq


