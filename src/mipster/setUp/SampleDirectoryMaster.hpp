#pragma once
/*
 * SampleDirectoryMaster.hpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */


#include "mipster/setUp/MipAnalysisDirectoryMaster.hpp"
#include "mipster/info/MipsSamplesNames.hpp"

namespace njhseq {
class SampleDirectoryMaster {
public:
	SampleDirectoryMaster(const MipAnalysisDirectoryMaster & masterDir, const MipFamSamp & mipSamp);
	void checkForAllDirectoriesThrow() const;

	void checkForExtractDirectoryThrow() const;
	void checkForBarCorDirectoryThrow() const;
	void checkForClusDirectoryThrow() const;

	void createExtractDirectory(bool overWrite) const;
	void ensureBarCorDirectoryExist() const;
	void ensureClusDirectoryExist() const;

	MipFamSamp mipSamp_;

	bfs::path masterSampleDir_;

	bfs::path extractDir_;
	bfs::path extractAlnCacheDir_;

	bfs::path barCorDir_;
	bfs::path barCorAlnCacheDir_;

	bfs::path clusDir_;
	bfs::path clusAlnCacheDir_;


	bfs::path getClusteredHapFnp(const std::string & mipName) const;

};


}  // namespace njhseq

