#pragma once
/*
 * SampleDirectoryMaster.hpp
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

#include "mipster/setUp/MipAnalysisDirectoryMaster.hpp"
#include "mipster/info/MipsSamplesNames.hpp"

namespace njhseq {
class SampleDirectoryMaster {
public:
	SampleDirectoryMaster(const MipAnalysisDirectoryMaster & masterDir, const MipFamSamp & mipSamp);
	void checkForAllDirectoriesThrow(bool keepCache) const;

	void checkForExtractDirectoryThrow(bool keepCache) const;
	void checkForBarCorDirectoryThrow(bool keepCache) const;
	void checkForClusDirectoryThrow(bool keepCache) const;

	void createExtractDirectory(bool overWrite) const;
	void ensureBarCorDirectoryExist(bool cacheAln) const;
	void ensureClusDirectoryExist(bool cacheAln) const;

	MipFamSamp mipSamp_;

	bfs::path masterSampleDir_;

	bfs::path extractDir_;
	bfs::path extractAlnCacheDir_;

	bfs::path barCorDir_;
	bfs::path barCorAlnCacheDir_;

	bfs::path clusDir_;
	bfs::path clusAlnCacheDir_;


	bfs::path getClusteredHapFnp(const std::string & mipName) const;
	bfs::path getClusteredInfoFnp(const std::string & mipName) const;

};


}  // namespace njhseq

