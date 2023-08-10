/*
 * SampleDirectoryMaster.cpp
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
#include "SampleDirectoryMaster.hpp"

namespace njhseq {

SampleDirectoryMaster::SampleDirectoryMaster(
		const MipAnalysisDirectoryMaster & masterDir, const MipFamSamp & mipSamp) :mipSamp_(mipSamp),
		masterSampleDir_(
				njh::files::join(masterDir.masterDir_.string(), mipSamp.samp_ + "/")),
				extractDir_(njh::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipExtraction/")),
				extractAlnCacheDir_(njh::files::join(masterSampleDir_.string(), "alnCache_mipExtraction/")),
				barCorDir_(njh::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipBarcodeCorrection/")),
				barCorAlnCacheDir_(njh::files::join(masterSampleDir_.string(), "alnCache_mipBarcodeCorrection/")),
				clusDir_(njh::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipClustering/")),
				clusAlnCacheDir_(njh::files::join(masterSampleDir_.string(), "alnCache_mipClustering/")){
}

void SampleDirectoryMaster::checkForAllDirectoriesThrow(bool keepCache) const {
	try {
		checkExistenceThrow(masterSampleDir_.string());
		checkForExtractDirectoryThrow(keepCache);
		checkForBarCorDirectoryThrow(keepCache);
		checkForClusDirectoryThrow(keepCache);
	} catch (std::exception &e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;

		ss << e.what() << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

void SampleDirectoryMaster::checkForExtractDirectoryThrow(bool keepCache) const{
	try {
		checkExistenceThrow(extractDir_.string());
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}
void SampleDirectoryMaster::checkForBarCorDirectoryThrow(bool keepCache) const{
	try {
		checkExistenceThrow(barCorDir_.string());
		if(keepCache){
			checkExistenceThrow(barCorAlnCacheDir_.string());
		}
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}
void SampleDirectoryMaster::checkForClusDirectoryThrow(bool keepCache) const{
	try {
		checkExistenceThrow(clusDir_.string());
		if(keepCache){
			checkExistenceThrow(clusAlnCacheDir_.string());
		}
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}



void SampleDirectoryMaster::createExtractDirectory(bool overWrite) const{
	njh::files::makeDir(njh::files::MkdirPar(extractDir_.string(), overWrite));
}

void SampleDirectoryMaster::ensureBarCorDirectoryExist(bool cacheAln) const{
	njh::files::makeDirP(njh::files::MkdirPar(barCorDir_.string()));
	if(cacheAln){
		njh::files::makeDirP(njh::files::MkdirPar(barCorAlnCacheDir_.string()));
	}
}

void SampleDirectoryMaster::ensureClusDirectoryExist(bool cacheAln) const{
	njh::files::makeDirP(njh::files::MkdirPar(clusDir_.string()));
	if(cacheAln){
		njh::files::makeDirP(njh::files::MkdirPar(clusAlnCacheDir_.string()));
	}
}

bfs::path SampleDirectoryMaster::getClusteredHapFnp(const std::string & mipName)const{
	return njh::files::make_path(clusDir_.string(), mipName,
			mipName + "_clustered.fastq.gz");
}

bfs::path SampleDirectoryMaster::getClusteredInfoFnp(const std::string & mipName)const{
	return njh::files::make_path(clusDir_.string(), mipName,"info.tab.txt");
}



}  // namespace njhseq

