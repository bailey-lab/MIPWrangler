/*
 * SampleDirectoryMaster.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */


#include "SampleDirectoryMaster.hpp"

namespace bibseq {

SampleDirectoryMaster::SampleDirectoryMaster(
		const MipAnalysisDirectoryMaster & masterDir, const MipFamSamp & mipSamp) :mipSamp_(mipSamp),
		masterSampleDir_(
				bib::files::join(masterDir.masterDir_.string(), mipSamp.samp_ + "/")),
				extractDir_(bib::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipExtraction/")),
				extractAlnCacheDir_(bib::files::join(masterSampleDir_.string(), "alnCache_mipExtraction/")),
				barCorDir_(bib::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipBarcodeCorrection/")),
				barCorAlnCacheDir_(bib::files::join(masterSampleDir_.string(), "alnCache_mipBarcodeCorrection/")),
				clusDir_(bib::files::join(masterSampleDir_.string(), mipSamp.samp_ + "_mipClustering/")),
				clusAlnCacheDir_(bib::files::join(masterSampleDir_.string(), "alnCache_mipClustering/")){
}

void SampleDirectoryMaster::checkForAllDirectoriesThrow() const{
	try {
		checkExistenceThrow(masterSampleDir_.string());
		checkForExtractDirectoryThrow();
		checkForBarCorDirectoryThrow();
		checkForClusDirectoryThrow();
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;

		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}

void SampleDirectoryMaster::checkForExtractDirectoryThrow() const{
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
void SampleDirectoryMaster::checkForBarCorDirectoryThrow() const{
	try {
		checkExistenceThrow(barCorDir_.string());
		checkExistenceThrow(barCorAlnCacheDir_.string());
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}
void SampleDirectoryMaster::checkForClusDirectoryThrow() const{
	try {
		checkExistenceThrow(clusDir_.string());
		checkExistenceThrow(clusAlnCacheDir_.string());
	} catch (std::exception & e) {
		std::stringstream ss;
		ss << "Error for mip analysis for sample:" << mipSamp_.samp_ << " , not all required directories were found,"
			 << " check to make sure you are in the correct directory" << std::endl;
		ss << e.what() << std::endl;
		throw std::runtime_error { ss.str() };
	}
}



void SampleDirectoryMaster::createExtractDirectory(bool overWrite) const{
	bib::files::makeDir(bib::files::MkdirPar(extractDir_.string(), overWrite));
}

void SampleDirectoryMaster::ensureBarCorDirectoryExist() const{
	bib::files::makeDirP(bib::files::MkdirPar(barCorDir_.string()));
	bib::files::makeDirP(bib::files::MkdirPar(barCorAlnCacheDir_.string()));
}

void SampleDirectoryMaster::ensureClusDirectoryExist() const{
	bib::files::makeDirP(bib::files::MkdirPar(clusDir_.string()));
	bib::files::makeDirP(bib::files::MkdirPar(clusAlnCacheDir_.string()));
}

bfs::path SampleDirectoryMaster::getClusteredHapFnp(const std::string & mipName)const{
	return bib::files::make_path(clusDir_.string(), mipName,
			mipName + "_clustered.fastq");
}





}  // namespace bibseq

