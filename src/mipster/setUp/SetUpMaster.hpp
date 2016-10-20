#pragma once
/*
 * SetUpMaster.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */

#include "mipster/setUp/setUpUtils.hpp"
#include "mipster/setUp/MipAnalysisDirectoryMaster.hpp"
#include "mipster/setUp/SampleDirectoryMaster.hpp"
#include "mipster/objects/MipCollection.hpp"
#include "mipster/info/MipsSamplesNames.hpp"

namespace bibseq {

class SetUpMaster {
public:
	SetUpMaster(const std::string & masterDir);

	MipAnalysisDirectoryMaster directoryMaster_;

	std::string mipArmFnp_;
	std::string mipsSampsNamesFnp_;

	std::string rawDataSuffix_;

	std::string mipServerName_;

	std::shared_ptr<MipCollection> mips_;
	std::shared_ptr<MipsSamplesNames> names_;

	void setServerName(const std::string & mipServerName);
	bfs::path getMipSerDir() const;

	void setMipArmFnp(const std::string & mipArmFnp);
	void setMipsSampsNamesFnp(const std::string & mipsSampsNamesFnp);
	void setRawDataSuffix(const std::string & rawDataSuffix);

	void loadMipsSampsInfo(uint32_t allowableArmErrors);

	void createDirStructSkeleton() const;

	void createDirStructSkeleton(const std::string & mipSampleFile,
			const std::string & mipArms);

	void createTopSampleDirs() const;
	void createPopClusMipDirs(uint32_t numThreads) const;

	void makeBarcodeCorDirs() const;
	void makeClusteringDirs() const;

	//get back a vector of warnings, if empty nothing went wrong
	VecStr checkDirStruct() const;

	bool checkForRawDataForSamp(const MipFamSamp & mipSampName) const;
	bool checkForExtractedMipFamForSamp(const MipFamSamp & mipSampName) const;
	bool checkForBarCorMipFamForSamp(const MipFamSamp & mipSampName) const;
	bool checkForClusteredMipFamForSamp(const MipFamSamp & mipSampName) const;

	bool checkForExtractedMipFamForSampThrow(
			const MipFamSamp & mipSampName) const;
	bool checkForBarCorMipFamForSampThrow(const MipFamSamp & mipSampName) const;

	std::vector<MipFamSamp> getSamplesWithRawData(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithExtracted(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithBarCor(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithClustered(uint32_t numThreads) const;
	std::vector<MipFamSamp> getMipFamsWithPopClustered(uint32_t numThreads) const;

	VecStr getMipGroupings() const;
	VecStr getAllMipTargets() const;
	VecStr getAllMipFamilies() const;

	std::string getGroupForMipFam(const std::string & mipFamName) const;
	VecStr getMipFamiliesForMipGroup(const std::string & groupName) const;

	bfs::path pathMipPopClusSampInfo(const MipFamSamp & mipSampName) const;
	bfs::path pathSampPopClusSampInfo(const MipFamSamp & mipSampName) const;

	bfs::path pathMipPopClusPopInfo(const MipFamSamp & mipSampName) const;
	bfs::path pathMipPopClusHaplo(const MipFamSamp & mipSampName) const;

	bfs::path pathSampleExtractInfo(const MipFamSamp & mipSampName) const;

	bfs::path pathMipExtractInfo(const std::string & mipTar) const;
	bfs::path pathSampleRawData(const MipFamSamp & mipSampName) const;

	bfs::path pathPopClusFinalHaplo(const MipFamSamp & mipSampName) const;
	//bfs::path pathPopClusOriginalHaplo(const MipFamSamp & mipSampName) const;

	bfs::path pathMipSampBarCorHap(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampBarCorBars(const MipFamSamp & mipSampName,
			const std::string & mipTarName) const;

	bfs::path pathSampDir(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampExtractDir(const MipTarFamSamp & mipTarSampName) const;
	bfs::path pathMipSampClusDir(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampBarCorDir(const MipFamSamp & mipSampName) const;

	void prepareMipAnalysisServer(uint32_t numThreads) const;
	table gatherExtractStats(const std::vector<MipFamSamp> & samplesExtracted,
			uint32_t numThreads) const;

};

}  // namespace bibseq
