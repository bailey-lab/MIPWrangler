#pragma once
/*
 * SetUpMaster.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */// MIPWrangler - A library for analyzing sequence data from molecular inversion probe analysis
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

#include "mipster/setUp/setUpUtils.hpp"
#include "mipster/setUp/MipAnalysisDirectoryMaster.hpp"
#include "mipster/setUp/SampleDirectoryMaster.hpp"
#include "mipster/objects/MipCollection.hpp"
#include "mipster/info/MipsSamplesNames.hpp"

#include <njhcpp/IO/BlockGzip.hpp>


namespace njhseq {

class SetUpMaster {
public:
	SetUpMaster(const bfs::path & masterDir);

	MipAnalysisDirectoryMaster directoryMaster_;

	bool development_{false};

	bfs::path mipArmFnp_;
	bfs::path mipsSampsNamesFnp_;

	std::string rawDataSuffix_;

	std::string firstReadSuffix_ = "_R1.fastq";
	std::string secondReadSuffix_ = "_R2.fastq";

	std::string mipServerName_;

	std::shared_ptr<MipCollection> mips_;
	std::shared_ptr<MipsSamplesNames> names_;

	std::unique_ptr<MultipleGroupMetaData> meta_;

	void setServerName(const std::string & mipServerName);
	bfs::path getMipSerDir() const;

	void setMipArmFnp(const bfs::path & mipArmFnp);
	void setMipsSampsNamesFnp(const bfs::path & mipsSampsNamesFnp);
	void setRawDataSuffix(const std::string & rawDataSuffix);
	void setMetaData(const bfs::path & metaFnp);

	void loadMipsSampsInfo(uint32_t allowableArmErrors);

	void createDirStructSkeleton() const;

	void createDirStructSkeleton(const bfs::path & mipSampleFile,
			const bfs::path & mipArms);

	void createTopSampleDirs() const;
	void createPopClusMipDirs(uint32_t numThreads) const;

	void makeBarcodeCorDirs(bool cacheAln) const;
	void makeClusteringDirs(bool cacheAln) const;

	//get back a vector of warnings, if empty nothing went wrong
	VecStr checkDirStruct() const;
	void checkDirStructThrow(const std::string & funcName) const;

	bool checkForRawDataForSamp(const MipFamSamp & mipSampName) const;
	bool checkForExtractedMipFamForSamp(const MipFamSamp & mipSampName) const;
	bool checkForBarCorMipFamForSamp(const MipFamSamp & mipSampName) const;
	bool checkForClusteredMipFamForSamp(const MipFamSamp & mipSampName) const;

	bool checkForExtractedMipFamForSampThrow(const MipFamSamp & mipSampName, bool cacheAln) const;
	bool checkForBarCorMipFamForSampThrow(const MipFamSamp & mipSampName, bool cacheAln) const;

	std::vector<MipFamSamp> getSamplesWithRawData(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithExtracted(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithBarCor(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithClustered(uint32_t numThreads) const;
	std::vector<MipFamSamp> getPairsWithPopClustered(uint32_t numThreads) const;
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

	//extract information
	bfs::path pathSampleExtractInfo(const MipFamSamp & mipSampName) const;
	bfs::path pathSampleExtractInfoByTarget(const MipFamSamp & mipSampName) const;
	bfs::path pathSampleExtractInfoSummary(const MipFamSamp & mipSampName) const;
	bfs::path pathSampleExtractStitchingInfo(const MipFamSamp & mipSampName) const;

	bfs::path pathMipExtractInfo(const std::string & mipTar) const;
	bfs::path pathSampleRawData(const MipFamSamp & mipSampName) const;

	bfs::path pathSampleRawDataFirstRead(const MipFamSamp & mipSampName) const;
	bfs::path pathSampleRawDataSecondRead(const MipFamSamp & mipSampName) const;

	bfs::path pathPopClusFinalHaplo(const MipFamSamp & mipSampName) const;
	//bfs::path pathPopClusOriginalHaplo(const MipFamSamp & mipSampName) const;

	/**@brief the mip family's final haplotypes for all samples, one gzip member per sample indexed by sample name in <file>.idx
	 *
	 * @param mipFam the mip family
	 * @return the path to the file
	 */
	bfs::path pathPopClusFinalHaplos(const std::string & mipFam) const;

	/**@brief whether a sample has final haplotypes for a mip family, in pathPopClusFinalHaplos() or else in the per sample file pathPopClusFinalHaplo()
	 *
	 * @param mipSampName the mip family and sample
	 * @return true if there are final haplotypes
	 */
	bool hasPopClusFinalHaplo(const MipFamSamp & mipSampName) const;

	/**@brief the fastq text of a sample's final haplotypes, from pathPopClusFinalHaplos() or else pathPopClusFinalHaplo()
	 *
	 * @param mipSampName the mip family and sample
	 * @return the fastq text
	 */
	std::string getPopClusFinalHaploText(const MipFamSamp & mipSampName) const;

	/**@brief read a sample's final haplotypes, the same as reading pathPopClusFinalHaplo() with SeqIOOptions::genFastqIn(fnp, processed)
	 *
	 * @param mipSampName the mip family and sample
	 * @param processed whether to process the read names for count/fraction
	 * @return the haplotypes
	 */
	std::vector<seqInfo> readPopClusFinalHaplo(const MipFamSamp & mipSampName, bool processed) const;

	/**@brief write a sample's final haplotypes to a gzipped fastq file if it doesn't exist yet, for code that needs them in a file
	 *
	 * @param mipSampName the mip family and sample
	 * @param fnp the file to write
	 */
	void writePopClusFinalHaploFile(const MipFamSamp & mipSampName, const bfs::path & fnp) const;

	bfs::path pathMipSampBarCorHap(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampBarCorBars(const MipFamSamp & mipSampName,
			const std::string & mipTarName) const;

	bfs::path pathSampDir(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampExtractDir(const MipTarFamSamp & mipTarSampName) const;
	bfs::path pathMipSampClusDir(const MipFamSamp & mipSampName) const;
	bfs::path pathMipSampBarCorDir(const MipFamSamp & mipSampName) const;

	bfs::path pathToAllPopInfo() const;

	void prepareMipAnalysisServer(uint32_t numThreads) const;

	table gatherExtractStats(const std::vector<MipFamSamp> & samplesExtracted,
			uint32_t numThreads) const;


	void writeAllExtractStatsFromSummary(
			const std::vector<MipFamSamp> & samplesExtracted,
			uint32_t numThreads,
			const OutOptions & outOpts) const;
	void writeAllExtractStatsFromInfoByTarget(
			const std::vector<MipFamSamp> & samplesExtracted,
			uint32_t numThreads,
			const OutOptions & outOpts) const;
	void writeAllExtractStitchStats(
			const std::vector<MipFamSamp> & samplesExtracted,
			uint32_t numThreads,
			const OutOptions & outOpts) const;

private:
	/**@brief the index of pathPopClusFinalHaplos(), read once per mip family
	 *
	 * @param mipFam the mip family
	 * @return the index, nullptr if the mip family has no such file
	 */
	std::shared_ptr<const njh::BlockGzipIndex> getPopClusFinalHaplosIndex(const std::string & mipFam) const;

	mutable std::mutex popClusFinalHaplosIndexesMut_;
	mutable std::unordered_map<std::string, std::shared_ptr<const njh::BlockGzipIndex>> popClusFinalHaplosIndexes_;





};

}  // namespace njhseq
