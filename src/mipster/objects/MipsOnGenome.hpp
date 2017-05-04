#pragma once
/*
 * MipsOnGenome.hpp
 *
 *  Created on: Jan 28, 2017
 *      Author: nick
 */

#include "mipster/objects/MipCollection.hpp"

namespace bibseq {

class MipsOnGenome {
public:
	MipsOnGenome(const bfs::path & mainDir, const bfs::path & mainInputDir,
			uint32_t numThreads);

	bfs::path mainDir_;
	bfs::path mainInputDir_;

	bfs::path genomeDir_;/**< directory with genome files*/
	bfs::path infoDir_; /**< directory with info files, specifically the mip arms file*/
	bfs::path armsDir_; /**< directory with fasta files of the just the arms of the mips*/
	bfs::path mapDir_; /**< directory with the bowtie2 hit results of the arm seqeuences to genomes*/
	bfs::path bedsDir_; /**< directory with the bowtie2 hits converted to bed files*/
	bfs::path fastaDir_; /**< directory with the fasta files of bed file regions extract from the genomes and collapse to similar sequences*/
	bfs::path tablesDir_; /**< directory to keep some summary tables*/

	bfs::path logDir_;

	bfs::path mipArmsFnp_;
	std::unique_ptr<MipCollection> mipArms_;

	uint32_t numThreads_ = 1;

private:
	std::string primaryGenome_ = "";
	std::set<std::string> selectedGenomes_;
public:

	std::string getPrimaryGenome();

	class Genome {
	public:
		Genome(const bfs::path & fnp);
		bfs::path fnp_;
		bfs::path fnpTwoBit_;

		void createTwoBit();

		void buildBowtie2Index() const;

		Json::Value chromosomeLengths() const;
	};

	std::unordered_map<std::string, std::unique_ptr<Genome>> genomes_;

	void checkInputThrow() const;

	void loadInGenomes();
	void setUpGenomes();
	void loadInArms();
	void createArmFiles();

	struct GenomeMip {
		std::string genome_;
		std::string mip_;

		std::string uid(const std::string & sep = "_") const;
	};

	void setMipArmsFnp(const bfs::path & mipArmsFnp);

	void mapArmsToGenomes();
	void genBeds(const comparison & allowableError);
	void genFastas();

	void genTables() const;

	void setPrimaryGenome(const std::string & genome);

	void setSelectedGenomes(const std::set<std::string> & genomes);
	void setSelectedGenomes(const VecStr & genomes);

	bfs::path pathToMipFasta(const std::string & mipName) const;
	bfs::path pathToMipBed(const std::string & mipName,
			const std::string & genomeName) const;
	bfs::path pathToAllInfoPrimaryGenome() const;
	bfs::path pathToAllInfoAllGenomes() const;
	bfs::path pathToExtractionCounts() const;

	VecStr getMips() const;
	VecStr getGenomes() const;

	std::vector<GenomeMip> genGenomeMipPairs() const;

	table getGenomeLocsForMipTar(const std::string & tar) const;
	table getGenomeLocsForAllMipTars() const;
	table getGenomeLocsForGenome(const std::string & genome) const;

	table getMipRegionStatsForGenome(const std::string & genome) const;

	table getMipTarStatsForGenome(const std::string & genome,
			const VecStr & mipTars, bool allRecords = false) const;
	table getMipTarStatsForGenomes(const VecStr & genomes, const VecStr & mipTars,
			bool allRecords = false) const;

	table genExtractionNumberTable() const;

};

}  // namespace bibseq

