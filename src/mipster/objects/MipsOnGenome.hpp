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
	MipsOnGenome(const bfs::path & mainDir, uint32_t numThreads);



	bfs::path mainDir_;

	bfs::path genomeDir_;/**< directory with genome files*/
	bfs::path infoDir_;  /**< directory with info files, specifically the mip arms file*/
	bfs::path armsDir_;  /**< directory with fasta files of the just the arms of the mips*/
	bfs::path mapDir_;   /**< directory with the bowtie2 hit results of the arm seqeuences to genomes*/
	bfs::path bedsDir_;  /**< directory with the bowtie2 hits converted to bed files*/
	bfs::path fastaDir_;  /**< directory with the fasta files of bed file regions extract from the genomes and collapse to similar sequences*/

	bfs::path logDir_;

	bfs::path mipArmsFnp_;
	std::unique_ptr<MipCollection> mipArms_;

	uint32_t numThreads_ = 1;
private:
	std::string primaryGenome_ = "";
public:

	class Genome {
	public:
		Genome(const bfs::path & fnp);
		bfs::path fnp_;
		bfs::path fnpTwoBit_;

		void createTwoBit();

		void buildBowtie2Index() const;
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

	void mapArmsToGenomes();

	void genBeds();

	void genFastas();

	void setPrimaryGenome(const std::string & genome);


	bfs::path pathToMipFasta(const std::string & mipName)const;

	VecStr getMips() const;
	VecStr getGenomes() const;

	std::vector<GenomeMip> genGenomeMipPairs()const;

};

}  // namespace bibseq
