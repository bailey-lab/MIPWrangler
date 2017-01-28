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

	bfs::path genomeDir_;
	bfs::path infoDir_;
	bfs::path mapDir_;
	bfs::path bedsDir_;
	bfs::path armsDir_;
	bfs::path logDir_;

	bfs::path mipArmsFnp_;
	std::unique_ptr<MipCollection> mipArms_;

	uint32_t numThreads_ = 1;

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

};

}  // namespace bibseq

