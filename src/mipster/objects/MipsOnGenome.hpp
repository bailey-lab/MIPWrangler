#pragma once
/*
 * MipsOnGenome.hpp
 *
 *  Created on: Jan 28, 2017
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
#include "mipster/objects/MipCollection.hpp"
#include <njhseq/objects/BioDataObject/reading.hpp>
#include <SeekDeep/utils.h>

namespace njhseq {

class MipsOnGenome {
public:

	struct pars{

		bfs::path mainDir = "";
		bfs::path inputDir = "";
		bfs::path mipArmsFnp = "";
		comparison allowableError;
		bool removeBeds = false;

		MultiGenomeMapper::inputParameters gMapperPars_;


	};

	MipsOnGenome(const pars & inputParameters);

	pars inputParameters_;


	bfs::path armsDir_; /**< directory with fasta files of the just the arms of the mips*/
	bfs::path mapDir_; /**< directory with the bowtie2 hit results of the arm seqeuences to genomes*/
	bfs::path bedsDir_; /**< directory with the bowtie2 hits converted to bed files*/
	bfs::path bedsPerGenomeDir_;
	bfs::path fastaDir_; /**< directory with the fast files of bed file regions extract from the genomes and collapse to similar sequences*/
	bfs::path fastaByFamilyDir_;
	bfs::path tablesDir_; /**< directory to keep some summary tables*/

	bfs::path logDir_;

	std::unique_ptr<MipCollection> mipArms_;

	MultiGenomeMapper gMapper_;


private:
public:


	void checkInputThrow() const;


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


	bfs::path pathToMipFasta(const std::string & mipName) const;
	bfs::path pathToMipFastaWithoutArms(const std::string & mipName) const;
	bfs::path pathToMipFamilyFasta(const std::string & mipFamily) const;
	bfs::path pathToMipFamilyFastaWithoutArms(const std::string & mipFamily) const;


	bfs::path pathToMipExtArmFasta(const std::string & mipName) const;
	bfs::path pathToMipLigArmFasta(const std::string & mipName) const;

	bfs::path pathToMipBed(const std::string & mipName,
			const std::string & genomeName) const;
	bfs::path pathToAllInfoPrimaryGenome() const;
	bfs::path pathToAllInfoForGenome(const std::string & genome) const;
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

	const static VecStr getMipTarStatsForGenomeHeader_;



	table genExtractionNumberTable() const;

};

}  // namespace njhseq

