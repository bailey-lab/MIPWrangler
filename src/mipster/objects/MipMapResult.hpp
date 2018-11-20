#pragma once
/*
 * MipMapResult.hpp
 *
 *  Created on: Jan 27, 2017
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


#include "mipster/common.h"
#include <njhseq/objects/BioDataObject/GenomicRegion.hpp>

namespace njhseq {

class MipMapResult {

public:

	MipMapResult(const std::string & mipName, const std::string & genomeName,
			const BamTools::BamAlignment & extAln,
			const BamTools::BamAlignment & ligAln);

	std::string mipName_;
	std::string genomeName_;
	BamTools::BamAlignment extAln_;
	BamTools::BamAlignment ligAln_;

	GenomicRegion region_;

	GenomicRegion extArmRegion_;
	GenomicRegion ligArmRegion_;

	bool isMapped() const;
	bool isConcordant() const;
	void setRegion(const BamTools::RefVector & refData);

};

/**@brief Get the results of mapping mip arms with bowtie to a genome, assumes the name of the genome is everything in filenmae up to the first underscore
 *
 * @param fnp the sorted bam file
 * @return mip map results
 */
std::vector<MipMapResult> getMipMapResults(const bfs::path & fnp, uint32_t insertSizeCutOff = 1000);

}  // namespace njhseq



