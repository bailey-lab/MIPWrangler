#pragma once
/*
 * MipMapResult.hpp
 *
 *  Created on: Jan 27, 2017
 *      Author: nick
 */



#include "mipster/common.h"
#include <elucidator/objects/BioDataObject/GenomicRegion.hpp>

namespace bibseq {

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

	bool isMapped() const;
	bool isConcordant() const;
	void setRegion(const BamTools::RefVector & refData);

};

/**@brief Get the results of mapping mip arms with bowtie to a genome, assumes the name of the genome is everything in filenmae up to the first underscore
 *
 * @param fnp the sorted bam file
 * @return mip map results
 */
std::vector<MipMapResult> getMipMapResults(const bfs::path & fnp);

}  // namespace bibseq



