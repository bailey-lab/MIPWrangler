#pragma once
/*
 * MipExtractionStats.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */





#include "mipster/common.h"
#include "mipster/objects/MipCollection.hpp"
#include "mipster/info/filterStats/SinlgeMipExtractInfo.hpp"

namespace bibseq {


class MipExtractionStats {

	std::unordered_map<std::string, SinlgeMipExtractInfo> stats_;
	uint32_t totalUnmatched_ = 0;
	uint32_t indeterminate_ = 0;

	uint32_t smallFragmentCount_ = 0;

public:


	MipExtractionStats(const std::string & sampName);

	std::string sampName_;

	void increaseCount(const std::string & name,
			SinlgeMipExtractInfo::extractCase eCase);
	void increaseIndeterminate();
	void increaseUnmatched();
	void increaseSmallFragment();
	uint32_t getGoodAmount(const std::string & name);
	VecStr getNames() const;
	bool haveStatFor(const std::string & name) const;
	std::vector<VecStr> outputContents(const MipCollection & mips, const std::string & delim);
	static std::string getNameForCase(SinlgeMipExtractInfo::extractCase eCase);
};


}
