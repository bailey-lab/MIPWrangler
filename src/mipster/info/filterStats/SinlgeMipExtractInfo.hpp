#pragma once
/*
 * SinlgeMipExtractInfo.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */




#include "mipster/common.h"


namespace bibseq {

class SinlgeMipExtractInfo {
	uint32_t good_ = 0;
	uint32_t failedMinLen_ = 0;
	uint32_t failedQual_ = 0;
	uint32_t failedLig_ = 0;
	uint32_t containsNs_ = 0;
public:
	friend class MipExtractionStats;
	enum class extractCase {
		GOOD,
		BADREVERSE,
		CONTAINSNS,
		MINLENBAD,
		MAXLENBAD,
		QUALITYFAILED,
		CONTAMINATION,
		NONE
	};

	void increaseCount(extractCase eCase);
	static std::string getNameForCase(extractCase eCase);
	std::string toStr(const std::string & delim) const;
	VecStr toVecStr() const;
	uint32_t getTotal() const;
};

}



