#pragma once
/*
 * BarcodeFilterStats.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */

#include "mipster/common.h"

namespace njhseq {

class BarcodeFilterStats {
public:
	class BarcodeFilterStat {
	public:
		BarcodeFilterStat(const std::string & mipTarget,
				const std::string & mipFamily);
		std::string mipTarget_;
		std::string mipFamily_;
		uint32_t initial_ = 0;
		uint32_t final_ = 0;
		std::vector<uint32_t> barCoverage_;
		uint32_t ligBarFilter_ = 0;
		uint32_t barFilter_ = 0;

		uint32_t totalFilter() const;
		static std::string printInfoHeader(const std::string & delim);
		void printInfo(std::ostream & out,const std::string & sampName, const std::string & delim) const;
	};
	BarcodeFilterStats(const std::string & sampName);
	std::string sampName_;
	std::unordered_map<std::string, BarcodeFilterStat> allStats_;
	void addFilterStat(const BarcodeFilterStat & stat);
	void printInfo(std::ostream & out, const std::string & delim) const;
};

}  // namespace njhseq

