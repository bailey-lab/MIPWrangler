/*
 * BarcodeFilterStats.cpp
 *
 *  Created on: Feb 8, 2016
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

#include "BarcodeFilterStats.hpp"

namespace njhseq {

BarcodeFilterStats::BarcodeFilterStat::BarcodeFilterStat(const std::string & mipTarget,
		const std::string & mipFamily) :
		mipTarget_(mipTarget), mipFamily_(mipFamily) {
}

uint32_t BarcodeFilterStats::BarcodeFilterStat::totalFilter() const {
	return ligBarFilter_ + barFilter_;
}

std::string BarcodeFilterStats::BarcodeFilterStat::printInfoHeader(const std::string & delim){
	return njh::conToStr(
			toVecStr("sample", "mipTarget", "mipFamily", "initialReadCnt", "finalReadCnt",
					"uniqueBarCnt", "avgBarCov", "medianBarCov", "totalReadCntFilter",
					"ligationBarFilterCnt", "barFilterCnt"), delim);
}

void BarcodeFilterStats::BarcodeFilterStat::printInfo(std::ostream & out,
		const std::string & sampName, const std::string & delim) const {
	out
			<< njh::conToStr(
					toVecStr(sampName, mipTarget_, mipFamily_, initial_,
							getPercentageString(final_, initial_), barCoverage_.size(),
							vectorMean(barCoverage_), vectorMedianCopy(barCoverage_),
							getPercentageString(totalFilter(), initial_),
							getPercentageString(ligBarFilter_, totalFilter()),
							getPercentageString(barFilter_, totalFilter())), delim)
			<< std::endl;
}

void BarcodeFilterStats::addFilterStat(const BarcodeFilterStat & stat) {
	if (njh::in(stat.mipTarget_, allStats_)) {
		throw std::runtime_error { "Error in : " + std::string(__PRETTY_FUNCTION__)
				+ " already contain stat for " + stat.mipTarget_ };
	}
	allStats_.emplace(stat.mipTarget_, stat);
}

BarcodeFilterStats::BarcodeFilterStats(const std::string & sampName) :
		sampName_(sampName) {

}

void BarcodeFilterStats::printInfo(std::ostream & out,
		const std::string & delim) const {
	out << BarcodeFilterStat::printInfoHeader(delim) << std::endl;
	auto keys = getVectorOfMapKeys(allStats_);
	njh::sort(keys);
	for (const auto & key : keys) {
		allStats_.at(key).printInfo(out,sampName_, delim);
	}
}

}  // namespace njhseq

