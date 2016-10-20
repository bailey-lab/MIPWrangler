/*
 * BarcodeFilterStats.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */

#include "BarcodeFilterStats.hpp"

namespace bibseq {

BarcodeFilterStats::BarcodeFilterStat::BarcodeFilterStat(const std::string & mipTarget,
		const std::string & mipFamily) :
		mipTarget_(mipTarget), mipFamily_(mipFamily) {
}

uint32_t BarcodeFilterStats::BarcodeFilterStat::totalFilter() const {
	return ligBarFilter_ + barFilter_;
}

std::string BarcodeFilterStats::BarcodeFilterStat::printInfoHeader(const std::string & delim){
	return bib::conToStr(
			toVecStr("sample", "mipTarget", "mipFamily", "initialReadCnt", "finalReadCnt",
					"uniqueBarCnt", "avgBarCov", "medianBarCov", "totalReadCntFilter",
					"ligationBarFilterCnt", "barFilterCnt"), delim);
}

void BarcodeFilterStats::BarcodeFilterStat::printInfo(std::ostream & out,
		const std::string & sampName, const std::string & delim) const {
	out
			<< bib::conToStr(
					toVecStr(sampName, mipTarget_, mipFamily_, initial_,
							getPercentageString(final_, initial_), barCoverage_.size(),
							vectorMean(barCoverage_), vectorMedianCopy(barCoverage_),
							getPercentageString(totalFilter(), initial_),
							getPercentageString(ligBarFilter_, totalFilter()),
							getPercentageString(barFilter_, totalFilter())), delim)
			<< std::endl;
}

void BarcodeFilterStats::addFilterStat(const BarcodeFilterStat & stat) {
	if (bib::in(stat.mipTarget_, allStats_)) {
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
	bib::sort(keys);
	for (const auto & key : keys) {
		allStats_.at(key).printInfo(out,sampName_, delim);
	}
}

}  // namespace bibseq

