/*
 * MipExtractionStats.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */




#include "MipExtractionStats.hpp"
#include "mipster/mipUtils/MipNameSorter.hpp"

namespace bibseq {
MipExtractionStats::MipExtractionStats(const std::string & sampName):sampName_(sampName){

}


void MipExtractionStats::increaseCount(const std::string & name, SinlgeMipExtractInfo::extractCase eCase){
	stats_[name].increaseCount(eCase);
}

void MipExtractionStats::increaseIndeterminate(){
	++indeterminate_;
}

void MipExtractionStats::increaseUnmatched(){
	++totalUnmatched_;
}

void MipExtractionStats::increaseSmallFragment(){
	++smallFragmentCount_;
}

uint32_t MipExtractionStats::getGoodAmount(const std::string & name){
	return stats_[name].good_;
}

VecStr MipExtractionStats::getNames()const{
	auto ret = getVectorOfMapKeys(stats_);
	bib::sort(ret);
	return ret;
}
std::string MipExtractionStats::getNameForCase(SinlgeMipExtractInfo::extractCase eCase) {
	return SinlgeMipExtractInfo::getNameForCase(eCase);
}

bool MipExtractionStats::haveStatFor(const std::string & name)const{
	return stats_.end() != stats_.find(name);
}

std::vector<VecStr> MipExtractionStats::outputContents(const MipCollection & mips,
		const std::string & delim) {
	std::vector<VecStr> ret;
	auto statKeys = getVectorOfMapKeys(stats_);
	MipNameSorter::sort(statKeys);
	SinlgeMipExtractInfo total;
	for (const auto & k : statKeys) {
		total.good_ += stats_[k].good_;
		total.failedLig_ += stats_[k].failedLig_;
		total.failedQual_ += stats_[k].failedQual_;
		total.failedMinLen_ += stats_[k].failedMinLen_;
		total.containsNs_ += stats_[k].containsNs_;
		total.badStitch_ += stats_[k].badStitch_;
		ret.emplace_back(
				concatVecs(VecStr {sampName_, k, mips.getFamilyForTarget(k) },
						stats_[k].toVecStr()));
	}
	uint32_t grandTotal = total.getTotal() + totalUnmatched_
			+ smallFragmentCount_ + indeterminate_;
	ret.emplace_back(concatVecs(VecStr { sampName_,"indeterminate", "indeterminate" }, VecStr {
			getPercentageString(indeterminate_, grandTotal), "0", "0", "0", "0", "0", "0"  }));
	ret.emplace_back(concatVecs(VecStr { sampName_,"unmatched", "unmatched" }, VecStr {
			getPercentageString(totalUnmatched_, grandTotal), "0", "0", "0", "0", "0", "0"  }));
	ret.emplace_back(
			concatVecs(VecStr { sampName_,"smallFragment", "smallFragment" }, VecStr {
					getPercentageString(smallFragmentCount_, grandTotal), "0", "0", "0", "0",
					"0", "0"  }));
	ret.emplace_back(
			concatVecs(VecStr { sampName_,"total", "total" },
					toVecStr(grandTotal, getPercentageString(total.good_, grandTotal),
							getPercentageString(total.failedLig_, grandTotal),
							getPercentageString(total.failedMinLen_, grandTotal),
							getPercentageString(total.failedQual_, grandTotal),
							getPercentageString(total.containsNs_, grandTotal),
							getPercentageString(total.badStitch_, grandTotal))));
	return ret;
}

}

