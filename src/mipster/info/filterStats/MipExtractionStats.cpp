/*
 * MipExtractionStats.cpp
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

#include "MipExtractionStats.hpp"
#include "mipster/mipUtils/MipNameSorter.hpp"

namespace njhseq {
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
	njh::sort(ret);
	return ret;
}
std::string MipExtractionStats::getNameForCase(SinlgeMipExtractInfo::extractCase eCase) {
	return SinlgeMipExtractInfo::getNameForCase(eCase);
}

bool MipExtractionStats::haveStatFor(const std::string & name)const{
	return stats_.end() != stats_.find(name);
}


table MipExtractionStats::outputContentsJustTargets(const MipCollection & mips){
	table ret(toVecStr(VecStr{"Sample", "mipTarget", "mipFamily"}, SinlgeMipExtractInfo::toVecStrHeader(minlen_, qualCheckStr)));

	auto statKeys = getVectorOfMapKeys(stats_);
	MipNameSorter::sort(statKeys);
	for (const auto & k : statKeys) {
		ret.addRow(toVecStr(
				sampName_,
				k,
				mips.getFamilyForTarget(k),
				stats_[k].toVecStr()));
	}
	return ret;
}

table MipExtractionStats::outputContentsSummary(const MipCollection & mips){
	table ret(toVecStr(VecStr{"Sample", "total"},
			VecStr{"unmatched", "indeterminate", "smallFragment"},
			SinlgeMipExtractInfo::toVecStrHeader(minlen_, qualCheckStr)
			));
	auto statKeys = getVectorOfMapKeys(stats_);
	SinlgeMipExtractInfo total;
	for (const auto & k : statKeys) {
		total.good_ += stats_[k].good_;
		total.failedLig_ += stats_[k].failedLig_;
		total.failedQual_ += stats_[k].failedQual_;
		total.failedMinLen_ += stats_[k].failedMinLen_;
		total.containsNs_ += stats_[k].containsNs_;
		total.badStitch_ += stats_[k].badStitch_;
	}
	uint32_t grandTotal = total.getTotal() + totalUnmatched_
			+ smallFragmentCount_ + indeterminate_;
	ret.addRow(toVecStr(sampName_,
			grandTotal,
			totalUnmatched_,
			indeterminate_,
			smallFragmentCount_,
			total.toVecStr()));
	return ret;
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

