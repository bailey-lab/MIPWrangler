/*
 * MipCollection.cpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include "MipCollection.hpp"
#include "mipster/mipUtils.h"
#include <unordered_map>

namespace bibseq {

void MipCollection::setAllAllowableArmError(uint32_t allowableArmError){
	for(auto & mip : mips_){
		mip.second.setAllowableErrorInArm(allowableArmError);
	}
}

void MipCollection::setAllWiggleRoomInArm(uint32_t wiggleRoom){
	for(auto & mip : mips_){
		mip.second.setWiggleRoomInArm(wiggleRoom);
	}
}

void MipCollection::setAllMinimumExpectedLen(size_t minimumExpectedLen){
	for(auto & mip : mips_){
		mip.second.setMinimumExpectedLen(minimumExpectedLen);
	}
}


VecStr MipCollection::getMipsForFamily(const VecStr & families) const{
	VecStr ret;
	for(const auto family : families){
		addOtherVec(ret, getMipsForFamily(family));
	}
	MipNameSorter::sort(ret);
	return ret;
}

VecStr MipCollection::getMipFamsForRegion(const std::string & region) const{
	std::set<std::string> retSet;
	for(const auto & m : mips_){
		if(region  == m.second.locGrouping_){
			retSet.insert(m.second.familyName_);
		}
	}
	VecStr ret(retSet.begin(), retSet.end());
	MipNameSorter::sort(ret);
	return ret;
}


VecStr MipCollection::getMipTarsForRegion(const std::string & region) const{
	VecStr ret;
	for(const auto & m : mips_){
		if(region  == m.second.locGrouping_){
			ret.emplace_back(m.first);
		}
	}
	MipNameSorter::sort(ret);
	return ret;
}


VecStr MipCollection::getMipTarsForRegions(const VecStr & regions) const{
	VecStr ret;
	for(const auto & region : regions){
		addOtherVec(ret, getMipTarsForRegion(region));
	}
	MipNameSorter::sort(ret);
	return ret;
}

VecStr MipCollection::getMipFamsForRegions(const VecStr & regions) const{
	VecStr ret;
	for(const auto & region : regions){
		addOtherVec(ret, getMipFamsForRegion(region));
	}
	MipNameSorter::sort(ret);
	return ret;
}

VecStr MipCollection::getMipTars() const{
	auto ret = getVectorOfMapKeys(mips_);
	MipNameSorter::sort(ret);
	return ret;
}

VecStr MipCollection::getMipRegions() const {
	std::set<std::string> locSet;
	for (const auto & m : mips_) {
		locSet.insert(m.second.locGrouping_);
	}
	VecStr ret { locSet.begin(), locSet.end() };
	MipNameSorter::sortByRegion(ret);
	return ret;
}

VecStr MipCollection::getMipFamilies() const {
	VecStr ret = mipFamilies_;
	MipNameSorter::sortByRegion(ret);
	return ret;
}

VecStr MipCollection::getMipRegionsForFams(const VecStr & mipFams) const{
	std::set<std::string> locSet;
	for (const auto & m : mips_) {
		if(bib::in(m.second.familyName_, mipFams)){
			locSet.insert(m.second.locGrouping_);
		}
	}
	VecStr ret{ locSet.begin(), locSet.end() };
	MipNameSorter::sortByRegion(ret);
	return ret;
}


MipCollection::MipCollection(const bfs::path & mipArmIdFile,
		uint32_t allowableArmError) :mipArmIdFnp_(mipArmIdFile),
		allowableArmError_(allowableArmError) {
	table mipInfo(mipArmIdFile.string(), "whitespace", true);
	VecStr neededColumns = { "mip_id", "extension_arm", "ligation_arm",
			"mip_family", "extension_barcode_length", "ligation_barcode_length",
			"gene_name", "mipset"};
	VecStr columnsNotFound;
	for (const auto & col : neededColumns) {
		if (!bib::in(col, mipInfo.columnNames_)) {
			columnsNotFound.emplace_back(col);
		}
	}
	if (!columnsNotFound.empty()) {
		std::stringstream ss;
		ss << "Need to have " << vectorToString(neededColumns, ",") << std::endl;
		ss << "Did not find " << vectorToString(columnsNotFound, ",") << std::endl;
		ss << "Only have " << vectorToString(mipInfo.columnNames_, ",")
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	for (const auto & row : mipInfo.content_) {
		if(std::all_of(row.begin(), row.end(), [](const std::string & col){ return "" == col;})){
			continue;
		}
		auto currentMipName = row[mipInfo.getColPos("mip_id")];
		if(bib::in(currentMipName, mips_)){
			std::stringstream ss;
			ss << "Error in: " << __PRETTY_FUNCTION__ <<  ": Collection already contains " << currentMipName << std::endl;
			throw std::runtime_error{ss.str()};
		}
		mips_[row[mipInfo.getColPos("mip_id")]] = Mip(
				estd::stou(
						row[mipInfo.getColPos("extension_barcode_length")]),
						estd::stou(
						row[mipInfo.getColPos("ligation_barcode_length")]),
				row[mipInfo.getColPos("ligation_arm")],
				row[mipInfo.getColPos("extension_arm")],
				row[mipInfo.getColPos("mip_id")],
				row[mipInfo.getColPos("mip_family")],
				row[mipInfo.getColPos("gene_name")],
				row[mipInfo.getColPos("mipset")]);
		mips_[row[mipInfo.getColPos("mip_id")]].setAllowableErrorInArm(allowableArmError);
	}
	std::set<std::string> mipFamilies;
	for (const auto & mip : mips_) {
		mipFamilies.insert(mip.second.familyName_);
		mipNamesForFamily_[mip.second.familyName_].emplace_back(mip.second.name_);
	}
	mipFamilies_ = VecStr(mipFamilies.begin(), mipFamilies.end());
}


VecStr MipCollection::getMipsForFamily(const std::string & family)const{
	VecStr ret;
	auto search = mipNamesForFamily_.find(family);
	if(search != mipNamesForFamily_.end()){
		ret = search->second;
	}
	MipNameSorter::sort(ret);
	return ret;
}

std::string MipCollection::getFamilyForTarget(const std::string & mipTarget) const{
	std::string ret = "";
	auto search = mips_.find(mipTarget);
	if(search != mips_.end()){
		ret = search->second.familyName_;
	}else{
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", couldn't find mip name: " << mipTarget  << std::endl;
		throw std::runtime_error{ss.str()};
	}
	return ret;
}

Mip MipCollection::determineBestMipInFamily(const seqInfo & read, Mip mip,
		aligner & alignerObjForGroupDet) const{
	// if part of a mip family which could have multiple mips
	// with same extraction arm, find best fitting ligatiion arm;
	if (mipNamesForFamily_.at(mip.familyName_).size() > 1) {
		double bestScore = 0.0;
		seqInfo backSeq;
		uint32_t sizeOfBackSeq = 2
				* (mip.ligationArm_.size() + mip.ligBarcodeLen_ + mip.wiggleRoomArm_);
		if (sizeOfBackSeq < read.seq_.size()) {
			backSeq = read.getSubRead(
					read.seq_.size() - sizeOfBackSeq);
		} else {
			backSeq = read;
		}
		for (const auto & mipName : mipNamesForFamily_.at(mip.familyName_)) {
			const auto & otherMip = mips_.at(mipName);
			alignerObjForGroupDet.alignCacheLocal(backSeq, otherMip.ligationArmObj_);
			//normalize score to arm length
			double currentScore = alignerObjForGroupDet.parts_.score_
					/ static_cast<double>(len(otherMip.ligationArmObj_));
			auto armPosMotifCurrent = otherMip.getPossibleExtArmPos(read);
			auto ligArmPosMotifCurrent = otherMip.getPossibleLigArmPos(read);
			//ensure that even if the ligation arm is better that the extension arm still fits
			//this could happen if there was a chimeric event that switched arms
			if (!armPosMotifCurrent.empty()
					&& !ligArmPosMotifCurrent.empty()) {
				if (currentScore > bestScore) {
					bestScore = currentScore;
					mip = mips_.at(mipName);
				}
			}
		}
	}
	return mip;
}

Mip MipCollection::determineBestMipInFamily(const PairedRead & seq, Mip mip,
		aligner & alignerObjForGroupDet) const{
	// if part of a mip family which could have multiple mips
	// with same extraction arm, find best fitting ligatiion arm;
	if (mipNamesForFamily_.at(mip.familyName_).size() > 1) {
		double bestScore = 0.0;
		seqInfo backSeq;
		uint32_t sizeOfBackSeq = 2
				* (mip.ligationArm_.size() + mip.ligBarcodeLen_ + mip.wiggleRoomArm_);
		if (sizeOfBackSeq < seq.mateSeqBase_.seq_.size()) {
			backSeq = seq.mateSeqBase_.getSubRead(
					seq.mateSeqBase_.seq_.size() - sizeOfBackSeq);
		} else {
			backSeq = seq.mateSeqBase_;
		}
		for (const auto & mipName : mipNamesForFamily_.at(mip.familyName_)) {
			const auto & otherMip = mips_.at(mipName);
			alignerObjForGroupDet.alignCacheLocal(backSeq, otherMip.ligationArmObj_);
			//normalize score to arm length
			double currentScore = alignerObjForGroupDet.parts_.score_
					/ static_cast<double>(len(otherMip.ligationArmObj_));
			auto armPosMotifCurrent = otherMip.getPossibleExtArmPos(seq.seqBase_);
			auto ligArmPosMotifCurrent = otherMip.getPossibleLigArmPos(seq.mateSeqBase_);
			//ensure that even if the ligation arm is better that the extension arm still fits
			//this could happen if there was a chimeric event that switched arms
			if (!armPosMotifCurrent.empty()
					&& !ligArmPosMotifCurrent.empty()) {
				if (currentScore > bestScore) {
					bestScore = currentScore;
					mip = mips_.at(mipName);
				}
			}
		}
	}
	return mip;
}

bool MipCollection::hasMipFamily(const std::string & mipFam) const {
	return bib::in(mipFam, mipFamilies_);
}

bool MipCollection::hasMipTarget(const std::string & mipTar) const {
	return bib::in(mipTar, mips_);
}

}  // namespace bibseq

