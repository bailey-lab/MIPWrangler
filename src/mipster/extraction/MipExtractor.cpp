/*
 * MipExtractor.cpp
 *
 *  Created on: Sep 9, 2017
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

#include "MipExtractor.hpp"
#include "mipster/info.h"
#include "mipster/mipUtils/MipNameSorter.hpp"

#include <SeekDeep/objects/PairedReadProcessor.hpp>

namespace njhseq {


MipExtractor::MipExtractor(){};
MipExtractor::MipExtractor(bool verbose): verbose_(verbose){};


void MipExtractor::extractFilterSampleForMipsPairedStitch(const std::vector<SeqIOOptions> & sampleIOOpts,
		const SetUpMaster & mipMaster,
		aligner & alignerObjForFamilyDet,
		aligner & alignerObjForStitching,
		const mipIllumArmExtractionPars & pars){

	njh::stopWatch watch;
	//std::cout << "On Thread: " << std::this_thread::get_id() << std::endl;
	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
	sampDirMaster.createExtractDirectory(pars.overWriteDirs);
	alignerObjForStitching.numberOfAlingmentsDone_ = 0;
	alignerObjForFamilyDet.resetAlnCache();
	alignerObjForFamilyDet.processAlnInfoInputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
	//set up sub directories
	bfs::path filteredOffDir = njh::files::makeDir(sampDirMaster.extractDir_.string(), njh::files::MkdirPar("filteredOff"));

	//create read stitcher

	PairedReadProcessor pProcessor(pars.processPairPars_);

	//create out files
	MultiSeqOutCache<PairedRead> mipOuts;
	mipOuts.setOpenLimit(pars.fileOpenLimit_/2);
	MultiSeqOutCache<seqInfo> mipStitchedOuts;
	mipStitchedOuts.setOpenLimit(pars.fileOpenLimit_/2);

	mipOuts.addReader("indeterminate",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "indeterminate").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("unmatchedReads",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "unmatchedReads").string(), sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("smallFragment",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "smallFragment").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));

	//extraction counts
	MipExtractionStats allExtractStats(pars.sampleName);
	allExtractStats.minCaptureLength_ = pars.minCaptureLength;
	if (pars.qFilPars_.checkingQFrac_) {
		allExtractStats.qualCheckStr = "failed_q"
				+ estd::to_string(pars.qFilPars_.qualCheck_) + "<"
				+ estd::to_string(pars.qFilPars_.qualCheckCutOff_);
	} else if (pars.qFilPars_.checkingQWindow) {
		allExtractStats.qualCheckStr = "failed_qW" + pars.qFilPars_.qualWindow_;
	} else {
		allExtractStats.qualCheckStr = "failed_quality(noneDone)";
	}

	std::unordered_map<std::string, PairedReadProcessor::ProcessedResultsCounts> pairStitchingCounts;

//	VecStr filterOutNames = { "_failedQuality", "_failedLigation", "_failedMinLen", "_containsNs", "_badStitch" };
	auto forStitchedOut = sampleIOOpts.front().out_;
	forStitchedOut.outExtention_ = ".fastq.gz";

	VecStr allMipTargets = mipMaster.getAllMipTargets();
	for (const auto & mip : allMipTargets) {
		bfs::path mipDirectory = njh::files::makeDir(sampDirMaster.extractDir_.string(), njh::files::MkdirPar(mip));

		mipOuts.addReader(mip,
				SeqIOOptions(
						njh::files::make_path(sampDirMaster.extractDir_.string(), mip, mip),
						sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
		mipOuts.addReader(mip + "_filteredOff",
				SeqIOOptions(
						njh::files::make_path(mipDirectory, mip
								+ "_filteredOff" ), sampleIOOpts.front().outFormat_,
										sampleIOOpts.front().out_) );

		if(pars.writeOutInitialExtractedPairs){
			mipOuts.addReader(mip + "_initial",
					SeqIOOptions(
							njh::files::make_path(mipDirectory, mip
									+ "_initial" ), sampleIOOpts.front().outFormat_,
											sampleIOOpts.front().out_) );
		}

		mipStitchedOuts.addReader(mip,
				SeqIOOptions(
						njh::files::make_path(sampDirMaster.extractDir_.string(), mip, mip),
						SeqIOOptions::outFormats::FASTQGZ, forStitchedOut));

		mipStitchedOuts.addReader(mip + "_filteredOff",
				SeqIOOptions(
						njh::files::make_path(mipDirectory, mip
								+ "_filteredOff" ),
						SeqIOOptions::outFormats::FASTQGZ, forStitchedOut) );

		pairStitchingCounts[mip] = PairedReadProcessor::ProcessedResultsCounts{};

	}






	for(const auto & sampleIOOpt : sampleIOOpts){
		//read in reads
		SeqIO readerOpt(sampleIOOpt);
		readerOpt.openIn();
		if (verbose_) {
			std::cout << sampleIOOpt.firstName_ << std::endl;
		}
		PairedRead seq;
		uint32_t readCount = 1;
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

		while (readerOpt.readNextRead(seq)) {
			if (readCount % 100 == 0 && verbose_) {
				std::cout << "\r" << "currently on " << readCount;
				std::cout.flush();
			}
			++readCount;
			//get length and resize aligner vector if needed
			uint64_t maxLen = alignerObjForFamilyDet.parts_.maxSize_;
			readVec::getMaxLength(seq, maxLen);
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

			if (maxLen > alignerObjForFamilyDet.parts_.maxSize_) {
				alignerObjForFamilyDet.parts_.setMaxSize(maxLen);
			}
			if (maxLen > alignerObjForStitching.parts_.maxSize_) {
				alignerObjForStitching.parts_.setMaxSize(maxLen);
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			bool found = false;
			std::unordered_map<std::string, std::pair<std::vector<Mip::ArmPosScore>, std::vector<Mip::ArmPosScore>>> possibleArms;
			std::unordered_map<std::string, std::vector<Mip::ArmPosScore>> possibleExtArms;
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

			for (const auto & mKey : allMipTargets) {
				const auto & mip = mipMaster.mips_->mips_.at(mKey);
				//check arm
				auto armPosMotif = mip.getPossibleExtArmPos(seq.seqBase_);
				//if found arm and at the right position continue on;
				if (!armPosMotif.empty()) {
					found = true;
					possibleExtArms.emplace(mip.name_, armPosMotif);
					auto ligArmPosMotif = mip.getPossibleLigArmPos(seq.mateSeqBase_);
					if(!ligArmPosMotif.empty()){
						possibleArms.emplace(mip.name_, std::make_pair(armPosMotif, ligArmPosMotif));
					}
				}
			}
			//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
			if(possibleArms.empty()){
				if(possibleExtArms.empty()){
					//no matches found
					allExtractStats.increaseUnmatched();
					mipOuts.add("unmatchedReads", seq);
				}else if(possibleExtArms.size() == 1){
					const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					//stitching
					++pairStitchingCounts[mip.name_].total;
					auto stitchedRes = pProcessor.processPairedEnd(seq,pairStitchingCounts[mip.name_], alignerObjForStitching);
					//for now just accepting r1 ends in r2 (no overlaps or perfect overlaps)
					SinlgeMipExtractInfo::extractCase eCase{SinlgeMipExtractInfo::extractCase::NONE};
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					if(stitchedRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2){
//					if (stitchedRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NONE ||
//							 stitchedRes.status_== PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP) {
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
						eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						//log and write read
						MetaDataInName failedStitchingMeta;
						failedStitchingMeta.addMeta("failedStitchCase", PairedReadProcessor::getOverlapStatusStr(stitchedRes.status_));
						seq.seqBase_.name_.append(failedStitchingMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedStitchingMeta.createMetaName());
						mipOuts.add(mip.name_ + "_filteredOff", seq);
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					} else {
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
						//quality control
						std::string failedQaulifierName;
						if(mip.getPossibleExtArmPos(*stitchedRes.combinedSeq_).empty()){
							MetaDataInName failedQcMeta;
							failedQcMeta.addMeta("failed", "ArmsNoLongerMatch");
							stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
							failedQaulifierName = "_filteredOff";
							eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						} else {
							eCase = mip.checkRead(*stitchedRes.combinedSeq_, pars.qFilPars_);
							failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
							if("" != failedQaulifierName){
								MetaDataInName failedQcMeta;
								failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
								stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
								failedQaulifierName = "_filteredOff";
							} else {
								//stitchedRes.combinedSeq_->outPutSeqAnsi(std::cout);
								BarcodeInfo barInfo = mip.determineExtBarcode(*stitchedRes.combinedSeq_);
								mip.determineLigBarcode(*stitchedRes.combinedSeq_, barInfo);
								MetaDataInName barMeta;
								barMeta.addMeta("extBar", barInfo.extBar_);
								if ("" != barInfo.ligBar_) {
									barMeta.addMeta("ligBar", barInfo.ligBar_);
								}
								barMeta.addMeta("fullBar_", barInfo.fullBar_);
								stitchedRes.combinedSeq_->name_.append(barMeta.createMetaName());
							}
						}
						//log and write read
						mipStitchedOuts.add(mip.name_ + failedQaulifierName, *stitchedRes.combinedSeq_);
					}
					allExtractStats.increaseCount(mip.name_, eCase);
				}else{
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					const auto & initialMip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
					auto currentMip = mipMaster.mips_->determineBestMipInFamily(seq, initialMip, alignerObjForFamilyDet);
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					//stitching
					++pairStitchingCounts[currentMip.name_].total;
					auto stitchedRes = pProcessor.processPairedEnd(seq,pairStitchingCounts[currentMip.name_], alignerObjForStitching);
					//for now just accepting r1 ends in r2 (no overlaps or perfect overlaps)
					SinlgeMipExtractInfo::extractCase eCase{SinlgeMipExtractInfo::extractCase::NONE};
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					if(stitchedRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2){
//					if (stitchedRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NONE ||
//							 stitchedRes.status_== PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP) {
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
						eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						//log and write read
						MetaDataInName failedStitchingMeta;
						failedStitchingMeta.addMeta("failedStitchCase", PairedReadProcessor::getOverlapStatusStr(stitchedRes.status_));
						seq.seqBase_.name_.append(failedStitchingMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedStitchingMeta.createMetaName());
						mipOuts.add(currentMip.name_ + "_filteredOff", seq);
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					}else{
						//quality control
						std::string failedQaulifierName;
						if(currentMip.getPossibleExtArmPos(*stitchedRes.combinedSeq_).empty()){
							MetaDataInName failedQcMeta;
							failedQcMeta.addMeta("failed", "ArmsNoLongerMatch");
							stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
							failedQaulifierName = "_filteredOff";
							eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						} else {
							eCase = currentMip.checkRead(*stitchedRes.combinedSeq_, pars.qFilPars_);
							failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
							if("" != failedQaulifierName){
								MetaDataInName failedQcMeta;
								failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
								stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
								failedQaulifierName = "_filteredOff";
							} else {
								//stitchedRes.combinedSeq_->outPutSeqAnsi(std::cout);
								BarcodeInfo barInfo = currentMip.determineExtBarcode(*stitchedRes.combinedSeq_);
								currentMip.determineLigBarcode(*stitchedRes.combinedSeq_, barInfo);
								MetaDataInName barMeta;
								barMeta.addMeta("extBar", barInfo.extBar_);
								if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
								barMeta.addMeta("fullBar_", barInfo.fullBar_);
								stitchedRes.combinedSeq_->name_.append(barMeta.createMetaName());
							}
						}

						//log and write read
						mipStitchedOuts.add(currentMip.name_ + failedQaulifierName, *stitchedRes.combinedSeq_);
					}
					allExtractStats.increaseCount(currentMip.name_, eCase);
				}
			} else if (1 == possibleArms.size()) {
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				//only one possible match
				const auto & mip = mipMaster.mips_->mips_.at(possibleArms.begin()->first);
				//stitching
				auto stitchedRes = pProcessor.processPairedEnd(seq,pairStitchingCounts[mip.name_], alignerObjForStitching);
				++pairStitchingCounts[mip.name_].total;
				if(pars.writeOutInitialExtractedPairs){
					mipOuts.add(mip.name_ + "_initial", seq);
				}
				//for now just accepting r1 ends in r2 (no overlaps or perfect overlaps)
				SinlgeMipExtractInfo::extractCase eCase{SinlgeMipExtractInfo::extractCase::NONE};
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				if(stitchedRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2){
//				if (stitchedRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NONE ||
//						 stitchedRes.status_== PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP) {


					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
					//log and write read
					MetaDataInName failedStitchingMeta;
					failedStitchingMeta.addMeta("failedStitchCase", PairedReadProcessor::getOverlapStatusStr(stitchedRes.status_));
					seq.seqBase_.name_.append(failedStitchingMeta.createMetaName());
					seq.mateSeqBase_.name_.append(failedStitchingMeta.createMetaName());
					mipOuts.add(mip.name_ + "_filteredOff", seq);
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				}else{
					//quality control
					std::string failedQaulifierName;
					if(mip.getPossibleExtArmPos(*stitchedRes.combinedSeq_).empty()){
						MetaDataInName failedQcMeta;
						failedQcMeta.addMeta("failed", "ArmsNoLongerMatch");
						stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
						failedQaulifierName = "_filteredOff";
						eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
					} else {
						eCase = mip.checkRead(*stitchedRes.combinedSeq_, pars.qFilPars_);
						failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
						if("" != failedQaulifierName){
							MetaDataInName failedQcMeta;
							failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
							stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
							failedQaulifierName = "_filteredOff";
						} else {
							//stitchedRes.combinedSeq_->outPutSeqAnsi(std::cout);
							BarcodeInfo barInfo = mip.determineExtBarcode(*stitchedRes.combinedSeq_);
							mip.determineLigBarcode(*stitchedRes.combinedSeq_, barInfo);
							MetaDataInName barMeta;
							barMeta.addMeta("extBar", barInfo.extBar_);
							if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
							barMeta.addMeta("fullBar_", barInfo.fullBar_);
							stitchedRes.combinedSeq_->name_.append(barMeta.createMetaName());
						}
					}
					//log and write read
					mipStitchedOuts.add(mip.name_ + failedQaulifierName, *stitchedRes.combinedSeq_);
				}
				allExtractStats.increaseCount(mip.name_, eCase);
			} else {
				//multiple possible matches
				uint32_t bestScore = 0;
				std::vector<std::string> bestMips;
				for(const auto & possible : possibleArms){
					uint32_t bestExtScore = 0;
					uint32_t bestLigScore = 0;
					for(const auto & ext : possible.second.first){
						if(ext.score_ > bestExtScore){
							bestExtScore = ext.score_;
						}
					}
					for(const auto & lig : possible.second.second){
						if(lig.score_ > bestLigScore){
							bestLigScore = lig.score_;
						}
					}
					uint32_t currentScore = bestLigScore + bestExtScore;
					if (currentScore > bestScore) {
						bestScore = currentScore;
						bestMips.clear();
						bestMips.emplace_back(possible.first);
					} else if (0 != currentScore && currentScore == bestScore) {
						bestMips.emplace_back(possible.first);
					}
				}
				//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				if(bestMips.size() == 1){
					const auto & mip = mipMaster.mips_->mips_.at(bestMips.front());
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					//stitching
					++pairStitchingCounts[mip.name_].total;
					auto stitchedRes = pProcessor.processPairedEnd(seq,pairStitchingCounts[mip.name_], alignerObjForStitching);
					//for now just accepting r1 ends in r2 (no overlaps or perfect overlaps)
					SinlgeMipExtractInfo::extractCase eCase{SinlgeMipExtractInfo::extractCase::NONE};
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					if(stitchedRes.status_ != PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2){
//					if(stitchedRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NONE ||
//						 stitchedRes.status_ == PairedReadProcessor::ReadPairOverLapStatus::NOOVERLAP){
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
						eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						//log and write read
						MetaDataInName failedStitchingMeta;
						failedStitchingMeta.addMeta("failedStitchCase", PairedReadProcessor::getOverlapStatusStr(stitchedRes.status_));
						seq.seqBase_.name_.append(failedStitchingMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedStitchingMeta.createMetaName());
						mipOuts.add(mip.name_ + "_filteredOff", seq);
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					}else{
						//quality control
						std::string failedQaulifierName;
						if(mip.getPossibleExtArmPos(*stitchedRes.combinedSeq_).empty()){
							MetaDataInName failedQcMeta;
							failedQcMeta.addMeta("failed", "ArmsNoLongerMatch");
							stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
							failedQaulifierName = "_filteredOff";
							eCase = SinlgeMipExtractInfo::extractCase::BADSTITCH;
						} else {
							eCase = mip.checkRead(*stitchedRes.combinedSeq_, pars.qFilPars_);
							failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
							if("" != failedQaulifierName){
								MetaDataInName failedQcMeta;
								failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
								stitchedRes.combinedSeq_->name_.append(failedQcMeta.createMetaName());
								failedQaulifierName = "_filteredOff";
							} else {
								//stitchedRes.combinedSeq_->outPutSeqAnsi(std::cout);
								BarcodeInfo barInfo = mip.determineExtBarcode(*stitchedRes.combinedSeq_);
								mip.determineLigBarcode(*stitchedRes.combinedSeq_, barInfo);
								MetaDataInName barMeta;
								barMeta.addMeta("extBar", barInfo.extBar_);
								if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
								barMeta.addMeta("fullBar_", barInfo.fullBar_);
								stitchedRes.combinedSeq_->name_.append(barMeta.createMetaName());
							}
						}
						//log and write read
						mipStitchedOuts.add(mip.name_ + failedQaulifierName, *stitchedRes.combinedSeq_);
						//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					}
					allExtractStats.increaseCount(mip.name_, eCase);
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
				}else if(bestMips.size() > 1){
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					//too many matches found
					allExtractStats.increaseIndeterminate();
					MetaDataInName indeterminateMeta;
					uint32_t bestCount = 0;
					for(const auto & best : bestMips){
						const auto & possible = possibleArms.at(best);
						uint32_t bestExtScore = 0;
						uint32_t bestLigScore = 0;

						for(const auto & ext : possible.first){
							if(ext.score_ > bestExtScore){
								bestExtScore = ext.score_;
							}
						}
						for(const auto & lig : possible.second){
							if(lig.score_ > bestLigScore){
								bestLigScore = lig.score_;
							}
						}
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_name",     best);
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_extScore", bestExtScore);
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_ligScore", bestLigScore);
						++bestCount;
					}
					//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
					seq.seqBase_.name_.append(indeterminateMeta.createMetaName() );
					seq.mateSeqBase_.name_.append(indeterminateMeta.createMetaName() );
					mipOuts.add("indeterminate", seq);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << std::endl;
					ss << "Best Mips vector is empty, which shouldn't be able to happen" << std::endl;
					ss << "For : " << sampleIOOpt.firstName_ << std::endl;
					throw std::runtime_error{ss.str()};
				}
			}
		}
		//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
		if (verbose_) {
			std::cout << std::endl;
		}
	}

	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	mipOuts.closeOutAll();
	mipStitchedOuts.closeOutAll();


	table extractionInfoTabByTarget = allExtractStats.outputContentsJustTargets(*mipMaster.mips_);
	extractionInfoTabByTarget.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoByTarget.txt", ".txt"), "\t", extractionInfoTabByTarget.hasHeader_));

	table extractionInfoTabSummary = allExtractStats.outputContentsSummary(*mipMaster.mips_);
	extractionInfoTabSummary.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoSummary.txt", ".txt"), "\t", extractionInfoTabSummary.hasHeader_));



	VecStr stitchResultsColNames{"Sample", "mipTarget", "mipFamily", "total", "r1EndsInR2", "r1BeginsInR2", "OverlapFail", "PerfectOverlap"};
	table stitchInfoByTarget(stitchResultsColNames);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	auto mipNameKeys = getVectorOfMapKeys(pairStitchingCounts);
	MipNameSorter::sort(mipNameKeys);
	for(const auto & mipKey : mipNameKeys){
		if(0 == pairStitchingCounts[mipKey].r1EndsInR2Combined && !pars.keepIntermediateFiles){
			auto extractionDirForMipKey = njh::files::make_path(sampDirMaster.extractDir_, mipKey);
			if(bfs::exists(extractionDirForMipKey)){
				njh::files::rmDirForce(extractionDirForMipKey);
			}
		}
		stitchInfoByTarget.addRow(pars.sampleName,
				mipKey,
				mipMaster.mips_->getFamilyForTarget(mipKey),
				pairStitchingCounts[mipKey].total,
				getPercentageString(pairStitchingCounts[mipKey].r1EndsInR2Combined,pairStitchingCounts[mipKey].total),
				getPercentageString(pairStitchingCounts[mipKey].r1BeginsInR2Combined,pairStitchingCounts[mipKey].total),
				getPercentageString(pairStitchingCounts[mipKey].overlapFail,pairStitchingCounts[mipKey].total),
				getPercentageString(pairStitchingCounts[mipKey].perfectOverlapCombined,pairStitchingCounts[mipKey].total));
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	stitchInfoByTarget.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "stitchInfoByTarget.txt", ".txt"), "\t", stitchInfoByTarget.hasHeader_));

	std::ofstream sampLog;
	openTextFile(sampLog, sampDirMaster.extractDir_.string() + "log.txt", ".txt", false, true);
	sampLog << "Ran on: " << njh::getCurrentDate() << std::endl;
	sampLog << "Number of Alignments Done: "
			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	sampLog << "Number of Alignments Done For Stitching: "
			<< alignerObjForStitching.numberOfAlingmentsDone_ << "\n";
	sampLog << "Run Time: " << watch.timeLapFormatted(6) << std::endl;
	if(pars.cacheAlignments){
		alignerObjForFamilyDet.processAlnInfoOutputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
	}
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
}


void MipExtractor::extractFilterSampleForMipsPaired(const std::vector<SeqIOOptions> & sampleIOOpts,
		const SetUpMaster & mipMaster,
		aligner & alignerObjForFamilyDet,
		const mipIllumArmExtractionPars & pars) {

	njh::stopWatch watch;
	//std::cout << "On Thread: " << std::this_thread::get_id() << std::endl;
	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
	sampDirMaster.createExtractDirectory(pars.overWriteDirs);
	alignerObjForFamilyDet.resetAlnCache();
	alignerObjForFamilyDet.processAlnInfoInputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
	//set up sub directories
	bfs::path filteredOffDir = njh::files::makeDir(sampDirMaster.extractDir_.string(), njh::files::MkdirPar("filteredOff"));

	//create out files
	MultiSeqOutCache<PairedRead> mipOuts;
	mipOuts.setOpenLimit(pars.fileOpenLimit_);

	mipOuts.addReader("indeterminate",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "indeterminate").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("unmatchedReads",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "unmatchedReads").string(), sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("smallFragment",
			SeqIOOptions(njh::files::make_path(filteredOffDir, "smallFragment").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	VecStr filterOutNames = { "_failedQuality", "_failedLigation", "_failedMinLen", "_containsNs" };

	VecStr allMipTargets = mipMaster.getAllMipTargets();
	for (const auto & mip : allMipTargets) {
		mipOuts.addReader(mip,
				SeqIOOptions(
						njh::files::make_path(sampDirMaster.extractDir_.string(), mip),
						sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
		mipOuts.addReader(mip + "_filteredOff",
				SeqIOOptions(
						njh::files::make_path(filteredOffDir, mip
								+ "_filteredOff" ), sampleIOOpts.front().outFormat_,
										sampleIOOpts.front().out_) );
	}

	MipExtractionStats allExtractStats(pars.sampleName);
	allExtractStats.minCaptureLength_ = pars.minCaptureLength;
	if (pars.qFilPars_.checkingQFrac_) {
		allExtractStats.qualCheckStr = "failed_q"
				+ estd::to_string(pars.qFilPars_.qualCheck_) + "<"
				+ estd::to_string(pars.qFilPars_.qualCheckCutOff_);
	} else if (pars.qFilPars_.checkingQWindow) {
		allExtractStats.qualCheckStr = "failed_qW" + pars.qFilPars_.qualWindow_;
	} else {
		allExtractStats.qualCheckStr = "failed_quality(noneDone)";
	}

	for(const auto & sampleIOOpt : sampleIOOpts){
		//read in reads
		SeqIO readerOpt(sampleIOOpt);
		readerOpt.openIn();
		if (verbose_) {
			std::cout << sampleIOOpt.firstName_ << std::endl;
		}
		PairedRead seq;
		uint32_t readCount = 1;
		while (readerOpt.readNextRead(seq)) {
			if (readCount % 100 == 0 && verbose_) {
				std::cout << "\r" << "currently on " << readCount;
				std::cout.flush();
			}
			++readCount;
			//get length and resize aligner vector if needed
			uint64_t maxLen = alignerObjForFamilyDet.parts_.maxSize_;
			readVec::getMaxLength(seq, maxLen);
			if (maxLen > alignerObjForFamilyDet.parts_.maxSize_) {
				alignerObjForFamilyDet.parts_.setMaxSize(maxLen);
			}
			bool found = false;
			std::unordered_map<std::string, std::pair<std::vector<Mip::ArmPosScore>, std::vector<Mip::ArmPosScore>>> possibleArms;
			std::unordered_map<std::string, std::vector<Mip::ArmPosScore>> possibleExtArms;
			for (const auto & mKey : allMipTargets) {
				const auto & mip = mipMaster.mips_->mips_.at(mKey);
				//check arm
				auto armPosMotif = mip.getPossibleExtArmPos(seq.seqBase_);
				//if found arm and at the right position continue on;
				if (!armPosMotif.empty()) {
					found = true;
					possibleExtArms.emplace(mip.name_, armPosMotif);
					auto ligArmPosMotif = mip.getPossibleLigArmPos(seq.mateSeqBase_);
					if(!ligArmPosMotif.empty()){
						possibleArms.emplace(mip.name_, std::make_pair(armPosMotif, ligArmPosMotif));
					}
				}
			}
			if(possibleArms.empty()){
				if(possibleExtArms.empty()){
					//no matches found
					allExtractStats.increaseUnmatched();
					mipOuts.add("unmatchedReads", seq);
				}else if(possibleExtArms.size() == 1){
					const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
					//quality control
					SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, pars.qFilPars_);
					std::string failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
					if("" != failedQaulifierName){
						MetaDataInName failedQcMeta;
						failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
						seq.seqBase_.name_.append(failedQcMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedQcMeta.createMetaName());
						failedQaulifierName = "_filteredOff";
					}else{
						BarcodeInfo barInfo = mip.determineExtBarcode(seq.seqBase_);
						mip.determineLigBarcode(seq.mateSeqBase_, barInfo);
						MetaDataInName barMeta;
						barMeta.addMeta("extBar", barInfo.extBar_);
						if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
						barMeta.addMeta("fullBar_", barInfo.fullBar_);
						seq.seqBase_.name_.append(barMeta.createMetaName());
						seq.mateSeqBase_.name_.append(barMeta.createMetaName());
					}
					//log and write read
					mipOuts.add(mip.name_ + failedQaulifierName, seq);
					allExtractStats.increaseCount(mip.name_, eCase);
				}else{
					const auto & mip = mipMaster.mips_->mips_.at(possibleExtArms.begin()->first);
					auto currentMip = mipMaster.mips_->determineBestMipInFamily(seq, mip, alignerObjForFamilyDet);
					//quality control
					SinlgeMipExtractInfo::extractCase eCase = currentMip.checkRead(seq, pars.qFilPars_);

					std::string failedQaulifierName = MipExtractionStats::getNameForCase(
							eCase);
					if("" != failedQaulifierName){
						MetaDataInName failedQcMeta;
						failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
						seq.seqBase_.name_.append(failedQcMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedQcMeta.createMetaName());
						failedQaulifierName = "_filteredOff";
					}else{
						BarcodeInfo barInfo = mip.determineExtBarcode(seq.seqBase_);
						mip.determineLigBarcode(seq.mateSeqBase_, barInfo);
						MetaDataInName barMeta;
						barMeta.addMeta("extBar", barInfo.extBar_);
						if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
						barMeta.addMeta("fullBar", barInfo.fullBar_);
						seq.seqBase_.name_.append(barMeta.createMetaName());
						seq.mateSeqBase_.name_.append(barMeta.createMetaName());
					}
					//log and write read
					mipOuts.add(currentMip.name_ + failedQaulifierName, seq);
					allExtractStats.increaseCount(currentMip.name_, eCase);
				}
			} else if (1 == possibleArms.size()) {
				//only one possible match
				const auto & mip = mipMaster.mips_->mips_.at(possibleArms.begin()->first);


				//quality control
				SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, pars.qFilPars_);

				std::string failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
				if("" != failedQaulifierName){
					MetaDataInName failedQcMeta;
					failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
					seq.seqBase_.name_.append(failedQcMeta.createMetaName());
					seq.mateSeqBase_.name_.append(failedQcMeta.createMetaName());
					failedQaulifierName = "_filteredOff";
				} else {
					BarcodeInfo barInfo = mip.determineExtBarcode(seq.seqBase_);
					mip.determineLigBarcode(seq.mateSeqBase_, barInfo);
					MetaDataInName barMeta;
					barMeta.addMeta("extBar", barInfo.extBar_);
					if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
					barMeta.addMeta("fullBar_", barInfo.fullBar_);
					seq.seqBase_.name_.append(barMeta.createMetaName());
					seq.mateSeqBase_.name_.append(barMeta.createMetaName());
				}
				//log and write read
				mipOuts.add(mip.name_ + failedQaulifierName, seq);
				allExtractStats.increaseCount(mip.name_, eCase);
			} else {
				//multiple possible matches
				uint32_t bestScore = 0;
				std::vector<std::string> bestMips;
				for(const auto & possible : possibleArms){
					uint32_t bestExtScore = 0;
					uint32_t bestLigScore = 0;
					for(const auto & ext : possible.second.first){
						if(ext.score_ > bestExtScore){
							bestExtScore = ext.score_;
						}
					}
					for(const auto & lig : possible.second.second){
						if(lig.score_ > bestLigScore){
							bestLigScore = lig.score_;
						}
					}
					uint32_t currentScore = bestLigScore + bestExtScore;
					if (currentScore > bestScore) {
						bestScore = currentScore;
						bestMips.clear();
						bestMips.emplace_back(possible.first);
					} else if (0 != currentScore && currentScore == bestScore) {
						bestMips.emplace_back(possible.first);
					}
				}
				if(bestMips.size() == 1){
					const auto & mip = mipMaster.mips_->mips_.at(bestMips.front());

					//quality control
					SinlgeMipExtractInfo::extractCase eCase = mip.checkRead(seq, pars.qFilPars_);

					std::string failedQaulifierName = MipExtractionStats::getNameForCase(eCase);
					if("" != failedQaulifierName){
						MetaDataInName failedQcMeta;
						failedQcMeta.addMeta("failed", '_' == failedQaulifierName.front()? failedQaulifierName.substr(1) : failedQaulifierName);
						seq.seqBase_.name_.append(failedQcMeta.createMetaName());
						seq.mateSeqBase_.name_.append(failedQcMeta.createMetaName());
						failedQaulifierName = "_filteredOff";
					}else{
						BarcodeInfo barInfo = mip.determineExtBarcode(seq.seqBase_);
						mip.determineLigBarcode(seq.mateSeqBase_, barInfo);
						MetaDataInName barMeta;
						barMeta.addMeta("extBar", barInfo.extBar_);
						if ("" != barInfo.ligBar_) barMeta.addMeta("ligBar", barInfo.ligBar_);
						barMeta.addMeta("fullBar_", barInfo.fullBar_);
						seq.seqBase_.name_.append(barMeta.createMetaName());
						seq.mateSeqBase_.name_.append(barMeta.createMetaName());
					}
					//log and write read
					mipOuts.add(mip.name_ + failedQaulifierName, seq);
					allExtractStats.increaseCount(mip.name_, eCase);
				}else if(bestMips.size() > 1){
					//too many matches found
					allExtractStats.increaseIndeterminate();
					MetaDataInName indeterminateMeta;
					uint32_t bestCount = 0;
					for(const auto & best : bestMips){
						const auto & possible = possibleArms.at(best);
						uint32_t bestExtScore = 0;
						uint32_t bestLigScore = 0;

						for(const auto & ext : possible.first){
							if(ext.score_ > bestExtScore){
								bestExtScore = ext.score_;
							}
						}
						for(const auto & lig : possible.second){
							if(lig.score_ > bestLigScore){
								bestLigScore = lig.score_;
							}
						}
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_name",     best);
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_extScore", bestExtScore);
						indeterminateMeta.addMeta("match_" + estd::to_string(bestCount)+ "_ligScore", bestLigScore);
						++bestCount;
					}
					seq.seqBase_.name_.append(indeterminateMeta.createMetaName() );
					seq.mateSeqBase_.name_.append(indeterminateMeta.createMetaName() );
					mipOuts.add("indeterminate", seq);
				}else{
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << std::endl;
					ss << "Best Mips vector is empty, which shouldn't be able to happen" << std::endl;
					ss << "For : " << sampleIOOpt.firstName_ << std::endl;
					throw std::runtime_error{ss.str()};
				}
			}
		}
		if (verbose_) {
			std::cout << std::endl;
		}
	}


	mipOuts.closeOutAll();


	table extractionInfoTabByTarget = allExtractStats.outputContentsJustTargets(*mipMaster.mips_);
	extractionInfoTabByTarget.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoByTarget.txt", ".txt"), "\t", extractionInfoTabByTarget.hasHeader_));

	table extractionInfoTabSummary = allExtractStats.outputContentsSummary(*mipMaster.mips_);
	extractionInfoTabSummary.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoSummary.txt", ".txt"), "\t", extractionInfoTabSummary.hasHeader_));

	std::ofstream sampLog;
	openTextFile(sampLog, sampDirMaster.extractDir_.string() + "log.txt", ".txt", false, true);
	sampLog << "Ran on: " << njh::getCurrentDate() << std::endl;
	sampLog << "Number of Alignments Done: "
			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	sampLog << "Run Time: " << watch.timeLapFormatted(6) << std::endl;
	if(pars.cacheAlignments){
		alignerObjForFamilyDet.processAlnInfoOutputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
	}
}


} /* namespace njhseq */
