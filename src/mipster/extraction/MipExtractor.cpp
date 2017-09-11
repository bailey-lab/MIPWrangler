/*
 * MipExtractor.cpp
 *
 *  Created on: Sep 9, 2017
 *      Author: nick
 */

#include "MipExtractor.hpp"
#include "mipster/info.h"

namespace bibseq {


MipExtractor::MipExtractor(){};
MipExtractor::MipExtractor(bool verbose): verbose_(verbose){};



void MipExtractor::extractFilterSampleForMipsPaired(const std::vector<SeqIOOptions> & sampleIOOpts,
		const SetUpMaster & mipMaster,
		aligner & alignerObjForFamilyDet,
		const mipIllumArmExtractionPars & pars) {

	bib::stopWatch watch;
	//std::cout << "On Thread: " << std::this_thread::get_id() << std::endl;
	SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", pars.sampleName));
	sampDirMaster.createExtractDirectory(pars.overWriteDirs);
	alignerObjForFamilyDet.resetAlnCache();
	alignerObjForFamilyDet.processAlnInfoInputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
	//set up sub directories
	bfs::path filteredOffDir = bib::files::makeDir(sampDirMaster.extractDir_.string(), bib::files::MkdirPar("filteredOff"));
	//create out files
	MultiSeqOutCache<PairedRead> mipOuts;
	mipOuts.setOpenLimit(pars.fileOpenLimit_);

	mipOuts.addReader("indeterminate",
			SeqIOOptions(bib::files::make_path(filteredOffDir, "indeterminate").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("unmatchedReads",
			SeqIOOptions(bib::files::make_path(filteredOffDir, "unmatchedReads").string(), sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	mipOuts.addReader("smallFragment",
			SeqIOOptions(bib::files::make_path(filteredOffDir, "smallFragment").string(),  sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
	VecStr filterOutNames = { "_failedQuality", "_failedLigation", "_failedMinLen", "_containsNs" };

	VecStr allMipTargets = mipMaster.getAllMipTargets();
	for (const auto & mip : allMipTargets) {
		mipOuts.addReader(mip,
				SeqIOOptions(
						bib::files::make_path(sampDirMaster.extractDir_.string(), mip),
						sampleIOOpts.front().outFormat_, sampleIOOpts.front().out_));
		mipOuts.addReader(mip + "_filteredOff",
				SeqIOOptions(
						bib::files::make_path(filteredOffDir, mip
								+ "_filteredOff" ), sampleIOOpts.front().outFormat_,
										sampleIOOpts.front().out_) );
	}

	MipExtractionStats allExtractStats(pars.sampleName);
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
						barMeta.addMeta("ligBar", barInfo.ligBar_);
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
						barMeta.addMeta("ligBar", barInfo.ligBar_);
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
				}else{
					BarcodeInfo barInfo = mip.determineExtBarcode(seq.seqBase_);
					mip.determineLigBarcode(seq.mateSeqBase_, barInfo);
					MetaDataInName barMeta;
					barMeta.addMeta("extBar", barInfo.extBar_);
					barMeta.addMeta("ligBar", barInfo.ligBar_);
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
						barMeta.addMeta("ligBar", barInfo.ligBar_);
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
	VecStr extracOnlyColNames {"sampleName", "mipTarget", "mipFamily", "readNumber",
			"goodReads", "failedLigationArm", "failedMinLen(<"
					+ estd::to_string(pars.minLen) + ")" };
	if (pars.qFilPars_.checkingQFrac_) {
		extracOnlyColNames.emplace_back(
				"failed_q" + estd::to_string(pars.qFilPars_.qualCheck_) + "<"
						+ estd::to_string(pars.qFilPars_.qualCheckCutOff_));
	} else if (pars.qFilPars_.checkingQWindow) {
		extracOnlyColNames.emplace_back("failed_qW" + pars.qFilPars_.qualWindow_);
	} else {
		extracOnlyColNames.emplace_back("failed_quality(noneDone)");
	}
	extracOnlyColNames.emplace_back("containsNs");
	table infoTabByTarget(extracOnlyColNames);
	infoTabByTarget.content_ = allExtractStats.outputContents(*mipMaster.mips_, "\t");
	infoTabByTarget.outPutContents(
			TableIOOpts(OutOptions(sampDirMaster.extractDir_.string() + "extractInfoByTarget.txt", ".txt"), "\t", infoTabByTarget.hasHeader_));
	std::ofstream sampLog;
	openTextFile(sampLog, sampDirMaster.extractDir_.string() + "log.txt", ".txt", false, true);
	sampLog << "Ran on: " << bib::getCurrentDate() << std::endl;
	sampLog << "Number of Alignments Done: "
			<< alignerObjForFamilyDet.numberOfAlingmentsDone_ << "\n";
	sampLog << "Run Time: " << watch.timeLapFormatted(6) << std::endl;
	alignerObjForFamilyDet.processAlnInfoOutputNoCheck(sampDirMaster.extractAlnCacheDir_.string(), verbose_);
}


} /* namespace bibseq */
