/*
 * info.cpp
 *
 *  Created on: Jan 12, 2015
 *      Author: nickhathaway
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
#include "info.hpp"
#include <unordered_map>

namespace njhseq {
namespace bfs = boost::filesystem;

VecStr processMipExtractInfoFile(table info){
	//table ret{VecStr{"matchingExtArm", "UnmatachedExtArm", "readsFailing_LigationArm",
		//"readsUsed","readsFailing_BarcodeFiltering", "readsFailing_Minlen", "readsFailing_Quality"}};
	//trim off (%) so elements can be treated as numbers
	info.trimElementsAtFirstOccurenceOf("(");
	uint32_t unmatched = estd::stou(info.content_.back()[1]);
	info.content_.erase(info.content_.end()-1);
	uint32_t readsUsed = getSumOfVecStr<uint32_t>(info.getColumn("readsUsed"));
	uint32_t readsFailingBarCodeFiltering = getSumOfVecStr<uint32_t>(info.getColumn("readsNotUsed"));
	uint32_t failedMinLen = getSumOfVecStr<uint32_t>(info.getColumnsStartWith("failedMinLen").getColumn(0));
	uint32_t failedQuality = getSumOfVecStr<uint32_t>(info.getColumnsStartWith("failed_q").getColumn(0));

	uint32_t failedLigationArm = getSumOfVecStr<uint32_t>(info.getColumn("failedLigationArm"));
	uint32_t readsMatchingExtArm = getSumOfVecStr<uint32_t>(info.getColumn("readNumber"));
	return toVecStr(readsMatchingExtArm, unmatched, failedLigationArm,
			readsUsed, readsFailingBarCodeFiltering, failedMinLen, failedQuality);
}

table getSampleStats(const std::string & dirName, bool verbose){
	auto sampDirs = njh::files::listAllFiles(dirName, false, VecStr{"samp"});
	if(sampDirs.empty()){
		if(verbose){
			std::cout << njh::bashCT::boldBlack("Extracting on dirs with GP instead of samp")<< std::endl;
		}
		sampDirs = njh::files::listAllFiles(dirName, false, {std::regex{".*/GP.*"}});
	}
	table finalInfo{VecStr{"sampleName","totalInitial","assembled", "discarded", "unassembled",
		"matchingExtArm", "UnmatachedExtArm", "readsFailing_LigationArm",
		"readsUsed","readsFailing_BarcodeFiltering", "readsFailing_Minlen", "readsFailing_Quality"}};
	for(const auto & d : sampDirs){
		if(d.second){
			auto pathName = d.first.string();
			njh::appendAsNeeded(pathName, "/");
			auto sName = d.first.filename().string();
			if(verbose){
				std::cout << njh::bashCT::bold << njh::bashCT::addColor(210) << "Processing: "
						<< njh::bashCT::addColor(105) << sName
						<< njh::bashCT::reset << std::endl;
			}
			//njh::scopedStopWatch watch(sName, true);
			auto assembled = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".assembled.fastq"),verbose);
			auto discarded = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".discarded.fastq"),verbose);
			auto unassembled = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".unassembled.forward.fastq"),verbose);
			uint32_t total = assembled + discarded + unassembled;
			if(assembled > 0){
				auto extractionDirs = njh::files::listAllFiles(pathName, false, VecStr{sName + ".assembled_mip"});
				std::map<bfs::path, bool, std::greater<bfs::path>> eDirs;
				for(const auto & ed : extractionDirs){
					if(ed.second){
						eDirs.emplace(ed);
					}
				}
				if (!eDirs.empty()) {
					auto resultsDir = njh::appendAsNeededRet(eDirs.begin()->first.string(), "/");
					table info { resultsDir + "info.txt", "\t", true };
					if(info.content_.empty()){
						continue;
					}
					auto processedInfo = processMipExtractInfoFile(info);
					finalInfo.content_.emplace_back(
							concatVecs(
									toVecStr(sName, total, assembled, discarded, unassembled),
									processedInfo));
				}else{
					std::cout << "Couldn't find latest analysis file for sample " << sName << std::endl;
				}
			}
		}
	}
	return finalInfo;
}

table getSampleStats(const std::string & dirName, bool verbose, const VecStr & sampNames){
	auto sampDirs = njh::files::listAllFiles(dirName, false, VecStr{});
	table finalInfo{VecStr{"sampleName","totalInitial","assembled", "discarded", "unassembled",
		"matchingExtArm", "UnmatachedExtArm", "readsFailing_LigationArm",
		"readsUsed","readsFailing_BarcodeFiltering", "readsFailing_Minlen", "readsFailing_Quality"}};
	for(const auto & d : sampDirs){
		if(!njh::in(d.first.filename().string(), sampNames)){
			continue;
		}
		if(d.second){
			auto pathName = d.first.string();
			njh::appendAsNeeded(pathName, "/");
			auto sName = d.first.filename().string();
			if(verbose){
				std::cout << njh::bashCT::bold << njh::bashCT::addColor(210) << "Processing: "
						<< njh::bashCT::addColor(105) << sName
						<< njh::bashCT::reset << std::endl;
			}
			//njh::scopedStopWatch watch(sName, true);
			auto assembled = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".assembled.fastq"),verbose);
			auto discarded = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".discarded.fastq"),verbose);
			auto unassembled = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".unassembled.forward.fastq"),verbose);
			uint32_t total = assembled + discarded + unassembled;
			if(assembled > 0){
				auto extractionDirs = njh::files::listAllFiles(pathName, false, VecStr{sName + ".assembled_mip"});
				std::map<bfs::path, bool, std::greater<bfs::path>> eDirs;
				for(const auto & ed : extractionDirs){
					if(ed.second){
						eDirs.emplace(ed);
					}
				}
				if (!eDirs.empty()) {
					auto resultsDir = njh::appendAsNeededRet(eDirs.begin()->first.string(), "/");
					table info { resultsDir + "info.txt", "\t", true };
					if(info.content_.empty()){
						continue;
					}
					auto processedInfo = processMipExtractInfoFile(info);
					finalInfo.content_.emplace_back(
							concatVecs(
									toVecStr(sName, total, assembled, discarded, unassembled),
									processedInfo));
				}else{
					std::cout << "Couldn't find latest analysis file for sample " << sName << std::endl;
				}
			}
		}
	}
	return finalInfo;
}

table getSampleMipStats(const std::string & dirName, bool verbose, const VecStr & sampNames){
	auto sampDirs = njh::files::listAllFiles(dirName, false, VecStr{});
	table mipFinalInfo;
	for(const auto & d : sampDirs){
		if(!njh::in(d.first.filename().string(), sampNames)){
			continue;
		}
		if(d.second){
			auto pathName = d.first.string();
			njh::appendAsNeeded(pathName, "/");
			auto sName = d.first.filename().string();
			if(verbose){
				std::cout << njh::bashCT::bold << njh::bashCT::addColor(210) << "Processing: "
						<< njh::bashCT::addColor(105) << sName
						<< njh::bashCT::reset << std::endl;
			}
			//njh::scopedStopWatch watch(sName, true);
			auto assembled = countSeqs(SeqIOOptions::genFastqIn(pathName + sName + ".assembled.fastq"), verbose);
			if(assembled > 0){
				auto extractionDirs = njh::files::listAllFiles(pathName, false, VecStr{sName + ".assembled_mip"});
				std::map<bfs::path, bool, std::greater<bfs::path>> eDirs;
				for(const auto & ed : extractionDirs){
					if(ed.second){
						eDirs.emplace(ed);
					}
				}
				if (!eDirs.empty()) {
					auto resultsDir = njh::appendAsNeededRet(eDirs.begin()->first.string(), "/");
					table info { resultsDir + "info.txt", "\t", true };
					if(info.content_.empty()){
						continue;
					}
					info.addColumn({sName}, "SampleName");
					if(mipFinalInfo.content_.empty()){
						mipFinalInfo = info;
					}else{
						mipFinalInfo.rbind(info, false);
					}
				}else{
					std::cout << "Couldn't find latest analysis file for sample " << sName << std::endl;
				}
			}
		}
	}
	return mipFinalInfo;
}


void updateNameWithBarinfo(sampleCluster & clus){
	clus.seqBase_.name_ = clus.firstReadName_;
	auto tPos = clus.seqBase_.name_.rfind("_t");
	clus.seqBase_.name_.at(tPos + 1) = 'B';
	clus.updateName();
}

void genInfoWithBars::increaseCounts(const genClusInfoWithBars & info){
	++totalClusterCnt_;
	totalBarcodeCnt_ += info.barcodeCnt_;
	totalReadCnt_ += info.readCnt_;
}


void processNameForBarReadInfo(const std::string & name,
		uint32_t & barNum, uint32_t & readNum){
	auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	auto underPos = name.rfind("_");
	//std::cout << vectorToString(toVecStr(rPos, bPos, underPos, name.substr(rPos + 2, bPos - rPos - 2),name.substr(bPos + 2, underPos - bPos - 2) ), "\n") << std::endl;
	readNum = estd::stou(name.substr(rPos + 2, bPos - rPos - 2));
	barNum = estd::stou(name.substr(bPos + 2, underPos - bPos - 2));
}
genClusInfoWithBars::genClusInfoWithBars(const seqInfo & seqBase,
		const std::string & firstReadName, uint32_t clusterID,
		const std::string & expectsStr) :
		seqBase_(seqBase.name_, seqBase.seq_, seqBase.qual_), clusterID_(clusterID), expectsStr_(
				expectsStr) {
	//scopedMessage mess{"genClusInfoWithBars_begin", "genClusInfoWithBars_end",std::cout, true};
	processNameForBarReadInfo(seqBase_.name_, barcodeCnt_, readCnt_);

	barcodeFrac_ = 0;
	readFrac_ = 0;
}

void genClusInfoWithBars::setFracInfo(const genInfoWithBars & totalInfo){
	barcodeFrac_ = barcodeCnt_ /static_cast<double>(totalInfo.totalBarcodeCnt_);
	readFrac_ = readCnt_ /static_cast<double>(totalInfo.totalReadCnt_);
}


std::string genClusInfoWithBars::getInfoHeader(const std::string & delim,
		bool checkingExpected) {
	VecStr headers { "c_clusterID", "c_name", "c_readCnt", "c_readFrac",
			"c_barcodeCnt", "c_barcodeFrac", "c_seq", "c_qual", "c_length" };
	if (checkingExpected) {
		headers.emplace_back("c_bestExpected");
	}
	return vectorToString(headers, delim);
}



std::string genClusInfoWithBars::getInfo(const std::string & delim, bool header,
		bool checkingExpected) const {
	if (header) {
		return genClusInfoWithBars::getInfoHeader(delim, checkingExpected);
	} else {
		if (checkingExpected) {
			return vectorToString(
					toVecStr(clusterID_, seqBase_.name_, readCnt_, readFrac_, barcodeCnt_,
							barcodeFrac_, seqBase_.seq_,
							seqBase_.getFastqQualString(SangerQualOffset), len(seqBase_),
							expectsStr_), delim);
		} else {
			return vectorToString(
					toVecStr(clusterID_, seqBase_.name_, readCnt_, readFrac_, barcodeCnt_,
							barcodeFrac_, seqBase_.seq_,
							seqBase_.getFastqQualString(SangerQualOffset), len(seqBase_)),
					delim);
		}
	}
}



VecStr genInfoWithBars::getInfoVec() const {
	return toVecStr(totalClusterCnt_, totalReadCnt_, totalBarcodeCnt_);
}

std::string genInfoWithBars::getInfo(const std::string & delim) const {
	return vectorToString(getInfoVec(), delim);
}

VecStr genInfoWithBars::getInfoVecHeader() {
	return VecStr { "totalClusterCnt", "totalReadCnt", "totalBarcodeCnt" };
}

std::string genInfoWithBars::getInfoHeader(const std::string & delim) {
	return vectorToString(getInfoVecHeader(), delim);
}

genSampeInfoWithBars::genSampeInfoWithBars(const std::string & sampName) :
		sampName_(sampName) {

}

void genSampeInfoWithBars::addUsedInfo(const genClusInfoWithBars & usedClus){
	clusUsed_.emplace_back(usedClus);
	clusInput_.emplace_back(usedClus);
	used_.increaseCounts(usedClus);
	input_.increaseCounts(usedClus);
}
void genSampeInfoWithBars::addChiExclInfo(const genClusInfoWithBars & chiExclClus){
	clusChiExcl_.emplace_back(chiExclClus);
	clusInput_.emplace_back(chiExclClus);
	chiExcl_.increaseCounts(chiExclClus);
	input_.increaseCounts(chiExclClus);
}
void genSampeInfoWithBars::addFreqExclInfo(const genClusInfoWithBars & freqExclClus){
	clusFreqExcl_.emplace_back(freqExclClus);
	clusInput_.emplace_back(freqExclClus);
	freqExcl_.increaseCounts(freqExclClus);
	input_.increaseCounts(freqExclClus);
}

void genSampeInfoWithBars::setFractions(){
	//input
	for(auto & info : clusInput_){
		info.setFracInfo(input_);
	}
	//used
	for(auto & info : clusUsed_){
		info.setFracInfo(used_);
	}
	//chi, set frac with the input total to indicate original frac
	for(auto & info : clusFreqExcl_){
		info.setFracInfo(input_);
	}
	//freq, set frac with the total input to indicate original frac
	for(auto & info : clusFreqExcl_){
		info.setFracInfo(input_);
	}
}

void sortGenClusInfoWithBarsVec(std::vector<genClusInfoWithBars> & vec){
	njh::sort(vec, [](const genClusInfoWithBars & info1, genClusInfoWithBars & info2){
		return info1.clusterID_ > info2.clusterID_;
	});
}




VecStr genSampeInfoWithBars::getSampInfo(bool setFractionsFirst) {
	if (setFractionsFirst) {
		setFractions();
	}
	return toVecStr(sampName_, used_.getInfoVec(), input_.getInfoVec(),
			chiExcl_.getInfoVec(), freqExcl_.getInfoVec());
}



std::string genSampeInfoWithBars::getSampInfoHeader() {
	return vectorToString(VecStr { "s_sName", "s_usedTotalClusterCnt",
			"s_usedTotalReadCnt", "s_usedTotalBarcodeCnt", "s_inputTotalClusterCnt",
			"s_inputTotalReadCnt", "s_inputTotalBarcodeCnt", "s_chiTotalClusterCnt",
			"s_chiTotalReadCnt", "s_chiTotalBarcodeCnt", "s_lowFreqTotalClusterCnt",
			"s_lowFreqTotalReadCnt", "s_lowFreqTotalBarcodeCnt" }, "\t");
}


genInputPopClusInfoWithBars::genInputPopClusInfoWithBars(
		const seqInfo & seqBase) :
		seqBase_(seqBase) {
	processNameForBarReadInfo(seqBase.name_, barcodeCnt_, readCnt_);
	auto underPos = seqBase.name_.rfind("_f");
	barcodeFrac_ = njh::lexical_cast<double>(seqBase.name_.substr(underPos + 2));
}

genPopClusInfoWithBars::genPopClusInfoWithBars(const std::string & popUid,
		const seqInfo & seqBase) :
		popUid_(popUid), mipPopUid_(seqBase.getStubName(true)), seqBase_(seqBase) {
	if(MetaDataInName::nameHasMetaData(mipPopUid_)){
		MetaDataInName::removeMetaDataInName(mipPopUid_);
	}
}

void genPopClusInfoWithBars::addInfo(const seqInfo & seqBaseInfo) {
	infos_.emplace_back(seqBaseInfo);
}

void genPopClusInfoWithBars::setFractions(const genInfoWithBars & info,
		uint32_t sampTotal) {

	sampFrac_ = sampCnt_ / static_cast<double>(sampTotal);
	barcodeFrac_ = barcodeCnt_ / static_cast<double>(info.totalBarcodeCnt_);
	readFrac_ = readCnt_ / static_cast<double>(info.totalReadCnt_);
}

void genPopClusInfoWithBars::update() {
	std::vector<double> barFracs;

	for (const auto & i : infos_) {
		++sampCnt_;
		barcodeCnt_ += i.barcodeCnt_;
		readCnt_ += i.readCnt_;
		barFracs.emplace_back(i.barcodeFrac_);
	}
	medianBarcodeFrac_ = vectorMedianRef(barFracs);
	meanBarcodeFrac_ = vectorMean(barFracs);
}

std::string genPopClusInfoWithBars::getInfoHeader(const std::string & delim) {
	return vectorToString(VecStr { "h_popUID", "h_mipPopUID", "h_sampleCnt", "h_sampleFrac",
			"h_medianBarcodeFrac", "h_meanBarcodeFrac", "h_readCnt", "h_readFrac",
			"h_barcodeCnt", "h_barcodeFrac", "h_inputNames", "h_seq", "h_qual" },
			delim);
}

std::string genPopClusInfoWithBars::getInfo(const std::string & delim) const {
	VecStr names;
	names.reserve(infos_.size());
	for (const auto & i : infos_) {
		names.emplace_back(i.seqBase_.name_);
	}
	return vectorToString(
			toVecStr(popUid_, mipPopUid_, sampCnt_, sampFrac_, medianBarcodeFrac_,
					meanBarcodeFrac_, readCnt_, readFrac_, barcodeCnt_, barcodeFrac_,
					vectorToString(names, ","), seqBase_.seq_,
					seqBase_.getFastqQualString(SangerQualOffset)), delim);
}



genPopInfoWithBars::genPopInfoWithBars(const std::string & targetName) :
		targetName_(targetName) {
	auto underScorePos = targetName_.find("_");
	if(std::string::npos != underScorePos){
		geneName_ = targetName_.substr(0, underScorePos);
	}else{
		geneName_ = targetName_;
	}
}

void genPopInfoWithBars::updateSampCnt(){
	std::set<std::string> sampNames_;
	for(const auto & i : infos_){
		for(const auto & subI : i.infos_){
		  VecStr toks = tokenizeString(subI.seqBase_.name_, ".");
		  sampNames_.emplace(njh::replaceString(toks[0], "CHI_", ""));
		}
	}

	totalInputSamples_ = sampNames_.size();
}

void genPopInfoWithBars::addInfo(const genPopClusInfoWithBars & info) {
	subInfoPositions_[info.popUid_] = infos_.size();
	infos_.emplace_back(info);
	totals_.totalBarcodeCnt_ += info.barcodeCnt_;
	totals_.totalReadCnt_ += info.readCnt_;
	++totals_.totalClusterCnt_;
	totalInputClusters_ += info.infos_.size();
}

void genPopInfoWithBars::updateInfoFracs() {
	updateSampCnt();
	for (auto & i : infos_) {
		i.setFractions(totals_, totalInputSamples_);
	}
}

std::string genPopInfoWithBars::getInfoHeader(const std::string & delim) {
	return vectorToString(VecStr {"p_geneName", "p_targetName", "p_sampleTotal",
			"p_totalInputClusters", "p_readTotal", "p_barcodeTotal",
			"p_finalHaplotypeNumber" }, delim);
}

std::string genPopInfoWithBars::getGenInfo(const std::string & delim) const {
	return vectorToString(
			toVecStr(geneName_, targetName_, totalInputSamples_, totalInputClusters_,  totals_.totalReadCnt_,
					totals_.totalBarcodeCnt_, totals_.totalClusterCnt_), delim);
}
std::string genPopInfoWithBars::getInfoForPopUid(const std::string & popUid, const std::string & delim)const{
	return getGenInfo(delim) + delim + infos_[subInfoPositions_.at(popUid)].getInfo(delim);
}



genSampeInfoWithBars generateSampInfo(const collapse::sampleCollapse & collapse){
  genSampeInfoWithBars info(collapse.sampName_);
  //add used
	uint32_t pos = 0;
	for (const auto &clus : iter::reversed(collapse.collapsed_.clusters_)) {
		info.addUsedInfo(
				genClusInfoWithBars(clus.seqBase_, clus.firstReadName_,
						collapse.collapsed_.clusters_.size() - 1 - pos, clus.expectsString));
		++pos;
	}
  //add excluded
  for(const auto & clus : collapse.excluded_.clusters_){
  	if(clus.seqBase_.name_.find("CHI")!=std::string::npos){
  		info.addChiExclInfo(genClusInfoWithBars(clus.seqBase_, clus.firstReadName_,
  				pos, clus.expectsString));
  	}else{
  		info.addFreqExclInfo(genClusInfoWithBars(clus.seqBase_, clus.firstReadName_,
  		  				pos, clus.expectsString));
  	}
  	++pos;
  }
  return info;
}



PopClusTabs printMipSampleCollapseInfo(
		collapse::SampleCollapseCollection & sampCollapses, bool checkingExpected,
		std::string targetName){
	if(!sampCollapses.popCollapse_){
		sampCollapses.loadInPreviousPop();
	}
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	std::stringstream infoFile;
	infoFile << "s_Sample";
	infoFile << "\t" << genPopInfoWithBars::getInfoHeader("\t");
	infoFile << "\t" << genPopClusInfoWithBars::getInfoHeader("\t");
	infoFile << "\t" << genSampeInfoWithBars::getSampInfoHeader();
	if(sampCollapses.groupMetaLoaded()){
		for(const auto & group : sampCollapses.groupMetaData_->groupData_){
			infoFile << "\t" << "s_" << group.first;
		}
	}
	infoFile << "\t" << genClusInfoWithBars::getInfoHeader("\t", checkingExpected);
	infoFile << std::endl;
	genPopInfoWithBars popInfo(targetName);
	for (const auto & clus : sampCollapses.popCollapse_->collapsed_.clusters_) {
		genPopClusInfoWithBars currentInfo(clus.seqBase_.getStubName(true),
				clus.seqBase_);
		for (const auto & subClus : clus.reads_) {
			currentInfo.addInfo(subClus->seqBase_);
		}
		currentInfo.update();
		popInfo.addInfo(currentInfo);
	}
	popInfo.updateInfoFracs();
	for (auto& samp : sampCollapses.popNames_.samples_) {
		sampCollapses.setUpSampleFromPrevious(samp);
		const auto & sampCollapse = sampCollapses.sampleCollapses_.at(samp);
		auto sampInfo = generateSampInfo(*sampCollapse);
		auto genSampInfo = sampInfo.getSampInfo(true);
		if(sampCollapses.groupMetaLoaded()){
			for(const auto & group : sampCollapses.groupMetaData_->groupData_){
				genSampInfo.emplace_back(group.second->getGroupForSample(samp));
			}
		}

		for (const auto & clus : sampInfo.clusUsed_) {
			infoFile << samp << "\t"
					<< popInfo.getInfoForPopUid(
							sampCollapses.popCollapse_->collapsed_.clusters_[sampCollapses.popCollapse_->collapsed_.subClustersPositions_.at(
									clus.seqBase_.getStubName(true))].seqBase_.getStubName(true),
							"\t");
			infoFile << "\t";
			infoFile << njh::conToStr(genSampInfo, "\t");
			infoFile << "\t" << clus.getInfo("\t", false, checkingExpected)
					<< std::endl;

		}
	}
	std::stringstream popInfoOutFile;
	popInfoOutFile << genPopInfoWithBars::getInfoHeader("\t");
	popInfoOutFile << "\t" << genPopClusInfoWithBars::getInfoHeader("\t");
	if(sampCollapses.groupMetaLoaded()){
		for (const auto & group : sampCollapses.groupMetaData_->groupData_) {
			popInfoOutFile << "\t" << "g_" + group.first + "_subGroupCounts"
					<< "\t" << "g_" + group.first + "_subGroupFracs";
		}
	}
	popInfoOutFile << std::endl;
	for (const auto & i : popInfo.infos_) {
		popInfoOutFile << popInfo.getGenInfo("\t")
				<< "\t" << i.getInfo("\t");
		if(sampCollapses.groupMetaLoaded()){
			std::vector<std::string> currentSamples;
			std::set<std::string> sampNames;
			for(const auto & subI : i.infos_){
			  VecStr toks = tokenizeString(subI.seqBase_.name_, ".");
			  sampNames.emplace(njh::replaceString(toks[0], "CHI_", ""));
			}
			currentSamples = std::vector<std::string>{sampNames.begin(), sampNames.end()};
			auto groupPopInfos = sampCollapses.groupMetaData_->getGroupPopInfos(currentSamples);
			for (const auto & info : groupPopInfos) {
				popInfoOutFile << "\t" << info.groupCountsStr() << "\t"
						<< info.groupFracsStr();
			}
		}
		popInfoOutFile << std::endl;
	}
	PopClusTabs ret;
	ret.sampTab_ = table(infoFile, "\t", true);
	ret.popTab_ = table(popInfoOutFile, "\t", true);
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	return ret;
}



} /* namespace njhseq */
