#pragma once
/*
 * info.hpp
 *
 *  Created on: Jan 12, 2015
 *      Author: nickhathaway
 */


#include <bibseq.h>
#include "mipster/common.h"

namespace bibseq {

VecStr processMipExtractInfoFile(table info);
table getSampleStats(const std::string & dirName, bool verbose);
table getSampleStats(const std::string & dirName, bool verbose, const VecStr & sampNames);
table getSampleMipStats(const std::string & dirName, bool verbose, const VecStr & sampNames);

void updateNameWithBarinfo(sampleCluster & clus);

struct genClusInfoWithBars;

struct genInfoWithBars {
	uint32_t totalReadCnt_ = 0;
	uint32_t totalBarcodeCnt_ = 0;
	uint32_t totalClusterCnt_ = 0;

	void increaseCounts(const genClusInfoWithBars & info);

	VecStr getInfoVec() const;
	std::string getInfo(const std::string & delim) const;

	static VecStr getInfoVecHeader();
	static std::string getInfoHeader(const std::string & delim);
};

struct genClusInfoWithBars {
	genClusInfoWithBars(const seqInfo & seqBase,
			const std::string & firstReadName, uint32_t clusterID,
			const std::string & expectsStr);
	seqInfo seqBase_;
	uint32_t clusterID_;
	uint32_t readCnt_;
	double readFrac_;
	uint32_t barcodeCnt_;
	double barcodeFrac_;
	std::string expectsStr_;

	void setFracInfo(const genInfoWithBars & totalInfo);
	std::string getInfo(const std::string & delim, bool header, bool checkingExpected)const;
	static std::string getInfoHeader(const std::string & delim, bool checkingExpected);
};

void sortGenClusInfoWithBarsVec(std::vector<genClusInfoWithBars> & vec);
class genSampeInfoWithBars {
public:
	genSampeInfoWithBars(const std::string & sampName);

	std::string sampName_;

	genInfoWithBars input_;

	genInfoWithBars used_;
	genInfoWithBars chiExcl_;
	genInfoWithBars freqExcl_;

	std::vector<genClusInfoWithBars> clusInput_;

	std::vector<genClusInfoWithBars> clusUsed_;
	std::vector<genClusInfoWithBars> clusChiExcl_;
	std::vector<genClusInfoWithBars> clusFreqExcl_;

	void addUsedInfo(const genClusInfoWithBars & usedClus);
	void addChiExclInfo(const genClusInfoWithBars & chiExclClus);
	void addFreqExclInfo(const genClusInfoWithBars & freqExclClus);

	void setFractions();

	VecStr getSampInfo(bool setFractionsFirst);
	static std::string getSampInfoHeader();

};


struct genInputPopClusInfoWithBars {

	genInputPopClusInfoWithBars(const seqInfo & seqBase);

	seqInfo seqBase_;
	uint32_t barcodeCnt_;
	uint32_t readCnt_;
	double barcodeFrac_;

};




struct genPopClusInfoWithBars {

	genPopClusInfoWithBars(const std::string & popUid, const seqInfo & seqBase);

	std::string popUid_;
	std::string mipPopUid_;
	seqInfo seqBase_;

	uint32_t sampCnt_ = 0;
	double sampFrac_ = 0;
	uint32_t readCnt_ = 0;
	double readFrac_ = 0;
	uint32_t barcodeCnt_ = 0;
	double barcodeFrac_ = 0;

	double medianBarcodeFrac_ = 0;
	double meanBarcodeFrac_ = 0;

	std::vector<genInputPopClusInfoWithBars> infos_;

	void addInfo(const seqInfo & seqBaseInfo);
	void setFractions(const genInfoWithBars & info, uint32_t sampTotal);

	void update();

	static std::string getInfoHeader(const std::string & delim);
	std::string getInfo(const std::string & delim)const;

};

struct genPopInfoWithBars {

	genPopInfoWithBars(const std::string & targetName);

	std::string targetName_;
	std::string geneName_;

	uint32_t totalInputClusters_ = 0;
	uint32_t totalInputSamples_ = 0;
	//in this total clusterTotal means sample total hap count
	genInfoWithBars totals_;

	std::vector<genPopClusInfoWithBars> infos_;
	std::unordered_map<std::string, uint32_t> subInfoPositions_;

	void addInfo(const genPopClusInfoWithBars & info);
	void updateInfoFracs();
	void updateSampCnt();

	static std::string getInfoHeader(const std::string & delim);
	std::string getInfoForPopUid(const std::string & popUid, const std::string & delim)const;
	std::string getGenInfo(const std::string & delim)const;
};

genSampeInfoWithBars generateSampInfo(const collapse::sampleCollapse & collapse);

struct PopClusTabs{
	table popTab_;
	table sampTab_;
};

PopClusTabs printMipSampleCollapseInfo(
		collapse::SampleCollapseCollection & sampCollapses, bool checkingExpected,
		std::string targetName);

} /* namespace bibseq */


