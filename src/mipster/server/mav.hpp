#pragma once
/*
 * mipAnalysisServer.hpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */



#include <seqServer.h>
#include "mipster/common.h"
#include "mipster/setUp/SetUpMaster.hpp"


namespace bibseq {





class mav: public bibseq::SeqApp {

	bfs::path masterDir_;
	bfs::path mipArmsFileName_;
	bfs::path mipsSamplesFile_;
	bfs::path serverResourceDir_;

	std::shared_ptr<SetUpMaster> mipMaster_;

	//extraction info
	std::unordered_map<std::string, TableCache> extractInfosBySamp_;
	std::unordered_map<std::string, TableCache> extractInfosByTar_;
	std::shared_ptr<TableCache> masterExtractInfo_;

	//pop clus info
	std::unordered_map<std::string, TableCache> popClusInfoBySample_;
	std::unordered_map<std::string, TableCache> popClusInfoByTar_;
	std::unordered_map<std::string, TableCache> popClusPopInfoByTar_;
	/*
	std::unordered_map<std::string, SeqIOOptsWithTime> popHapsByTar_;
	std::unordered_map<std::string,
			std::unordered_map<std::string, SeqIOOptsWithTime>> finalHapsByTarSamp_;
	std::unordered_map<std::string,
			std::unordered_map<std::string, SeqIOOptsWithTime>> originalHapsByTarBySamp_;
			*/

	void redirect(std::shared_ptr<restbed::Session> session,
			std::string errorMessage);

	void getAllNamesHandler(std::shared_ptr<restbed::Session> session);
	void mainPageHandler(std::shared_ptr<restbed::Session> session);
	//one region
	void oneRegionPageHandler(std::shared_ptr<restbed::Session> session);
	void getNamesForRegionHandler(std::shared_ptr<restbed::Session> session);
	//one mip fam
	void oneMipFamPageHandler(std::shared_ptr<restbed::Session> session);
	void getNamesForMipFamHandler(std::shared_ptr<restbed::Session> session);

	//extraction
	//extraction per target
	void initialReadStatsPerMipTarPageHandler(std::shared_ptr<restbed::Session> session);
	void samplesForExtractedMipHandler(std::shared_ptr<restbed::Session> session);
	void getInitialReadStatsPerMipTarPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getInitialReadStatsPerMipTarHandler(std::shared_ptr<restbed::Session> session);
	//extraction per samples
	void initialReadStatsPerSamplePageHandler(std::shared_ptr<restbed::Session> session);
	void extractedMipsForSampleHandler(std::shared_ptr<restbed::Session> session);
	void getInitialReadStatsPerSamplePostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getInitialReadStatsPerSampleHandler(std::shared_ptr<restbed::Session> session);

	//overall extraction
	void initialReadStatsPageHandler(std::shared_ptr<restbed::Session> session);
	void getInitialReadStatsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getInitialReadStatsHandler(std::shared_ptr<restbed::Session> session);

	//one sample, all mips
	void oneSampAllMipDataPageHandler(std::shared_ptr<restbed::Session> session);
	void mipFamNamesForSampHandler(std::shared_ptr<restbed::Session> session);
	void getOneSampAllMipDataPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getOneSampAllMipDataHandler(std::shared_ptr<restbed::Session> session);

	//one region one sample
	void regionInfoForSampPageHandler(std::shared_ptr<restbed::Session> session);
	void getMipOverlapGraphDataHandler(std::shared_ptr<restbed::Session> session);

	//one mip all samps
	void oneMipAllSampsPageHandler(std::shared_ptr<restbed::Session> session);
	void getOneMipAllSampsDataPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getOneMipAllSampsDataHandler(std::shared_ptr<restbed::Session> session);
	void getOneMipPopSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getOneMipPopSeqsHandler(std::shared_ptr<restbed::Session> session);
	void getOneMipAllSampsPopDataPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getOneMipAllSampsPopDataHandler(std::shared_ptr<restbed::Session> session);

	//one mip one samp
	void oneMipOneSampDataPageHandler(std::shared_ptr<restbed::Session> session);
	void getMipOneSampOriginalSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getMipOneSampOriginalSeqsHandler(std::shared_ptr<restbed::Session> session);
	void getMipOneSampFinalSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);
	void getMipOneSampFinalSeqsHandler(std::shared_ptr<restbed::Session> session);
	void getOneMipOneSampsDataHandler(std::shared_ptr<restbed::Session> session);
public:

	mav(const Json::Value & config);

	virtual std::vector<std::shared_ptr<restbed::Resource>> getAllResources();
	virtual VecStr requiredOptions() const;

	typedef bibseq::SeqApp super;

	/**@brief get all the sample names, mip grouping names, mip family names, mip target names
	 *
	 * @return
	 */
	std::shared_ptr<restbed::Resource> getAllNames();
	std::shared_ptr<restbed::Resource> mainPage();

	//one region
	std::shared_ptr<restbed::Resource> oneRegionPage();
	std::shared_ptr<restbed::Resource> getNamesForRegion();
	//one mip fam
	std::shared_ptr<restbed::Resource> oneMipFamPage();
	std::shared_ptr<restbed::Resource> getNamesForMipFam();

	//extraction
	//extraction per target
	std::shared_ptr<restbed::Resource> initialReadStatsPerMipTarPage();
	std::shared_ptr<restbed::Resource> samplesForExtractedMip();
	std::shared_ptr<restbed::Resource> getInitialReadStatsPerMipTar();
	//extraction per samples
	std::shared_ptr<restbed::Resource> initialReadStatsPerSamplePage();
	std::shared_ptr<restbed::Resource> extractedMipsForSample();
	std::shared_ptr<restbed::Resource> getInitialReadStatsPerSample();
	//extraction overall
	std::shared_ptr<restbed::Resource> initialReadStatsPage();
	std::shared_ptr<restbed::Resource> getInitialReadStats();
	//one sample, all mips
	std::shared_ptr<restbed::Resource> oneSampAllMipDataPage();
	std::shared_ptr<restbed::Resource> mipFamNamesForSamp();
	std::shared_ptr<restbed::Resource> getOneSampAllMipData();
	//one region one sample
	std::shared_ptr<restbed::Resource> regionInfoForSampPage();
	std::shared_ptr<restbed::Resource> getMipOverlapGraphData();
	//one mip all samples
	std::shared_ptr<restbed::Resource> oneMipAllSampsPage();
	std::shared_ptr<restbed::Resource> getOneMipAllSampsData();
	std::shared_ptr<restbed::Resource> getOneMipPopSeqs();
	std::shared_ptr<restbed::Resource> getOneMipAllSampsPopData();
	//one mip one samp
	std::shared_ptr<restbed::Resource> oneMipOneSampDataPage();
	std::shared_ptr<restbed::Resource> getMipOneSampOriginalSeqs();
	std::shared_ptr<restbed::Resource> getMipOneSampFinalSeqs();
	std::shared_ptr<restbed::Resource> getOneMipOneSampsData();

	/*
	void mainPage();
	//all names


	//extraction info
	void showInitialReadStats();
	void getInitialReadStats();

	void showInitialReadStatsPerSample(std::string sampName);
	void getInitialReadStatsPerSample(std::string sampName);

	void showInitialReadStatsPerMipTar(std::string mipTar);
	void getInitialReadStatsPerMipTar(std::string mipTar);

	void extractedMipsForSample(std::string sampName);
	void samplesForExtractedMip(std::string mipTar);

	//samp pop info all mips
	void showOneSampAllMipData(std::string sampName);
	void getOneSampAllMipData(std::string sampName);
	void mipFamNamesForSamp(std::string sampName);

	//gene overview
	void showGeneInfo(std::string geneName);
	void mipFamNamesForGene(std::string geneName);
	void sampNamesForGene(std::string geneName);
	void mipOverlapGraphData(std::string geneName, std::string sampName);

	//one mip fam
	void showMipInfo(std::string mipName);
	void getGeneNameForMipName(std::string mipName);
	void getSampNamesForMipName(std::string mipName);

	//one mip all samps pop info
	void showOneMipAllSampsData(std::string mipName);
	void getOneMipAllSampsData(std::string mipName);
	void getOneMipAllSampsPopData(std::string mipName);
	void getMipPopSeqs(std::string mipName);

	void showOneMipOneSampsData(std::string mipName, std::string sampName);
	void getOneMipOneSampsData(std::string mipName, std::string sampName);
	void getOneMipOneSampsPopData(std::string mipName, std::string sampName);
	void getMipOneSampOriginalSeqs(std::string mipName, std::string sampName);
	void getMipOneSampFinalSeqs(std::string mipName, std::string sampName);

	void showGeneInfoForSampData(std::string geneName, std::string sampName);
	*/
};


int mavRunner(const bib::progutils::CmdArgs & inputCommands);

}  // namespace bibseq





