#pragma once
/*
 * mgv.hpp
 *
 *  Created on: Jan 30, 2017
 *      Author: nick
 */




#include "mipster/objects/MipsOnGenome.hpp"
#include <seqServer.h>

namespace bibseq {

class mgv: public bibseq::SeqApp {
	typedef bibseq::SeqApp super;

	bfs::path mainDir_;
	bfs::path serverResourceDir_;
	std::unique_ptr<MipsOnGenome> mipsInfo_;

public:

	mgv(const Json::Value & config);

	virtual std::vector<std::shared_ptr<restbed::Resource>> getAllResources();

	virtual VecStr requiredOptions() const;


	void getAllNamesHandler(std::shared_ptr<restbed::Session> session);
	void getMipsForRegionHandler(std::shared_ptr<restbed::Session> session);
	void getMipSeqsHandler(std::shared_ptr<restbed::Session> session);
	void getMipSeqsPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);

	void getGenomeLensHandler(std::shared_ptr<restbed::Session> session);
	void getGenomeLensPostHandler(std::shared_ptr<restbed::Session> session,
			const restbed::Bytes & body);



	//page handlers
	void mainPageHandler(std::shared_ptr<restbed::Session> session);
	void showMipSeqsHandler(std::shared_ptr<restbed::Session> session);
	void showRegionInfoHandler(std::shared_ptr<restbed::Session> session);

	//getting general info
	std::shared_ptr<restbed::Resource> getMipsForRegion();
	std::shared_ptr<restbed::Resource> getAllNames();
	std::shared_ptr<restbed::Resource> getGenomeLens();

	//get mip seqs
	std::shared_ptr<restbed::Resource> getMipSeqs();

	std::shared_ptr<restbed::Resource> mainPage();
	std::shared_ptr<restbed::Resource> showRegionInfo();
	std::shared_ptr<restbed::Resource> showMipSeqs();


	void redirect(std::shared_ptr<restbed::Session> session,
			std::string errorMessage);


	static void mgvErrorHandler(const int statusCode,
			const std::exception& exception,
			const std::shared_ptr<restbed::Session>& session);

};


}  // namespace bibseq




