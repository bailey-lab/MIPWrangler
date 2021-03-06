#pragma once
/*
 * mgv.hpp
 *
 *  Created on: Jan 30, 2017
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


#include "mipster/objects/MipsOnGenome.hpp"
#include <seqServer.h>

namespace njhseq {

class mgv: public njhseq::SeqApp {
	typedef njhseq::SeqApp super;

	bfs::path mainDir_;
	bfs::path serverResourceDir_;
	std::unique_ptr<MipsOnGenome> mipsInfo_;

	std::unique_ptr<TableCache> allTarInfo_;

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

	void getMipRegionInfoForGenomeHandler(
			std::shared_ptr<restbed::Session> session);
	void getMipRegionInfoForGenomePostHandler(
			std::shared_ptr<restbed::Session> session, const restbed::Bytes & body);


	void getMipTarInfoForGenomeHandler(
			std::shared_ptr<restbed::Session> session);
	void getMipTarInfoForGenomePostHandler(
			std::shared_ptr<restbed::Session> session, const restbed::Bytes & body);



	void getMipGenomeLocsHandler(std::shared_ptr<restbed::Session> session);

	//page handlers
	void mainPageHandler(std::shared_ptr<restbed::Session> session);
	void showMipSeqsHandler(std::shared_ptr<restbed::Session> session);
	void showRegionInfoHandler(std::shared_ptr<restbed::Session> session);

	void getInfoMipArmsInfoHandler(std::shared_ptr<restbed::Session> session);
	void getAllMipTarInfoAllGenomesHandler(std::shared_ptr<restbed::Session> session);

	//getting general info
	//void (std::shared_ptr<restbed::Session> session);
	std::shared_ptr<restbed::Resource> getMipGenomeLocs();
	std::shared_ptr<restbed::Resource> getMipsForRegion();
	std::shared_ptr<restbed::Resource> getAllNames();
	std::shared_ptr<restbed::Resource> getGenomeLens();

	std::shared_ptr<restbed::Resource> getMipRegionInfoForGenome();
	std::shared_ptr<restbed::Resource> getMipTarInfoForGenome();

	std::shared_ptr<restbed::Resource> getInfoMipArmsInfo();
	std::shared_ptr<restbed::Resource> getAllMipTarInfoAllGenomes();

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


}  // namespace njhseq




