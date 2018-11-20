/*
 * mav_regionOverView.cpp
 *
 *  Created on: Feb 16, 2016
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
#include "mav.hpp"

namespace njhseq {


void mav::oneRegionPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto regionName = request->get_path_parameter("regionName");
	auto regions = mipMaster_->getMipGroupings();
	if (njh::in(regionName, regions)) {
		auto body = genHtmlDoc(rootName_, pages_.at("oneRegionInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< regionName << "\n";
		ss << "Options are " << njh::conToStr(regions, ", ") << "\n";
		ss << "Redirecting...." << "\n";
		redirect(session, ss.str());
	}
}


void mav::getNamesForRegionHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	auto regionName = request->get_path_parameter("regionName");
	auto regions = mipMaster_->getMipGroupings();
	Json::Value ret;
	if (njh::in(regionName, regions)) {
		VecStr mipNames = mipMaster_->getMipFamiliesForMipGroup(regionName);
		std::set<std::string> sampNames;
		for (const auto & mipName : mipNames) {
			auto search = popClusInfoByTar_.find(mipName);
			if (search != popClusInfoByTar_.end()) {
				auto currentSamps = search->second.get().getColumnLevels("s_sName");
				std::copy(currentSamps.begin(), currentSamps.end(),
						std::inserter(sampNames, sampNames.end()));
			}
		}
		ret["mipFamilies"] = njh::json::toJson(mipNames);
		ret["samples"] = njh::json::toJson(sampNames);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< regionName << "\n";
		std::cerr << "Options are " << njh::conToStr(regions, ", ") << "\n";
	}
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

std::shared_ptr<restbed::Resource> mav::oneRegionPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "showRegionInfo" },
	{"regionName", UrlPathFactory::pat_wordNumsDash_}}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						oneRegionPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getNamesForRegion() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "getNamesForRegion" }, {
					"regionName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getNamesForRegionHandler(session);
					}));
	return resource;
}


} //namespace njhseq

