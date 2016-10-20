/*
 * mav_regionOverView.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: nick
 */


#include "mav.hpp"

namespace bibseq {


void mav::oneRegionPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto regionName = request->get_path_parameter("regionName");
	auto regions = mipMaster_->getMipGroupings();
	if (bib::in(regionName, regions)) {
		auto body = genHtmlDoc(rootName_, pages_.at("oneRegionInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< regionName << "\n";
		ss << "Options are " << bib::conToStr(regions, ", ") << "\n";
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
	if (bib::in(regionName, regions)) {
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
		ret["mipFamilies"] = bib::json::toJson(mipNames);
		ret["samples"] = bib::json::toJson(sampNames);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< regionName << "\n";
		std::cerr << "Options are " << bib::conToStr(regions, ", ") << "\n";
	}
	auto body = bib::json::writeAsOneLine(ret);
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


} //namespace bibseq

