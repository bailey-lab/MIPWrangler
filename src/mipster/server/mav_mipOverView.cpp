/*
 * mav_mipOverView.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: nick
 */


#include "mav.hpp"

namespace njhseq {


void mav::oneMipFamPageHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	if (njh::in(mipFam,popClusInfoByTar_)) {
		auto body = genHtmlDoc(rootName_, pages_.at("oneMipInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for mip "
				<< mipFam << "\n";
		ss << "Redirecting...." << "\n";
		redirect(session, ss.str());
	}
}

void mav::getNamesForMipFamHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	Json::Value ret;
	if (njh::in(mipFam,popClusInfoByTar_)) {
		std::set<std::string> sampNames;
		auto search = popClusInfoByTar_.find(mipFam);
		if (search != popClusInfoByTar_.end()) {
			auto currentSamps = search->second.get().getColumnLevels("s_sName");
			std::copy(currentSamps.begin(), currentSamps.end(),
					std::inserter(sampNames, sampNames.end()));
		}
		ret["regionName"] = mipMaster_->getGroupForMipFam(mipFam);
		ret["samples"] = njh::json::toJson(sampNames);
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": error, no information for mip "
				<< mipFam << "\n";
	}

	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}


std::shared_ptr<restbed::Resource> mav::oneMipFamPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "showMipInfo" }, {
					"mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						oneMipFamPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getNamesForMipFam() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "getNamesForMipFam" }, {
					"mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getNamesForMipFamHandler(session);
					}));
	return resource;
}



} //namesapce njhseq


