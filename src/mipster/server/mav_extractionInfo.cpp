/*
 * mav_extractionInfo.cpp
 *
 *  Created on: Feb 15, 2016
 *      Author: nick
 */



#include "mav.hpp"

namespace bibseq {

void mav::initialReadStatsPageHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto body = genHtmlDoc(rootName_, pages_.at("initialExtractStats.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

void mav::getInitialReadStatsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
	bib::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "samples" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	}else{
		VecStr samples = bib::json::jsonArrayToVec<std::string>(postData["samples"], [](const Json::Value & val){ return val.asString();});
		auto containsSampleName = [&samples](const std::string & str) {
			return bib::in(str, samples);
		};
		auto tab = masterExtractInfo_->get().extractByComp("Sample", containsSampleName);
		tab.trimElementsAtFirstOccurenceOf("(");
		ret = tableToJsonByRow(tab, "Sample", VecStr{"Sample", "raw", "matchingExtArm"});
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getInitialReadStatsHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getInitialReadStatsPostHandler(ses, body);
					}));
}

void mav::initialReadStatsPerMipTarPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mipTar");
	auto search = extractInfosByTar_.find(mipTar);
	if (search != extractInfosByTar_.end() && search->second.opts_.in_.inExists()) {
		auto body = genHtmlDoc(rootName_, pages_.at("initialExtractStatsMipTarget.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		if(search == extractInfosByTar_.end()){
			ss << __PRETTY_FUNCTION__ << ": error, no extraction information for mip target "
					<< mipTar << "\n";
			ss << "options are " << bib::conToStr(bib::getVecOfMapKeys(extractInfosByTar_), ", ") << "\n";
		}else if(!search->second.opts_.in_.inExists()){
			ss << __PRETTY_FUNCTION__ << ": error, extraction information for mip target "
					<< mipTar << " in file " << search->second.opts_.in_.inFilename_ << " doesn't exist" << "\n";
		}else{
			ss << __PRETTY_FUNCTION__ << ": error, in getting extraction information for mip target "
								<< mipTar << "\n";
		}
		ss << "Redirecting....";
		redirect(session, ss.str());
	}
}

void mav::initialReadStatsPerSamplePageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	auto search = extractInfosBySamp_.find(sample);
	if (search != extractInfosBySamp_.end() && search->second.opts_.in_.inExists()) {
		auto body = genHtmlDoc(rootName_, pages_.at("initialExtractStatsSample.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		if(search == extractInfosByTar_.end()){
			ss << __PRETTY_FUNCTION__ << ": error, no extraction information for sample "
					<< sample << "\n";
			ss << "options are " << bib::conToStr(bib::getVecOfMapKeys(extractInfosByTar_), ", ") << "\n";
		}else if(!search->second.opts_.in_.inExists()){
			ss << __PRETTY_FUNCTION__ << ": error, extraction information for sample "
					<< sample << " in file " << search->second.opts_.in_.inFilename_ << " doesn't exist" << "\n";
		}else{
			ss << __PRETTY_FUNCTION__ << ": error, in getting extraction information for sample "
								<< sample << "\n";
		}
		ss << "Redirecting....";
		redirect(session, ss.str());
	}
}

void mav::samplesForExtractedMipHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mipTar");
	Json::Value ret;
	auto search = extractInfosByTar_.find(mipTar);
	if(search != extractInfosByTar_.end() && search->second.opts_.in_.inExists()){
		ret["samples"] = bib::json::toJson(search->second.get().getColumn("sampleName"));
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "couldn't find mipTar: " << mipTar << std::endl;
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getInitialReadStatsPerMipTarPostHandler(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mipTar");
	const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
	bib::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "samples" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	}else{
		VecStr samples = bib::json::jsonArrayToVec<std::string>(postData["samples"], [](const Json::Value & val){ return val.asString();});
		auto search = extractInfosByTar_.find(mipTar);
		if(search != extractInfosByTar_.end() && search->second.opts_.in_.inExists()){
			auto containsSampleName = [&samples](const std::string & str) {
				return bib::in(str, samples);
			};
			auto tab = search->second.get().extractByComp("sampleName", containsSampleName);
			tab.trimElementsAtFirstOccurenceOf("(");
			ret = tableToJsonByRow(tab, "sampleName", VecStr{"sampleName", "readNumber", "goodReads"});
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "couldn't find mipTar: " << mipTar << std::endl;
		}
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getInitialReadStatsPerMipTarHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getInitialReadStatsPerMipTarPostHandler(ses, body);
					}));
}

void mav::extractedMipsForSampleHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	Json::Value ret;
	auto search = extractInfosBySamp_.find(sample);
	if(search != extractInfosByTar_.end() && search->second.opts_.in_.inExists()){
		ret["mipTargets"] = bib::json::toJson(search->second.get().getColumn("mipTarget"));
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "couldn't find sampName: " << sample << std::endl;
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getInitialReadStatsPerSamplePostHandler(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
	bib::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "mipTargets" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	}else{
		VecStr mipTargets = bib::json::jsonArrayToVec<std::string>(postData["mipTargets"], [](const Json::Value & val){ return val.asString();});
		auto search = extractInfosBySamp_.find(sample);
		if(search != extractInfosBySamp_.end() && search->second.opts_.in_.inExists()){
			auto containsMipTarget = [&mipTargets](const std::string & str) {
				return bib::in(str, mipTargets);
			};
			auto tab = search->second.get().extractByComp("mipTarget", containsMipTarget);
			tab.trimElementsAtFirstOccurenceOf("(");
			ret = tableToJsonByRow(tab, "mipTarget", VecStr{"mipTarget", "readNumber", "goodReads"});
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "couldn't find sample: " << sample << std::endl;
		}
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getInitialReadStatsPerSampleHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getInitialReadStatsPerSamplePostHandler(ses, body);
					}));
}


std::shared_ptr<restbed::Resource> mav::initialReadStatsPerMipTarPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showInitialReadStatsPerMipTar" }, { "mipTar",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						initialReadStatsPerMipTarPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::samplesForExtractedMip(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"samplesForExtractedMip" }, { "mipTar",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						samplesForExtractedMipHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getInitialReadStatsPerMipTar() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getInitialReadStatsPerMipTar" }, { "mipTar",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getInitialReadStatsPerMipTarHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::initialReadStatsPerSamplePage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showInitialReadStatsPerSample" }, { "sample",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						initialReadStatsPerSamplePageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::extractedMipsForSample(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"extractedMipsForSample" }, { "sample",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						extractedMipsForSampleHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getInitialReadStatsPerSample() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getInitialReadStatsPerSample" }, { "sample",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getInitialReadStatsPerSampleHandler(session);
					}));
	return resource;
}


std::shared_ptr<restbed::Resource> mav::initialReadStatsPage(){

	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showInitialReadStats" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						initialReadStatsPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getInitialReadStats() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getInitialReadStats" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getInitialReadStatsHandler(session);
					}));
	return resource;
}



}  // namespace bibseq
