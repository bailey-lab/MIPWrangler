/*
 * mgv.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: nick
 */

#include "mgv.hpp"


#include "mipster/mipUtils.h"

namespace bibseq {
mgv::mgv(const Json::Value & config) :
		SeqApp(config) {
	mainDir_ = config["masterDir"].asString();
	serverResourceDir_ = bib::appendAsNeededRet(config["resources"].asString(),
			"/");
	mipsInfo_ = std::make_unique<MipsOnGenome>(mainDir_,
			config["numThreads"].asUInt());
	mipsInfo_->loadInArms();
	mipsInfo_->loadInGenomes();
	if (config.isMember("primaryGenome")) {
		mipsInfo_->setPrimaryGenome(config["primaryGenome"].asString());
	}

	 jsFiles_->addFiles(
	 bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mgv/js"), ".js"));

	 cssFiles_->addFiles(
	 bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mgv/css"), ".css"));

	addScripts(bib::files::make_path(serverResourceDir_, "mgv"));

	//register mip seqs
	for (const auto & mipName : mipsInfo_->getMips()) {
		auto mipFaPath = mipsInfo_->pathToMipFasta(mipName);
		if (bfs::exists(mipFaPath)) {
			seqs_->updateAddCache(mipName,
					SeqIOOptions::genFastaIn(mipFaPath, false));
		}
	}
}

std::vector<std::shared_ptr<restbed::Resource>> mgv::getAllResources() {
	auto ret = super::getAllResources();
	//main page
	ret.emplace_back(mainPage());
	//basic info
	ret.emplace_back(getAllNames());
	ret.emplace_back(getGenomeLens());

	//show mip sequences
	ret.emplace_back(getMipSeqs());
	ret.emplace_back(showMipSeqs());

	ret.emplace_back(getMipsForRegion());
	ret.emplace_back(showRegionInfo());




	return ret;
}

VecStr mgv::requiredOptions() const {
	return concatVecs(super::requiredOptions(), VecStr { "masterDir",
			"numThreads", "primaryGenome", "resources" });
}
void mgv::mainPageHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto body = genHtmlDoc(rootName_, pages_.at("mainPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

void mgv::getAllNamesHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	Json::Value ret;

	auto genomes = mipsInfo_->getGenomes();
	bib::sort(genomes);
	ret["genomes"] = bib::json::toJson(genomes);
	auto mips = mipsInfo_->getMips();
	MipNameSorter::sort(mips);
	ret["mips"] = bib::json::toJson(mips);
	auto regions = mipsInfo_->mipArms_->getMipRegions();
	bib::sort(regions);
	ret["mipRegions"] = bib::json::toJson(regions);

	ret["primaryGenome"] = bib::json::toJson(mipsInfo_->getPrimaryGenome());

	auto body = bib::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}






void mgv::getMipsForRegionHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	auto regionName = request->get_path_parameter("regionName");
	auto regions = mipsInfo_->mipArms_->getMipRegions();
	Json::Value ret;
	if (bib::in(regionName, regions)) {
		VecStr mipNames = mipsInfo_->mipArms_->getMipTarsForRegion(regionName);
		MipNameSorter::sort(mipNames);
		ret["mips"] = bib::json::toJson(mipNames);
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




void mgv::redirect(std::shared_ptr<restbed::Session> session,
		std::string errorMessage) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	std::cerr << errorMessage << std::endl;
	auto body = genHtmlDoc(rootName_, pages_.at("redirectPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}

std::shared_ptr<restbed::Resource> mgv::getAllNames() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getAllNames" } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getAllNamesHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::getMipsForRegion() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getMipsForRegion" },
			{"regionName", UrlPathFactory::pat_wordNumsDash_} }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipsForRegionHandler(session);
					}));
	return resource;
}





std::shared_ptr<restbed::Resource> mgv::mainPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						mainPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::showMipSeqs() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "showMipSeqs" }, { "mip", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						showMipSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::showRegionInfo() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "showRegionInfo" },
			{ "regionName", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
							showRegionInfoHandler(session);
					}));
	return resource;
}



void mgv::showRegionInfoHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto regionName = request->get_path_parameter("regionName");
	auto regions = mipsInfo_->mipArms_->getMipRegions();
	Json::Value ret;
	if (bib::in(regionName, regions)) {
		auto body = genHtmlDoc(rootName_, pages_.at("showRegionInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << "Mip Name wasn't found: " << regionName << "\n";
		ss << "options are "
				<< bib::conToStr(regions, ", ")
				<< "\n";
		ss << "Redirecting....";
		redirect(session, ss.str());
	}

}

void mgv::showMipSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mip");
	if (!bib::in(mipTar, mipsInfo_->mipArms_->mips_)) {
		std::stringstream ss;
		ss << "Mip Name wasn't found: " << mipTar << "\n";
		ss << "options are "
				<< bib::conToStr(bib::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ")
				<< "\n";
		ss << "Redirecting....";
		redirect(session, ss.str());
	} else {
		auto body = genHtmlDoc(rootName_, pages_.at("showMipSeqs.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}
}

std::shared_ptr<restbed::Resource> mgv::getMipSeqs() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getMipSeqs" }, { "mip", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::getGenomeLens() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getGenomeLens" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getGenomeLensHandler(session);
					}));
	return resource;
}




void mgv::getMipSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = 0;
	request->get_header("Content-Length", content_length);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getMipSeqsPostHandler(ses, body);
					}));
}

void mgv::getGenomeLensHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = 0;
	request->get_header("Content-Length", content_length);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
							getGenomeLensPostHandler(ses, body);
					}));
}




void mgv::getGenomeLensPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
	bib::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "genomes" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto genomes = bib::json::jsonArrayToVec<std::string>(postData["genomes"], [](const Json::Value & val){ return val.asString();});
		VecStr missingGenomes;
		for(const auto & genome : genomes){
			if(!bib::in(genome, mipsInfo_->genomes_)){
				missingGenomes.emplace_back(genome);
			}else{
				ret[genome] = mipsInfo_->genomes_.at(genome)->chromosomeLengths();
			}
		}
		if(!missingGenomes.empty()){
			std::cerr << "Missing info for genomes " << bib::conToStr(missingGenomes, ", ") << std::endl;
		}
	}

	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mgv::getMipSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mip");
	const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
	Json::Value ret;
	if (!bib::in(mipTar, mipsInfo_->mipArms_->mips_)) {
		std::cerr << __PRETTY_FUNCTION__ << " error, no matching mip to " << mipTar
				<< std::endl;
		std::cerr << "options are "
				<< bib::conToStr(bib::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ")
				<< "\n";
	} else {
		uint32_t sesUid = std::numeric_limits<uint32_t>::max();
		//check to see if there is a session already started associated with this seq
		if (!postData.isMember("sessionUID")) {
			sesUid = startSeqCacheSession();
		} else {
			sesUid = postData["sessionUID"].asUInt();
		}
		/* could select just for certain genomes here
		 VecStr popUIDs = bib::json::jsonArrayToVec<std::string>(postData["popUIDs"],
		 [](const Json::Value & val) {return val.asString();});
		 seqsBySession_[sesUid]->cache_.at(mipFam).reload();
		 seqsBySession_[sesUid]->cache_.at(mipFam).toggleSeqs(
		 [&popUIDs](const readObject & seq) {
		 return bib::in(seq.seqBase_.getStubName(false), popUIDs);
		 });
		 */
		//make sure seqs aren't empty, viewer doesn't know how to handle that, if it isn't make sure to remove the placeholder seq if it is there
		seqsBySession_[sesUid]->cache_.at(mipTar).ensureNonEmptyReads();
		ret = seqsBySession_[sesUid]->getJson(mipTar);
		ret["sessionUID"] = bib::json::toJson(sesUid);
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mgv::mgvErrorHandler(const int statusCode, const std::exception& exception,
		const std::shared_ptr<restbed::Session>& session) {
	std::cerr << "statusCode: " << statusCode << std::endl;
	std::cerr << exception.what() << std::endl;
	if (session->is_open()) {
		session->close(statusCode, exception.what(), { { "Server", "Restbed" } });
	}
}

}  // namespace bibseq

