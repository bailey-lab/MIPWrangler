/*
 * mgv.cpp
 *
 *  Created on: Jan 30, 2017
 *      Author: nick
 */

#include "mgv.hpp"


#include "mipster/mipUtils.h"

#include <elucidator/objects/BioDataObject.h>

namespace njhseq {
mgv::mgv(const Json::Value & config) :
		SeqApp(config) {
	mainDir_ = config["masterDir"].asString();
	serverResourceDir_ = njh::appendAsNeededRet(config["resources"].asString(),
			"/");
	MipsOnGenome::pars mogPars;
	mogPars.mainDir = mainDir_;
	mogPars.inputDir = config["inputDir"].asString();
	mogPars.gMapperPars_.genomeDir_ = njh::files::make_path(mogPars.inputDir, "genomes");
	mogPars.gMapperPars_.gffDir_ = njh::files::make_path(mogPars.inputDir, "info/gff");

	mogPars.gMapperPars_.numThreads_= config["numThreads"].asUInt();
	if (config.isMember("selectedGenomes")) {
		auto selectedGenomes = njh::json::jsonArrayToSet<std::string>(
				config["selectedGenomes"],
				[](const Json::Value & val) {return val.asString();});

		mogPars.gMapperPars_.selectedGenomes_= selectedGenomes;

	}
	if ("" != config["mipArmsFnp"].asString()) {
		mogPars.mipArmsFnp = config["mipArmsFnp"].asString();
	}

	mogPars.gMapperPars_.primaryGenome_ = config["primaryGenome"].asString();
	mipsInfo_ = std::make_unique<MipsOnGenome>(mogPars);

	mipsInfo_->loadInArms();
	mipsInfo_->gMapper_.loadInGenomes();
	mipsInfo_->gMapper_.setUpGenomes();


	jsFiles_->addFiles(
			njh::files::gatherFiles(
					njh::files::make_path(serverResourceDir_, "mgv/js"), ".js"));

	cssFiles_->addFiles(
			njh::files::gatherFiles(
					njh::files::make_path(serverResourceDir_, "mgv/css"), ".css"));

	addScripts(njh::files::make_path(serverResourceDir_, "mgv"));

	//register mip seqs
	for (const auto & mipName : mipsInfo_->getMips()) {
		auto mipFaPath = mipsInfo_->pathToMipFasta(mipName);
		if (bfs::exists(mipFaPath)) {
			seqs_->updateAddCache(mipName,
					SeqIOOptions::genFastaIn(mipFaPath, false));
		}
	}
	//add all tar info table
	allTarInfo_ = std::make_unique<TableCache>(
			TableIOOpts(InOptions(mipsInfo_->pathToAllInfoPrimaryGenome()), "\t",
					true));
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

	ret.emplace_back(getMipGenomeLocs());

	ret.emplace_back(getMipRegionInfoForGenome());
	ret.emplace_back(getMipTarInfoForGenome());


	ret.emplace_back(getInfoMipArmsInfo());
	ret.emplace_back(getAllMipTarInfoAllGenomes());

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
	njh::sort(genomes);
	ret["genomes"] = njh::json::toJson(genomes);
	auto mips = mipsInfo_->getMips();
	MipNameSorter::sort(mips);
	ret["mips"] = njh::json::toJson(mips);
	auto regions = mipsInfo_->mipArms_->getMipRegions();
	MipNameSorter::sort(regions, MipNameSorter::regionNamePat);
	ret["mipRegions"] = njh::json::toJson(regions);

	ret["primaryGenome"] = njh::json::toJson(mipsInfo_->gMapper_.pars_.primaryGenome_);

	auto body = njh::json::writeAsOneLine(ret);
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
	if (njh::in(regionName, regions)) {
		VecStr mipNames = mipsInfo_->mipArms_->getMipTarsForRegion(regionName);
		MipNameSorter::sort(mipNames);
		ret["mips"] = njh::json::toJson(mipNames);
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

std::shared_ptr<restbed::Resource> mgv::getMipGenomeLocs() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getMipGenomeLocs" }, { "mipTar", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipGenomeLocsHandler(session);
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
	if (njh::in(regionName, regions)) {
		auto body = genHtmlDoc(rootName_, pages_.at("showRegionInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << "Mip Name wasn't found: " << regionName << "\n";
		ss << "options are "
				<< njh::conToStr(regions, ", ")
				<< "\n";
		ss << "Redirecting....";
		redirect(session, ss.str());
	}

}

void mgv::showMipSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto mipTar = request->get_path_parameter("mip");
	if (!njh::in(mipTar, mipsInfo_->mipArms_->mips_)) {
		std::stringstream ss;
		ss << "Mip Name wasn't found: " << mipTar << "\n";
		ss << "options are "
				<< njh::conToStr(njh::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ")
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

std::shared_ptr<restbed::Resource> mgv::getMipRegionInfoForGenome() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getMipRegionInfoForGenome" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getMipRegionInfoForGenomeHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::getMipTarInfoForGenome() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getMipTarInfoForGenome" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getMipTarInfoForGenomeHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::getInfoMipArmsInfo() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getInfoMipArmsInfo" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getInfoMipArmsInfoHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mgv::getAllMipTarInfoAllGenomes() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ },
			{ "getAllMipTarInfoAllGenomes" } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
		getAllMipTarInfoAllGenomesHandler(session);
					}));
	return resource;
}







void mgv::getMipSeqsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getMipSeqsPostHandler(ses, body);
					}));
}





void mgv::getMipGenomeLocsHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	auto mipTar = request->get_path_parameter("mipTar");

	Json::Value ret;
	if (njh::in(mipTar, mipsInfo_->mipArms_->mips_)) {
		auto locs = mipsInfo_->getGenomeLocsForMipTar(mipTar);
		//locs.sortTable("genome", false);
		ret = tableToJsonByRow(locs, "genome");

	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< mipTar << "\n";
		std::cerr << "Options are " << njh::conToStr(mipsInfo_->mipArms_->getMipTars(), ", ") << "\n";
	}
	auto body = njh::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void mgv::getInfoMipArmsInfoHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	table mipArms(mipsInfo_->inputParameters_.mipArmsFnp, "\t", true);

	auto body = njh::json::writeAsOneLine(tableToJsonByRow(mipArms, "mip_family"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void mgv::getAllMipTarInfoAllGenomesHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();

	table allInfo(mipsInfo_->pathToAllInfoAllGenomes(), "\t", true);

	auto body = njh::json::writeAsOneLine(tableToJsonByRow(allInfo, "region"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}







void mgv::getGenomeLensHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
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
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "genomes" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		auto genomes = njh::json::jsonArrayToVec<std::string>(postData["genomes"], [](const Json::Value & val){ return val.asString();});
		VecStr missingGenomes;
		for(const auto & genome : genomes){
			if(!njh::in(genome, mipsInfo_->gMapper_.genomes_)){
				missingGenomes.emplace_back(genome);
			}else{
				ret[genome] = mipsInfo_->gMapper_.genomes_.at(genome)->chromosomeLengths();
			}
		}
		if(!missingGenomes.empty()){
			std::cerr << "Missing info for genomes " << njh::conToStr(missingGenomes, ", ") << std::endl;
		}
	}
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}


void mgv::getMipRegionInfoForGenomeHandler(std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
		getMipRegionInfoForGenomePostHandler(ses, body);
					}));
}




void mgv::getMipRegionInfoForGenomePostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "genome" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		std::string genome = postData["genome"].asString();
		if (!njh::in(genome, mipsInfo_->gMapper_.genomes_)) {
			std::cerr << "Missing info for genomes " << genome << std::endl;
		} else {
			auto retTab = mipsInfo_->getMipRegionStatsForGenome(genome);
			ret = tableToJsonByRow(retTab, "region");
		}
	}
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mgv::getMipTarInfoForGenomeHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getMipTarInfoForGenomePostHandler(ses, body);
					}));
}




void mgv::getMipTarInfoForGenomePostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "genome", "mipTars" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		std::string genome = postData["genome"].asString();
		auto mipTars = njh::json::jsonArrayToVec<std::string>(postData["mipTars"], [](const Json::Value & val){return val.asString();});
		if (!njh::in(genome, mipsInfo_->gMapper_.genomes_)) {
			std::cerr << "Missing info for genomes " << genome << std::endl;
		} else {
			VecStr missingTars;
			VecStr presentMipTars;
			for(const auto & mipTar : mipTars){
				if(!njh::in(mipTar, mipsInfo_->mipArms_->mips_)){
					missingTars.emplace_back(mipTar);
				}else{
					presentMipTars.emplace_back(mipTar);
				}
			}
			if(!missingTars.empty()){
				std::cerr << "Missing info for mip targets " << njh::conToStr(getVectorOfMapKeys(mipsInfo_->mipArms_->mips_), ", ") << std::endl;
			}
			MipNameSorter::sort(presentMipTars);
			//auto retTab = mipsInfo_->getMipTarStatsForGenome(genome, presentMipTars);
			auto retTab =
					allTarInfo_->get().extractByComp("target",
							[&presentMipTars](const std::string & row) {return njh::in(row, presentMipTars);});
			ret = tableToJsonByRow(retTab, "region");
		}
	}
	auto retBody = njh::json::writeAsOneLine(ret);
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
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	Json::Value ret;
	if (!njh::in(mipTar, mipsInfo_->mipArms_->mips_)) {
		std::cerr << __PRETTY_FUNCTION__ << " error, no matching mip to " << mipTar
				<< std::endl;
		std::cerr << "options are "
				<< njh::conToStr(njh::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ")
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
		 VecStr popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
		 [](const Json::Value & val) {return val.asString();});
		 seqsBySession_[sesUid]->cache_.at(mipFam).reload();
		 seqsBySession_[sesUid]->cache_.at(mipFam).toggleSeqs(
		 [&popUIDs](const readObject & seq) {
		 return njh::in(seq.seqBase_.getStubName(false), popUIDs);
		 });
		 */
		//make sure seqs aren't empty, viewer doesn't know how to handle that, if it isn't make sure to remove the placeholder seq if it is there
		seqsBySession_[sesUid]->cache_.at(mipTar).ensureNonEmptyReads();
		ret = seqsBySession_[sesUid]->getJson(mipTar);
		ret["sessionUID"] = njh::json::toJson(sesUid);
	}
	auto retBody = njh::json::writeAsOneLine(ret);
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

}  // namespace njhseq

