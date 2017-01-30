
//  mipsterMipExplorerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterMipExplorerRunner.hpp"


namespace bibseq {

mipsterMipExplorerRunner::mipsterMipExplorerRunner() :
		bib::progutils::programRunner(
				{ addFunc("viewMipsOnGenome", viewMipsOnGenome, false),
					addFunc("setUpViewMipsOnGenome", setUpViewMipsOnGenome, false),
				},
				"mipsterMipExplorer") {
}//



int mipsterMipExplorerRunner::setUpViewMipsOnGenome(
		const bib::progutils::CmdArgs & inputCommands) {
	bfs::path mainDir = "";
	uint32_t numThreads = 1;
	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(mainDir, "--masterDir", "The master directory", true);
	setUp.finishSetUp(std::cout);

	MipsOnGenome mips(mainDir, numThreads);

	mips.loadInArms();
	mips.loadInGenomes();

	mips.setUpGenomes();
	mips.createArmFiles();
	mips.mapArmsToGenomes();
	mips.genBeds();
	mips.genFastas();



	return 0;
}

class mgv: public bibseq::SeqApp {
	typedef bibseq::SeqApp super;

	bfs::path mainDir_;
	bfs::path serverResourceDir_;
	std::unique_ptr<MipsOnGenome> mipsInfo_;

public:

	mgv(const Json::Value & config) :
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
		/*
		 jsFiles_->addFiles(
		 bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mav/js"), ".js"));
		 cssFiles_->addFiles(
		 bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mav/css"), ".css"));
		 */
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

	virtual std::vector<std::shared_ptr<restbed::Resource>> getAllResources(){
		auto ret = super::getAllResources();
		//main page
		ret.emplace_back(mainPage());
		//basic info
		ret.emplace_back(getAllNames());

		//show mip sequences
		ret.emplace_back(getMipSeqs());
		ret.emplace_back(showMipSeqs());

		return ret;
	}

	virtual VecStr requiredOptions() const{
		return concatVecs(super::requiredOptions(), VecStr { "masterDir",
				"numThreads", "primaryGenome", "resources"});
	}
	void mainPageHandler(std::shared_ptr<restbed::Session> session){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto body = genHtmlDoc(rootName_, pages_.at("mainPage.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}


	void getAllNamesHandler(std::shared_ptr<restbed::Session> session){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto request = session->get_request();
		Json::Value ret;

		auto genomes = mipsInfo_->getGenomes();
		bib::sort(genomes);
		ret["genomes"] = bib::json::toJson(genomes);
		auto mips = mipsInfo_->getMips();
		bib::sort(mips);
		ret["mips"] = bib::json::toJson(mips);
		auto body = bib::json::writeAsOneLine(ret);
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateAppJsonHeader(body);
		session->close(restbed::OK, body, headers);
	}

	void redirect(std::shared_ptr<restbed::Session> session,
			std::string errorMessage) {
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		std::cerr << errorMessage << std::endl;
		auto body = genHtmlDoc(rootName_, pages_.at("redirectPage.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}

	std::shared_ptr<restbed::Resource> getAllNames(){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto resource = std::make_shared<restbed::Resource>();
		resource->set_path(
				UrlPathFactory::createUrl( { { rootName_ }, { "getAllNames" }}));
		resource->set_method_handler("GET",
				std::function<void(std::shared_ptr<restbed::Session>)>(
						[this](std::shared_ptr<restbed::Session> session) {
							getAllNamesHandler(session);
						}));
		return resource;
	}

	std::shared_ptr<restbed::Resource> mainPage(){
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

	std::shared_ptr<restbed::Resource> showMipSeqs(){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto resource = std::make_shared<restbed::Resource>();
		auto path =  UrlPathFactory::createUrl( { { rootName_ },{"showMipSeqs"}, {"mip", UrlPathFactory::pat_wordNumsDash_} });
		std::cout << path << std::endl;
		resource->set_path(UrlPathFactory::createUrl( { { rootName_ },{"showMipSeqs"} ,{"mip", UrlPathFactory::pat_wordNumsDash_} }));
		resource->set_method_handler("GET",
				std::function<void(std::shared_ptr<restbed::Session>)>(
						[this](std::shared_ptr<restbed::Session> session) {
							showMipSeqsHandler(session);
						}));
		return resource;
	}

	void showMipSeqsHandler(std::shared_ptr<restbed::Session> session){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		const auto request = session->get_request();
		auto mipTar = request->get_path_parameter("mip");
		if (!bib::in(mipTar, mipsInfo_->mipArms_->mips_)) {
			std::stringstream ss;
			ss << "Mip Name wasn't found: " << mipTar << "\n";
			ss << "options are " << bib::conToStr(bib::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ") << "\n";
			ss << "Redirecting....";
			redirect(session, ss.str());
		} else {
			auto body = genHtmlDoc(rootName_, pages_.at("showMipSeqs.js"));
			const std::multimap<std::string, std::string> headers =
					HeaderFactory::initiateTxtHtmlHeader(body);
			session->close(restbed::OK, body, headers);
		}
	}

	std::shared_ptr<restbed::Resource> getMipSeqs(){
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto resource = std::make_shared<restbed::Resource>();
		resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
				"getMipSeqs" },
				{ "mip", UrlPathFactory::pat_wordNumsDash_ } }));
		resource->set_method_handler("POST",
				std::function<void(std::shared_ptr<restbed::Session>)>(
						[this](std::shared_ptr<restbed::Session> session) {
							getMipSeqsHandler(session);
						}));
		return resource;
	}

	void getMipSeqsHandler(
			std::shared_ptr<restbed::Session> session) {
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

	void getMipSeqsPostHandler(
			std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
		auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
		auto request = session->get_request();
		auto mipTar = request->get_path_parameter("mip");
		const auto postData = bib::json::parse(std::string(body.begin(), body.end()));
		Json::Value ret;
		if (!bib::in(mipTar, mipsInfo_->mipArms_->mips_)) {
			std::cerr << __PRETTY_FUNCTION__ << " error, no matching mip to " << mipTar << std::endl;
			std::cerr << "options are " << bib::conToStr(bib::getVecOfMapKeys(mipsInfo_->mipArms_->mips_), ", ") << "\n";
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


};

void mgvErrorHandler(const int statusCode, const std::exception& exception,
		const std::shared_ptr<restbed::Session>& session) {
	std::cerr << "statusCode: " << statusCode << std::endl;
	std::cerr << exception.what() << std::endl;
	if(session->is_open()){
		session->close(statusCode, exception.what(), { { "Server", "Restbed" } });
	}
}

int mipsterMipExplorerRunner::viewMipsOnGenome(
		const bib::progutils::CmdArgs & inputCommands) {
	bfs::path mainDir = "";
	std::string primaryGenome = "";
	uint32_t numThreads = 1;

	SeqAppCorePars seqServerCorePars;
	seqServerCorePars.name_ = "mgv0";
	seqServerCorePars.port_ = 10000;
	bfs::path resourceDirName = bib::files::make_path(MIPWrangler_INSTALLDIR,
			"etc/serverResources");
	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(mainDir, "--masterDir", "The master directory", true);
	setUp.setOption(primaryGenome, "--primaryGenome", "The primary genome", true);
	setUp.setOption(resourceDirName, "--resourceDirName",
			"Name of the resource Directory where the js and hmtl is located",
			!bfs::exists(resourceDirName));
	resourceDirName = bib::appendAsNeededRet(resourceDirName.string(), "/");
	seqServerCorePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

  //check for html/js/css resource files
  checkExistenceThrow(resourceDirName);
  checkExistenceThrow(bib::files::make_path(resourceDirName,"mgv/"));


  Json::Value appConfig;
  seqServerCorePars.addCoreOpts(appConfig);
  appConfig["resources"] = bib::json::toJson(resourceDirName);
  appConfig["primaryGenome"] = bib::json::toJson(primaryGenome);
  appConfig["masterDir"] = bib::json::toJson(mainDir);
  if(setUp.pars_.verbose_){
  	std::cout << seqServerCorePars.getAddress() << std::endl;
  }

  mgv server(appConfig);
	auto resources = server.getAllResources();

	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(seqServerCorePars.port_);
	settings->set_default_header("Connection", "close");

	restbed::Service service;
	service.set_error_handler(mgvErrorHandler);
	for(const auto & resource : resources){
		service.publish(resource);
	}
	try {
		service.start(settings);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
	}


	return 0;

}



                    
} // namespace bibseq
