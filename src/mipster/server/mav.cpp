/*
 * mipAnalysisServer.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */




#include "mav.hpp"
#include "mipster/parameters.h"


namespace bibseq {





mav::mav(const Json::Value & config) :
				bibseq::SeqApp(config){

	serverResourceDir_ = bib::appendAsNeededRet(config["resources"].asString(),"/");
	masterDir_ = config["masterDir"].asString();
	mipArmsFileName_ = config["mipArmsFileName"].asString();
	mipsSamplesFile_ = config["mipsSamplesFile"].asString();
	std::string seqFileSuffix = config["seqFileSuffix"].asString();
	mipMaster_ = std::make_shared<SetUpMaster>(masterDir_.string());
	mipMaster_->setMipArmFnp(mipArmsFileName_.string());
	mipMaster_->setMipsSampsNamesFnp(mipsSamplesFile_.string());
	mipMaster_->loadMipsSampsInfo(2);
	mipMaster_->setRawDataSuffix(seqFileSuffix);
	mipMaster_->setServerName(rootName_.substr(1));

	bib::stopWatch watch;
	watch.setLapName("Extraction Info");
	{
		//set up extraction info
		masterExtractInfo_ = std::make_shared<TableCache>(
				TableIOOpts(
						InOptions(bfs::path(
								mipMaster_->getMipSerDir().string()
										+ "extractionInfo/allExtractInfo.tab.txt")), "\t", true));
		auto samplesExtracted = parseJsonForMipSamps(config["samplesExtracted"]);
		for (const auto & sampleExtracted : samplesExtracted) {
			extractInfosBySamp_.emplace(sampleExtracted.samp_,
					TableCache(TableIOOpts(InOptions(mipMaster_->pathSampleExtractInfo(sampleExtracted)),
							"\t", true)));
		}

		for (const auto & mipTar : mipMaster_->getAllMipTargets()) {
			TableIOOpts mipExtractOpt(InOptions(mipMaster_->pathMipExtractInfo(mipTar)),
									"\t", true);
			extractInfosByTar_.emplace(mipTar,
					TableCache(mipExtractOpt));
		}
	}
	watch.startNewLap("Pop Clus Info");
	{
		auto mipFamsPopClustered = parseJsonForMipSamps(
				config["mipFamsPopClustered"]);
		//pop info by mip family
		for (const auto & mipFam : mipFamsPopClustered) {
			popClusInfoByTar_.emplace(mipFam.mipFam_,
					TableCache(
							TableIOOpts(
									InOptions(
											mipMaster_->pathMipPopClusSampInfo(mipFam)),
									"\t", true)));
			popClusPopInfoByTar_.emplace(mipFam.mipFam_,
					TableCache(
							TableIOOpts(
									InOptions(mipMaster_->pathMipPopClusPopInfo(mipFam)),
									"\t", true)));
			seqs_->updateAddCache(mipFam.mipFam_,
					SeqIOOptions::genFastqIn(
							mipMaster_->pathMipPopClusHaplo(mipFam).string(), true));
		}
		//pop info by sample
		for (const auto & samp : mipMaster_->names_->samples_) {
			TableIOOpts sampPopClusOpt(
					InOptions(
							mipMaster_->pathSampPopClusSampInfo(MipFamSamp("", samp))),
					"\t", true);
			if (sampPopClusOpt.in_.inExists()) {
				popClusInfoBySample_.emplace(samp, TableCache(sampPopClusOpt));
			}
		}

		auto mipFamsSampsPopClustered = parseJsonForMipSamps(
						config["mipFamsSampsPopClustered"]);
		for (const auto & mipSamp : mipFamsSampsPopClustered) {
			seqs_->updateAddCache(mipSamp.mipFam_ + "_" + mipSamp.samp_ + "_final",
					SeqIOOptions::genFastqIn(
							mipMaster_->pathPopClusFinalHaplo(mipSamp).string(), true));
			/*seqs_->updateAddCache(mipSamp.mipFam_ + "_" + mipSamp.samp_ + "_original",
					SeqIOOptions::genFastqIn(
							mipMaster_->pathPopClusOriginalHaplo(mipSamp).string(), true));*/
		}
	}


	watch.startNewLap("URL mapping set up");
	//add html pages
	jsFiles_->addFiles(
			bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mav/js"), ".js"));
	cssFiles_->addFiles(
			bib::files::gatherFiles(bib::files::make_path(serverResourceDir_, "mav/css"), ".css"));



	addScripts(bib::files::make_path(serverResourceDir_, "mav"));


	if(debug_){
		std::cout << "Finished set up" << std::endl;
		watch.logLapTimes(std::cout, true, 6, true);
	}
}


VecStr mav::requiredOptions() const {
	return catenateVectors(super::requiredOptions(), VecStr { "masterDir",
			"samplesExtracted", "masterDir", "mipsSamplesFile", "mipArmsFileName",
			"seqFileSuffix", "samplesExtracted", "mipFamsPopClustered",
			"mipFamsSampsPopClustered" });
}

void mav::mainPageHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto body = genHtmlDoc(rootName_, pages_.at("mainPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}


void mav::getAllNamesHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto projectName = request->get_path_parameter("projectName");
	Json::Value ret;

	auto samps = mipMaster_->names_->samples_;
	bib::sort(samps);
	ret["samples"] = bib::json::toJson(samps);
	auto groupings = mipMaster_->getMipGroupings();
	bib::sort(groupings);
	ret["mipRegions"] = bib::json::toJson(groupings);
	auto famlies = mipMaster_->getAllMipFamilies();
	bib::sort(famlies);
	ret["mipFamilies"] = bib::json::toJson(famlies);
	auto targets = mipMaster_->getAllMipTargets();
	bib::sort(targets);
	ret["mipTargets"] = bib::json::toJson(targets);

	auto body = bib::json::writeAsOneLine(ret);
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(body);
	session->close(restbed::OK, body, headers);
}

void mav::redirect(std::shared_ptr<restbed::Session> session,
		std::string errorMessage) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	std::cerr << errorMessage << std::endl;
	auto body = genHtmlDoc(rootName_, pages_.at("redirectPage.js"));
	const std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateTxtHtmlHeader(body);
	session->close(restbed::OK, body, headers);
}





std::vector<std::shared_ptr<restbed::Resource>> mav::getAllResources() {
	auto ret = super::getAllResources();

	//main page
	ret.emplace_back(mainPage());
	//basic info
	ret.emplace_back(getAllNames());

	//one region
	ret.emplace_back(oneRegionPage());
	ret.emplace_back(getNamesForRegion());

	//one mip fam
	ret.emplace_back(oneMipFamPage());
	ret.emplace_back(getNamesForMipFam());

	//extraction
	//extraction per target
	ret.emplace_back(initialReadStatsPerMipTarPage());
	ret.emplace_back(samplesForExtractedMip());
	ret.emplace_back(getInitialReadStatsPerMipTar());
	//extraction per samples
	ret.emplace_back(initialReadStatsPerSamplePage());
	ret.emplace_back(extractedMipsForSample());
	ret.emplace_back(getInitialReadStatsPerSample());
	//extraction overall
	ret.emplace_back(initialReadStatsPage());
	ret.emplace_back(getInitialReadStats());

	//one sample, all mips
	ret.emplace_back(oneSampAllMipDataPage());
	ret.emplace_back(mipFamNamesForSamp());
	ret.emplace_back(getOneSampAllMipData());

	//one region, one samp
	ret.emplace_back(regionInfoForSampPage());
	ret.emplace_back(getMipOverlapGraphData());

	//one region, one samp
	ret.emplace_back(oneMipAllSampsPage());
	ret.emplace_back(getOneMipAllSampsData());
	ret.emplace_back(getOneMipPopSeqs());
	ret.emplace_back(getOneMipAllSampsPopData());

	//one mip one samp
	ret.emplace_back(oneMipOneSampDataPage());
	ret.emplace_back(getMipOneSampOriginalSeqs());
	ret.emplace_back(getMipOneSampFinalSeqs());
	ret.emplace_back(getOneMipOneSampsData());
	return ret;
}

std::shared_ptr<restbed::Resource> mav::mainPage() {
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

std::shared_ptr<restbed::Resource> mav::getAllNames(){
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





void errorHandler(const int statusCode, const std::exception& exception,
		const std::shared_ptr<restbed::Session>& session) {
	std::cerr << "statusCode: " << statusCode << std::endl;
	std::cerr << exception.what() << std::endl;
	if(session->is_open()){
		session->close(statusCode, exception.what(), { { "Server", "Restbed" } });
	}
}

int mavRunner(const bib::progutils::CmdArgs & inputCommands){
	bibseq::seqSetUp setUp(inputCommands);
	SeqAppCorePars seqServerCorePars;
	seqServerCorePars.name_ = "mip0";
	seqServerCorePars.port_ = 10000;
	bfs::path resourceDirName = bib::files::make_path(MIPWrangler_INSTALLDIR,
			"etc/serverResources");
	mipCorePars mipCorepars;
	mipCorepars.processDefaults(setUp);
	setUp.setOption(mipCorepars.numThreads, "--numThreads",
			"Number of threads to use");
	setUp.setOption(resourceDirName, "--resourceDirName",
			"Name of the resource Directory where the js and hmtl is located",
			!bfs::exists(resourceDirName));
	setUp.setOption(mipCorepars.seqFileSuffix, "--seqFileSuffix",
			"The ending of the sequence append to sample name");
	resourceDirName = bib::appendAsNeededRet(resourceDirName.string(), "/");
	setUp.processVerbose();
	setUp.processDebug();
	seqServerCorePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(mipCorepars.masterDir);
	mipMaster.setMipArmFnp(mipCorepars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(mipCorepars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(mipCorepars.allowableErrors);
	mipMaster.setServerName(seqServerCorePars.name_);
	mipMaster.setRawDataSuffix(mipCorepars.seqFileSuffix);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

  //check for html/js/css resource files
  checkExistenceThrow(resourceDirName);
  checkExistenceThrow(bib::files::make_path(resourceDirName,"mav/"));

  //prepare servers
  mipMaster.prepareMipAnalysisServer(mipCorepars.numThreads);
  //collect master table info for all extraction info
  auto samplesExtracted = mipMaster.getPairsWithExtracted(mipCorepars.numThreads);
  auto mipFamsPopClustered = mipMaster.getMipFamsWithPopClustered(mipCorepars.numThreads);
  auto mipFamsSampsPopClustered = mipMaster.getPairsWithClustered(mipCorepars.numThreads);


  Json::Value appConfig;
  seqServerCorePars.addCoreOpts(appConfig);
  mipCorepars.addCorePathsToConfig(appConfig);
  appConfig["resources"] = bib::json::toJson(resourceDirName);
  appConfig["samplesExtracted"] = bib::json::toJson(samplesExtracted);
  appConfig["mipFamsPopClustered"] = bib::json::toJson(mipFamsPopClustered);
  appConfig["mipFamsSampsPopClustered"] = bib::json::toJson(mipFamsSampsPopClustered);
  if(setUp.pars_.verbose_){
  	std::cout << seqServerCorePars.getAddress() << std::endl;
  }

  mav server(appConfig);
	auto resources = server.getAllResources();

	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(seqServerCorePars.port_);
	settings->set_default_header("Connection", "close");

	restbed::Service service;
	service.set_error_handler(errorHandler);
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

}  // namespace bibseq
