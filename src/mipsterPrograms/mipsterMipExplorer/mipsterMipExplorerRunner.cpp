
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
	bfs::path inputDir = "";
	std::string primaryGenome = "";
	std::string selectGenomes = "";
	uint32_t numThreads = 1;
	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.setOption(primaryGenome, "--primaryGenome", "The primary genome", true);
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(mainDir, "--masterDir", "The master output directory", true);
	setUp.setOption(inputDir, "--inputDir", "The master input directory, with arm sequences and genomes", true);
	setUp.setOption(selectGenomes, "--selectGenomes", "Extract info from only these genomes, default is all genomes found");
	setUp.finishSetUp(std::cout);

	MipsOnGenome mips(mainDir, inputDir, numThreads);
	if("" != selectGenomes){
		auto genomes = tokenizeString(selectGenomes, ",");
		genomes.emplace_back(primaryGenome);
		std::set<std::string> genomeSet{genomes.begin(), genomes.end()};
		mips.setSelectedGenomes(genomeSet);
	}
	mips.loadInArms();
	mips.loadInGenomes();
	if ("" != primaryGenome) {
		mips.setPrimaryGenome(primaryGenome);
	}
	mips.setUpGenomes();
	mips.createArmFiles();
	mips.mapArmsToGenomes();
	mips.genBeds();
	mips.genFastas();

	mips.genTables();
	return 0;
}




int mipsterMipExplorerRunner::viewMipsOnGenome(
		const bib::progutils::CmdArgs & inputCommands) {
	bfs::path mainDir = "";
	bfs::path inputDir = "";
	std::string primaryGenome = "";
	uint32_t numThreads = 1;
	std::string selectGenomes = "";
	SeqAppCorePars seqServerCorePars;
	seqServerCorePars.name_ = "mgv0";
	seqServerCorePars.port_ = 10000;
	bfs::path resourceDirName = bib::files::make_path(MIPWrangler_INSTALLDIR,
			"etc/serverResources");
	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(mainDir, "--masterDir", "The master output directory directory with arm extraction results", true);
	setUp.setOption(inputDir, "--inputDir", "The master input directory, with arm sequences and genomes", true);
	setUp.setOption(primaryGenome, "--primaryGenome", "The primary genome", true);
	setUp.setOption(selectGenomes, "--selectGenomes", "Extract info from only these genomes, default is all genomes found");

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
  appConfig["inputDir"] = bib::json::toJson(inputDir);
	if("" != selectGenomes){
		auto genomes = tokenizeString(selectGenomes, ",");
		genomes.emplace_back(primaryGenome);
		std::set<std::string> genomeSet{genomes.begin(), genomes.end()};
		appConfig["selectedGenomes"] = bib::json::toJson(genomeSet);
	}

  if(setUp.pars_.verbose_){
  	std::cout << seqServerCorePars.getAddress() << std::endl;
  }

  mgv server(appConfig);
	auto resources = server.getAllResources();

	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(seqServerCorePars.port_);
	settings->set_default_header("Connection", "close");
	settings->set_worker_limit(4);
	restbed::Service service;
	service.set_error_handler(mgv::mgvErrorHandler);
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
