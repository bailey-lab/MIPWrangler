
//  mipsterMipExplorerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterMipExplorerRunner.hpp"


namespace njhseq {

mipsterMipExplorerRunner::mipsterMipExplorerRunner() :
		njh::progutils::ProgramRunner(
				{ addFunc("viewMipsOnGenome", viewMipsOnGenome, false),
					addFunc("setUpViewMipsOnGenome", setUpViewMipsOnGenome, false),
					addFunc("mipsAgainstHaplotypes", mipsAgainstHaplotypes, false),
				},//
				"mipsterMipExplorer") {
}//


int mipsterMipExplorerRunner::mipsAgainstHaplotypes(
		const njh::progutils::CmdArgs & inputCommands) {
	MipsOnGenome::pars inputPars;
	bfs::path mipArmsFnp = "";
	std::string selectedMips = "";
	mipsterMipExplorerSetUp setUp(inputCommands);

	setUp.processDebug();
	setUp.processVerbose();
	setUp.processComparison(inputPars.allowableError);
	setUp.processReadInNames({"--fasta"}, true);
	setUp.setOption(inputPars.gMapperPars_.numThreads_, "--numThreads", "Number of Threads");
	setUp.setOption(selectedMips, "--selectedMips", "Selected Mips");
	setUp.setOption(mipArmsFnp, "--mipArmsFnp", "Mip Arms Fnp", true);
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);
	njh::stopWatch watch;
	watch.setLapName("Initial");
	setUp.startARunLog(setUp.pars_.directoryName_);
	auto genomeDir = njh::files::make_path(setUp.pars_.directoryName_,"genomes");
	auto infoDir = njh::files::make_path(setUp.pars_.directoryName_,"info");
	njh::files::makeDir(njh::files::MkdirPar{genomeDir});
	njh::files::makeDir(njh::files::MkdirPar{infoDir});
	auto seqOutOpts = SeqIOOptions::genFastaOut(njh::files::make_path(genomeDir, "inputSeqs.fasta"));

	auto mipArms = std::make_unique<MipCollection>(mipArmsFnp, 6);
	auto selectedMipsCon = getInputValues(selectedMips, ",");
	OutOptions mipArmsOpts(njh::files::make_path(infoDir, "mip_arms.tab.txt"));
	{
		OutputStream mipArmsOut(mipArmsOpts);
		mipArmsOut << njh::conToStr(Mip::writeInfoLineHeader(), "\t") << std::endl;
		for(const auto & m : mipArms->mips_){
			if(selectedMipsCon.empty() || selectedMipsCon.front() == "" || njh::in(m.first, selectedMipsCon)){
				m.second.writeInfoLine(mipArmsOut);
			}
		}
	}

	inputPars.gMapperPars_.primaryGenome_ = "inputSeqs";
	inputPars.mainDir = setUp.pars_.directoryName_;
	inputPars.inputDir = setUp.pars_.directoryName_;
	inputPars.gMapperPars_.genomeDir_ = njh::files::make_path(inputPars.inputDir, "genomes");
	inputPars.gMapperPars_.gffDir_ = njh::files::make_path(inputPars.inputDir, "info/gff");


	SeqInput reader(setUp.pars_.ioOptions_);
	auto seqs = reader.readAllReads<seqInfo>();
	SeqOutput::write(seqs, seqOutOpts);


	MipsOnGenome mips(inputPars);
	watch.startNewLap("loadInArms");
	mips.loadInArms();
	watch.startNewLap("createArmFiles");
	mips.createArmFiles();
	watch.startNewLap("loadInGenomes");
	mips.gMapper_.loadInGenomes();
	watch.startNewLap("setUpGenomes");
	mips.gMapper_.setUpGenomes();
	watch.startNewLap("mapArmsToGenomes");
	mips.mapArmsToGenomes();
	watch.startNewLap("genBeds");
	mips.genBeds(inputPars.allowableError);
	if(setUp.pars_.debug_){
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return 0;
}

int mipsterMipExplorerRunner::setUpViewMipsOnGenome(
		const njh::progutils::CmdArgs & inputCommands) {

	MipsOnGenome::pars inputPars;

	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.processComparison(inputPars.allowableError);
	setUp.setOption(inputPars.gMapperPars_.primaryGenome_, "--primaryGenome", "The primary genome", true);
	setUp.setOption(inputPars.gMapperPars_.numThreads_, "--numThreads", "Number of Threads");
	setUp.setOption(inputPars.mipArmsFnp, "--mipArmsFnp", "Mip Arms Fnp", true);
	setUp.setOption(inputPars.removeBeds, "--removeBeds", "Whether to remove the bed file to do a re-extraction or not, usefull if changing number of errors allowed");
	setUp.setOption(inputPars.mainDir, "--masterDir", "The master output directory", true);
	setUp.setOption(inputPars.inputDir, "--inputDir", "The master input directory, with arm sequences and genomes", true);
	std::string selectGenomes = "";
	setUp.setOption(selectGenomes, "--selectGenomes", "Extract info from only these genomes, default is all genomes found");
	if("" != selectGenomes){
		auto selectGenomesToks = tokenizeString(selectGenomes, ",");
		inputPars.gMapperPars_.selectedGenomes_ = std::set<std::string>(selectGenomesToks.begin(), selectGenomesToks.end());
	}
	inputPars.gMapperPars_.genomeDir_ = njh::files::make_path(inputPars.inputDir, "genomes");
	inputPars.gMapperPars_.gffDir_ = njh::files::make_path(inputPars.inputDir, "info/gff");
	setUp.finishSetUp(std::cout);

	njh::stopWatch watch;
	watch.setLapName("Initial");
	MipsOnGenome mips(inputPars);
	setUp.startARunLog(mips.logDir_.string());

	watch.startNewLap("loadInArms");
	mips.loadInArms();
	watch.startNewLap("createArmFiles");
	mips.createArmFiles();
	watch.startNewLap("loadInGenomes");
	mips.gMapper_.loadInGenomes();
	watch.startNewLap("setUpGenomes");
	mips.gMapper_.setUpGenomes();
	watch.startNewLap("mapArmsToGenomes");
	mips.mapArmsToGenomes();
	if(inputPars.removeBeds){
		njh::files::rmDirForce(mips.bedsDir_);
		njh::files::makeDir(njh::files::MkdirPar(mips.bedsDir_));
		njh::files::makeDir(njh::files::MkdirPar(mips.bedsPerGenomeDir_));

	}
	watch.startNewLap("genBeds");
	mips.genBeds(inputPars.allowableError);
	watch.startNewLap("genFastas");
	mips.genFastas();
	watch.startNewLap("genTables");
	mips.genTables();

	if(setUp.pars_.debug_){
		watch.logLapTimes(std::cout, true, 6, true);
	}
	return 0;
}




int mipsterMipExplorerRunner::viewMipsOnGenome(
		const njh::progutils::CmdArgs & inputCommands) {
	bfs::path mainDir = "";
	bfs::path inputDir = "";
	std::string primaryGenome = "";
	uint32_t numThreads = 1;
	std::string selectGenomes = "";bfs::path mipArmsFnp = "";
	SeqAppCorePars seqServerCorePars;
	seqServerCorePars.name_ = "mgv0";
	seqServerCorePars.port_ = 10000;
	bfs::path resourceDirName = njh::files::make_path(MIPWrangler_INSTALLDIR,
			"etc/serverResources");
	mipsterMipExplorerSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(numThreads, "--numThreads", "Number of Threads");
	setUp.setOption(mainDir, "--masterDir", "The master output directory directory with arm extraction results", true);
	setUp.setOption(inputDir, "--inputDir", "The master input directory, with arm sequences and genomes", true);
	setUp.setOption(primaryGenome, "--primaryGenome", "The primary genome", true);
	setUp.setOption(selectGenomes, "--selectGenomes", "Extract info from only these genomes, default is all genomes found");
	setUp.setOption(mipArmsFnp, "--mipArmsFnp", "Mip Arms Fnp", true);

	setUp.setOption(resourceDirName, "--resourceDirName",
			"Name of the resource Directory where the js and hmtl is located",
			!bfs::exists(resourceDirName));
	resourceDirName = njh::appendAsNeededRet(resourceDirName.string(), "/");
	seqServerCorePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

  //check for html/js/css resource files
  checkExistenceThrow(resourceDirName);
  checkExistenceThrow(njh::files::make_path(resourceDirName,"mgv/"));


  Json::Value appConfig;
  seqServerCorePars.addCoreOpts(appConfig);
  appConfig["resources"] = njh::json::toJson(resourceDirName);
  appConfig["primaryGenome"] = njh::json::toJson(primaryGenome);
  appConfig["masterDir"] = njh::json::toJson(mainDir);
  appConfig["inputDir"] = njh::json::toJson(inputDir);
  appConfig["mipArmsFnp"] = njh::json::toJson(mipArmsFnp);

	if("" != selectGenomes){
		auto genomes = tokenizeString(selectGenomes, ",");
		genomes.emplace_back(primaryGenome);
		std::set<std::string> genomeSet{genomes.begin(), genomes.end()};
		appConfig["selectedGenomes"] = njh::json::toJson(genomeSet);
	}

  if(setUp.pars_.verbose_){
  	std::cout << seqServerCorePars.getAddress() << std::endl;
  }

  mgv server(appConfig);
	auto resources = server.getAllResources();

	auto settings = std::make_shared<restbed::Settings>();
	settings->set_port(seqServerCorePars.port_);
	settings->set_default_header("Connection", "close");
	settings->set_bind_address(seqServerCorePars.bindAddress_);
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



                    
} // namespace njhseq
