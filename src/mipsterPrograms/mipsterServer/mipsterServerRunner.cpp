
//  mipsterServerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterServerRunner.hpp"


namespace bibseq {

mipsterServerRunner::mipsterServerRunner()
    : bib::progutils::ProgramRunner({
																		 addFunc("mipAnalysisServerSetUp", mipAnalysisServerSetUp, false),
																		 addFunc("mav", mavRunner, false)},
                    "mipsterServer") {}//



int mipsterServerRunner::mipAnalysisServerSetUp(const bib::progutils::CmdArgs & inputCommands) {
	SeqAppCorePars seqServerCorePars;
	mipCorePars mipCorepars;

	mipsterServerSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();

	seqServerCorePars.name_ = "mip0";
	seqServerCorePars.port_ = 10000;

	mipCorepars.processDefaults(setUp);
	setUp.setOption(seqServerCorePars.name_, "--name", "Name of the server", true);
	setUp.setOption(mipCorepars.numThreads, "--numThreads",
			"Number of threads to use");
	setUp.setOption(mipCorepars.seqFileSuffix, "--seqFileSuffix",
			"The ending of the sequence append to sample name");

	//seqServerCorePars.setCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(mipCorepars.masterDir);
	mipMaster.setMipArmFnp(mipCorepars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(mipCorepars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(mipCorepars.allowableErrors);
	mipMaster.setServerName(seqServerCorePars.name_);
	mipMaster.setRawDataSuffix(mipCorepars.seqFileSuffix);
	if("" != mipCorepars.sampleMetaFnp){
		mipMaster.setMetaData(mipCorepars.sampleMetaFnp);
	}
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}


  //prepare servers
  mipMaster.prepareMipAnalysisServer(mipCorepars.numThreads);

	return 0;
}

                    
} // namespace bibseq
