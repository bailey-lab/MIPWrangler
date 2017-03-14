
//  mipsterMipTesterRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterMipTesterRunner.hpp"


namespace bibseq {

mipsterMipTesterRunner::mipsterMipTesterRunner() :
		bib::progutils::programRunner(
				{
					addFunc("testingVariationCalling", testingVariationCalling, false)
				},
				"mipsterMipTester") {
}//







int mipsterMipTesterRunner::testingVariationCalling(
		const bib::progutils::CmdArgs & inputCommands) {
	mipCorePars pars;
	mipsterMipTesterSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.processDefaults(setUp);
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to utilize");
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	auto finalResultPairs = mipMaster.getPairsWithPopClustered(pars.numThreads);



	return 0;

}



                    
} // namespace bibseq
