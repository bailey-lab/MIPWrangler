
//  mipsterSetUpRunner.cpp
//
//  Created by Nick Hathaway on 2016/02/13.
//  Copyright (c) 2016 Nick Hathaway. All rights reserved.
//

    
#include "mipsterSetUpRunner.hpp"
    
namespace njhseq {

mipsterSetUpRunner::mipsterSetUpRunner()
    : njh::progutils::ProgramRunner({
																		 addFunc("makeMipPopClusDirectories", makeMipPopClusDirectories, true),
																		 addFunc("checkDirectoryStructure", checkDirectoryStructure, false),
																		 addFunc("createSkeletonDirectoryStructure", createSkeletonDirectoryStructure, true),
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 addFunc("checkForRawData", checkForRawData, true),
																		 addFunc("checkForExtracted", checkForExtracted, true),
																		 addFunc("checkForBarCor", checkForBarCor, true),
																		 addFunc("checkForClustered", checkForClustered, true),
																		 addFunc("cloneAnalysisDirectory", cloneAnalysisDirectory, false)},
                    "MIPWranglerSetUp") {}//



int mipsterSetUpRunner::cloneAnalysisDirectory(const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	std::string dir = "";
	setUp.setOption(dir, "--newMasterDir",
			"The new master directory for a clone of the current set up", true);
	setUp.setOption(pars.seqFileSuffix, "--seqFileSuffix",
			"The ending of the sequence append to sample name");
	pars.processDefaults(setUp);
	setUp.finishSetUp(std::cout);

	SetUpMaster newMipMaster(dir);
	newMipMaster.createDirStructSkeleton(pars.mipsSamplesFile,
			pars.mipArmsFileName);
	newMipMaster.setRawDataSuffix(pars.seqFileSuffix);

	SetUpMaster oldMipMaster(pars.masterDir);
	oldMipMaster.setRawDataSuffix(pars.seqFileSuffix);

	for (const auto & samp : newMipMaster.names_->samples_) {
		bfs::create_symlink(bfs::canonical(oldMipMaster.pathSampleRawData(MipFamSamp("", samp))),
				newMipMaster.pathSampleRawData(MipFamSamp("", samp)));
	}
	return 0;
}

int mipsterSetUpRunner::checkDirectoryStructure(const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	setUp.setOption(pars.masterDir, "--masterDir", "The master directory", true);
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(pars.masterDir);
	auto warnings = mipMaster.checkDirStruct();
	if(warnings.empty()){
		std::cout << njh::bashCT::boldGreen("PASS!") << std::endl;
	}else{
		std::cout << njh::bashCT::boldRed("FAIL!") << std::endl;
		std::cout << "Warnings:" << std::endl;
		std::cout << njh::conToStr(warnings, "\n") << std::endl;
	}
	return 0;
}


int mipsterSetUpRunner::createSkeletonDirectoryStructure(const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	setUp.setOption(pars.masterDir, "--masterDir", "The master directory", true);
	setUp.setOption(pars.mipArmsFileName, "--mipArmsFilename",
			"Name of the mip arms file");
	setUp.setOption(pars.mipsSamplesFile, "--mipSampleFile",
			"Mip sample filename, two columns, one column is mips, other is samples");
	setUp.finishSetUp(std::cout);

	SetUpMaster mipMaster(pars.masterDir);
	if ("" == pars.mipArmsFileName && "" == pars.mipsSamplesFile) {
		mipMaster.createDirStructSkeleton();
	} else {
		mipMaster.createDirStructSkeleton(pars.mipsSamplesFile, pars.mipArmsFileName);
	}
	return 0;
}
                    

int mipsterSetUpRunner::makeMipPopClusDirectories(const njh::progutils::CmdArgs & inputCommands) {
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	pars.numThreads = 1;
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use to more quickly make the directories");
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.createPopClusMipDirs(pars.numThreads);
  return 0;
}



int mipsterSetUpRunner::checkForRawData(
		const njh::progutils::CmdArgs & inputCommands) {
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	pars.processDefaults(setUp);
	setUp.setOption(pars.seqFileSuffix, "--seqFileSuffix",
			"The ending of the sequence append to sample name");
	pars.numThreads = 1;
	setUp.setOption(pars.numThreads, "--numThreads",
			"Number of threads to use to more quickly make the directories");
	setUp.processDebug();
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setRawDataSuffix(pars.seqFileSuffix);
	auto pairs = mipMaster.getSamplesWithRawData(pars.numThreads);
	printMipSampVec(pairs);
	return 0;
}

int mipsterSetUpRunner::checkForExtracted(const njh::progutils::CmdArgs & inputCommands) {
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	pars.processDefaults(setUp);
	pars.numThreads = 1;
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use to more quickly make the directories");
	setUp.processDebug();
	setUp.processVerbose();
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto pairs = mipMaster.getPairsWithExtracted(pars.numThreads);
	printMipSampVec(pairs);
  return 0;
}

int mipsterSetUpRunner::checkForBarCor(const njh::progutils::CmdArgs & inputCommands) {
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();
	pars.numThreads = 1;
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use to more quickly make the directories");
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto pairs = mipMaster.getPairsWithBarCor(pars.numThreads);
	printMipSampVec(pairs);
  return 0;
}

int mipsterSetUpRunner::checkForClustered(const njh::progutils::CmdArgs & inputCommands) {
	mipsterSetUpSetUp setUp(inputCommands);
	mipCorePars pars;
	pars.processDefaults(setUp);
	setUp.processDebug();
	setUp.processVerbose();
	pars.numThreads = 1;
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to use to more quickly make the directories");
	setUp.finishSetUp(std::cout);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);

	auto pairs = mipMaster.getPairsWithClustered(pars.numThreads);
	printMipSampVec(pairs);
  return 0;
}


                    
} // namespace njhseq
