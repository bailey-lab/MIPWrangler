
//  mipsterServerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
// MIPWrangler - A library for analyzing sequence data from molecular inversion probe analysis
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of MIPWrangler.
//
// MIPWrangler is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MIPWrangler is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MIPWrangler.  If not, see <http://www.gnu.org/licenses/>.
//
    
#include "mipsterServerRunner.hpp"


namespace njhseq {

mipsterServerRunner::mipsterServerRunner()
    : njh::progutils::ProgramRunner({
																		 addFunc("mipAnalysisServerSetUp", mipAnalysisServerSetUp, false),
																		 addFunc("mav", mavRunner, false)},
                    "mipsterServer") {}//



int mipsterServerRunner::mipAnalysisServerSetUp(const njh::progutils::CmdArgs & inputCommands) {
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
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}


  //prepare servers
  mipMaster.prepareMipAnalysisServer(mipCorepars.numThreads);

	return 0;
}

                    
} // namespace njhseq
