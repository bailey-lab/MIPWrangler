
//  mipsterMipExplorerRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterMipExplorerRunner.hpp"


namespace bibseq {

mipsterMipExplorerRunner::mipsterMipExplorerRunner() :
		bib::progutils::programRunner(
				{ addFunc("viewMipsOnGenome", viewMipsOnGenome, false) },
				"mipsterMipExplorer") {
}



int mipsterMipExplorerRunner::viewMipsOnGenome(
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



	return 0;
}



                    
} // namespace bibseq
