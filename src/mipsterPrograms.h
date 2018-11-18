#pragma once

// Created on 2014/12/30
// Including headers in programs


#include "mipsterPrograms/mipsterAnalysis.h"
#include "mipsterPrograms/mipsterServer.h"
#include "mipsterPrograms/mipsterUtils.h"
#include "mipsterPrograms/mipsterSim.h"
#include "mipsterPrograms/mipsterSetUp.h"
#include "mipsterPrograms/mipsterMipExplorer.h"
#include "mipsterPrograms/mipsterMipTester.h"


namespace njhseq {

class mipsterRunner: public njh::progutils::OneRing {
public:
	mipsterRunner();
};
mipsterRunner::mipsterRunner() :
		njh::progutils::OneRing(
				{ addRing<mipsterAnalysisRunner>(),
					addRing<mipsterServerRunner>(),
					addRing<mipsterUtilsRunner>(),
					addRing<mipsterSetUpRunner>(),
				  addRing<mipsterSimRunner>(),
					addRing<mipsterMipExplorerRunner>(),
					addRing<mipsterMipTesterRunner>()}, { }, "MIPWrangler",
					"1", "0", "0-dev") {
}
} // namespace njhseq

