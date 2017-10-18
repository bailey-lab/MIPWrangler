// Created on 2014/12/30
// main.cpp

#include "mipsterPrograms.h"

namespace bibseq {

class mipsterRunner: public bib::progutils::OneRing {
public:
	mipsterRunner();
};
mipsterRunner::mipsterRunner() :
		bib::progutils::OneRing(
				{ addRing<mipsterAnalysisRunner>(),
					addRing<mipsterServerRunner>(),
					addRing<mipsterUtilsRunner>(),
					addRing<mipsterSetUpRunner>(),
				  addRing<mipsterSimRunner>(),
					addRing<mipsterMipExplorerRunner>(),
					addRing<mipsterMipTesterRunner>()}, { }, "MIPWrangler",
					"1", "0", "0-dev") {
}
} // namespace bibseq

int main(int argc, char* argv[]) {
	try {
		bibseq::mipsterRunner mRunner;
		return mRunner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
