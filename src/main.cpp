// Created on 2014/12/30
// main.cpp

#include "mipsterPrograms.h"

namespace bibseq {

class mipsterRunner: public bib::progutils::oneRing {
public:
	mipsterRunner();
};
mipsterRunner::mipsterRunner() :
		bib::progutils::oneRing(
				{ addRing<mipsterAnalysisRunner>(),
					addRing<mipsterServerRunner>(),
					addRing<mipsterUtilsRunner>(),
					addRing<mipsterSetUpRunner>(),
				  addRing<mipsterSimRunner>(),
					addRing<mipsterMipExplorerRunner>()}, { }, "mipster") {
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
