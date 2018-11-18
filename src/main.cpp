// Created on 2014/12/30
// main.cpp

#include "mipsterPrograms.h"



int main(int argc, char* argv[]) {
	try {
		njhseq::mipsterRunner mRunner;
		return mRunner.run(argc, argv);
	} catch (std::exception & e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
}
