/*
 * setUpUtils.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */

#include "setUpUtils.hpp"

namespace bibseq {



bool checkForFileDirExistence(const bfs::path & fnp, VecStr & warnings) {
	if (!bfs::exists(fnp)) {
		warnings.emplace_back(
				"Directory or file " + bib::bashCT::boldRed(fnp.string())
						+ " doesn't exist");
		return false;
	}
	return true;
}

}  // namespace bibseq

