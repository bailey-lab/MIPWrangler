/*
 * setUpUtils.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */

#include "setUpUtils.hpp"

namespace njhseq {



bool checkForFileDirExistence(const bfs::path & fnp, VecStr & warnings) {
	if (!bfs::exists(fnp)) {
		warnings.emplace_back(
				"Directory or file " + njh::bashCT::boldRed(fnp.string())
						+ " doesn't exist");
		return false;
	}
	return true;
}

}  // namespace njhseq

