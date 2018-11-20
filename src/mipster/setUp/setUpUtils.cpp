/*
 * setUpUtils.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: nick
 */
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

