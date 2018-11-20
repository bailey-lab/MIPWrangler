/*
 * utils.cpp
 *
 *  Created on: Feb 8, 2016
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

#include "utils.hpp"

namespace njhseq {

void checkExistenceThrow(const bfs::path & dirName, const std::string & funcName){
	if(!bfs::exists(dirName)){
		std::stringstream ss;
		ss << "Error in : " << funcName << std::endl;
		ss << dirName << " needs to be created before running " << funcName << std::endl;
		throw std::runtime_error{ss.str()};
	}
}

void checkExistenceThrow(const bfs::path & dirName){
	if(!bfs::exists(dirName)){
		std::stringstream ss;
		ss << "Error" << std::endl;
		ss << dirName << " needs to exist "<< std::endl;
		throw std::runtime_error{ss.str()};
	}
}

bool requireExternalProgramThrow(const std::string & program){
	auto hasProgram = njh::sys::hasSysCommand(program);
	if (!hasProgram) {
		std::stringstream ss;
		ss << njh::bashCT::boldBlack(program)
				<< njh::bashCT::boldRed(
						" is not in path or may not be executable, cannot be used")
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	return hasProgram;
}




}  // namespace njhseq


