/*
 * utils.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


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


