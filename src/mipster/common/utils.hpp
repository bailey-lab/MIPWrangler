#pragma once
/*
 * utils.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */



#include "mipster/common/allExtLibraryIncludes.h"
#include "mipster/common/allSystemIncludes.h"


namespace bibseq {

namespace bfs = boost::filesystem;
void checkExistenceThrow(const bfs::path & dirName, const std::string & funcName);
void checkExistenceThrow(const bfs::path & dirName);


bool requireExternalProgramThrow(const std::string & program);


}  // namespace bibseq

