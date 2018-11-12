/*
 * mipUtils.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */


#include "mipUtils.hpp"

namespace njhseq {

void processNameForBarReadCounts(const std::string & name,
		uint32_t & barNum, uint32_t & readNum){
	auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	auto underPos = name.rfind("_");
	//std::cout << vectorToString(toVecStr(rPos, bPos, underPos, name.substr(rPos + 2, bPos - rPos - 2),name.substr(bPos + 2, underPos - bPos - 2) ), "\n") << std::endl;
	readNum = estd::stou(name.substr(rPos + 2, bPos - rPos - 2));
	barNum = estd::stou(name.substr(bPos + 2, underPos - bPos - 2));
}

uint32_t getReadCntFromMipName(const std::string & name) {
	auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	//auto underPos = name.rfind("_");
	return estd::stou(name.substr(rPos + 2, bPos - rPos - 2));
}

uint32_t getBarcodeCntFromMipName(const std::string & name) {
	//auto rPos = name.rfind("_R");
	auto bPos = name.rfind("_B");
	auto underPos = name.rfind("_");
	return estd::stou(name.substr(bPos + 2, underPos - bPos - 2));
}

}  // namespace njhseq
