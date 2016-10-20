#pragma once
/*
 * BarcodeInfo.hpp
 *
 *  Created on: Dec 10, 2015
 *      Author: nick
 */




#include "mipster/common.h"


namespace bibseq {

class BarcodeInfo {
public:
	BarcodeInfo();
	BarcodeInfo(const std::string & mipId, const std::string & extBar,
			const std::string & ligBar);
	BarcodeInfo(const std::string & mipId, const std::string & extBar);
	void setFullBar();

	std::string mipId_;
	std::string extBar_;
	std::string ligBar_;
	std::string fullBar_;
};





}  // namespace bibseq

