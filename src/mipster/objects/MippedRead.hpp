#pragma once
/*
 * MippedRead.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


#include "mipster/common.h"
#include "mipster/objects/BarcodeInfo.hpp"

namespace bibseq {

class MippedRead : public readObject{
public:

	using readObject::readObject;

	std::shared_ptr<BarcodeInfo> barInfo_;
};

}  // namespace bibseq




