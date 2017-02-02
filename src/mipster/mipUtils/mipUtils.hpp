#pragma once
/*
 * mipUtils.hpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */


#include "mipster/common.h"

namespace bibseq {


class MipNameSorter{
public:

	static void sort(VecStr & names);

	static void sort(VecStr & names, const std::regex & namePat);

};

}  // namespace bibseq


