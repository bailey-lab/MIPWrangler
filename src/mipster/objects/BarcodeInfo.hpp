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

	std::string mipId_;   /**< the mip id the barcode belongs to */
	std::string extBar_;  /**< the extension barcode */
	std::string ligBar_;  /**< the ligation barcode */
	std::string fullBar_; /** the ext and lig barcode pasted together to create just one object to reference as a barcode*/

	Json::Value toJson() const;

};





}  // namespace bibseq

