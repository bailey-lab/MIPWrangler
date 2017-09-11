/*
 * BarcodeInfo.cpp
 *
 *  Created on: Dec 10, 2015
 *      Author: nick
 */


#include "BarcodeInfo.hpp"

namespace bibseq {

BarcodeInfo::BarcodeInfo() {
}

BarcodeInfo::BarcodeInfo(const std::string & mipId, const std::string & extBar,
		const std::string & ligBar) :
		mipId_(mipId), extBar_(extBar), ligBar_(ligBar) {
	setFullBar();
}

BarcodeInfo::BarcodeInfo(const std::string & mipId, const std::string & extBar) :
		mipId_(mipId), extBar_(extBar), ligBar_("") {
	setFullBar();
}
void BarcodeInfo::setFullBar() {
	fullBar_ = extBar_ + ligBar_;
}

Json::Value BarcodeInfo::toJson() const{
	Json::Value ret;
	ret["class"] = bib::json::toJson(bib::getTypeName(*this));

	ret["fullBar_"] = bib::json::toJson(fullBar_);
	ret["extBar_"] = bib::json::toJson(extBar_);
	ret["ligBar_"] = bib::json::toJson(ligBar_);
	ret["mipId_"] = bib::json::toJson(mipId_);

	return ret;
}



}  // namespace bibseq

