/*
 * BarcodeInfo.cpp
 *
 *  Created on: Dec 10, 2015
 *      Author: nick
 */


#include "BarcodeInfo.hpp"

namespace njhseq {

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
	ret["class"] = njh::json::toJson(njh::getTypeName(*this));

	ret["fullBar_"] = njh::json::toJson(fullBar_);
	ret["extBar_"] = njh::json::toJson(extBar_);
	ret["ligBar_"] = njh::json::toJson(ligBar_);
	ret["mipId_"] = njh::json::toJson(mipId_);

	return ret;
}



}  // namespace njhseq

