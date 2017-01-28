/*
 * MetaDataInName.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: nick
 */


#include "MetaDataInName.hpp"

namespace bibseq {


MetaDataInName::MetaDataInName() {

}

MetaDataInName::MetaDataInName(const std::string & str) {
	processNameForMeta(str, false);
}

void MetaDataInName::processNameForMeta(const std::string & name, bool replace){
	auto firstBracket = name.find("[");
	if(std::string::npos == firstBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find [ in " << name << std::endl;
		throw std::runtime_error{ss.str()};
	}
	auto secondBracket = name.find("]", firstBracket);
	if(std::string::npos == secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", could not find ] in " << name  << " after " << firstBracket
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	if(firstBracket > secondBracket){
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__ << ", [ must come before ] "
				<< "\n";
		ss << "[ pos: " << firstBracket << "; ] pos: " << secondBracket << "\n"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto toks = tokenizeString(name.substr(firstBracket + 1, secondBracket - firstBracket - 1), ";");
	for(const auto & tok : toks){
		auto subToks = tokenizeString(tok, "=");
		if(2 != subToks.size()){
			std::stringstream ss;
			ss << "Error in : " << __PRETTY_FUNCTION__
					<< "values should be separated by one =, not " << tok
					<< std::endl;
			throw std::runtime_error{ss.str()};
		}else{
			addMeta(subToks[0], subToks[1], replace);
		}
	}
}

bool MetaDataInName::containsMeta(const std::string & key) const {
	return meta_.find(key) != meta_.end();
}

std::string MetaDataInName::getMeta(const std::string & key) const {
	auto search = meta_.find(key);
	if (search != meta_.end()) {
		return search->second;
	}else{
		std::stringstream ss;
		ss << __FILE__ << " - " << __LINE__ << " : " << __PRETTY_FUNCTION__ << ", error no meta field " << key << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(meta_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return "";
}

std::string MetaDataInName::createMetaName() const {
	std::string newMeta = "[";
	auto metaKeys = bib::getVecOfMapKeys(meta_);
	bib::sort(metaKeys);
	for (const auto & metaKey : metaKeys) {
		const auto & meta = meta_.at(metaKey);
		if ("[" != newMeta) {
			newMeta.append(";" + metaKey + "=" + meta);
		} else {
			newMeta.append(metaKey + "=" + meta);
		}
	}
	newMeta.append("]");
	return newMeta;
}

void MetaDataInName::resetMetaInName(std::string & name,
		size_t pos) {
	auto firstBracket = name.find("[");
	auto secondBracket = name.find("]", firstBracket);
	std::string newMeta = createMetaName();
	if (std::string::npos != firstBracket
			&& std::string::npos != secondBracket) {
		name = name.substr(0, firstBracket) + newMeta
				+ name.substr(secondBracket + 1);
	} else {
		if (std::numeric_limits<size_t>::max() != pos && pos < name.size()) {
			name.insert(name.begin() + pos, newMeta.begin(), newMeta.end());
		} else {
			name += newMeta;
		}
	}
}

bool MetaDataInName::nameHasMetaData(const std::string & name) {
	auto firstBracket = name.find("[");
	if (std::string::npos == firstBracket) {
		return false;
	}
	auto secondBracket = name.find("]", firstBracket);
	if (std::string::npos == secondBracket) {
		return false;
	}
	if (firstBracket > secondBracket) {
		return false;
	}
	return true;
}

}  // namespace bibseq


