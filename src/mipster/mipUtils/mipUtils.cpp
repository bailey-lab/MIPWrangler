/*
 * mipUtils.cpp
 *
 *  Created on: Jul 26, 2016
 *      Author: nick
 */


#include "mipUtils.hpp"

namespace bibseq {

void MipNameSorter::sort(VecStr & names){
	std::regex namePat{"(.*)mip([0-9]+)$"};
	sort(names, namePat);
}

void MipNameSorter::sort(VecStr & names, const std::regex & namePat){
	if(namePat.mark_count() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error regex pattern must have at least 2 sub patterns capture to sort by" << std::endl;
	}
	bib::sort(names, [&namePat](const std::string & name1,
			const std::string & name2){
		std::smatch match1;
		std::smatch match2;
		if(std::regex_match(name1, match1, namePat) &&
				std::regex_match(name2, match2, namePat)){
			return match1[1] == match2[1] ? bib::lexical_cast<uint32_t>(match1[2]) < bib::lexical_cast<uint32_t>(match2[2]) : name1 < name2;
		}else{
			return name1 < name2;
		}
	});
}

}  // namespace bibseq
