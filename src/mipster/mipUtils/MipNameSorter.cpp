/*
 * MipNameSorter.cpp
 *
 *  Created on: Feb 21, 2017
 *      Author: nick
 */


// MIPWrangler - A library for analyzing sequence data from molecular inversion probe analysis
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of MIPWrangler.
//
// MIPWrangler is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MIPWrangler is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MIPWrangler.  If not, see <http://www.gnu.org/licenses/>.
//
#include "MipNameSorter.hpp"

namespace njhseq {

const std::regex MipNameSorter::mipNamePat{"(.*)mip([0-9]+)[_]*.*$"};
const std::regex MipNameSorter::regionNamePat{R"(^([A-Za-z]+)(\d+)(.*))"};


bool MipNameSorter::compareNames(const std::string & name1,
		const std::string & name2, const std::regex & namePat) {
	//assumes pattern has been checked for match count of at least 2
	std::smatch match1;
	std::smatch match2;
	if (std::regex_match(name1, match1, namePat)
			&& std::regex_match(name2, match2, namePat)) {
		return
				match1[1] == match2[1] ?
						estd::stou(match1[2])
								< estd::stou(match2[2]) :
						name1 < name2;
	} else {
		return name1 < name2;
	}
}

bool MipNameSorter::compareNames(const std::string & name1,
		const std::string & name2, const std::regex & namePat,
		const std::regex & secondaryPat) {
	std::smatch match1;
	std::smatch match2;
	if (std::regex_match(name1, match1, namePat)
			&& std::regex_match(name2, match2, namePat)) {
		if (match1[1] == match2[1]) {
			return estd::stou(match1[2])
					< estd::stou(match2[2]);
		} else {
			return compareNames(name1, name2, secondaryPat);
		}
	} else {
		return name1 < name2;
	}
}


void MipNameSorter::sort(VecStr & names, const std::regex & namePat) {
	sort(names, namePat, std::function<std::string(const std::string &)>([](const std::string & str1){ return str1;}));
}

void MipNameSorter::sort(VecStr & names, const std::regex & namePat, const std::regex & secondaryPat){
	sort(names, namePat,secondaryPat, std::function<std::string(const std::string &)>([](const std::string & str1){ return str1;}));
}

}  // namespace njhseq

