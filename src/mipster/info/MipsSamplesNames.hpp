#pragma once
/*
 * MipsSamplesNames.hpp
 *
 *  Created on: Feb 5, 2016
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

#include "mipster/common.h"

namespace njhseq {
struct MipFamSamp {
	MipFamSamp(std::string mipFam, std::string samp) :
			mipFam_(mipFam), samp_(samp) {
	}
	std::string mipFam_;
	std::string samp_;

	Json::Value toJson() const;
};

struct MipTarFamSamp {
	MipTarFamSamp(std::string mipTar, MipFamSamp mipFamSamp) :
			mipTar_(mipTar), mipFamSamp_(mipFamSamp) {
	}
	std::string mipTar_;
	MipFamSamp mipFamSamp_;

	Json::Value toJson() const;
};

class MipsSamplesNames {
public:
	MipsSamplesNames(const bfs::path & mipSampleFilename);
	MipsSamplesNames(const VecStr & mips, const VecStr & samples);
	VecStr mips_;
	VecStr samples_;

	void setSamples(const VecStr & samples);
	void setMips(const VecStr & mips);
	void write(std::ostream & out) const;
	bool hasSample(const std::string & samp) const;
	bool hasMip(const std::string & mip) const;
	std::vector<MipFamSamp> createAllPairings() const;

	std::set<std::string> getSetSampNames() const;

};

void printMipSampVec(const std::vector<MipFamSamp> & mipSamps,
		std::ostream & out = std::cout);

std::vector<MipFamSamp> parseJsonForMipSamps(const Json::Value & val);
std::vector<MipFamSamp> parseJsonForMipSamps(const std::string & str);
std::vector<MipFamSamp> parseJsonForMipSamps(std::istream & is);
std::string toJsonStr(const std::vector<MipFamSamp> & mipSamps);



}  // namespace njhseq




