#pragma once
/*
 * MipsSamplesNames.hpp
 *
 *  Created on: Feb 5, 2016
 *      Author: nick
 */


#include "mipster/common.h"

namespace bibseq {
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
	MipsSamplesNames(const std::string & mipSampleFilename);
	MipsSamplesNames(const VecStr & mips, const VecStr & samples);
	VecStr mips_;
	VecStr samples_;

	void setSamples(const VecStr & samples);
	void setMips(const VecStr & mips);
	void write(std::ostream & out) const;
	bool hasSample(const std::string & samp) const;
	bool hasMip(const std::string & mip) const;
	std::vector<MipFamSamp> createAllPairings() const;
};

void printMipSampVec(const std::vector<MipFamSamp> & mipSamps,
		std::ostream & out = std::cout);

std::vector<MipFamSamp> parseJsonForMipSamps(const Json::Value & val);
std::vector<MipFamSamp> parseJsonForMipSamps(const std::string & str);
std::vector<MipFamSamp> parseJsonForMipSamps(std::istream & is);
std::string toJsonStr(const std::vector<MipFamSamp> & mipSamps);



}  // namespace bibseq




