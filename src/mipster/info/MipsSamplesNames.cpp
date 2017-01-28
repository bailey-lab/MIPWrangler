/*
 * MipsSamplesNames.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: nick
 */

#include "MipsSamplesNames.hpp"

namespace bibseq {

Json::Value MipFamSamp::toJson() const {
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["mipFam_"] = bib::json::toJson(mipFam_);
	ret["samp_"] = bib::json::toJson(samp_);
	return ret;
}

Json::Value MipTarFamSamp::toJson() const{
	Json::Value ret;
	ret["class"] = bib::getTypeName(*this);
	ret["mipTar_"] = bib::json::toJson(mipTar_);
	ret["mipFamSamp_"] = bib::json::toJson(mipFamSamp_);
	return ret;
}
MipsSamplesNames::MipsSamplesNames(const VecStr & mips, const VecStr & samples) :
		mips_(mips), samples_(samples) {
	bib::sort(samples_);
	bib::sort(mips_);
	//remove blanks, this often happens because the columns are different lengths
	removeElement(samples_, std::string(""));
	removeElement(mips_, std::string(""));
	//remove duplicates
	removeDuplicates(samples_);
	removeDuplicates(mips_);
}
MipsSamplesNames::MipsSamplesNames(const bfs::path & mipSampleFilename) {
	table mipSampInfo(mipSampleFilename, "\t", true);
	VecStr missingCols;
	VecStr neededCols { "mips", "samples" };
	for (const auto & col : neededCols) {
		if (!bib::in(col, mipSampInfo.columnNames_)) {
			missingCols.emplace_back(col);
		}
	}
	if (!missingCols.empty()) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", missing the following columns from " << mipSampleFilename << ", "
				<< bib::conToStr(missingCols, ",") << std::endl;
		ss << "Need the following columns, " << bib::conToStr(neededCols)
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	mips_ = mipSampInfo.getColumn("mips");
	for(auto & m : mips_){
		bib::trim(m);
	}
	samples_ = mipSampInfo.getColumn("samples");
	for(auto & s : samples_){
		bib::trim(s);
	}
	bib::sort(samples_);
	bib::sort(mips_);
	//remove blanks, this often happens because the columns are different lengths
	removeElement(samples_, std::string(""));
	removeElement(mips_, std::string(""));
	//remove duplicates
	removeDuplicates(samples_);
	removeDuplicates(mips_);
}

void MipsSamplesNames::setSamples(const VecStr & samples){
	samples_ = samples;
	bib::sort(samples_);
}

void MipsSamplesNames::setMips(const VecStr & mips){
	mips_ = mips;
	bib::sort(mips_);
}

void MipsSamplesNames::write(std::ostream & out)const{
	std::vector<VecStr> output;
	VecStr header;
	if(samples_.size() > mips_.size()){
		header = VecStr{"samples", "mips"};
		for(const auto & samp : samples_){
			output.emplace_back(VecStr{samp});
		}
		for(const auto pos : iter::range(mips_.size())){
			output[pos].emplace_back(mips_[pos]);
		}
	}else{
		header = VecStr{"mips", "samples"};
		for(const auto & mip : mips_ ){
			output.emplace_back(VecStr{mip});
		}
		for(const auto pos : iter::range(samples_.size())){
			output[pos].emplace_back(samples_[pos]);
		}
	}
	out << bib::conToStr(header, "\t") << std::endl;
	for(const auto & outline : output){
		out << bib::conToStr(outline, "\t") << std::endl;
	}
}



std::vector<MipFamSamp> MipsSamplesNames::createAllPairings()const{
	std::vector<MipFamSamp> allPairings;
	allPairings.reserve(mips_.size() * samples_.size());
	for (const auto & mipFam : mips_) {
		for (const auto & samp : samples_) {
			allPairings.emplace_back(mipFam, samp);
		}
	}
	return allPairings;
}

bool MipsSamplesNames::hasSample(const std::string & samp) const {
	return bib::in(samp, samples_);
}
bool MipsSamplesNames::hasMip(const std::string & mip) const {
	return bib::in(mip, mips_);
}


void printMipSampVec(const std::vector<MipFamSamp> & mipSamps, std::ostream & out){
	Json::FastWriter fWriter;
	for(const auto & mipSamp : mipSamps){
		out << fWriter.write(bib::json::toJson(mipSamp));
	}
}


std::vector<MipFamSamp> parseJsonForMipSamps(const Json::Value & val){
	std::vector<MipFamSamp> ret;
	for(const auto & obj : val){
		ret.emplace_back(obj["mipFam_"].asString(), obj["samp_"].asString());
	}
	return ret;
}

std::vector<MipFamSamp> parseJsonForMipSamps(const std::string & str){
	std::vector<MipFamSamp> ret;
	Json::Value root = bib::json::parse(str);
	for(const auto & obj : root){
		ret.emplace_back(obj["mipFam_"].asString(), obj["samp_"].asString());
	}
	return ret;
}

std::vector<MipFamSamp> parseJsonForMipSamps(std::istream & is){
	std::vector<MipFamSamp> ret;
	Json::Reader jReader;
	Json::Value root;
	jReader.parse(is, root);
	return parseJsonForMipSamps(root);
}

std::string toJsonStr(const std::vector<MipFamSamp> & mipSamps){
	Json::FastWriter fWriter;
	return fWriter.write(bib::json::toJson(mipSamps));
}

}  // namespace bibseq


