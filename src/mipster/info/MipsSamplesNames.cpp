/*
 * MipsSamplesNames.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: nick
 */

#include "MipsSamplesNames.hpp"
#include "mipster/mipUtils/MipNameSorter.hpp"

namespace njhseq {

Json::Value MipFamSamp::toJson() const {
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["mipFam_"] = njh::json::toJson(mipFam_);
	ret["samp_"] = njh::json::toJson(samp_);
	return ret;
}

Json::Value MipTarFamSamp::toJson() const{
	Json::Value ret;
	ret["class"] = njh::getTypeName(*this);
	ret["mipTar_"] = njh::json::toJson(mipTar_);
	ret["mipFamSamp_"] = njh::json::toJson(mipFamSamp_);
	return ret;
}
MipsSamplesNames::MipsSamplesNames(const VecStr & mips, const VecStr & samples) :
		mips_(mips), samples_(samples) {
	njh::sort(samples_);
	MipNameSorter::sort(mips_);
	//remove blanks, this often happens because the columns are different lengths
	removeElement(samples_, std::string(""));
	removeElement(mips_, std::string(""));
	//remove duplicates
	removeDuplicates(samples_);
	removeDuplicates(mips_);
	MipNameSorter::sort(mips_);
}
MipsSamplesNames::MipsSamplesNames(const bfs::path & mipSampleFilename) {
	table mipSampInfo(mipSampleFilename, "\t", true);
	VecStr missingCols;
	VecStr neededCols { "mips", "samples" };
	for (const auto & col : neededCols) {
		if (!njh::in(col, mipSampInfo.columnNames_)) {
			missingCols.emplace_back(col);
		}
	}
	if (!missingCols.empty()) {
		std::stringstream ss;
		ss << "Error in : " << __PRETTY_FUNCTION__
				<< ", missing the following columns from " << mipSampleFilename << ", "
				<< njh::conToStr(missingCols, ",") << std::endl;
		ss << "Need the following columns, " << njh::conToStr(neededCols)
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	mips_ = mipSampInfo.getColumn("mips");
	for(auto & m : mips_){
		njh::trim(m);
	}
	samples_ = mipSampInfo.getColumn("samples");
	for(auto & s : samples_){
		njh::trim(s);
	}
	njh::sort(samples_);
	MipNameSorter::sort(mips_);
	//remove blanks, this often happens because the columns are different lengths
	removeElement(samples_, std::string(""));
	removeElement(mips_, std::string(""));
	//remove duplicates
	removeDuplicates(samples_);
	removeDuplicates(mips_);
	MipNameSorter::sort(mips_);
}

void MipsSamplesNames::setSamples(const VecStr & samples){
	samples_ = samples;
	njh::sort(samples_);
}

void MipsSamplesNames::setMips(const VecStr & mips){
	mips_ = mips;
	MipNameSorter::sort(mips_);
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
	out << njh::conToStr(header, "\t") << std::endl;
	for(const auto & outline : output){
		out << njh::conToStr(outline, "\t") << std::endl;
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


std::set<std::string> MipsSamplesNames::getSetSampNames() const {
	return std::set<std::string> { samples_.begin(), samples_.end() };
}

bool MipsSamplesNames::hasSample(const std::string & samp) const {
	return njh::in(samp, samples_);
}
bool MipsSamplesNames::hasMip(const std::string & mip) const {
	return njh::in(mip, mips_);
}


void printMipSampVec(const std::vector<MipFamSamp> & mipSamps, std::ostream & out){
	for(const auto & mipSamp : mipSamps){
		out << njh::json::writeAsOneLine(njh::json::toJson(mipSamp));
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
	Json::Value root = njh::json::parse(str);
	for(const auto & obj : root){
		ret.emplace_back(obj["mipFam_"].asString(), obj["samp_"].asString());
	}
	return ret;
}

std::vector<MipFamSamp> parseJsonForMipSamps(std::istream & is){
	std::vector<MipFamSamp> ret;
	Json::Value root = njh::json::parseStream(is);
	return parseJsonForMipSamps(root);
}

std::string toJsonStr(const std::vector<MipFamSamp> & mipSamps){
	return njh::json::writeAsOneLine(njh::json::toJson(mipSamps));
}

}  // namespace njhseq


