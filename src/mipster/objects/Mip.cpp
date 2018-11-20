/*
 * mipObj.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
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
#include "Mip.hpp"

namespace njhseq {

Mip::Mip() :
		extBarcodeLen_(0), ligationArmMotObj_(""), igationArmMotObj5_to_3prime_(""), extentionArmMotObj_("") {
}


Mip::Mip(uint32_t extMbLen, uint32_t ligMbLen, const std::string & ligationArm,
		const std::string & extentionArm, const std::string & name,
		const std::string & familyName, const std::string & locGrouping,
		const std::string & mipSet) :
		name_(name),
		familyName_(familyName),
		regionGroup_(locGrouping),
		mipSet_(mipSet),
		extBarcodeLen_(extMbLen),
		ligBarcodeLen_(ligMbLen),
		ligationArm_(seqUtil::reverseComplement(stringToUpperReturn(ligationArm), "DNA")),
		ligationArm5_to_3prime_(stringToUpperReturn(ligationArm)),
		extentionArm_(stringToUpperReturn(extentionArm)),
		ligationArmObj_(seqInfo("ligationArm", ligationArm_)),
		extentionArmObj_(seqInfo("extentionArm", extentionArm_)),
		ligationArmMotObj_(ligationArm_),
		igationArmMotObj5_to_3prime_(stringToUpperReturn(ligationArm)),
		extentionArmMotObj_(extentionArm_) {
	//check for names

	std::regex namePat{".*_mip([0-9]+)[_]*.*$"};
	//std::regex namePat{".*mip([0-9]+)[_.*]*$"};
	std::smatch match;

	bool error = false;
	std::stringstream errorOutput;
	if(!std::regex_match(name_, match, namePat)){
		error = true;
		errorOutput << "Error in name format for target name " << name_ << ", must end with _mip[0-9]" << "\n";
	}
	std::smatch familyMatch;
	if(!std::regex_match(familyName_, familyMatch, namePat)){
		error = true;
		errorOutput << "Error in name format for family name " << familyName_ << ", must end with _mip[0-9]" << "\n";
	}
	if(error){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", errors found in constructing mip" << "\n";
		ss << errorOutput.str();
		throw std::runtime_error{ss.str()};
	}

}



BarcodeInfo Mip::determineExtBarcodeTrim(seqInfo & read) const {
	//assumes these arms are already determined to be there and isn't completely safe to run if
	//they aren't
	auto armPosMotif = getPossibleExtArmPos(read);
	//log  barcode
	std::string extBarcode = read.seq_.substr(
			armPosMotif.front().pos_ - extBarcodeLen_, extBarcodeLen_);
	read.trimFront(armPosMotif.front().pos_ + extentionArm_.size());
	return BarcodeInfo(name_, extBarcode);
}

BarcodeInfo Mip::determineExtBarcode(const seqInfo & read) const {
	//assumes these arms are already determined to be there and isn't completely safe to run if
	//they aren't
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	auto armPosMotif = getPossibleExtArmPos(read);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;
	//std::cout << "\t" <<"armPosMotif.size(): " << armPosMotif.size() << std::endl;
	//log  barcode
	std::string extBarcode = read.seq_.substr(
			armPosMotif.front().pos_ - extBarcodeLen_, extBarcodeLen_);
	//std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << std::endl;

	return BarcodeInfo(name_, extBarcode);
}

void Mip::determineLigBarcodeTrim(seqInfo & read, BarcodeInfo & info) const {
	//assumes these arms are already determined to be there and isn't completely safe to run if
	//they aren't
	// find ligation arm
	std::string ligBarcode = "";
	auto ligArmPos = getPossibleLigArmPos(read);
	if (ligBarcodeLen_ > 0) {
		ligBarcode = read.seq_.substr(ligArmPos.back().pos_ + len(ligationArmObj_),
				ligBarcodeLen_);
	}

	read.trimBack(ligArmPos.back().pos_);
	info.ligBar_ = ligBarcode;
	info.setFullBar();
}

void Mip::determineLigBarcode(const seqInfo & read, BarcodeInfo & info) const {
	//assumes these arms are already determined to be there and isn't completely safe to run if
	//they aren't
	// find ligation arm
	std::string ligBarcode = "";
	if (ligBarcodeLen_ > 0) {
		auto ligPosMotiff = getPossibleLigArmPos(read);
		ligBarcode = read.seq_.substr(ligPosMotiff.back().pos_ + len(ligationArmObj_),
				ligBarcodeLen_);
	}
	info.ligBar_ = ligBarcode;
	info.setFullBar();
}

BarcodeInfo Mip::determineBarcodesTrim(seqInfo & read) const {
	auto info = determineExtBarcodeTrim(read);
	determineLigBarcodeTrim(read, info);
	return info;
}

BarcodeInfo Mip::determineBarcodes(const seqInfo & read) const {
	auto info = determineExtBarcode(read);
	determineLigBarcode(read, info);
	return info;
}

void Mip::setAllowableErrorInArm(uint32_t allowableErrors) {
	allowableErrors_ = allowableErrors;
}

void Mip::setWiggleRoomInArm(uint32_t wiggleRoom) {
	wiggleRoomArm_ = wiggleRoom;
}

void Mip::setMinimumExpectedLen(size_t minimumExpectedLen) {
	minimumExpectedLen_ = minimumExpectedLen;
}


SinlgeMipExtractInfo::extractCase Mip::checkRead(seqInfo & read,
		const QualFilteringPars & qFilPars) const {
	//assumes good until a failure is reached;
	SinlgeMipExtractInfo::extractCase eCase =
			SinlgeMipExtractInfo::extractCase::GOOD;
	//quality control
	bool passQC = true;
	//check ligation arm
	if(passQC){
		auto ligArmPos = getPossibleLigArmPos(read);
		if (ligArmPos.empty()) {
			passQC = false;
			//if fails ligation arm and is below minimum expected length
			//add failure to bad min length as it probably just doesn't have the arm
			if (len(read) < minimumExpectedLen_) {
				eCase = SinlgeMipExtractInfo::extractCase::MINLENBAD;
				read.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
			}else{
				eCase = SinlgeMipExtractInfo::extractCase::BADREVERSE;
				read.name_.append("_failedLigArm");
			}
		}
	}

	//check minimum length
	if (passQC) {
		if (len(read) < minimumExpectedLen_) {
			eCase = SinlgeMipExtractInfo::extractCase::MINLENBAD;
			read.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
			passQC = false;
		}
	}

	//check quality scores
	if (passQC) {
		bool failedQuality = false;
		if (qFilPars.checkingQFrac_) {
			if (read.getQualCheck(qFilPars.qualCheck_) < qFilPars.qualCheckCutOff_) {
				failedQuality = true;
				read.name_.append("_failedQC");
			}
		} else if (qFilPars.checkingQWindow) {
			if (!seqUtil::checkQualityWindow(qFilPars.qualityWindowLength_,
					qFilPars.qualityWindowThres_, qFilPars.qualityWindowStep_,
					read.qual_)) {
				failedQuality = true;
				read.name_.append("_failedQCWindow");
			}
		}
		if (failedQuality) {
			eCase = SinlgeMipExtractInfo::extractCase::QUALITYFAILED;
			passQC = false;
		}
	}

	//check for Ns
	if (passQC) {
		if (std::string::npos != read.seq_.find('N')) {
			eCase = SinlgeMipExtractInfo::extractCase::CONTAINSNS;
			read.name_.append("_Contains_Ns");
			passQC = false;
		}
	}
	return eCase;
}

SinlgeMipExtractInfo::extractCase Mip::checkRead(PairedRead & seq,
		const QualFilteringPars & qFilPars) const {
	//assumes good until a failure is reached;
	SinlgeMipExtractInfo::extractCase eCase =
			SinlgeMipExtractInfo::extractCase::GOOD;
	//quality control
	//check ligation arm
	auto ligArmPos = getPossibleLigArmPos(seq.mateSeqBase_);
	if (ligArmPos.empty()) {
		//if fails ligation arm and is below minimum expected length
		//add failure to bad min length as it probably just doesn't have the arm
		if (len(seq.seqBase_) < minimumExpectedLen_ || len(seq.mateSeqBase_) < minimumExpectedLen_) {
			eCase = SinlgeMipExtractInfo::extractCase::MINLENBAD;
			seq.seqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
			seq.mateSeqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
		}else{
			eCase = SinlgeMipExtractInfo::extractCase::BADREVERSE;
			seq.seqBase_.name_.append("_failedLigArm");
			seq.mateSeqBase_.name_.append("_failedLigArm");
		}
		return eCase;
	}

	//check minimum length
	if (len(seq.seqBase_) < minimumExpectedLen_ || len(seq.mateSeqBase_) < minimumExpectedLen_) {
		eCase = SinlgeMipExtractInfo::extractCase::MINLENBAD;
		seq.seqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
		seq.mateSeqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
		return eCase;
	}

	//check if trimming to bad quality
	if(qFilPars.trimAtQual_){
		readVecTrimmer::trimAtFirstQualScore(seq.seqBase_, qFilPars.trimAtQualCutOff_);
		readVecTrimmer::trimToLastQualScore(seq.mateSeqBase_, qFilPars.trimAtQualCutOff_);
		//check minimum length again after trimming
		if (len(seq.seqBase_) < minimumExpectedLen_ || len(seq.mateSeqBase_) < minimumExpectedLen_) {
			eCase = SinlgeMipExtractInfo::extractCase::MINLENBAD;
			seq.seqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
			seq.mateSeqBase_.name_.append("_len<" + estd::to_string(minimumExpectedLen_));
			return eCase;
		}
	}



	//check quality scores
	bool failedQuality = false;
	if (qFilPars.checkingQFrac_) {
		if (seq.getQualCheck(qFilPars.qualCheck_) < qFilPars.qualCheckCutOff_ ) {
			failedQuality = true;
			seq.seqBase_.name_.append("_failedQC");
			seq.mateSeqBase_.name_.append("_failedQC");
		}
	} else if (qFilPars.checkingQWindow) {
		if (!seqUtil::checkQualityWindow(qFilPars.qualityWindowLength_,
				qFilPars.qualityWindowThres_, qFilPars.qualityWindowStep_,
				seq.seqBase_.qual_) ||
				!seqUtil::checkQualityWindow(qFilPars.qualityWindowLength_,
								qFilPars.qualityWindowThres_, qFilPars.qualityWindowStep_,
								seq.mateSeqBase_.qual_)) {
			failedQuality = true;
			seq.seqBase_.name_.append("_failedQCWindow");
			seq.mateSeqBase_.name_.append("_failedQCWindow");
		}
	}
	if (failedQuality) {
		eCase = SinlgeMipExtractInfo::extractCase::QUALITYFAILED;
		return eCase;
	}

	//check for Ns
	if (std::string::npos != seq.seqBase_.seq_.find('N') || std::string::npos != seq.mateSeqBase_.seq_.find('N')) {
		eCase = SinlgeMipExtractInfo::extractCase::CONTAINSNS;
		seq.seqBase_.name_.append("_Contains_Ns");
		seq.mateSeqBase_.name_.append("_Contains_Ns");
		return eCase;
	}

	return eCase;
}


std::vector<Mip::ArmPosScore> Mip::getPossibleExtArmPos(const seqInfo & read) const {
	std::vector<Mip::ArmPosScore> ret;
	auto positions = extentionArmMotObj_.findPositionsFull(read.seq_, allowableErrors_, extBarcodeLen_ + wiggleRoomArm_,
			extBarcodeLen_ + wiggleRoomArm_ + extentionArm_.size());
	for(const auto pos : positions){
		auto score = extentionArmMotObj_.scoreMotif(read.seq_.begin() + pos, read.seq_.begin() + pos + extentionArmMotObj_.size());
		ret.emplace_back(Mip::ArmPosScore{pos, score});
	}
	return ret;
}

std::vector<Mip::ArmPosScore> Mip::getPossibleLigArmPos(const seqInfo & read) const {
	std::vector<Mip::ArmPosScore> ret;
	std::vector<size_t> positions;
	if (wiggleRoomArm_ + ligBarcodeLen_ + ligationArm_.length() < len(read)) {
		positions = ligationArmMotObj_.findPositionsFull(read.seq_, allowableErrors_,
				len(read) - (wiggleRoomArm_ + ligBarcodeLen_ + ligationArm_.length()),
				len(read) - (wiggleRoomArm_ + ligBarcodeLen_));
		for(const auto pos : positions){
			auto score = ligationArmMotObj_.scoreMotif(read.seq_.begin() + pos, read.seq_.begin() + pos + ligationArmMotObj_.size());
			ret.emplace_back(Mip::ArmPosScore{pos, score});
		}
	}
	return ret;
}

std::vector<Mip::ArmPosScore> Mip::getPossibleLigArmPosFront(const seqInfo & read) const {

	std::vector<Mip::ArmPosScore> ret;
	auto positions = igationArmMotObj5_to_3prime_.findPositionsFull(read.seq_, allowableErrors_, ligBarcodeLen_ + wiggleRoomArm_,
			ligBarcodeLen_ + wiggleRoomArm_ + ligationArm5_to_3prime_.size());
	for(const auto pos : positions){
		auto score = igationArmMotObj5_to_3prime_.scoreMotif(read.seq_.begin() + pos, read.seq_.begin() + pos + igationArmMotObj5_to_3prime_.size());
		ret.emplace_back(Mip::ArmPosScore{pos, score});
	}
	return ret;
}




VecStr Mip::writeInfoLineHeader(){
	return VecStr{"mip_family","mip_id",
		"extension_arm","ligation_arm",
		"extension_barcode_length",
		"ligation_barcode_length","gene_name",
		"mipset"};
}

void Mip::writeInfoLine(std::ostream & out) const{
	out << familyName_
			<< "\t" << name_
			<< "\t" << extentionArm_
			<< "\t" << seqUtil::reverseComplement(ligationArm_, "DNA")
			<< "\t" << extBarcodeLen_
			<< "\t" << ligBarcodeLen_
			<< "\t" << regionGroup_
			<< "\t" << mipSet_ << "\n";
}

void Mip::writeOutArms(const OutOptions & opts) const{
	auto extOpts = SeqIOOptions::genFastaOut(njh::files::make_path(opts.outFilename_, name_ + "_ext-arm").string());
	auto ligOpts = SeqIOOptions::genFastaOut(njh::files::make_path(opts.outFilename_, name_ + "_lig-arm").string());
	extOpts.out_.overWriteFile_ = opts.overWriteFile_;
	extOpts.out_.append_ = opts.append_;
	ligOpts.out_.overWriteFile_ = opts.overWriteFile_;
	ligOpts.out_.append_ = opts.append_;
	seqInfo extArm("[mipTar=" + name_ + ";mipFam=" + familyName_ +";]", extentionArm_);
	seqInfo ligArm("[mipTar=" + name_ + ";mipFam=" + familyName_ +";]", seqUtil::reverseComplement(ligationArm_, "DNA"));
	SeqOutput::write(std::vector<seqInfo>{extArm}, extOpts);
	SeqOutput::write(std::vector<seqInfo>{ligArm}, ligOpts);
}


uint32_t Mip::getMipNumFromName(const std::string & name)  {
	std::regex namePat { "(.*)mip([0-9]+)$" };

	std::smatch match;
	uint32_t mipNumber = std::numeric_limits<uint32_t>::max();
	if (std::regex_match(name, match, namePat)) {
		mipNumber = estd::stou(match[2]);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error in processing name, didn't match expected pattern = "
				<< "(.*)mip([0-9]+)$" << "\n";
		throw std::runtime_error { ss.str() };
	}
	return mipNumber;
}

} /* namespace njhseq */
