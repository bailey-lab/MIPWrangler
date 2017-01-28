/*
 * mipObj.cpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include "Mip.hpp"

namespace bibseq {

Mip::Mip() :
		extBarcodeLen_(0), ligationArmMotObj_(""), extentionArmMotObj_("") {
}


Mip::Mip(uint32_t extMbLen, uint32_t ligMbLen, const std::string & ligationArm,
		const std::string & extentionArm, const std::string & name,
		const std::string & familyName, const std::string & locGrouping) :
		name_(name), familyName_(familyName),locGrouping_(locGrouping),
			extBarcodeLen_(extMbLen), ligBarcodeLen_(
				ligMbLen), ligationArm_(
				seqUtil::reverseComplement(stringToUpperReturn(ligationArm), "DNA")), extentionArm_(
				stringToUpperReturn(extentionArm)), ligationArmObj_(
				seqInfo("ligationArm", ligationArm_)), extentionArmObj_(
				seqInfo("extentionArm", extentionArm_)), ligationArmMotObj_(
				ligationArm_), extentionArmMotObj_(extentionArm_) {

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
	auto armPosMotif = getPossibleExtArmPos(read);
	//log  barcode
	std::string extBarcode = read.seq_.substr(
			armPosMotif.front().pos_ - extBarcodeLen_, extBarcodeLen_);
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

void Mip::writeOutArms(const OutOptions & opts) const{
	auto extOpts = SeqIOOptions::genFastaOut(bib::files::make_path(opts.outFilename_, name_ + "_ext-arm").string());
	auto ligOpts = SeqIOOptions::genFastaOut(bib::files::make_path(opts.outFilename_, name_ + "_lig-arm").string());
	extOpts.out_.overWriteFile_ = opts.overWriteFile_;
	extOpts.out_.append_ = opts.append_;
	ligOpts.out_.overWriteFile_ = opts.overWriteFile_;
	ligOpts.out_.append_ = opts.append_;
	seqInfo extArm("[mipTar=" + name_ + ";mipFam=" + familyName_ +";]", extentionArm_);
	seqInfo ligArm("[mipTar=" + name_ + ";mipFam=" + familyName_ +";]", seqUtil::reverseComplement(ligationArm_, "DNA"));
	SeqOutput::write(std::vector<seqInfo>{extArm}, extOpts);
	SeqOutput::write(std::vector<seqInfo>{ligArm}, ligOpts);
}

} /* namespace bibseq */
