#pragma once
/*
 * mipObj.hpp
 *
 *  Created on: Dec 30, 2014
 *      Author: nickhathaway
 */

#include <njhseq/objects/helperObjects/motif.hpp>
#include <njhseq/objects/seqObjects/readObject.hpp>

#include "mipster/common.h"
#include "mipster/objects/BarcodeInfo.hpp"
#include "mipster/info/filterStats/SinlgeMipExtractInfo.hpp"

namespace njhseq {

class Mip {
public:
	Mip();

	Mip(uint32_t extMbLen, uint32_t ligBmLen, const std::string & ligationArm,
			const std::string & extentionArm, const std::string & name,
			const std::string & familyName, const std::string & locGrouping,
			const std::string & mipSet);

	//members
	std::string name_;
	std::string familyName_ = "";
	std::string regionGroup_ = "";
	std::string mipSet_ = "";

	uint32_t extBarcodeLen_;
	uint32_t ligBarcodeLen_ = 0;

	std::string ligationArm_;
	std::string ligationArm5_to_3prime_;
	std::string extentionArm_;

	seqInfo ligationArmObj_;
	seqInfo extentionArmObj_;

	motif ligationArmMotObj_;
	motif igationArmMotObj5_to_3prime_;
	motif extentionArmMotObj_;

	struct ArmPosScore {
	public:
		size_t pos_ = std::numeric_limits<size_t>::max();
		uint32_t score_ = std::numeric_limits<uint32_t>::max();
	};

	uint32_t allowableErrors_ = 0;
	uint32_t wiggleRoomArm_ = 0;
	size_t minimumExpectedLen_ = std::numeric_limits<size_t>::max();

	void setAllowableErrorInArm(uint32_t allowableErrors);
	void setWiggleRoomInArm(uint32_t wiggleRoom);
	void setMinimumExpectedLen(size_t minimumExpectedLen);

	std::vector<ArmPosScore> getPossibleExtArmPos(const seqInfo & read) const;
	std::vector<ArmPosScore> getPossibleLigArmPos(const seqInfo & read) const;
	std::vector<ArmPosScore> getPossibleLigArmPosFront(const seqInfo & read) const;

	BarcodeInfo determineExtBarcodeTrim(seqInfo & read) const;
	BarcodeInfo determineExtBarcode(const seqInfo & read) const;

	void determineLigBarcodeTrim(seqInfo & read, BarcodeInfo & info) const;
	void determineLigBarcode(const seqInfo & read, BarcodeInfo & info) const;

	BarcodeInfo determineBarcodesTrim(seqInfo & read) const;
	BarcodeInfo determineBarcodes(const seqInfo & read) const;

	SinlgeMipExtractInfo::extractCase checkRead(seqInfo & read,
			const QualFilteringPars & qFilPars) const;

	SinlgeMipExtractInfo::extractCase checkRead(PairedRead & read,
			const QualFilteringPars & qFilPars) const;


	void writeOutArms(const OutOptions & opts) const;

	static VecStr writeInfoLineHeader();
	void writeInfoLine(std::ostream & out) const;

	static uint32_t getMipNumFromName(const std::string & name);

};

} /* namespace njhseq */

