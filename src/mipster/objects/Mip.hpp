#pragma once
/*
 * mipObj.hpp
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
#include <njhseq/objects/helperObjects/motif.hpp>
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>

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
	std::string familyName_;
	std::string regionGroup_;
	std::string mipSet_;

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

	std::vector<PairedReadProcessor::ReadPairOverLapStatus> allowableStatuses{PairedReadProcessor::ReadPairOverLapStatus::R1ENDSINR2};

	struct ArmPosScore {
	public:
		size_t pos_ = std::numeric_limits<size_t>::max();
		uint32_t score_ = std::numeric_limits<uint32_t>::max();
	};

	uint32_t allowableErrors_ = 0;
	uint32_t wiggleRoomArm_ = 0;
	uint32_t min_capture_length_ = std::numeric_limits<uint32_t>::max();

	void setAllowableErrorInArm(uint32_t allowableErrors);
	void setWiggleRoomInArm(uint32_t wiggleRoom);
	void setMinCaptureLength(uint32_t min_capture_length);

	[[nodiscard]] std::vector<ArmPosScore> getPossibleExtArmPos(const seqInfo & read) const;
	[[nodiscard]] std::vector<ArmPosScore> getPossibleLigArmPos(const seqInfo & read) const;
	[[nodiscard]] std::vector<ArmPosScore> getPossibleLigArmPosFront(const seqInfo & read) const;

	BarcodeInfo determineExtBarcodeTrim(seqInfo & read) const;
	[[nodiscard]] BarcodeInfo determineExtBarcode(const seqInfo & read) const;

	void determineLigBarcodeTrim(seqInfo & read, BarcodeInfo & info) const;
	void determineLigBarcode(const seqInfo & read, BarcodeInfo & info) const;

	BarcodeInfo determineBarcodesTrim(seqInfo & read) const;
	[[nodiscard]] BarcodeInfo determineBarcodes(const seqInfo & read) const;

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

