#pragma once
/*
 * MipCollection.hpp
 *
 *  Created on: Jan 17, 2016
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
#include "mipster/objects/Mip.hpp"


namespace njhseq {

/**@brief Class to hold mip information about molecular barcodes and arms
 *
 */
class MipCollection {
public:

	/**@brief Construct with arm file and how many errors to allow in arms
	 *
	 * @param mipArmIdFile The id file that holds all the info for the mips
	 * @param allowableArmError The number of errors to allow in arms
	 */
	MipCollection(const bfs::path & mipArmIdFile, uint32_t allowableArmError);

	bfs::path mipArmIdFnp_;/**< The file info was read from*/
	uint32_t allowableArmError_; /**< allowable error in the arms*/
	std::unordered_map<std::string, Mip> mips_; /**< Map to hold mip info, key is the mip target name*/
	std::unordered_map<std::string, VecStr> mipNamesForFamily_; /**< Map to hold mip target names for mip family*/
	VecStr mipFamilies_; /**< All mip family names*/

	VecStr getMipsForFamily(const std::string & family) const;
	VecStr getMipsForFamily(const VecStr & families) const;

	VecStr getMipTarsForRegion(const std::string & region) const;
	VecStr getMipTarsForRegions(const VecStr & regions) const;
	VecStr getMipFamsForRegion(const std::string & region) const;
	VecStr getMipFamsForRegions(const VecStr & regions) const;
	VecStr getMipRegions() const;
	VecStr getMipFamilies() const;
	VecStr getMipTars() const;
	VecStr getMipRegionsForFams(const VecStr & mipFams) const;

	std::string getFamilyForTarget(const std::string & mipTarget) const;

	bool hasMipFamily(const std::string & mipFam) const;
	bool hasMipTarget(const std::string & mipTar) const;

	Mip determineBestMipInFamily(const seqInfo & read, Mip mip,
			aligner & alignerObjForGroupDet) const;

	Mip determineBestMipInFamily(const PairedRead & read, Mip mip,
			aligner & alignerObjForGroupDet) const;


	void setAllAllowableArmError(uint32_t allowableArmError);
	void setAllWiggleRoomInArm(uint32_t wiigleRoom);
	void setAllMinCaptureLength(uint32_t min_capture_length, bool force = false);

	void writeMipArmsFile(const OutOptions & outOpts) const;
};

}  // namespace njhseq

