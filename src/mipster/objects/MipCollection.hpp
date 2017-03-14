#pragma once
/*
 * MipCollection.hpp
 *
 *  Created on: Jan 17, 2016
 *      Author: nick
 */

#include "mipster/objects/Mip.hpp"


namespace bibseq {

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


	void setAllAllowableArmError(uint32_t allowableArmError);
	void setAllWiggleRoomInArm(uint32_t wiigleRoom);
	void setAllMinimumExpectedLen(size_t minimumExpectedLen);
};

}  // namespace bibseq

