/*
 * MipMapResult.cpp
 *
 *  Created on: Jan 27, 2017
 *      Author: nick
 */




#include "MipMapResult.hpp"
#include <elucidator/BamToolsUtils.h>




namespace njhseq {
MipMapResult::MipMapResult(const std::string & mipName, const std::string & genomeName,
		const BamTools::BamAlignment & extAln,
		const BamTools::BamAlignment & ligAln) :
		mipName_(mipName), genomeName_(genomeName), extAln_(extAln), ligAln_(
				ligAln) {

}

bool MipMapResult::isMapped() const {
	return extAln_.IsMapped() && ligAln_.IsMapped();
}

bool MipMapResult::isConcordant() const {
	return extAln_.RefID == ligAln_.RefID;
}

void MipMapResult::setRegion(const BamTools::RefVector & refData) {
	MetaDataInName meta;
	meta.addMeta("genome", genomeName_, false);
	meta.addMeta("mipTar", mipName_, false);
	//full region
	region_.uid_ = meta.createMetaName();
	region_.chrom_ = refData[extAln_.RefID].RefName;
	region_.reverseSrand_ = extAln_.IsReverseStrand();
	region_.start_ = std::min(extAln_.Position, ligAln_.Position);
	region_.end_ = std::max(extAln_.GetEndPosition(), ligAln_.GetEndPosition());
	//extension arm
	extArmRegion_.uid_ = meta.createMetaName() + "-ext";
	extArmRegion_.chrom_ = refData[extAln_.RefID].RefName;
	extArmRegion_.reverseSrand_ = extAln_.IsReverseStrand();
	extArmRegion_.start_ = extAln_.Position;
	extArmRegion_.end_ = extAln_.GetEndPosition();
	//ligation arm
	ligArmRegion_.uid_ = meta.createMetaName() + "-lig";
	ligArmRegion_.chrom_ = refData[extAln_.RefID].RefName;
	ligArmRegion_.reverseSrand_ = extAln_.IsReverseStrand();
	ligArmRegion_.start_ = ligAln_.Position;
	ligArmRegion_.end_ = ligAln_.GetEndPosition();
}

std::vector<MipMapResult> getMipMapResults(const bfs::path & fnp, uint32_t insertSizeCutOff){
	if(bfs::exists(fnp)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << fnp << " doesn't exist" << "\n";
	}
	BamTools::BamReader bReader;
	BamTools::BamAlignment bAln;
	bReader.Open(fnp.string());
	checkBamOpenThrow(bReader, fnp);
	loadBamIndexThrow(bReader);
	auto refs = bReader.GetReferenceData();
	BamAlnsCache alnCache;
	std::vector<MipMapResult> results;

	std::string genomeFile = fnp.filename().string();
	std::string genome = genomeFile.substr(0, genomeFile.find("_"));
	while (bReader.GetNextAlignment(bAln)) {
		if(!bAln.IsPrimaryAlignment()){
			continue;
		}

		if(std::abs(bAln.InsertSize) > insertSizeCutOff){
			//skipp really large insert sizes;
			continue;
		}
		if (bAln.IsPaired()) {
			if(bAln.RefID != bAln.MateRefID){
				//discordant mapping, continue;
				continue;
			}
			if (bAln.MatePosition == bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//if mapped to the same place and the mate is yet to be encountered
					//enter into cache for until mate is encountered
					alnCache.add(bAln);
					continue;
				}
			}
			if (bAln.MatePosition <= bAln.Position) {
				if (!alnCache.has(bAln.Name)) {
					//since input should be sorted if matePosition is less than this positiion
					//it should be in the cache therefore program was unable to find mate

					//do orphaned operation
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ": Error, all input should be paired" << "\n";
					ss << "From: " << fnp << "\n";
					throw std::runtime_error{ss.str()};
				} else {
					auto search = alnCache.get(bAln.Name);
					MetaDataInName meta(bAln.Name);
					if (!meta.containsMeta("mipTar")) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__
								<< ": Error, name should contain meta data filed "
								<< njh::bashCT::boldRed("mipTar") << std::endl;
						throw std::runtime_error { ss.str() };
					}
					std::unique_ptr<MipMapResult> result;
					if (bAln.IsFirstMate()) {
						result = std::make_unique<MipMapResult>(meta.getMeta("mipTar"),
								genome, bAln, *search);
					} else {
						result = std::make_unique<MipMapResult>(meta.getMeta("mipTar"),
								genome, *search, bAln);
					}
					if (result->isConcordant() && result->isMapped()) {
						result->setRegion(refs);
					}
					results.emplace_back(*result);
					// now that operations have been computed, remove first mate found from cache
					alnCache.remove(search->Name);
				}
			} else {
				//enter into cache for until mate is encountered
				alnCache.add(bAln);
			}
		}else{
			//unpaired read
			//do unpaired read operation
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, all input should be paired " << "\n";
			ss << "From: " << fnp << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	return results;
}


}  // namespace njhseq

