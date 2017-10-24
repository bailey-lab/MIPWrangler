/*
 * filtering.cpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


#include "filtering.hpp"


namespace bibseq {





//filtering on barcodes
void filterOnMultipleLigBar(
		const std::unordered_map<std::string,
				std::vector<std::shared_ptr<MippedRead>>>& ligBarReads,
				std::vector<std::vector<std::shared_ptr<MippedRead>>>& readsPerBarcode,
				BarcodeFilterStats::BarcodeFilterStat& tarStat) {
	if (ligBarReads.size() == 1) {
		readsPerBarcode.push_back(ligBarReads.begin()->second);
	} else {
		// look to see if there is a majority for the second barcodes
		// if there is a tie take all ties
		/**@todo improve what to do with the second barcode disagreement */
		bool tiedMax = false;
		auto ligKeys = getVectorOfMapKeys(ligBarReads);
		std::vector<uint32_t> cov;
		for(const auto & ligBar : ligKeys) {
			cov.emplace_back(ligBarReads.at(ligBar).size());
		}
		auto maxEl = std::max_element(cov.begin(), cov.end());
		auto positions = getPositionsOfTarget(cov, *maxEl);
		if(positions.size() > 1) {
			tiedMax = true;
		}
		//add the reads with the most second barcode including ties
		for(const auto & pos : positions) {
			readsPerBarcode.push_back(ligBarReads.at(ligKeys[pos]));
		}
		for(const auto & pos : iter::range<uint32_t>(ligKeys.size())){
			if(!bib::in(pos, positions)){
				/**@todo report/write out the reads not used */
				tarStat.ligBarFilter_+= ligBarReads.at(ligKeys[pos]).size();
			}
		}
	}
}

std::shared_ptr<MippedRead> filterWithBarcodeCoverage(
		const std::vector<std::shared_ptr<MippedRead>> & barReads,aligner & alignerObj,
		const SeqSetUpPars & setUpPars, const mipBarcodeCorrectionPars & pars,
		BarcodeFilterStats::BarcodeFilterStat& tarStat) {
	std::shared_ptr<MippedRead> correctedRead = barReads.front();
	BarcodeInfo currentBarInfo = *(correctedRead->barInfo_);
	if (setUpPars.debug_) {
		std::cout << "Currently processing barcode: " << currentBarInfo.fullBar_
				<< std::endl;
	}

	if (barReads.size() > 1) {
		//barcode coverage is greater than 1 read
		//this can be taken advantage of to correct for errors by clustering these reads
		if (setUpPars.debug_) {
			std::cout << "\tCurrently clustering barcode: " << currentBarInfo.fullBar_
					<< std::endl;
		}
		//first collapse on simply string comparison, just an identical match collapse
		auto identicalClusters = clusterCollapser::collapseIdenticalReads(
				barReads, pars.qualRep); //median
		std::vector<cluster> clusters =
				baseCluster::convertVectorToClusterVector<cluster>(
						identicalClusters);
		if (setUpPars.debug_) {
			std::cout << "\tCollapsed to unique clusters for barcode: "
					<< currentBarInfo.fullBar_ << std::endl;
			std::cout << "\tFound " << clusters.size()
					<< " clusters for barcode: " << currentBarInfo.fullBar_ << std::endl;
		}
		if (clusters.size() > 1) {
			if (setUpPars.debug_) {
				std::cout << "\tSorting clusters for barcode: " << currentBarInfo.fullBar_
						<< std::endl;
			}
			//sort reads so more abundant and higher quality reads are higher up
			readVecSorter::sort(clusters);
			if (setUpPars.debug_) {
				std::cout << "\tDone Sorting clusters for barcode: " << currentBarInfo.fullBar_
						<< std::endl;
			}
			std::vector<uint32_t> positions(clusters.size());
			bib::iota<uint32_t>(positions, 0);
			//iterate over reads with smallest to highest comparison
			if (setUpPars.debug_) {
				std::cout << "\tStarting the clustering on id for barcode: "
						<< currentBarInfo.fullBar_ << std::endl;
			}
			for (const auto & readPos : iter::reversed(positions)) {
				if (setUpPars.debug_) {
					std::cout << "\t\tOn first cluster " << readPos
							<< " for barcode: " << currentBarInfo.fullBar_ << std::endl;
				}
				if (clusters[readPos].remove) {
					continue;
				}
				for (const auto & subReadPos : iter::range(readPos)) {
					if (setUpPars.debug_) {
						std::cout << "\t\tOn sub cluster " << subReadPos
								<< " for barcode: " << currentBarInfo.fullBar_ << std::endl;
					}
					if (clusters[subReadPos].remove) {
						continue;
					}
					if (pars.useReadLen) {
						if (setUpPars.debug_) {
							std::cout << "\t\tChecking lengths" << std::endl;
						}
						if (uAbsdiff(clusters[readPos].seqBase_.seq_.length(),
								clusters[subReadPos].seqBase_.seq_.length())
								> pars.readlenDiff) {
							if (setUpPars.debug_) {
								std::cout << "\t\tLengths differ by more than "
										<< pars.readlenDiff << " differ by "
										<< uAbsdiff(clusters[readPos].seqBase_.seq_.length(),
												clusters[subReadPos].seqBase_.seq_.length())
										<< std::endl;
							}
							continue;
						}
					}
					alignerObj.alignCacheGlobal(clusters[subReadPos], clusters[readPos]);
					if (setUpPars.debug_) {
						std::cout << "\t\tComparing" << std::endl;
					}
					alignerObj.profilePrimerAlignment(clusters[subReadPos],
							clusters[readPos]);
					alignerObj.comp_.recalcMismatchQuality(setUpPars.qScorePars_);
					alignerObj.comp_.setEventBaseIdentityHq();
					//cluster together reads that are .98 identity close to each other in a greedy fashion
					if (alignerObj.comp_.distances_.eventBasedIdentityHq_
							>= pars.barcodeIdentity) {
						if (setUpPars.debug_) {
							std::cout << "\t\tPassed identity check collapsing"
									<< std::endl;
						}
						clusters[subReadPos].addRead(clusters[readPos]);
						clusters[readPos].remove = true;
						break;
					}
				}
			}
			if (setUpPars.debug_) {
				std::cout << "\tDone the clustering on id for barcode: "
						<< currentBarInfo.fullBar_ << std::endl;
			}
			//remove the reads that had been added
			clusters = readVecSplitter::splitVectorOnRemove(clusters).first;
			//calculate a consensus for the reads
			clusterVec::allCalculateConsensus(clusters, alignerObj, true);
		}

		//choose a representative sequence for the barcode by using the top most abundant,highest quality cluster
		readVecSorter::sort(clusters);
		correctedRead = std::make_shared<MippedRead>(
				clusters.front().seqBase_);
		correctedRead->barInfo_ = std::make_shared<BarcodeInfo>(currentBarInfo);
		if(clusters.size() > 1){
			for(const auto pos : iter::range<uint32_t>(1,clusters.size())){
				alignerObj.alignCacheGlobal(clusters.front(), clusters[pos]);
				alignerObj.profilePrimerAlignment(clusters.front(),
						clusters[pos]);
				alignerObj.comp_.recalcMismatchQuality(setUpPars.qScorePars_);
				alignerObj.comp_.setEventBaseIdentityHq();
				/**@todo may want to print out the number of clusters left,
				 * so they can be examined for why they might not have been clustered in
				 * and how far off they are and the average size of the filtered clusters */
				tarStat.barFilter_ += clusters[pos].seqBase_.cnt_;
			}
		}
		if (setUpPars.debug_) {
			std::cout << "Done clustering barcode: " << currentBarInfo.fullBar_ << std::endl;
		}
	}
	return correctedRead;
}

}  // namespace bibseq

