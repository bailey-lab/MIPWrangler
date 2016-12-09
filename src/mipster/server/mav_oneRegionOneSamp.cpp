/*
 * mav_oneGeneOneSamp.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: nick
 */


#include "mav.hpp"
#include <experimental/objects/dataContainers/graphs.h>


namespace bibseq {


void mav::regionInfoForSampPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	auto region = request->get_path_parameter("region");
	std::stringstream ss;
	bool redirecting = false;
		//check region
	if (!bib::in(region, mipMaster_->getMipGroupings())) {
		ss << "Region name wasn't found :" << region << std::endl;
		redirecting = true;
	}
	//check for sample @todo check to see if region specifically has info on sample
	if (popClusInfoBySample_.find(sample) == popClusInfoBySample_.end()) {
		ss << "Sample name wasn't found :" << sample << std::endl;
		redirecting = true;
	}

	if (redirecting) {
		redirect(session, ss.str());
	}else{
		auto body = genHtmlDoc(rootName_, pages_.at("oneRegionOneSamp.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}

}

void mav::getMipOverlapGraphDataHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	auto region = request->get_path_parameter("region");
	Json::Value ret;
	if (bib::in(region, mipMaster_->getMipGroupings())) {
		//check for sample @todo check to see if region specifically has info on sample
		if (popClusInfoBySample_.find(sample) != popClusInfoBySample_.end()) {
			auto mipFams = mipMaster_->mips_->getMipFamsForRegion(region);
			std::vector<std::shared_ptr<seqInfo>> finalSeqs;
			for (const auto & mipName : mipFams) {
				std::string searchTerm = mipName + "_" + sample + "_final";
				if (bib::in(searchTerm, seqs_->cache_)) {
					addOtherVec(finalSeqs,
							seqs_->cache_.at(searchTerm).ioOpts_.getPtrs<seqInfo>());
				} else {
					std::cerr << "Sample and Mip fam seqs not found :"
							<< sample << " - " << mipName << std::endl;
				}
			}
			uint32_t overlapMin = 5;
			uint64_t maxSize = 0;
			std::map<uint32_t, std::vector<std::shared_ptr<seqInfo>>>seqsByMipNum;
			SeqOverlapGraph graph;

			for (const auto & seq : finalSeqs) {
				readVec::getMaxLength(seq, maxSize);

				std::unordered_map<std::string, std::string> meta;
				seq->processNameForMeta(meta);
				if (bib::has(meta, "mipFam")) {
					auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
					uint32_t mipNum = bib::lexical_cast<uint32_t>(
							bib::replaceString(mipFamToks.back(), "mip", ""));
					seqsByMipNum[mipNum].push_back(seq);
				} else {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__
							<< ": Error, should have meta data for mipFam" << std::endl;
					throw std::runtime_error { ss.str() };
				}
			}

			aligner alignerObj(maxSize, gapScoringParameters(5, 1, 0, 0, 0, 0),
					substituteMatrix(2, -2), false);
			comparison noErrors;
			//check to see if target should be reversed
			for (const auto & mips : seqsByMipNum) {
				auto search = seqsByMipNum.find(mips.first + 1);
				if (seqsByMipNum.end() != search) {
					for (const auto & seq1 : mips.second) {
						uint32_t forCount = 0;
						uint32_t revCount = 0;
						for (const auto & seq2 : search->second) {
							alignerObj.alignCacheGlobal(*seq1, *seq2);
							alignerObj.profilePrimerAlignment(*seq1, *seq2);
							auto forComp = alignerObj.comp_;
							auto seq2Rev = *seq2;
							seq2Rev.reverseComplementRead(true, true);
							alignerObj.alignCacheGlobal(*seq1, seq2Rev);
							alignerObj.profilePrimerAlignment(*seq1, seq2Rev);
							bool passFor = noErrors.passErrorProfile(forComp);
							bool passRev = noErrors.passErrorProfile(alignerObj.comp_);
							if (passFor && passRev) {
								if (forComp.distances_.query_.covered_
										>= alignerObj.comp_.distances_.query_.covered_) {
									++forCount;
								} else {
									++revCount;
								}
							} else if (passFor) {
								++forCount;
							} else if (passRev) {
								++revCount;
							}
						}
						if (revCount > forCount) {
							for (auto & seq : search->second) {
								seq->reverseComplementRead(true, true);
							}
						}
					}
				}
				bool notProcessed = false;
				if (notProcessed) {
					for (auto & seq : mips.second) {
						seq->frac_ = 1.0 / len(mips.second);
					}
				}
				//add to graph
				for (const auto & seq : mips.second) {
					graph.addNode(seq);
				}
			}
			for (const auto & mips : seqsByMipNum) {
				auto search = seqsByMipNum.find(mips.first + 1);
				if (seqsByMipNum.end() != search) {
					for (const auto & seq1 : mips.second) {
						for (const auto & seq2 : search->second) {
							alignerObj.alignCacheGlobal(*seq1, *seq2);
							alignerObj.profilePrimerAlignment(*seq1, *seq2);
							if (noErrors.passErrorProfile(alignerObj.comp_)
									&& alignerObj.comp_.distances_.query_.covered_
											>= overlapMin) {
								bool seq1GapFront = '-'
										== alignerObj.alignObjectA_.seqBase_.seq_.front();
								bool seq1GapBack = '-'
										== alignerObj.alignObjectA_.seqBase_.seq_.back();
								bool seq2GapFront = '-'
										== alignerObj.alignObjectB_.seqBase_.seq_.front();
								bool seq2GapBack = '-'
										== alignerObj.alignObjectB_.seqBase_.seq_.back();
								if ((seq1GapFront && seq1GapBack)
										|| (seq2GapFront && seq2GapBack)) {
									//no overlap
								} else {
									if (seq1GapBack || seq2GapFront) {
										graph.addEdge(seq1->name_, seq2->name_);
									} else {
										graph.addEdge(seq2->name_, seq1->name_);
									}
								}
							}
						}
					}
				}
			}
			ret = graph.createSankeyOutput();
			uint32_t maxNum = 0;
			for (const auto & seqs : seqsByMipNum) {
				if (len(seqs.second) > maxNum) {
					maxNum = len(seqs.second);
				}
			}
			auto outColors = bib::njhColors(maxNum);
			bibseq::VecStr outColorsStrs;
			outColorsStrs.reserve(outColors.size());
			for (const auto & c : outColors) {
				outColorsStrs.emplace_back("#" + c.hexStr_);
			}

			for (auto & node : ret["nodes"]) {
				std::unordered_map<std::string, std::string> meta;
				seqInfo seq(node["name"].asString(), "A");
				seq.processNameForMeta(meta);
				if (bib::has(meta, "mipFam")) {
					auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
					uint32_t mipNum = bib::lexical_cast<uint32_t>(
							bib::replaceString(mipFamToks.back(), "mip", ""));
					node["mipNum"] = mipNum;
					std::string clusterId = seq.name_.substr(seq.name_.find("].") + 2,
							seq.name_.rfind("_R") - (seq.name_.find("].") + 2));
					node["color"] = outColorsStrs[bib::lexical_cast<uint32_t>(clusterId)];
				} else {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__
							<< ": Error, should have meta data for mipFam" << std::endl;
					throw std::runtime_error { ss.str() };
				}
			}
		} else {
			std::cerr << "Sample name wasn't found :" << sample << std::endl;
		}
	} else {
		std::cerr << "Region name wasn't found :" << region << std::endl;
	}

	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);

}

std::shared_ptr<restbed::Resource> mav::regionInfoForSampPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showRegionInfoForSamp" },
			{"region", UrlPathFactory::pat_wordNumsDash_},
			{"sample", UrlPathFactory::pat_wordNumsDash_}  }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						regionInfoForSampPageHandler(session);
					}));
	return resource;
}
std::shared_ptr<restbed::Resource> mav::getMipOverlapGraphData(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"mipOverlapGraphData" },
			{"region", UrlPathFactory::pat_wordNumsDash_},
			{"sample", UrlPathFactory::pat_wordNumsDash_}  }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipOverlapGraphDataHandler(session);
					}));
	return resource;
}

/*
void mav::showGeneInfoForSampData(std::string geneName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"geneName", geneName},{"sampName", sampName}}), std::cout, debug_);
	if(popClusInfoBySample_.find(sampName) == popClusInfoBySample_.end()){
		std::cerr << __PRETTY_FUNCTION__ << ": Sample Name wasn't found: " << sampName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get();
	}else{
		auto search = pages_.find("oneGeneOneSamp");
		response().out() << search->second.get();
	}
}

void mav::mipOverlapGraphData(std::string geneName, std::string sampleName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, { { "geneName",
			geneName }, {"sampleName", sampleName} }), std::cout, debug_);
	ret_json();
	if (!bib::in(geneName, mipMaster_->getMipGroupings())) {
		std::cerr << "Gene name wasn't found :" << geneName << std::endl;
		response().out() << "";
	} else {

		auto mipFams = mipMaster_->mips_->getMipFamsForRegion(geneName);
		std::vector<std::shared_ptr<seqInfo>> finalSeqs;
		for(const auto & mipName : mipFams ){
			if(bib::in(mipName, finalHapsByTarSamp_)){
				if(bib::in(sampleName, finalHapsByTarSamp_.at(mipName))){
					addOtherVec(finalSeqs, finalHapsByTarSamp_.at(mipName).at(sampleName).getPtrs<seqInfo>());
				}else{
					std::cerr << "Mip fam name wasn't found :" << mipName << std::endl;
				}
			}else{
				std::cerr << "Mip fam name wasn't found :" << mipName << std::endl;
			}
		}
		uint32_t overlapMin = 5;
		uint64_t maxSize = 0;
		std::map<uint32_t, std::vector<std::shared_ptr<seqInfo>>> seqsByMipNum;
		SeqOverlapGraph graph;


		for(const auto & seq : finalSeqs){
			readVec::getMaxLength(seq, maxSize);

			std::unordered_map<std::string, std::string> meta;
			seq->processNameForMeta(meta);
			if(bib::has(meta, "mipFam")){
				auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
				uint32_t mipNum = bib::lexical_cast<uint32_t>(bib::replaceString(mipFamToks.back(), "mip", ""));
				seqsByMipNum[mipNum].push_back(seq);
			}else{
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam" << std::endl;
				throw std::runtime_error{ss.str()};
			}
		}

		aligner alignerObj(maxSize, gapScoringParameters(5,1,0,0,0,0), substituteMatrix(2,-2),
				false);
		comparison noErrors;
		//check to see if target should be reversed
		for (const auto & mips : seqsByMipNum) {
			auto search = seqsByMipNum.find(mips.first + 1);
			if (seqsByMipNum.end() != search) {
				for(const auto & seq1 : mips.second){
					uint32_t forCount = 0;
					uint32_t revCount = 0;
					for(const auto & seq2 : search->second){
						alignerObj.alignCacheGlobal(*seq1, *seq2);
						alignerObj.profilePrimerAlignment(*seq1, *seq2);
						auto forComp = alignerObj.comp_;
						auto seq2Rev = *seq2;
						seq2Rev.reverseComplementRead(true, true);
						alignerObj.alignCacheGlobal(*seq1, seq2Rev);
						alignerObj.profilePrimerAlignment(*seq1, seq2Rev);
						bool passFor = noErrors.passErrorProfile(forComp);
						bool passRev = noErrors.passErrorProfile(alignerObj.comp_);
						if(passFor && passRev){
							if(forComp.distances_.query_.covered_ >= alignerObj.comp_.distances_.query_.covered_){
								++forCount;
							}else{
								++revCount;
							}
						}else if(passFor){
							++forCount;
						}else if(passRev){
							++revCount;
						}
					}
					if(revCount > forCount){
						for(auto & seq : search->second){
							seq->reverseComplementRead(true, true);
						}
					}
				}
			}
			bool notProcessed = false;
			if(notProcessed){
				for(auto & seq : mips.second){
					seq->frac_ = 1.0/len(mips.second);
				}
			}
			//add to graph
			for(const auto & seq : mips.second){
				graph.addNode(seq);
			}
		}
		for (const auto & mips : seqsByMipNum) {
			auto search = seqsByMipNum.find(mips.first + 1);
			if (seqsByMipNum.end() != search) {
				for(const auto & seq1 : mips.second){
					for(const auto & seq2 : search->second){
						alignerObj.alignCacheGlobal(*seq1, *seq2);
						alignerObj.profilePrimerAlignment(*seq1, *seq2);
						if(noErrors.passErrorProfile(alignerObj.comp_)
								&& alignerObj.comp_.distances_.query_.covered_ >= overlapMin){
							bool seq1GapFront = '-' == alignerObj.alignObjectA_.seqBase_.seq_.front();
							bool seq1GapBack  = '-' == alignerObj.alignObjectA_.seqBase_.seq_.back();
							bool seq2GapFront = '-' == alignerObj.alignObjectB_.seqBase_.seq_.front();
							bool seq2GapBack  = '-' == alignerObj.alignObjectB_.seqBase_.seq_.back();
							if ((seq1GapFront && seq1GapBack)
									|| (seq2GapFront && seq2GapBack)) {
								//no overlap
							} else {
								if (seq1GapBack || seq2GapFront) {
									graph.addEdge(seq1->name_, seq2->name_);
								} else {
									graph.addEdge(seq2->name_, seq1->name_);
								}
							}
						}
					}
				}
			}
		}
		auto sankeyOutput = graph.createSankeyOutput();
		uint32_t maxNum = 0;
		for(const auto & seqs : seqsByMipNum){
			if(len(seqs.second) > maxNum){
				maxNum = len(seqs.second);
			}
		}
		auto outColors = bib::njhColors(maxNum);
		bibseq::VecStr outColorsStrs;
		outColorsStrs.reserve(outColors.size());
		for(const auto & c : outColors) {
			outColorsStrs.emplace_back("#" + c.hexStr_);
		}


		for (auto & node : sankeyOutput["nodes"]) {
			std::unordered_map<std::string, std::string> meta;
			seqInfo seq(node["name"].asString(), "A");
			seq.processNameForMeta(meta);
			if (bib::has(meta, "mipFam")) {
				auto mipFamToks = bib::tokenizeString(meta.at("mipFam"), "_");
				uint32_t mipNum = bib::lexical_cast<uint32_t>(
						bib::replaceString(mipFamToks.back(), "mip", ""));
				node["mipNum"] = mipNum;
				std::string clusterId = seq.name_.substr(seq.name_.find("].") + 2, seq.name_.rfind("_R") - (seq.name_.find("].") + 2));
				node["color"] = outColorsStrs[bib::lexical_cast<uint32_t>(clusterId)];
			} else {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, should have meta data for mipFam"
						<< std::endl;
				throw std::runtime_error { ss.str() };
			}
		}
		response().out() << sankeyOutput;
	}
}
*/

}  // namespace bibseq
