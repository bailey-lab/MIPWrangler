/*
 * mav_oneMapAllSamps.cpp
 *
 *  Created on: Feb 16, 2016
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

#include "mav.hpp"

namespace njhseq {



void mav::oneMipAllSampsPageHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	if (njh::in(mipFam,popClusInfoByTar_)) {
		auto body = genHtmlDoc(rootName_, pages_.at("oneMipAllSampsInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	} else {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for mip "
				<< mipFam << "\n";
		ss << "Redirecting...." << "\n";
		redirect(session, ss.str());
	}
}

VecStr genColorsStrsForNamesMultiples(VecStr popNames){
	removeDuplicates(popNames);
	auto popColors = getColorsForNames(popNames);
	VecStr popColorsStrs;
	auto keys = getVectorOfMapKeys(popColors);
	njh::sort(keys);
	for(const auto & pColorKey : keys){
		popColorsStrs.emplace_back(popColors[pColorKey].getHexStr());
	}
	return popColorsStrs;
}

void mav::getOneMipAllSampsDataPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "samples" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		VecStr samples = njh::json::jsonArrayToVec<std::string>(postData["samples"],
				[](const Json::Value & val) {return val.asString();});
		auto search = popClusInfoByTar_.find(mipFam);
		if (search != popClusInfoByTar_.end()) {
			auto containsSampName = [&samples](const std::string & str) {
				return njh::in(str, samples);
			};
			auto trimedTab = search->second.get().extractByComp("s_sName",
					containsSampName);
			trimedTab.sortTable("s_sName", "c_clusterID", false);
			//,,,
			ret = tableToJsonByRow(trimedTab, "s_Sample", VecStr { "s_Sample",
					"h_popUID", "c_readCnt", "c_barcodeCnt" });
			/*
			ret = tableToJsonByRow(trimedTab, "s_Sample", VecStr { "s_Sample",
					"p_sampleTotal", "p_barcodeTotal", "p_finalHaplotypeNumber",
					"h_popUID", "h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac",
					"h_barcodeCnt", "h_barcodeFrac", "s_sName", "s_usedTotalClusterCnt",
					"s_usedTotalBarcodeCnt", "c_clusterID", "c_name", "c_readCnt",
					"c_readFrac", "c_barcodeCnt", "c_barcodeFrac" });
			*/
			auto popUIDs = trimedTab.getColumn("h_popUID");
			removeDuplicates(popUIDs);
			njh::sort(popUIDs);
			auto popColorsStrs = genColorsStrsForNamesMultiples(popUIDs);
			ret["popColors"] = njh::json::toJson(popColorsStrs);
			ret["popUIDs"] = njh::json::toJson(popUIDs);
			auto popCounts = njhseq::countVec(trimedTab.getColumn("h_popUID"));
			auto popColors = njh::njhColors(popCounts.size());
		}
	}
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getOneMipAllSampsDataHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getOneMipAllSampsDataPostHandler(ses, body);
					}));
}

void mav::getOneMipPopSeqsPostHandler(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		uint32_t sesUid = std::numeric_limits<uint32_t>::max();
		//check to see if there is a session already started associated with this seq
		if (!postData.isMember("sessionUID")) {
			sesUid = startSeqCacheSession();
		} else {
			sesUid = postData["sessionUID"].asUInt();
		}
		VecStr popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"],
				[](const Json::Value & val) {return val.asString();});
		seqsBySession_[sesUid]->cache_.at(mipFam).reload();
		seqsBySession_[sesUid]->cache_.at(mipFam).toggleSeqs(
				[&popUIDs](const readObject & seq) {
					return njh::in(seq.seqBase_.getStubName(false), popUIDs);
				});
		//make sure seqs aren't empty, viewer doesn't know how to handle that, if it isn't make sure to remove the placeholder seq if it is there
		seqsBySession_[sesUid]->cache_.at(mipFam).ensureNonEmptyReads();
		ret = seqsBySession_[sesUid]->getJson(mipFam);
		ret["sessionUID"] = njh::json::toJson(sesUid);
	}
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getOneMipPopSeqsHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getOneMipPopSeqsPostHandler(ses, body);
					}));
}

void mav::getOneMipAllSampsPopDataPostHandler(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "popUIDs" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	}else{
		VecStr popUIDs = njh::json::jsonArrayToVec<std::string>(postData["popUIDs"], [](const Json::Value & val){ return val.asString();});
		auto search = popClusPopInfoByTar_.find(mipFam);
		if (search == popClusPopInfoByTar_.end()) {
			std::cerr << __PRETTY_FUNCTION__ << ": Couldn't find mipName: " << mipFam
					<< std::endl;
		} else {
			auto tab = search->second.get();
			auto containsPopUID = [&popUIDs](const std::string & str) {
				return njh::in(str, popUIDs);
			};
			auto trimedTab = search->second.get().extractByComp("h_popUID",
					containsPopUID);
			ret = tableToJsonByRow(trimedTab, "h_popUID", VecStr {
					"p_geneName", "p_targetName", "p_sampleTotal", "p_totalInputClusters",
					"p_readTotal", "p_barcodeTotal", "p_finalHaplotypeNumber", "h_popUID",
					"h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac",
					"h_meanBarcodeFrac", "h_readFrac", "h_barcodeCnt", "h_barcodeFrac",
					"h_inputNames", "h_seq", "h_qual" }, VecStr { "p_sampleTotal",
					"p_totalInputClusters", "p_readTotal", "p_barcodeTotal",
					"p_finalHaplotypeNumber" });
		}
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getOneMipAllSampsPopDataHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	//get_header("Content-Length", content_length)
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes & body)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getOneMipAllSampsPopDataPostHandler(ses, body);
					}));
}


std::shared_ptr<restbed::Resource> mav::oneMipAllSampsPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(
			UrlPathFactory::createUrl( { { rootName_ }, { "showOneMipAllSampsData" }, {
					"mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						oneMipAllSampsPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getOneMipAllSampsData() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getOneMipAllSampsData" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getOneMipAllSampsDataHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getOneMipPopSeqs(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getOneMipPopSeqs" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getOneMipPopSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getOneMipAllSampsPopData(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getOneMipAllSampsPopData" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getOneMipAllSampsPopDataHandler(session);
					}));
	return resource;
}




}  // namespace njhseq
