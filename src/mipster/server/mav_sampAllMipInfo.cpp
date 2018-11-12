/*
 * mav_sampAllMipInfo.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: nick
 */


#include "mav.hpp"

namespace njhseq {


void mav::oneSampAllMipDataPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	if (popClusInfoBySample_.find(sample) == popClusInfoBySample_.end()) {
		std::stringstream ss;
		ss << "Sample Name wasn't found: " << sample << "\n";
		ss << "Redirecting....";
		redirect(session, ss.str());
	} else {
		auto body = genHtmlDoc(rootName_, pages_.at("oneSampAllMipInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}
}

void mav::mipFamNamesForSampHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	Json::Value ret;
	auto search = popClusInfoBySample_.find(sample);
	if (search == popClusInfoBySample_.end()) {
		std::cerr << __PRETTY_FUNCTION__ << ": Sample Name wasn't found: " << sample << "\n";
	} else {
		ret["mipFamilies"] = njh::json::toJson(
								search->second.get().getColumnLevels("p_targetName"));
	}

	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getOneSampAllMipDataPostHandler(
		std::shared_ptr<restbed::Session> session, const restbed::Bytes & body) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto sample = request->get_path_parameter("sample");
	const auto postData = njh::json::parse(std::string(body.begin(), body.end()));
	njh::json::MemberChecker checker(postData);
	Json::Value ret;
	if (checker.failMemberCheck( { "mipFamilies" }, __PRETTY_FUNCTION__)) {
		std::cerr << checker.message_.str() << std::endl;
	} else {
		VecStr mips = njh::json::jsonArrayToVec<std::string>(postData["mipFamilies"], [](const Json::Value & val){return val.asString();});
		//printVector(mips);
		auto search = popClusInfoBySample_.find(sample);
		if (popClusInfoBySample_.find(sample) == popClusInfoBySample_.end()) {
			std::cerr << "Sample Name wasn't found: " << sample << std::endl;
		} else {
			auto containsMipName = [&mips](const std::string & str) {
				return njh::in(str, mips);
			};
			//printVector(mips);
			auto trimedTab = search->second.get().extractByComp("p_targetName",
					containsMipName);
			//std::cout << "samp: " << sampName << std::endl;
			//std::cout << "size: " << trimedTab.content_.size() << std::endl;
			//trimedTab.outPutContentOrganized(std::cout);
			trimedTab.sortTable("p_geneName", "p_targetName", "c_clusterID", false);
			ret = tableToJsonByRow(trimedTab, "p_targetName", VecStr { "p_targetName",
					"h_popUID", "c_readCnt", "c_barcodeCnt" });
			/*
			ret = tableToJsonByRow(trimedTab, "p_targetName", VecStr {
					"p_geneName", "p_targetName", "p_sampleTotal", "p_totalInputClusters",
					"p_readTotal", "p_barcodeTotal", "p_finalHaplotypeNumber", "h_popUID",
					"h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac", "h_readFrac",
					"h_barcodeCnt", "h_barcodeFrac s_sName", "s_usedTotalClusterCnt",
					"s_usedTotalBarcodeCnt", "s_inputTotalReadCnt",
					"s_inputTotalBarcodeCnt", "c_clusterID", "c_name", "c_readCnt",
					"c_readFrac", "c_barcodeCnt", "c_barcodeFrac" });
			*/
			std::unordered_map<std::string, uint32_t> mipClusIdCounts;
			for (const auto & m : trimedTab.getColumn("c_clusterID")) {
				++mipClusIdCounts[m];
			}
			auto outColors = njh::njhColors(mipClusIdCounts.size());
			njhseq::VecStr outColorsStrs;
			outColorsStrs.reserve(outColors.size());
			for (const auto & c : outColors) {
				outColorsStrs.emplace_back("#" + c.hexStr_);
			}
			ret["popColors"] = njh::json::toJson(outColorsStrs);
		}
	}
	auto retBody = njh::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getOneSampAllMipDataHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = request->get_header("Content-Length", 0);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes &)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getOneSampAllMipDataPostHandler(ses, body);
					}));
}

std::shared_ptr<restbed::Resource> mav::oneSampAllMipDataPage() {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showOneSampAllMipData" }, { "sample",
			UrlPathFactory::pat_wordNumsDash_ } }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						oneSampAllMipDataPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::mipFamNamesForSamp(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"mipFamNamesForSamp" }, {"sample", UrlPathFactory::pat_wordNumsDash_}  }));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						mipFamNamesForSampHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getOneSampAllMipData(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getOneSampAllMipData" }, {"sample", UrlPathFactory::pat_wordNumsDash_}  }));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getOneSampAllMipDataHandler(session);
					}));
	return resource;
}

} // namespace njhseq


