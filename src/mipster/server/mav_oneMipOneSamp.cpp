/*
 * mav_oneMipOneSamp.cpp
 *
 *  Created on: Feb 16, 2016
 *      Author: nick
 */

#include "mav.hpp"

namespace bibseq {

void mav::oneMipOneSampDataPageHandler(
		std::shared_ptr<restbed::Session> session) {
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	auto sample = request->get_path_parameter("sample");
	auto search = popClusInfoByTar_.find(mipFam);
	bool redirecting = false;
	std::stringstream ss;
	if (search == popClusInfoByTar_.end()) {
		redirecting = true;
		ss << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipFam
				<< std::endl;
	} else {
		if (!bib::in(sample, search->second.get().getColumn("s_sName"))) {
			redirecting = true;
			ss << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sample
					<< " for mipName: " << mipFam << std::endl;
		}
	}
	if (redirecting) {
		redirect(session, ss.str());
	} else {
		auto body = genHtmlDoc(rootName_, pages_.at("oneMipOneSampInfo.js"));
		const std::multimap<std::string, std::string> headers =
				HeaderFactory::initiateTxtHtmlHeader(body);
		session->close(restbed::OK, body, headers);
	}
}

void mav::getMipOneSampOriginalSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	auto sample = request->get_path_parameter("sample");
	Json::Value postData;
	if(!body.empty()){
		postData = bib::json::parse(std::string(body.begin(), body.end()));
	}
	Json::Value ret;

	std::string searchTerm = mipFam + "_" + sample + "_original";
	if (bib::in(searchTerm, seqs_->cache_)) {
		uint32_t sesUid = std::numeric_limits<uint32_t>::max();
		//check to see if there is a session already started associated with this seq
		if (!postData.isMember("sessionUID")) {
			sesUid = startSeqCacheSession();
		} else {
			sesUid = postData["sessionUID"].asUInt();
		}
		ret = seqsBySession_[sesUid]->getJson(searchTerm);
		ret["sessionUID"] = bib::json::toJson(sesUid);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": couldn't find " << searchTerm
				<< " Sample and Mip fam seqs not found :" << sample << " - " << mipFam
				<< std::endl;
	}

	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getMipOneSampOriginalSeqsHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = 0;
	request->get_header("Content-Length", content_length);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes &)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getMipOneSampOriginalSeqsPostHandler(ses, body);
					}));
}

void mav::getMipOneSampFinalSeqsPostHandler(std::shared_ptr<restbed::Session> session,
		const restbed::Bytes & body){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	auto sample = request->get_path_parameter("sample");
	Json::Value postData;
	if(!body.empty()){
		postData = bib::json::parse(std::string(body.begin(), body.end()));
	}
	Json::Value ret;

	std::string searchTerm = mipFam + "_" + sample + "_final";
	if (bib::in(searchTerm, seqs_->cache_)) {
		uint32_t sesUid = std::numeric_limits<uint32_t>::max();
		//check to see if there is a session already started associated with this seq
		if (!postData.isMember("sessionUID")) {
			sesUid = startSeqCacheSession();
		} else {
			sesUid = postData["sessionUID"].asUInt();
		}
		ret = seqsBySession_[sesUid]->getJson(searchTerm );
		ret["sessionUID"] = bib::json::toJson(sesUid);
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": couldn't find " << searchTerm << " Sample and Mip fam seqs not found :"
				<< sample << " - " << mipFam << std::endl;
	}

	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}

void mav::getMipOneSampFinalSeqsHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	const auto request = session->get_request();
	auto heads = request->get_headers();
	size_t content_length = 0;
	request->get_header("Content-Length", content_length);
	session->fetch(content_length,
			std::function<
					void(std::shared_ptr<restbed::Session>, const restbed::Bytes &)>(
					[this](std::shared_ptr<restbed::Session> ses, const restbed::Bytes & body) {
						getMipOneSampFinalSeqsPostHandler(ses, body);
					}));
}

void mav::getOneMipOneSampsDataHandler(std::shared_ptr<restbed::Session> session){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto request = session->get_request();
	auto mipFam = request->get_path_parameter("mipFam");
	auto sample = request->get_path_parameter("sample");
	Json::Value ret;
	auto search = popClusInfoByTar_.find(mipFam);
	if (search != popClusInfoByTar_.end()) {
		if(bib::in(sample,search->second.get().getColumn("s_sName"))){
			auto containsSampName = [&sample](const std::string & str) {
				return str == sample;
			};
			auto trimedTab = search->second.get().extractByComp("s_sName", containsSampName);
			auto popUids = trimedTab.getColumnLevels("h_popUID");
			auto containsPopUID = [&popUids](const std::string & str) {
				return bib::in(str, popUids);
			};
			ret = tableToJsonByRow(
					popClusPopInfoByTar_.find(mipFam)->second.get().extractByComp("h_popUID",
							containsPopUID), "h_popUID", VecStr {
				"p_geneName", "p_targetName", "p_sampleTotal", "p_totalInputClusters",
				"p_readTotal", "p_barcodeTotal", "p_finalHaplotypeNumber", "h_popUID",
				"h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac",
				"h_meanBarcodeFrac", "h_readFrac", "h_barcodeCnt", "h_barcodeFrac",
				"h_inputNames", "h_seq", "h_qual" }, VecStr {
							"p_sampleTotal", "p_totalInputClusters", "p_readTotal",
							"p_barcodeTotal", "p_finalHaplotypeNumber" });
		} else {
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :"
					<< sample << " for mipName: " << mipFam << std::endl;
		}
	} else {
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :"
				<< mipFam << std::endl;
	}
	auto retBody = bib::json::writeAsOneLine(ret);
	std::multimap<std::string, std::string> headers =
			HeaderFactory::initiateAppJsonHeader(retBody);
	headers.emplace("Connection", "close");
	session->close(restbed::OK, retBody, headers);
}


std::shared_ptr<restbed::Resource> mav::oneMipOneSampDataPage(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"showOneMipOneSampData" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ },
			{ "sample", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						oneMipOneSampDataPageHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getMipOneSampOriginalSeqs(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getMipOneSampOriginalSeqs" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ },
			{ "sample", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipOneSampOriginalSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getMipOneSampFinalSeqs(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getMipOneSampFinalSeqs" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ },
			{ "sample", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("POST",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getMipOneSampFinalSeqsHandler(session);
					}));
	return resource;
}

std::shared_ptr<restbed::Resource> mav::getOneMipOneSampsData(){
	auto mess = messFac_->genLogMessage(__PRETTY_FUNCTION__);
	auto resource = std::make_shared<restbed::Resource>();
	resource->set_path(UrlPathFactory::createUrl( { { rootName_ }, {
			"getOneMipOneSampsData" },
			{ "mipFam", UrlPathFactory::pat_wordNumsDash_ },
			{ "sample", UrlPathFactory::pat_wordNumsDash_ }}));
	resource->set_method_handler("GET",
			std::function<void(std::shared_ptr<restbed::Session>)>(
					[this](std::shared_ptr<restbed::Session> session) {
						getOneMipOneSampsDataHandler(session);
					}));
	return resource;
}

/*
void mav::showOneMipOneSampsData(std::string mipName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"mipName", mipName},{"sampName", sampName}}), std::cout, debug_);
	auto search = popClusInfoByTar_.find(mipName);
	if (search != popClusInfoByTar_.end()) {
		if(bib::in(sampName,search->second.get().getColumn("s_sName"))){
			auto pageSearch = pages_.find("oneMipOneSampInfo");
			response().out() << pageSearch->second.get();
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sampName << " for mipName: " << mipName << std::endl;
			auto search = pages_.find("redirectPage");
			response().out() << search->second.get();
		}
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipName << std::endl;
		auto search = pages_.find("redirectPage");
		response().out() << search->second.get();
	}
}

void mav::getOneMipOneSampsData(std::string mipName, std::string sampName){
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, {{"mipName", mipName},{"sampName", sampName}}), std::cout, debug_);
	ret_json();
	Json::Value ret;
	std::map<std::string, std::string> testMap;


	auto search = popClusInfoByTar_.find(mipName);
	if (search != popClusInfoByTar_.end()) {
		if(bib::in(sampName,search->second.get().getColumn("s_sName"))){
			auto containsSampName = [&sampName](const std::string & str) {
				return str == sampName;
			};
			auto trimedTab = search->second.get().extractByComp("s_sName", containsSampName);
			trimedTab.sortTable("c_clusterID", false);
			ret = tableToJsonByRow(trimedTab, "h_popUID", VecStr { "s_Sample",
				"p_sampleTotal", "p_barcodeTotal", "p_finalHaplotypeNumber", "h_popUID",
				"h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac", "h_barcodeCnt",
				"h_barcodeFrac", "s_sName", "s_usedTotalClusterCnt",
				"s_usedTotalBarcodeCnt", "c_clusterID", "c_name", "c_readCnt",
				"c_readFrac", "c_barcodeCnt", "c_barcodeFrac" });
			auto popCounts = bibseq::countVec(trimedTab.getColumn("h_popUID"));
			auto popColors = bib::njhColors(popCounts.size());
			VecStr popColorsStrs(popColors.size(), "");
			uint32_t count = 0;
			uint32_t halfCount = 0;
			for(const auto & cPos : iter::range(popColors.size())) {
				uint32_t pos = 0;
				if(cPos %2 == 0) {
					pos = popColors.size()/2 + halfCount;
					++halfCount;
				} else {
					pos = count;
					++count;
				}
				popColorsStrs[cPos] = "#" + popColors[pos].hexStr_;
			}
			ret["popColors"] = bib::json::toJson(popColorsStrs);
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sampName << " for mipName: " << mipName << std::endl;
		}
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipName << std::endl;
	}
	response().out() << ret;
}

void mav::getOneMipOneSampsPopData(std::string mipName, std::string sampName) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, { { "mipName",
			mipName }, { "sampName", sampName } }), std::cout, debug_);
	ret_json();
	Json::Value ret;
	auto search = popClusInfoByTar_.find(mipName);
	if (search != popClusInfoByTar_.end()) {
		if(bib::in(sampName,search->second.get().getColumn("s_sName"))){
			auto containsSampName = [&sampName](const std::string & str) {
				return str == sampName;
			};
			auto trimedTab = search->second.get().extractByComp("s_sName", containsSampName);
			auto popUids = trimedTab.getColumnLevels("h_popUID");
			auto containsPopUID = [&popUids](const std::string & str) {
				return bib::in(str, popUids);
			};
			ret = tableToJsonByRow(
					popClusPopInfoByTar_.find(mipName)->second.get().extractByComp("h_popUID",
							containsPopUID), "h_popUID", VecStr {
				"p_geneName", "p_targetName", "p_sampleTotal", "p_totalInputClusters",
				"p_readTotal", "p_barcodeTotal", "p_finalHaplotypeNumber", "h_popUID",
				"h_sampleCnt", "h_sampleFrac", "h_medianBarcodeFrac",
				"h_meanBarcodeFrac", "h_readFrac", "h_barcodeCnt", "h_barcodeFrac",
				"h_inputNames", "h_seq", "h_qual" }, VecStr {
							"p_sampleTotal", "p_totalInputClusters", "p_readTotal",
							"p_barcodeTotal", "p_finalHaplotypeNumber" });
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sampName << " for mipName: " << mipName << std::endl;
		}
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipName << std::endl;
	}
	response().out() << ret;
}

void mav::getMipOneSampOriginalSeqs(std::string mipName, std::string sampName) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, { { "mipName",
			mipName }, { "sampName", sampName } }), std::cout, debug_);
	auto search = originalHapsByTarBySamp_.find(mipName);
	ret_json();
	Json::Value ret;
	if(search != originalHapsByTarBySamp_.end() ){
		auto sampSearch = search->second.find(sampName);
		if(sampSearch != search->second.end()){
			std::string mipUid = mipName +"_" +sampName + "_OriginalSeqs";
			if (!seqs_->containsRecord(mipUid) || !seqs_->recordValid(mipUid)
					|| sampSearch->second.outDated()) {
				seqs_->updateAddCache(mipUid,
						std::make_shared<std::vector<readObject>>(sampSearch->second.get<readObject>()));
			}
			ret = seqs_->getJson(mipUid);
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sampName << " for mipName: " << mipName << std::endl;
		}
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipName << std::endl;
	}
	response().out() << ret;
}

void mav::getMipOneSampFinalSeqs(std::string mipName, std::string sampName) {
	bib::scopedMessage mess(messStrFactory(__PRETTY_FUNCTION__, { { "mipName",
			mipName }, { "sampName", sampName } }), std::cout, debug_);
	auto search = finalHapsByTarSamp_.find(mipName);
	ret_json();
	Json::Value ret;
	if(search != finalHapsByTarSamp_.end() ){
		auto sampSearch = search->second.find(sampName);
		if(sampSearch != search->second.end()){
			std::string mipUid = mipName +"_" +sampName + "_FinalSeqs";
			if (!seqs_->containsRecord(mipUid) || !seqs_->recordValid(mipUid)
					|| sampSearch->second.outDated()) {
				seqs_->updateAddCache(mipUid,
						std::make_shared<std::vector<readObject>>(sampSearch->second.get<readObject>()));
			}

			ret = seqs_->getJson(mipUid);
		}else{
			std::cerr << __PRETTY_FUNCTION__ << ": " << "Samp name wasn't found :" << sampName << " for mipName: " << mipName << std::endl;
		}
	}else{
		std::cerr << __PRETTY_FUNCTION__ << ": " << "Mip name wasn't found :" << mipName << std::endl;
	}
	response().out() << ret;
}
*/
}  // namespace bibseq
