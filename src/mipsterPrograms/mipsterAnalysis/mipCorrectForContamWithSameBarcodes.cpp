/*
 * mipCorrectForContamWithSameBarcodes.cpp
 *
 *  Created on: Feb 18, 2016
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
#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"

namespace njhseq {

void correctForContamForMipFam(const mipCorrectForContamWithSameBarcodesPars & pars,
		const SetUpMaster & mipMaster){
	//std::cout << __FILE__ << " " << __LINE__ << std::endl;
	auto mipFamTars = mipMaster.mips_->getMipsForFamily(pars.mipName);
	//key1 = samp, key2 = mipFam, val = vec of barcodes to remove
	std::unordered_map<std::string, std::unordered_map<std::string, VecStr>> barcodesToRemove;
	for(const auto & mipName : mipFamTars){
		std::vector<MipFamSamp> barcodesExists;
		for (const auto & samp : mipMaster.names_->samples_) {
			MipFamSamp mipSamp(mipName, samp);
			if (bfs::exists(mipMaster.pathMipSampBarCorBars(MipFamSamp(pars.mipName, samp), mipSamp.mipFam_))) {
				barcodesExists.emplace_back(mipSamp);
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::unordered_map<std::string, strCounter> barCounts;
		for (const auto & bars : barcodesExists) {
			std::ifstream barFile(mipMaster.pathMipSampBarCorBars(MipFamSamp(pars.mipName, bars.samp_), bars.mipFam_).string());
			if (!barFile) {
				std::stringstream ss;
				ss << "Error in opening " << mipMaster.pathMipSampBarCorBars(MipFamSamp(pars.mipName, bars.samp_), bars.mipFam_)
						<< std::endl;
				throw std::runtime_error { ss.str() };
			}
			std::string line = "";
			while (njh::files::crossPlatGetline(barFile, line)) {
				if (!njh::beginsWith(line, "Barcode")) {
					auto toks = tokenizeString(line, "\t");
					if (2 != toks.size()) {
						std::stringstream ss;
						ss << "Error in reading barcode line: " << std::endl;
						ss << line << std::endl;
						ss << "lines should have two columns, this line had " << toks.size()
								<< std::endl;
						throw std::runtime_error { ss.str() };
					}
					barCounts[bars.samp_].increaseCountByString(toks[0],
							estd::stou(toks[1]));
				}
			}
		}
		//std::cout << __FILE__ << " " << __LINE__ << std::endl;
		std::unordered_map<std::string, VecStr> sharedBars;
		for (const auto & sampCount : barCounts) {
			for (const auto & bars : sampCount.second.counts_) {
				sharedBars[bars.first].emplace_back(sampCount.first);
			}
		}
		table sharedTab(VecStr { "SharedBar", "sample", "readCnt" });

		for (const auto & shared : sharedBars) {
			if (shared.second.size() > 1) {
				uint32_t bestCount = 0;
				VecStr bestSamp;
				for (const auto & samp : shared.second) {
					if(barCounts[samp].counts_[shared.first] > bestCount){
						bestCount = barCounts[samp].counts_[shared.first];
						bestSamp = {samp};
					}else if(barCounts[samp].counts_[shared.first] == bestCount){
						bestSamp.emplace_back(samp);
					}
				}
				VecStr samplesToRemove = shared.second;
				//if the best sample has no ties get rid of all others
				//if a tie for best remove from all samples
				if(bestSamp.size() == 1){
					removeElement(samplesToRemove, bestSamp.front());
				}
				for(const auto & samp : samplesToRemove){
					if(barCounts[samp].counts_[shared.first] <= pars.readCutOff){
						barcodesToRemove[samp][mipName].emplace_back(shared.first);
					}
				}
			}
		}
	}
//	std::cout << __FILE__ << " " << __LINE__ << std::endl;

	for (const auto & samps : barcodesToRemove) {
		auto hapFile = mipMaster.pathMipSampBarCorHap(MipFamSamp(pars.mipName, samps.first));
		bfs::path tempHapFile = mipMaster.pathMipSampBarCorHap(MipFamSamp(pars.mipName, samps.first)).string() + "_temp.gz";
		bfs::rename(hapFile, tempHapFile);
//		std::cout << "tempHapFile: " << tempHapFile << std::endl;
//		std::cout << "hapFile :" << hapFile << std::endl;
		SeqIOOptions opts = SeqIOOptions::genFastqInOutGz(tempHapFile,hapFile,true);
//		std::cout << njh::json::toJson(opts) << std::endl;
		SeqIO reader(opts);
		reader.openIn();
		reader.openOut();
		readObject read;
		while(reader.readNextRead(read)){
			//std::cout << __FILE__ << " " << __LINE__ << std::endl;
			read.processNameForMeta();
			if (read.containsMeta("bar") && read.containsMeta("mipTar")) {
				auto search = samps.second.find(read.getMeta("mipTar"));
				if(search != samps.second.end()){
					if(njh::in(read.getMeta("bar"), search->second)){
						read.addMeta("remove", true, true);
						read.resetMetaInName();
					}else{
						read.addMeta("remove", false, true);
						read.resetMetaInName();
					}
				}else{
					read.addMeta("remove", false, true);
					read.resetMetaInName();
				}
				reader.write(read);
			} else {
				std::stringstream ss;
				ss << "Error in " << __PRETTY_FUNCTION__
						<< ", error in parsing read for meta data: bar and mipTar" << std::endl;
				ss << "Current meta data is: "
						<< njh::conToStr(getVectorOfMapKeys(read.meta_.meta_), ",") << std::endl;
				throw std::runtime_error { ss.str() };
			}
		}
		reader.closeIn();
		reader.closeOut();
		bfs::remove(tempHapFile);
	}
}

int mipsterAnalysisRunner::mipCorrectForContamWithSameBarcodes(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterAnalysisSetUp setUp(inputCommands);
	mipCorrectForContamWithSameBarcodesPars pars;
	setUp.setUpMipCorrectForContamWithSameBarcodes(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);

	correctForContamForMipFam(pars, mipMaster);

	return 0;
}

int mipsterAnalysisRunner::mipCorrectForContamWithSameBarcodesMultiple(
		const njh::progutils::CmdArgs & inputCommands) {
	// parameters
	mipsterAnalysisSetUp setUp(inputCommands);
	mipCorrectForContamWithSameBarcodesParsMultiple pars;
	setUp.setUpMipCorrectForContamWithSameBarcodesMultiple(pars);
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	auto warnings = mipMaster.checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipMaster.mips_->setAllWiggleRoomInArm(pars.wiggleRoom);
	Json::Value logInfo;
	std::ofstream logFile;
	openTextFile(logFile,njh::files::make_path(
			mipMaster.directoryMaster_.logsDir_, pars.logFilename), ".json",
			pars.overWriteLog, true);
	logInfo["date"] = getCurrentDate();
	logInfo["workingDir"] = inputCommands.workingDir_;
	logInfo["command"] = inputCommands.commandLine_;
	std::vector<std::thread> threads;
	njh::concurrent::LockableQueue<std::string> mipsQueue(
			mipMaster.names_->mips_);
	std::mutex logMut;
	auto & pairingLog = logInfo["PerMip"];
	auto runCorrectForContamWithSameBarcodes = [&logMut,&pairingLog](njh::concurrent::LockableQueue<std::string>& mipsQueue,
			const mipCorrectForContamWithSameBarcodesParsMultiple& pars,
			const SetUpMaster & mipMaster){
		std::string mipName = "";
		while (mipsQueue.getVal(mipName)) {
			njh::stopWatch watch;
			auto currentSampPars = pars.createForMip(mipName);
			correctForContamForMipFam(currentSampPars,mipMaster);
			{
				std::lock_guard<std::mutex> logLock(logMut);
				auto & currentPair = pairingLog[mipName];
				currentPair["mip"] = mipName;
				currentPair["runTime"] = watch.totalTimeFormatted(6);
			}
		}
	};

	for (uint32_t t = 0; t < pars.numThreads; ++t) {
		threads.emplace_back(
				std::thread(runCorrectForContamWithSameBarcodes, std::ref(mipsQueue), std::cref(pars),
						std::cref(mipMaster)));
	}

	for (auto & t : threads) {
		t.join();
	}
	logInfo["totalRunTime"] = setUp.timer_.totalTimeFormatted(6);
	logFile << logInfo;

	return 0;
}

}  // namespace njhseq
