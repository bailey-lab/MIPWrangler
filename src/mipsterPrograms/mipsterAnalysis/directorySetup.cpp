/*
 * directorySetup.cpp
 *
 *  Created on: Feb 13, 2016
 *      Author: nick
 */

#include "mipsterAnalysisRunner.hpp"
#include "mipsterAnalysisSetUp.hpp"

namespace njhseq {

int mipsterAnalysisRunner::runGzExtractStitch(const njh::progutils::CmdArgs & inputCommands) {
	mipsterAnalysisSetUp setUp(inputCommands);
	runGzExtractStitchPars pars;
	setUp.setUpRunGzExtractStitch(pars);
	//check for required external programs
	if(pars.usePear){
		requireExternalProgramThrow("pear");
	}else if(pars.usePanda){
		requireExternalProgramThrow("pandaseq");
	}else{
		requireExternalProgramThrow("flash");
	}
#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	//apple zcat is stupid and requires files end with .Z because why not
	//using brew install gnutls instead
	std::string zcatCmd = "gzcat";
#else
	std::string zcatCmd = "zcat";
#endif
	requireExternalProgramThrow(zcatCmd);
	if(pars.trim){
		requireExternalProgramThrow("sickle");
	}
	//set up the mip analysis directory structure
	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	if("" != pars.sampleMetaFnp){
		mipMaster.setMetaData(pars.sampleMetaFnp);
	}
	mipMaster.createDirStructSkeleton();
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << njh::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	pars.dir = njh::appendAsNeededRet(pars.dir.string(), "/");

	if(!bfs::exists(pars.dir)){
		std::stringstream ss;
		ss << njh::bashCT::boldRed(pars.dir.string()) << " doesn't exist";
		throw std::runtime_error {ss.str()};
	}

	std::ofstream logFile;
	openTextFile(logFile, njh::files::make_path( mipMaster.directoryMaster_.logsDir_,
			pars.logFilename), ".json", setUp.pars_.ioOptions_.out_);


	auto files = njh::files::listAllFiles(pars.dir, false, {std::regex{R"(.*.fastq.gz)"}});
	if(setUp.pars_.debug_){
		std::cout << "Files: " << std::endl;
		printOutMapContents(files, "\t", std::cout);
	}

	std::unordered_map<std::string, VecStr> readPairs;
	std::unordered_map<std::string, VecStr> readPairsUnrecognized;
	for(const auto & f : files){
		auto filename = f.first.filename().string();
		std::string sampName = filename.substr(0, filename.find("_"));
		if(njh::in(sampName, mipMaster.names_->samples_)){
			readPairs[sampName].emplace_back(filename);
		}else{
			readPairsUnrecognized[sampName].emplace_back(filename);
		}
	}
	auto keys = getVectorOfMapKeys(readPairs);
	njh::sort(keys);
	VecStr samplesExtracted;
	VecStr samplesEmpty;
	njh::concurrent::LockableQueue<std::string> filesKeys(keys);
	Json::Value logs;
	logs["mainCommand"] = setUp.commands_.commandLine_;
	Json::Value & extractLog = logs["extractLog"];
	std::mutex logsLock;
	const std::string zcatTempCmd = zcatCmd + " " + pars.dir.string() + "REPLACETHIS.gz > "
			+ mipMaster.directoryMaster_.masterDir_.string() + "DIRNAME/REPLACETHIS";
	auto extractFiles =
			[&pars,&samplesExtracted,&samplesEmpty](njh::concurrent::LockableQueue<std::string> & filesKeys,
					const std::string & zcatTempCmd, const std::string & outputDirectory,
					const std::unordered_map<std::string, VecStr>& readPairs,
					Json::Value & logs,std::mutex & logsMut) {
				std::string key = "";
				VecStr currentEmptySamps;
				VecStr currentSamplesExtracted;
				std::unordered_map<std::string, Json::Value> currentLogs;
				while(filesKeys.getVal(key)){
					Json::Value log;
					njh::stopWatch watch;
					std::stringstream out;
					std::stringstream err;
					bool success = true;
					bool allEmpty = true;
					out << key << std::endl;
					auto sampDir = njh::files::makeDir(outputDirectory, njh::files::MkdirPar(key));
					Json::Value zCatLogs;

					for(const auto & f : readPairs.at(key)){
						auto fqName = njh::files::bfs::path(f);
						auto inputName = njh::files::make_path(pars.dir,fqName);
						//if the file is empty copy to an empty directory to avoid concatenating errors
						if(bfs::file_size(inputName) == 0){
							std::string emptyDir = njh::files::makeDirP(sampDir, njh::files::MkdirPar("emptyFiles")).string();
							bfs::copy(inputName, njh::files::join(emptyDir, bfs::path(inputName).filename().string()));
							continue;
						}
						allEmpty = false;
						auto currentZcatTempCmd = njh::replaceString(zcatTempCmd, "DIRNAME", key);
						currentZcatTempCmd = njh::replaceString(currentZcatTempCmd, "REPLACETHIS", fqName.replace_extension("").string());
						out << currentZcatTempCmd << std::endl;
						auto runOut = njh::sys::run({currentZcatTempCmd});
						if(!runOut){
							err << "Ran into error\n";
							err << runOut.stdErr_ << "\n";
							success = false;
						}
						zCatLogs[f] = runOut.toJson();
					}
					if(allEmpty){
						currentEmptySamps.emplace_back(key);
						success = false;
					}
					auto r1_files = njh::files::listAllFiles(sampDir, false, {std::regex{R"(.*_R1_.*.fastq)"}});

					if(r1_files.size() > 1 && success && !allEmpty){
						Json::Value laneCatLogs;
						auto sampRawDir = njh::files::makeDir(sampDir, njh::files::MkdirPar("raw"));
						std::string r1Cmd = "cat";
						for(const auto & r1f : r1_files){
							r1Cmd.append(" " + r1f.first.string());
						}
						r1Cmd.append(" >" + outputDirectory + key + "/" + key + "_R1.fastq");
						auto r2_files = njh::files::listAllFiles(sampDir, false, {std::regex{R"(.*_R2_.*.fastq)"}});
						std::string r2Cmd = "cat";
						for(const auto & r2f : r2_files){
							r2Cmd.append(" " + r2f.first.string());
						}
						r2Cmd.append(" >" + outputDirectory + key + "/" + key + "_R2.fastq");
						out << r1Cmd << std::endl;
						auto r1_runOut = njh::sys::run({r1Cmd});
						out << r2Cmd << std::endl;
						auto r2_runOut = njh::sys::run({r2Cmd});
						if(!r1_runOut || !r2_runOut){
							err << "Error in concatenating files for " << key << std::endl;
							success = false;
						}
						laneCatLogs["R1"] = r1_runOut.toJson();
						laneCatLogs["R2"] = r2_runOut.toJson();

						for(const auto & r1f : r1_files){
							bfs::rename(r1f.first, sampRawDir.string() + r1f.first.filename().string());
						}
						for(const auto & r2f : r2_files){
							bfs::rename(r2f.first, sampRawDir.string() + r2f.first.filename().string());
						}
						log["laneCatLogs"] = laneCatLogs;
					} else if (r1_files.size() == 1 && success && !allEmpty) {
						auto r2_files = njh::files::listAllFiles(sampDir, false, {std::regex{R"(.*_R2_.*.fastq)"}});
						bfs::rename(r1_files.begin()->first,  outputDirectory + key + "/" + key + "_R1.fastq");
						bfs::rename(r2_files.begin()->first,  outputDirectory + key + "/" + key + "_R2.fastq");
					}
					if(success){
						currentSamplesExtracted.emplace_back(key);
					}
					log["zCatLogs"] = zCatLogs;
					log["sample"] = key;
					log["time"] = watch.totalTime();
					log["out"] = out.str();
					log["err"] = err.str();
					log["success"] = njh::json::toJson(success);
					currentLogs[key] = log;
				}
				{
					std::lock_guard<std::mutex> lock(logsMut);
					addOtherVec(samplesEmpty, currentEmptySamps);
					addOtherVec(samplesExtracted, currentSamplesExtracted);
					for(const auto & keyLog : currentLogs){
						logs[keyLog.first] = keyLog.second;
					}
				}
			};
	std::vector<std::thread> threads;
	for(uint32_t tNum = 0; tNum < pars.numThreads; ++tNum){
		threads.emplace_back(std::thread(extractFiles, std::ref(filesKeys), std::cref(zcatTempCmd),
				std::cref(mipMaster.directoryMaster_.masterDir_.string()), std::cref(readPairs), std::ref(extractLog), std::ref(logsLock)));
	}
	for(auto & th : threads){
		th.join();
	}
	logFile << logs << std::endl;
	std::unordered_map<std::string,njh::sys::RunOutput> stitchOutputs;
	//pear
	std::string pearCmdTemp;
	if (pars.usePear) {
		if(pars.trim){
			pearCmdTemp = "pear -f " + mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/trimmed_REPLACE_*R1*.fastq -r "
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/trimmed_REPLACE_*R2*.fastq -o "
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/REPLACE --min-overlap "
					+ njh::lexical_cast<std::string>(pars.minOverlap) + " -j ";
		}else{
			pearCmdTemp = "pear -f " + mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/REPLACE_*R1*.fastq -r "
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/REPLACE_*R2*.fastq -o "
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/REPLACE --min-overlap "
					+ njh::lexical_cast<std::string>(pars.minOverlap) + " -j ";
		}

		if (pars.numThreads > 10) {
			pearCmdTemp += estd::to_string(std::ceil(pars.numThreads / 10.0)) + " >"
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/pearOut.txt";
		} else {
			pearCmdTemp += estd::to_string(pars.numThreads) + " >"
					+ mipMaster.directoryMaster_.masterDir_.string()
					+ "REPLACE/pearOut.txt";
		}
	}
	//flash
	std::string flashCmdTemplate;
	if(pars.trim){
		flashCmdTemplate = "flash "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R1.fastq "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R2.fastq -o "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE"
				" --min-overlap " + estd::to_string(pars.minOverlap) +
				" --max-overlap " + estd::to_string(pars.maxOverlap) +
				" -x " + estd::to_string(pars.mismatchDensity) +
				" -t " + estd::to_string(pars.numThreads) +
				" 2>&1 | tee "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/flash.log";
	}else{
		flashCmdTemplate = "flash "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R1*.fastq "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R2*.fastq -o "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE"
				" --min-overlap " + estd::to_string(pars.minOverlap) +
				" --max-overlap " + estd::to_string(pars.maxOverlap) +
				" -x " + estd::to_string(pars.mismatchDensity) +
				" -t " + estd::to_string(pars.numThreads) +
				" 2>&1 | tee "
				+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/flash.log";
	}

	//pandaseq
	//
	//-g pandaLog.txt
	std::string pandaseqCmdTemplate;
	if(pars.trim){
		pandaseqCmdTemplate = "pandaseq "
				"-f " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R1.fastq "
				"-r " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R2.fastq "
				"-w " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE.pandaseq.fastq "
				" -F "
				" -o " + estd::to_string(pars.minOverlap) +
				" -O " + estd::to_string(pars.maxOverlap) +
				" -T " + estd::to_string(pars.numThreads) +
				" -g " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/pandaseqLog.txt";
	}else{
		pandaseqCmdTemplate = "pandaseq "
				"-f " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R1*.fastq "
				"-r " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R2*.fastq "
				"-w " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE.pandaseq.fastq "
				" -F "
				" -o " + estd::to_string(pars.minOverlap) +
				" -O " + estd::to_string(pars.maxOverlap) +
				" -T " + estd::to_string(pars.numThreads) +
				" -g " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/pandaseqLog.txt";
	}

	std::string sickleCmdTemplate = "sickle pe "
			"-f " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R1*.fastq "
			"-r " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/REPLACE_*R2*.fastq -t sanger "
			"-o " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R1.fastq "
			"-p " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/trimmed_REPLACE_R2.fastq "
			"-s " + mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/singles_REPLACE.fastq"
			" 2>&1 | tee "
			+ mipMaster.directoryMaster_.masterDir_.string() + "REPLACE/sickle.log";



	VecStr samplesStitched;
	if(pars.usePear){
		if(pars.numThreads > 10){
			njh::concurrent::LockableQueue<std::string> keyQueue(keys);
			std::mutex pearOutMut;
			auto runPear = [&keyQueue,&pearOutMut,&stitchOutputs,&pearCmdTemp,&samplesEmpty](){
				std::string k = "";
				std::unordered_map<std::string,njh::sys::RunOutput> currentPearOutputs;
				while(keyQueue.getVal(k)){
					if(njh::in(k, samplesEmpty)){
						continue;
					}
					auto currentPearCmdTemp = njh::replaceString(pearCmdTemp, "REPLACE", k);
					auto runOut = njh::sys::run( { currentPearCmdTemp });
					currentPearOutputs.emplace(k,runOut);
				}
				{
					std::lock_guard<std::mutex> pearOutLock(pearOutMut);
					for(const auto & pearOut : currentPearOutputs){
						stitchOutputs[pearOut.first] = pearOut.second;
					}
				}
			};
			std::vector<std::thread> threads;
			for(uint32_t t = 0; t < 10; ++t){
				threads.emplace_back(std::thread(runPear));
			}
			for(auto & t : threads){
				t.join();
			}
		}else{
			for (const auto & k : keys) {
				if(njh::in(k, samplesEmpty)){
					continue;
				}
				auto currentPearCmdTemp = njh::replaceString(pearCmdTemp, "REPLACE", k);
				auto runOut = njh::sys::run( { currentPearCmdTemp });
				stitchOutputs.emplace(k,runOut);
			}
		}
	}else if(pars.usePanda){
		if(pars.numThreads > 10){
			njh::concurrent::LockableQueue<std::string> keyQueue(keys);
			std::mutex stitchOutputsMut;
			auto runPear = [&pars,&keyQueue,&stitchOutputsMut,&stitchOutputs,&pandaseqCmdTemplate,&sickleCmdTemplate,&samplesEmpty](){
				std::string k = "";
				std::unordered_map<std::string,njh::sys::RunOutput> currentOutputs;
				while(keyQueue.getVal(k)){
					if(njh::in(k, samplesEmpty)){
						continue;
					}
					if(njh::in(k, samplesEmpty)){
						continue;
					}
					auto currentPandaCmdTemp = njh::replaceString(pandaseqCmdTemplate, "REPLACE", k);
					if (pars.trim) {
						auto currentSickleCmdTemp = njh::replaceString(sickleCmdTemplate, "REPLACE",
								k);
						auto runOut = njh::sys::run( { currentSickleCmdTemp });
						currentOutputs.emplace(k + "_" + "sickle", runOut);
					}
					auto runOut = njh::sys::run( { currentPandaCmdTemp });
					currentOutputs.emplace(k,runOut);
				}
				{
					std::lock_guard<std::mutex> pearOutLock(stitchOutputsMut);
					for(const auto & pearOut : currentOutputs){
						stitchOutputs[pearOut.first] = pearOut.second;
					}
				}
			};
			std::vector<std::thread> threads;
			for(uint32_t t = 0; t < 10; ++t){
				threads.emplace_back(std::thread(runPear));
			}
			for(auto & t : threads){
				t.join();
			}
		}else{
			for (const auto & k : keys) {
				if(njh::in(k, samplesEmpty)){
					continue;
				}
				auto currentPandaCmdTemp = njh::replaceString(pandaseqCmdTemplate, "REPLACE", k);
				if (pars.trim) {
					auto currentSickleCmdTemp = njh::replaceString(sickleCmdTemplate, "REPLACE",
							k);
					auto runOut = njh::sys::run( { currentSickleCmdTemp });
					stitchOutputs.emplace(k + "_" + "sickle", runOut);
				}
				auto runOut = njh::sys::run( { currentPandaCmdTemp });
				stitchOutputs.emplace(k,runOut);
			}
		}
	}else{
		if(pars.numThreads > 10){
			njh::concurrent::LockableQueue<std::string> keyQueue(keys);
			std::mutex stitchOutputsMut;
			auto runPear = [&pars,&keyQueue,&stitchOutputsMut,&stitchOutputs,&flashCmdTemplate,&sickleCmdTemplate,&samplesEmpty](){
				std::string k = "";
				std::unordered_map<std::string,njh::sys::RunOutput> currentOutputs;
				while(keyQueue.getVal(k)){
					if(njh::in(k, samplesEmpty)){
						continue;
					}
					if(njh::in(k, samplesEmpty)){
						continue;
					}
					auto currentFlashCmdTemp = njh::replaceString(flashCmdTemplate, "REPLACE", k);
					if (pars.trim) {
						auto currentSickleCmdTemp = njh::replaceString(sickleCmdTemplate, "REPLACE",
								k);
						auto runOut = njh::sys::run( { currentSickleCmdTemp });
						currentOutputs.emplace(k + "_" + "sickle", runOut);
					}
					auto runOut = njh::sys::run( { currentFlashCmdTemp });
					currentOutputs.emplace(k,runOut);
				}
				{
					std::lock_guard<std::mutex> pearOutLock(stitchOutputsMut);
					for(const auto & pearOut : currentOutputs){
						stitchOutputs[pearOut.first] = pearOut.second;
					}
				}
			};
			std::vector<std::thread> threads;
			for(uint32_t t = 0; t < 10; ++t){
				threads.emplace_back(std::thread(runPear));
			}
			for(auto & t : threads){
				t.join();
			}
		}else{
			for (const auto & k : keys) {
				if(njh::in(k, samplesEmpty)){
					continue;
				}
				auto currentFlashCmdTemp = njh::replaceString(flashCmdTemplate, "REPLACE", k);
				if (pars.trim) {
					auto currentSickleCmdTemp = njh::replaceString(sickleCmdTemplate, "REPLACE",
							k);
					auto runOut = njh::sys::run( { currentSickleCmdTemp });
					stitchOutputs.emplace(k + "_" + "sickle", runOut);
				}
				auto runOut = njh::sys::run( { currentFlashCmdTemp });
				stitchOutputs.emplace(k,runOut);
			}
		}
	}
	for(const auto & output : stitchOutputs){
		logFile << output.second.toJson() << std::endl;
		if(!output.second.success_){
			std::cerr << "Error in running pear with " << output.second.cmd_ << std::endl;;
			std::cerr << njh::bashCT::boldRed(output.second.stdErr_) << std::endl;
		}else{
			if(!njh::containsSubString(output.first, "sickle")){
				samplesStitched.emplace_back(output.first);
			}
		}
	}
	//copy over resources;

	bfs::copy(pars.mipArmsFileName,njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "mip_arm_id.tab.txt"));
	njh::files::makeDirP(njh::files::MkdirPar(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo")));
	std::ofstream outSamplesFoundFile;
	openTextFile(outSamplesFoundFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesFound.tab.txt")));
	std::ofstream outSamplesEmptyfile;
	openTextFile(outSamplesEmptyfile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/outSamplesEmpty.tab.txt")));
	std::ofstream gzCatSamplesFile;
	openTextFile(gzCatSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/gzCatSamples.tab.txt")));
	std::ofstream pearExtractSamplesFile;
	openTextFile(pearExtractSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/samplesStitched.tab.txt")));
	std::ofstream samplesMissingFile;
	openTextFile(samplesMissingFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "sampleExtractInfo/samplesMissing.tab.txt")));
	auto allFoundSamps = concatVecs(samplesExtracted, samplesEmpty);
	njh::sort(allFoundSamps);
	njh::sort(samplesEmpty);
	njh::sort(samplesExtracted);
	njh::sort(samplesStitched);
	if(!allFoundSamps.empty()) outSamplesFoundFile << njh::conToStr(allFoundSamps, "\n") << std::endl;
	if(!samplesEmpty.empty()) outSamplesEmptyfile << njh::conToStr(samplesEmpty, "\n") << std::endl;
	if(!samplesExtracted.empty()) gzCatSamplesFile << njh::conToStr(samplesExtracted, "\n") << std::endl;
	if(!samplesStitched.empty()) pearExtractSamplesFile << njh::conToStr(samplesStitched, "\n") << std::endl;
	VecStr samplesMissing;
	for(const auto & samp : mipMaster.names_->samples_){
		if(!njh::in(samp, samplesStitched)){
			samplesMissing.push_back(samp);
		}
	}
	if(!samplesMissing.empty()) samplesMissingFile << njh::conToStr(samplesMissing, "\n") << std::endl;
	std::ofstream allMipsSamplesFile;
	openTextFile(allMipsSamplesFile,OutOptions(njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "allMipsSamplesNames.tab.txt")));
	MipsSamplesNames goodSamples = *mipMaster.names_;
	goodSamples.setSamples(samplesStitched);
	goodSamples.write(allMipsSamplesFile);
	bfs::copy(pars.mipsSamplesFile,njh::files::join(mipMaster.directoryMaster_.resourceDir_.string(), "original_allMipsSamplesNames.tab.txt"));
	mipMaster.createPopClusMipDirs(pars.numThreads);

	if (nullptr != mipMaster.meta_) {
		bfs::copy_file(mipMaster.meta_->groupingsFile_,
				njh::files::make_path(mipMaster.directoryMaster_.masterDir_,
						"resources", "samplesMeta.tab.txt"));
	}

	return 0;
}


} // namespace njhseq

