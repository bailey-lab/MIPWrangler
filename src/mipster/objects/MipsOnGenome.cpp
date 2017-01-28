/*
 * MipsOnGenome.cpp
 *
 *  Created on: Jan 28, 2017
 *      Author: nick
 */




#include "MipsOnGenome.hpp"
#include "mipster/objects/MipMapResult.hpp"
#include <TwoBit.h>

namespace bibseq {

MipsOnGenome::MipsOnGenome(const bfs::path & mainDir, uint32_t numThreads) :
		mainDir_(mainDir), numThreads_(numThreads) {
	genomeDir_ = bib::files::make_path(mainDir, "genomes");
	infoDir_ = bib::files::make_path(mainDir, "info");
	mipArmsFnp_ = bib::files::make_path(mainDir, "info", "mip_arms.tab.txt");
	mapDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("mapped"));
	bedsDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("beds"));
	armsDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("arms"));
	logDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("logs"));
	checkInputThrow();
	requireExternalProgramThrow("bowtie2");

}

void MipsOnGenome::checkInputThrow()const{
	std::stringstream ss;
	bool failed = false;
	auto checkForPath = [&failed,&ss](const bfs::path & fnp ){
		if(!bfs::exists(fnp)){
			failed = true;
			ss << bib::bashCT::boldRed(fnp.string())<< " needs to exist "<< "\n";
		}
	};
	checkForPath(mainDir_);
	checkForPath(genomeDir_);
	checkForPath(infoDir_);
	checkForPath(mipArmsFnp_);
	if(failed){
		std::stringstream outSS;
		outSS << __PRETTY_FUNCTION__ << ", error in checking directory set up" << "\n";
		outSS << ss.str();
		throw std::runtime_error{outSS.str()};
	}
}


MipsOnGenome::Genome::Genome(const bfs::path & fnp) :
		fnp_(fnp) {
	checkExistenceThrow(fnp_, __PRETTY_FUNCTION__);
}


void MipsOnGenome::Genome::createTwoBit() {
	auto twoBitFnp = fnp_;
	twoBitFnp.replace_extension(".2bit");
	TwoBit::faToTwoBitPars pars;
	pars.inputFilename = fnp_.string();
	pars.outFilename = twoBitFnp.string();
	bool buildTwoBit = false;
	if(!bfs::exists(twoBitFnp)){
		buildTwoBit = true;
	}else if(bib::files::firstFileIsOlder(twoBitFnp, fnp_) ){
		buildTwoBit = true;
		pars.overWrite = true;
	}
	if(buildTwoBit){
		TwoBit::fastasToTwoBit(pars);
	}
	fnpTwoBit_ = pars.outFilename;
}

void MipsOnGenome::Genome::buildBowtie2Index() const {
	std::string indexCmd = "bowtie2-build {GENOMEPREFIX}.fasta {GENOMEPREFIX}";
	auto genomePrefix = bib::files::removeExtension(fnp_.string());
	indexCmd = bib::replaceString(indexCmd, "{GENOMEPREFIX}", genomePrefix);
	bool buildIndex = false;
	auto testPath = genomePrefix + ".1.bt2";
	if (!bfs::exists(testPath)
			|| bib::files::firstFileIsOlder(testPath, fnp_)) {
		buildIndex = true;
	}
	if (buildIndex) {
		auto runOutput = bib::sys::run( { indexCmd });
		if (!runOutput.success_) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << " failed to build index for " << fnp_
					<< "\n";
			ss << runOutput.toJson() << "\n";
			throw std::runtime_error { ss.str() };
		}
	}
}

void MipsOnGenome::loadInGenomes(){
	auto fastaFiles = bib::files::gatherFiles(genomeDir_, ".fasta", false);
	for(const auto & f : fastaFiles){
		genomes_[bib::files::removeExtension(f.filename().string())] = std::make_unique<Genome>(f);
	}
}

void MipsOnGenome::setUpGenomes(){
	for(auto & gen : genomes_){
		gen.second->createTwoBit();
		gen.second->buildBowtie2Index();
	}
}

void MipsOnGenome::loadInArms(){
	mipArms_ = std::make_unique<MipCollection>(mipArmsFnp_, 6);
}

void MipsOnGenome::createArmFiles(){
	OutOptions opts(armsDir_);
	opts.overWriteFile_ = true;
	for(const auto & m : mipArms_->mips_){
		m.second.writeOutArms(opts);
	}
}



std::string MipsOnGenome::GenomeMip::uid(const std::string & sep) const {
	return genome_ + sep + mip_;
}


void MipsOnGenome::mapArmsToGenomes() {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip{gen.first, m.second.name_});
		}
	}
	bib::concurrent::LockableQueue<GenomeMip> pairsQueue(pairs);
	std::mutex logMut;
	Json::Value log;

	auto mapArm = [this, &pairsQueue,&logMut,&log](){
		GenomeMip pair;
		while(pairsQueue.getVal(pair)){
			std::stringstream ss;
			bool succes = false;
			std::string outStub = bib::files::make_path(mapDir_,
					pair.uid()).string();
			std::string outCheck = outStub + ".sorted.bam";
			if (!bfs::exists(outCheck)
					|| bib::files::firstFileIsOlder(outCheck, genomes_.at(pair.genome_)->fnp_)
					|| bib::files::firstFileIsOlder(outCheck, mipArms_->mipArmIdFnp_)) {
				std::stringstream bowtie2Cmd;
				bowtie2Cmd
						<< " bowtie2 -p 1 -D 20 -R 3 -N 1 -L 18 -i S,1,0.5 --gbar 1 -k 1 --end-to-end "
						<< "-x " << bib::files::make_path(genomeDir_, pair.genome_)
						<< " -f -1 " << bib::files::make_path(armsDir_, pair.mip_)
						<< "_ext-arm.fasta" << " -2 "
						<< bib::files::make_path(armsDir_, pair.mip_)
						<< "_lig-arm.fasta" << " -S " << outStub + ".sam";
				std::stringstream samtoolsCmds;
				samtoolsCmds << "samtools view -Sb " << outStub << ".sam | "
						<< "samtools sort - -o " << outStub << ".sorted.bam && "
						<< "samtools index " << outStub << ".sorted.bam";
				auto bowtie2RunOutput = bib::sys::run( { bowtie2Cmd.str() });
				if (bowtie2RunOutput.success_) {
					auto samtoolsRunOutput = bib::sys::run( { samtoolsCmds.str() });
					succes = true;
					if (!samtoolsRunOutput.success_) {
						ss << "Failed to sort " << outStub << ".sam" << std::endl;
						ss << samtoolsRunOutput.stdErr_ << std::endl;
					}
				} else {
					ss << "Failed to map " << pair.mip_ << " to "
							<< pair.genome_ << std::endl;
					ss << bowtie2RunOutput.stdErr_ << std::endl;
				}
			}
			{
				std::lock_guard<std::mutex> lock(logMut);
				log[pair.uid()]["succes"] = succes;
				log[pair.uid()]["message"] = ss.str();
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads_; ++t){
		threads.emplace_back(std::thread(mapArm));
	}
	for(auto & t : threads){
		t.join();
	}

	std::ofstream outFile;
	OutOptions logOpts(bib::files::make_path(logDir_, "mapLog-" + bib::getCurrentDate() + ".json"));
	logOpts.overWriteFile_ = true;
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;
}

void MipsOnGenome::genBeds() {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip{gen.first, m.second.name_});
		}
	}
	bib::concurrent::LockableQueue<GenomeMip> pairsQueue(pairs);
	std::mutex logMut;
	Json::Value log;
	auto genBeds = [this, &pairsQueue,&logMut,&log](){
		GenomeMip pair;
		while(pairsQueue.getVal(pair)){
			std::string outStub = bib::files::make_path(mapDir_,
					pair.genome_ + "_" + pair.mip_).string();
			std::string outCheck = outStub + ".sorted.bam";
			std::stringstream ss;
			bool succes = false;
			if (!bfs::exists(outCheck)) {
				std::lock_guard<std::mutex> lock(logMut);
				ss << "Failed to find " << outCheck << " for " << pair.mip_ << " to "
						<< pair.genome_ << std::endl;
			}else{
				auto results = getMipMapResults(outCheck);
				if(results.empty()){
					ss << "Failed to get results from " << outCheck << " for " << pair.mip_ << " to "
							<< pair.genome_ << std::endl;
				}else if(results.size() > 1){
					ss << "Was expecting only to get 1 result from " << outCheck << " for " << pair.mip_ << " to "
													<< pair.genome_ << " but got " << results.size( ) << " instead" << std::endl;
				}else{
					std::ofstream outFile;
					OutOptions bedOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + ".bed"));

					bedOpts.overWriteFile_ = true;
					openTextFile(outFile, bedOpts);

					if (results.front().isConcordant() && results.front().isMapped()) {
						outFile << results.front().region_.genBedRecordCore().toDelimStr() << std::endl;
						succes =  true;
					} else {
						if (!results.front().isMapped()) {
							if (!results.front().extAln_.IsMapped()) {
								ss << "Result for " << results.front().mipName_ << " mapped to "
										<< results.front().genomeName_ << " extension arm didn't map "
										<< std::endl;
							}
							if (!results.front().ligAln_.IsMapped()) {
								ss << "Result for " << results.front().mipName_ << " mapped to "
										<< results.front().genomeName_ << " ligation arm didn't map" << std::endl;
							}
						} else if (!results.front().isConcordant()) {
							ss << "Result for " << results.front().mipName_ << " mapped to "
									<< results.front().genomeName_ << " had discordant arms" << std::endl;
						}
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(logMut);
				log[pair.uid()]["succes"] = succes;
				log[pair.uid()]["message"] = ss.str();
			}

		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads_; ++t){
		threads.emplace_back(std::thread(genBeds));
	}
	for(auto & t : threads){
		t.join();
	}

	std::ofstream outFile;
	OutOptions logOpts(bib::files::make_path(logDir_, "bedLog-" + bib::getCurrentDate() + ".json"));
	logOpts.overWriteFile_ = true;
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;

}

}  // namespace bibseq

