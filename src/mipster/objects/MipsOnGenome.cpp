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
	fastaDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("fastas"));
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
	if(fastaFiles.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, there has to be at least one genome in " << genomeDir_ << "\n";
		ss << "Found none ending with .fasta" << "\n";
		throw std::runtime_error{ss.str()};
	}
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
	log["date"] = bib::getCurrentDateFull();
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
			}else{
				succes = true;
				ss << outCheck << " is up to date";
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
	logOpts.outFilename_ = bib::files::findNonexitantFile(logOpts.outFilename_);
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;
}
void MipsOnGenome::genFastas() {
	const VecStr mips = bib::getVecOfMapKeys( mipArms_->mips_);
	const VecStr genomes = bib::getVecOfMapKeys( genomes_	);
	bib::concurrent::LockableQueue<std::string> mipQueue(mips);
	std::mutex logMut;
	Json::Value log;
	log["date"] = bib::getCurrentDateFull();
	auto genFastasFunc = [this, &mipQueue,&genomes,&logMut,&log](){
		std::string mipName;
		while(mipQueue.getVal(mipName)){
			std::stringstream ss;
			bool succes = false;
			auto outOpts = SeqIOOptions::genFastaOut(bib::files::make_path(fastaDir_, mipName));
			outOpts.out_.overWriteFile_ = true;
			std::unordered_map<std::string, std::shared_ptr<InOptions>> bedOpts;
			bool needsUpdate = false;
			for(const auto & genome : genomes){
				std::shared_ptr<InOptions> bedOpt = std::make_shared<InOptions>(bib::files::make_path(bedsDir_, genome + "_" + mipName + ".bed"));
				//there is a possibility that the bed creatin failed due to bad mapping but this also doesn't take into account if the beds haven't been created yet
				if(bedOpt->inExists()){
					if(outOpts.out_.outExists() && bib::files::firstFileIsOlder(outOpts.out_.outName(), bedOpt->inFilename_)){
						needsUpdate = true;
					}
					bedOpts[genome] = bedOpt;
				}
			}
			if(bedOpts.empty()){
				succes = false;
				ss << "No bed files found for " << mipName << "\n";
			}else{
				if(outOpts.outExists() && !needsUpdate){
					succes = true;
					ss << outOpts.out_.outName() << " already up to date";
				}else{
					std::vector<seqInfo> seqs;
					for(const auto & bedOpt : bedOpts){
						std::string genome = bedOpt.first;
						auto regions = gatherRegions(bedOpt.second->inFilename_.string(), "", false);
						if(regions.empty()){
							succes = false;
							ss << "Error in parsing " << bedOpt.second->inFilename_ << "\n";
						}else{
							TwoBit::TwoBitFile twoBitFile(bib::files::make_path(genomeDir_, genome + ".2bit"));
							std::string seq = "";
							twoBitFile[regions.front().chrom_]->getSequence(seq, regions.front().start_, regions.front().end_);
							//consider leaving lower case
							bib::strToUpper(seq);
							if(regions.front().reverseSrand_){
								seq = seqUtil::reverseComplement(seq, "DNA");
							}
							seqs.emplace_back(seqInfo(genome, seq));
						}
					}
					if(seqs.empty()){
						succes = false;
						ss << "Failed to extract any sequences from bed files " << "\n";
					}else{
						std::vector<seqInfo> outputSeqs;
						for(const auto & seq : seqs){
							if(outputSeqs.empty()){
								outputSeqs.emplace_back(seq);
							}else{
								bool foundSame = false;
								for( auto & outSeq : outputSeqs){
									if(seq.seq_ == outSeq.seq_ ){
										outSeq.name_ += "-" + seq.name_;
										foundSame = true;
										break;
									}
								}
								if(!foundSame){
									outputSeqs.emplace_back(seq);
								}
							}
						}
						for(auto & outSeq : outputSeqs){
							if(std::string::npos != outSeq.name_.find('-')){
								auto gs = bib::tokenizeString(outSeq.name_, "-");
								bib::sort(gs);
								outSeq.name_ = bib::conToStr(gs, "-");
							}
						}
						succes = true;
						SeqOutput::write(outputSeqs, outOpts);
					}
				}
			}
			{
				std::lock_guard<std::mutex> lock(logMut);
				log[mipName]["succes"] = succes;
				log[mipName]["message"] = ss.str();
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < numThreads_; ++t){
		threads.emplace_back(std::thread(genFastasFunc));
	}
	for(auto & t : threads){
		t.join();
	}
	std::ofstream outFile;
	OutOptions logOpts(bib::files::make_path(logDir_, "fastaLog-" + bib::getCurrentDate() + ".json"));
	logOpts.outFilename_ = bib::files::findNonexitantFile(logOpts.outFilename_);
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;
}


std::vector<MipsOnGenome::GenomeMip> MipsOnGenome::genGenomeMipPairs()const{
	std::vector<GenomeMip> pairs;
	for (const auto & gen : genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip{gen.first, m.second.name_});
		}
	}
	return pairs;
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
	log["date"] = bib::getCurrentDateFull();
	auto genBedsFunc = [this, &pairsQueue,&logMut,&log](){
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
					OutOptions bedOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + ".bed"));
					if(bedOpts.outExists() && bib::files::firstFileIsOlder(outCheck ,bedOpts.outName())){
						succes = true;
						ss << bedOpts.outName() << " already up to date" << "\n";
					}else{
						if (results.front().isConcordant() && results.front().isMapped()) {
							std::ofstream outFile;
							bedOpts.overWriteFile_ = true;
													openTextFile(outFile, bedOpts);
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
		threads.emplace_back(std::thread(genBedsFunc));
	}
	for(auto & t : threads){
		t.join();
	}
	std::ofstream outFile;
	OutOptions logOpts(bib::files::make_path(logDir_, "bedLog-" + bib::getCurrentDate() + ".json"));
	logOpts.outFilename_ = bib::files::findNonexitantFile(logOpts.outFilename_);
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;
}

void MipsOnGenome::setPrimaryGenome(const std::string & genome){
	if(!bib::in(genome, genomes_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no genome found matching" << genome << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(genomes_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	primaryGenome_ = genome;
}

bfs::path MipsOnGenome::pathToMipFasta(const std::string & mipName)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching" << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaDir_, mipName + ".fasta");
}

VecStr MipsOnGenome::getMips() const{
	auto ret = bib::getVecOfMapKeys(mipArms_->mips_);
	bib::sort(ret);
	return ret;
}

VecStr MipsOnGenome::getGenomes() const{
	VecStr ret;
	VecStr genomes = bib::getVecOfMapKeys(genomes_);
	if("" != primaryGenome_){
		ret.emplace_back(primaryGenome_);
		removeElement(genomes, primaryGenome_);

	}
	bib::sort(genomes);
	addOtherVec(ret, genomes);
	return ret;
}

}  // namespace bibseq
