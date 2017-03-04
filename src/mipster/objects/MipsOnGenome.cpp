/*
 * MipsOnGenome.cpp
 *
 *  Created on: Jan 28, 2017
 *      Author: nick
 */

#include "MipsOnGenome.hpp"
#include "mipster/objects/MipMapResult.hpp"
#include "mipster/mipUtils.h"
#include <elucidator/objects/BioDataObject.h>
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
	tablesDir_ = bib::files::makeDirP(mainDir, bib::files::MkdirPar("tables"));
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
	fnpTwoBit_ = fnp_;
	fnpTwoBit_.replace_extension(".2bit");
}
void MipsOnGenome::Genome::createTwoBit() {

	TwoBit::faToTwoBitPars pars;
	pars.inputFilename = fnp_.string();
	pars.outFilename = fnpTwoBit_.string();
	bool buildTwoBit = false;
	if(!bfs::exists(fnpTwoBit_)){
		buildTwoBit = true;
	}else if(bib::files::firstFileIsOlder(fnpTwoBit_, fnp_) ){
		buildTwoBit = true;
		pars.overWrite = true;
	}
	if(buildTwoBit){
		TwoBit::fastasToTwoBit(pars);
	}
	fnpTwoBit_ = pars.outFilename;
}


Json::Value MipsOnGenome::Genome::chromosomeLengths() const{
	Json::Value ret;

	TwoBit::TwoBitFile genFile(fnpTwoBit_);
	auto lens = genFile.getSeqLens();
	auto lenKeys = bib::getVecOfMapKeys(lens);
	bib::sort(lenKeys);
	for(const auto & lenKey : lenKeys	){
		Json::Value lenObj;
		lenObj["name"] = lenKey;
		lenObj["len"] = lens[lenKey];
		ret.append(lenObj);
	}
	return ret;
}

std::string MipsOnGenome::getPrimaryGenome(){
	return primaryGenome_;
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
			auto trimmedOutOpts = SeqIOOptions::genFastaOut(bib::files::make_path(fastaDir_, "trimmed_" + mipName));
			trimmedOutOpts.out_.overWriteFile_ = true;
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
					std::vector<seqInfo> trimmedSeqs;
					for(const auto & bedOpt : bedOpts){
						std::string genome = bedOpt.first;
						auto regions =    gatherRegions(bedOpt.second->inFilename_.string(), "", false);
						auto extRegions = gatherRegions(bib::replaceString(bedOpt.second->inFilename_.string(), ".bed", "-ext.bed"), "", false);
						auto ligRegions = gatherRegions(bib::replaceString(bedOpt.second->inFilename_.string(), ".bed", "-lig.bed"), "", false);
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
							seqInfo trimmedSeq(genome, seq);
							readVecTrimmer::trimOffForwardBases(trimmedSeq, extRegions.front().getLen());
							readVecTrimmer::trimOffEndBases(trimmedSeq, ligRegions.front().getLen());
							seqs.emplace_back(seqInfo(genome, seq));
							trimmedSeqs.emplace_back(trimmedSeq);
						}
					}
					if(seqs.empty()){
						succes = false;
						ss << "Failed to extract any sequences from bed files " << "\n";
					}else{
						auto collapseSimSeqs = [](std::vector<seqInfo> & seqs){
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
							return outputSeqs;
						};

						succes = true;
						auto outputSeqs = collapseSimSeqs(seqs);
						SeqOutput::write(outputSeqs, outOpts);
						auto trimedOutputSeqs = collapseSimSeqs(trimmedSeqs);
						SeqOutput::write(trimedOutputSeqs, trimmedOutOpts);
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


table MipsOnGenome::getGenomeLocsForMipTar(const std::string & tar) const{
	if(!bib::in(tar, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for regionName "
				<< tar << "\n";
		ss << "Options are " << bib::conToStr(mipArms_->getMipTars(), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	SeqInput seqReader(SeqIOOptions::genFastaIn(pathToMipFasta(tar)));
	auto seqs = seqReader.readAllReadsPtrs<readObject>();
	std::unordered_map<std::string, std::shared_ptr<readObject>> seqsByName;
	for (const auto & seq : seqs) {
		auto toks = tokenizeString(seq->seqBase_.name_, "-");
		for (const auto & tok : toks) {
			seqsByName[tok] = seq;
		}
	}
	table locs(VecStr{"genome","chrom","start", "stop", "strand", "length",  "GCContent", "longestHomopolymer" });
	for(const auto & genome : genomes_){
		auto bedFnp = bib::files::make_path(pathToMipBed(tar, genome.first));
		if(bfs::exists(bedFnp)){
			BedRecordCore bedCore;
			BioDataFileIO<BedRecordCore> bedReader((IoOptions(InOptions(bedFnp))));
			bedReader.openIn();
			uint32_t longestHomopolymer = 0;
			double gcContent = 0;
			if (bib::in(genome.first, seqsByName)) {
				seqsByName[genome.first]->setLetterCount();
				seqsByName[genome.first]->counter_.resetAlphabet(true);
				seqsByName[genome.first]->counter_.setFractions();
				seqsByName[genome.first]->counter_.calcGcContent();
				gcContent = seqsByName[genome.first]->counter_.gcContent_;
				seqsByName[genome.first]->createCondensedSeq();
				longestHomopolymer = vectorMaximum(
						seqsByName[genome.first]->condensedSeqCount);
			}
			while(bedReader.readNextRecord(bedCore)){
				locs.addRow(genome.first, bedCore.chrom_, bedCore.chromStart_, bedCore.chromEnd_, bedCore.strand_,bedCore.length(),
						gcContent, longestHomopolymer);
			}
		}
	}
	return locs;
}

table MipsOnGenome::getGenomeLocsForAllMipTars() const {
	table locs(VecStr { "genome", "target", "chrom", "start", "stop", "strand",
			"length", "GCContent", "longestHomopolymer" });
	for (const auto & tar : getMips()) {
		SeqInput seqReader(SeqIOOptions::genFastaIn(pathToMipFasta(tar)));
		auto seqs = seqReader.readAllReadsPtrs<readObject>();
		std::unordered_map<std::string, std::shared_ptr<readObject>> seqsByName;
		for (const auto & seq : seqs) {
			auto toks = tokenizeString(seq->seqBase_.name_, "-");
			for (const auto & tok : toks) {
				seqsByName[tok] = seq;
			}
		}
		for(const auto & genome : getGenomes()){
			auto bedFnp = bib::files::make_path(bedsDir_, genome + "_" + tar + ".bed");
			if (bfs::exists(bedFnp)) {
				BedRecordCore bedCore;
				BioDataFileIO<BedRecordCore> bedReader((IoOptions(InOptions(bedFnp))));
				bedReader.openIn();
				uint32_t longestHomopolymer = 0;
				double gcContent = 0;
				if (bib::in(genome, seqsByName)) {
					seqsByName[genome]->setLetterCount();
					seqsByName[genome]->counter_.resetAlphabet(true);
					seqsByName[genome]->counter_.setFractions();
					seqsByName[genome]->counter_.calcGcContent();
					gcContent = seqsByName[genome]->counter_.gcContent_;
					seqsByName[genome]->createCondensedSeq();
					longestHomopolymer = vectorMaximum(
							seqsByName[genome]->condensedSeqCount);
				}
				while (bedReader.readNextRecord(bedCore)) {
					locs.addRow(genome, tar, bedCore.chrom_, bedCore.chromStart_,
							bedCore.chromEnd_, bedCore.strand_,
							bedCore.length(),
							gcContent, longestHomopolymer);
				}
			}
		}
	}
	return locs;
}

table MipsOnGenome::getGenomeLocsForGenome(const std::string & genome) const {
	if (!bib::in(genome, genomes_)) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": error, no information for genome " << genome
				<< "\n";
		ss << "Options are " << bib::conToStr(getGenomes(), ", ") << "\n";
		throw std::runtime_error { ss.str() };
	}
	table locs(VecStr { "genome", "region", "target", "chrom", "start", "end",
			"strand", "length", "GCContent", "longestHomopolymer" });
	for (const auto & tar : getMips()) {
		auto bedFnp = bib::files::make_path(bedsDir_, genome + "_" + tar + ".bed");
		if (bfs::exists(bedFnp)) {
			SeqInput seqReader(SeqIOOptions::genFastaIn(pathToMipFasta(tar)));
			auto seqs = seqReader.readAllReadsPtrs<readObject>();
			std::unordered_map<std::string, std::shared_ptr<readObject>> seqsByName;
			for (const auto & seq : seqs) {
				auto toks = tokenizeString(seq->seqBase_.name_, "-");
				for (const auto & tok : toks) {
					seqsByName[tok] = seq;
				}
			}
			BedRecordCore bedCore;
			BioDataFileIO<BedRecordCore> bedReader((IoOptions(InOptions(bedFnp))));
			bedReader.openIn();
			uint32_t longestHomopolymer = 0;
			double gcContent = 0;
			if (bib::in(genome, seqsByName)) {
				seqsByName[genome]->setLetterCount();
				seqsByName[genome]->counter_.resetAlphabet(true);
				seqsByName[genome]->counter_.setFractions();
				seqsByName[genome]->counter_.calcGcContent();
				gcContent = seqsByName[genome]->counter_.gcContent_;
				seqsByName[genome]->createCondensedSeq();
				longestHomopolymer = vectorMaximum(
						seqsByName[genome]->condensedSeqCount);
			}
			while (bedReader.readNextRecord(bedCore)) {
				locs.addRow(genome, mipArms_->mips_.at(tar).locGrouping_, tar,
						bedCore.chrom_, bedCore.chromStart_, bedCore.chromEnd_,
						bedCore.strand_, bedCore.length(), gcContent, longestHomopolymer);
			}
		}
	}
	return locs;
}

table MipsOnGenome::getMipRegionStatsForGenome(
		const std::string & genome) const {
	auto locs = getGenomeLocsForGenome(genome);
	auto locsSplit = locs.splitTableOnColumn("region");
	table ret(VecStr { "region", "numOfTargets", "approxAreaCovered", "chrom",
			"start", "end" });
	auto locKeys = getVectorOfMapKeys(locsSplit);
	MipNameSorter::sortByRegion(locKeys);
	for (const auto & locKey : locKeys) {
		std::string reg = locKey;
		uint32_t numOfTargets = locsSplit[locKey].nRow();
		std::string chrom = locsSplit[locKey].content_.front()[locsSplit[locKey].getColPos(
				"chrom")];
		uint32_t minStart = vectorMinimum(
				vecStrToVecNum<uint32_t>(locsSplit[locKey].getColumn("start")));
		uint32_t maxStop = vectorMaximum(
				vecStrToVecNum<uint32_t>(locsSplit[locKey].getColumn("end")));
		ret.addRow(reg, numOfTargets, maxStop - minStart, chrom, minStart, maxStop);
	}
	return ret;
}

table MipsOnGenome::getMipTarStatsForGenome(const std::string & genome,
		const VecStr & mipTars) const{
	table ret(VecStr { "region", "target", "genome", "chrom", "start", "end", "strand",
			"length", "containsTandemRepeat", "PossibleLengthVariation", "variantNum",
			"totalVariantsPossible", "variantRatio", "GCContent", "longestHomopolymer" });
	for (const auto & mipTar : mipTars) {
		auto bedFnp = bib::files::make_path(pathToMipBed(mipTar, genome));
		if(!bfs::exists(bedFnp)){
			continue;
		}
		//std::cout << mipTar << std::endl;
		auto seqsOpts = SeqIOOptions::genFastaIn(pathToMipFasta(mipTar));
		uint32_t totalHapsPossible = 0;
		uint32_t hapNum = 0;
		bool containsTandems = false;
		bool lengthVariation = false;
		std::vector<uint32_t> readLens;
		readObject refSeq;
		seqInfo seq;
		SeqInput reader(seqsOpts);
		reader.openIn();
		while(reader.readNextRead(seq)){
			auto tandems = aligner::findTandemRepeatsInSequence(seq.seq_, 2, -2, -7, 20);
			for(const auto & tandem : tandems){
				if(tandem.numberOfRepeats_ > 2){
					containsTandems = true;
					break;
				}
			}
			readLens.emplace_back(len(seq));
			auto toks = tokenizeString(seq.name_, "-");
			totalHapsPossible+= toks.size();
			++hapNum;
			if(bib::containsSubString(seq.name_, genome)){
				refSeq = readObject(seq);
			}
		}
		auto minLen = vectorMinimum(readLens);
		auto maxLen = vectorMaximum(readLens);
		if(maxLen - minLen > 4){
			lengthVariation = true;
		}
		BedRecordCore bedCore;
		BioDataFileIO<BedRecordCore> bedReader((IoOptions(InOptions(bedFnp))));
		bedReader.openIn();
		bedReader.readNextRecord(bedCore);
		double gcContent = 0;
		uint32_t longestHomopolymer = 0;
		if("" != refSeq.seqBase_.name_){
			refSeq.setLetterCount();
			refSeq.counter_.resetAlphabet(true);
			refSeq.counter_.setFractions();
			refSeq.counter_.calcGcContent();
			gcContent = refSeq.counter_.gcContent_;
			refSeq.createCondensedSeq();
			longestHomopolymer = vectorMaximum(refSeq.condensedSeqCount);
		}
		ret.addRow(mipArms_->mips_[mipTar].locGrouping_,
							mipTar, genome,
							bedCore.chrom_, bedCore.chromStart_, bedCore.chromEnd_,
							bedCore.strand_, bedCore.length(),
							containsTandems ? "yes": "no",
							lengthVariation ? "yes": "no",
							hapNum, totalHapsPossible,
							static_cast<double>(hapNum)/totalHapsPossible,
							gcContent, longestHomopolymer);
	}
	return ret;
}

std::vector<MipsOnGenome::GenomeMip> MipsOnGenome::genGenomeMipPairs() const {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip { gen.first, m.second.name_ });
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
				uint32_t insertSizeCutoff = 1000;
				//temporary fix
				if(bib::containsSubString(pair.mip_, "full")){
					insertSizeCutoff = 2500;
				}
				auto results = getMipMapResults(outCheck, insertSizeCutoff);
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
							OutOptions bedExtOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + "-ext.bed"));
							OutOptions bedLigOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + "-lig.bed"));
							//full region
							std::ofstream outFile;
							bedOpts.overWriteFile_ = true;
							bedOpts.openFile(outFile);
							outFile << results.front().region_.genBedRecordCore().toDelimStr() << std::endl;
							//ext
							std::ofstream outExtFile;
							bedExtOpts.overWriteFile_ = true;
							bedExtOpts.openFile(outExtFile);
							outExtFile << results.front().extArmRegion_.genBedRecordCore().toDelimStr() << std::endl;
							//lig
							std::ofstream outLigFile;
							bedLigOpts.overWriteFile_ = true;
							bedLigOpts.openFile(outLigFile);
							outLigFile << results.front().ligArmRegion_.genBedRecordCore().toDelimStr() << std::endl;
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
		ss << __PRETTY_FUNCTION__ << ": Error, no genome found matching = " << genome << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(genomes_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	primaryGenome_ = genome;
}


bfs::path MipsOnGenome::pathToMipBed(const std::string & mipName, const std::string & genome)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(!bib::in(genome, genomes_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no genome found matching " << genome << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(genomes_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(bedsDir_, genome + "_" + mipName + ".bed");
}


bfs::path MipsOnGenome::pathToAllInfoPrimaryGenome() const {
	return bib::files::make_path(tablesDir_,
			"allTarInfo_" + primaryGenome_ + ".tab.txt");
}


bfs::path MipsOnGenome::pathToMipFasta(const std::string & mipName)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaDir_, mipName + ".fasta");
}

VecStr MipsOnGenome::getMips() const {
	auto ret = bib::getVecOfMapKeys(mipArms_->mips_);
	MipNameSorter::sort(ret);
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




void MipsOnGenome::genTables() const{
	if("" != primaryGenome_){
		auto allTarInfo = getMipTarStatsForGenome(primaryGenome_, getMips());
		auto tabOpts = TableIOOpts::genTabFileOut(pathToAllInfoPrimaryGenome());
		tabOpts.out_.overWriteFile_ = true;
		allTarInfo.outPutContents(tabOpts);
	}
}




}  // namespace bibseq
