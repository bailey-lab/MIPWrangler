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
#include <elucidator/seqToolsUtils.h>
#include <TwoBit.h>




namespace bibseq {



MipsOnGenome::MipsOnGenome(const pars & inputParameters) :
		inputParameters_(inputParameters),
		gMapper_(inputParameters.gMapperPars_) {

	//std::cout << bib::json::toJson(inputParameters.gMapperPars_) << std::endl;

	bib::files::makeDirP(bib::files::MkdirPar { inputParameters_.mainDir });
	mapDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("mapped"));
	bedsDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("beds"));
	bedsPerGenomeDir_ = bib::files::makeDirP(bedsDir_, bib::files::MkdirPar("perGenome"));

	fastaDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("fastas"));
	fastaByFamilyDir_ = bib::files::makeDirP(fastaDir_, bib::files::MkdirPar("byFamily"));
	armsDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("arms"));
	logDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("logs"));
	tablesDir_ = bib::files::makeDirP(inputParameters_.mainDir, bib::files::MkdirPar("tables"));

	checkInputThrow();

	requireExternalProgramThrow("bowtie2");

}

void MipsOnGenome::checkInputThrow() const {
	std::stringstream ss;
	bool failed = false;
	auto checkForPath = [&failed,&ss](const bfs::path & fnp ) {
		if(!bfs::exists(fnp)) {
			failed = true;
			ss << bib::bashCT::boldRed(fnp.string())<< " needs to exist "<< "\n";
		}
	};
	checkForPath(inputParameters_.inputDir);
	checkForPath(inputParameters_.mipArmsFnp);
	checkForPath(inputParameters_.gMapperPars_.genomeDir_);

	if (failed) {
		std::stringstream outSS;
		outSS << __PRETTY_FUNCTION__ << ", error in checking directory set up"
				<< "\n";
		outSS << ss.str();
		throw std::runtime_error { outSS.str() };
	}
}


void MipsOnGenome::loadInArms(){
	mipArms_ = std::make_unique<MipCollection>(inputParameters_.mipArmsFnp, 6);
}

void MipsOnGenome::createArmFiles(){
	OutOptions opts(armsDir_);
	bib::files::makeDirP(bib::files::MkdirPar(bib::files::make_path(armsDir_, "md5s")));
	opts.overWriteFile_ = true;

	for(const auto & m : mipArms_->mips_){
		if(!opts.outExists()){
			m.second.writeOutArms(opts);
		}else{
			OutOptions md5sumFile(bib::files::make_path(armsDir_, "md5s", m.first));
			if(!md5sumFile.outExists()){
				auto md5Sum = bib::md5(m.second.extentionArm_ + m.second.ligationArm_);
				std::ofstream outMd5File;
				md5sumFile.openFile(outMd5File);
				outMd5File << md5Sum;
				m.second.writeOutArms(opts);
			}else{
				auto md5Sum = bib::md5(m.second.extentionArm_ + m.second.ligationArm_);
				auto md5SumFromFile = bib::files::get_file_contents(md5sumFile.outName(), false);
//				std::cout << md5Sum << std::endl;
//				std::cout << md5SumFromFile << std::endl;
				if(md5Sum != md5SumFromFile){
					m.second.writeOutArms(opts);
				}
			}
		}
	}
}

std::string MipsOnGenome::GenomeMip::uid(const std::string & sep) const {
	return genome_ + sep + mip_;
}


void MipsOnGenome::setMipArmsFnp(const bfs::path & mipArmsFnp){
	bib::files::checkExistenceThrow(mipArmsFnp, __PRETTY_FUNCTION__);
	inputParameters_.mipArmsFnp = mipArmsFnp;
}



void MipsOnGenome::mapArmsToGenomes() {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : gMapper_.genomes_) {
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
			std::string outCheckExt = outStub + "_ext.sorted.bam";
			std::string outCheckLig = outStub + "_lig.sorted.bam";
			if (!  bfs::exists(outCheckExt)
					|| bib::files::firstFileIsOlder(outCheckExt, gMapper_.genomes_.at(pair.genome_)->fnp_)
					|| bib::files::firstFileIsOlder(outCheckExt, pathToMipExtArmFasta(pair.mip_)) ||
					!  bfs::exists(outCheckLig)
					|| bib::files::firstFileIsOlder(outCheckLig, gMapper_.genomes_.at(pair.genome_)->fnp_)
					|| bib::files::firstFileIsOlder(outCheckLig, pathToMipLigArmFasta(pair.mip_))) {
				std::stringstream bowtie2CmdExt;
				bowtie2CmdExt
						<< " bowtie2 -D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end "
						<< "-x " << bib::files::make_path(inputParameters_.gMapperPars_.genomeDir_, pair.genome_)
						<< " -f -U " << bib::files::make_path(armsDir_, pair.mip_)
						<< "_ext-arm.fasta" << " | samtools view - -b | samtools sort - -o " << outStub << "_ext.sorted.bam"
						<< " " << "&& samtools index " << outStub << "_ext.sorted.bam";
				std::stringstream bowtie2CmdLig;
				bowtie2CmdLig
						<< " bowtie2 -D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -a --end-to-end "
						<< "-x " << bib::files::make_path(inputParameters_.gMapperPars_.genomeDir_, pair.genome_)
						<< " -f -U " << bib::files::make_path(armsDir_, pair.mip_)
						<< "_lig-arm.fasta" << " | samtools view - -b | samtools sort - -o " << outStub << "_lig.sorted.bam"
						<< " " << "&& samtools index " << outStub << "_lig.sorted.bam";

				auto bowtie2CmdExtRunOutput = bib::sys::run( { bowtie2CmdExt.str() });
				auto bowtie2CmdLigRunOutput = bib::sys::run( { bowtie2CmdLig.str() });
				if (!bowtie2CmdExtRunOutput.success_) {
					ss << "Failed to map extension arm of " << pair.mip_ << " to "
							<< pair.genome_ << std::endl;
					ss << bowtie2CmdExtRunOutput.stdErr_ << std::endl;
				}
				if (!bowtie2CmdLigRunOutput.success_) {
					ss << "Failed to map ligation arm of " << pair.mip_ << " to "
							<< pair.genome_ << std::endl;
					ss << bowtie2CmdLigRunOutput.stdErr_ << std::endl;
				}
				succes = true;
			}else{
				succes = true;
				ss << outCheckExt << " and " <<  outCheckLig << " is up to date";
			}
			{
				std::lock_guard<std::mutex> lock(logMut);
				log[pair.uid()]["succes"] = succes;
				log[pair.uid()]["message"] = ss.str();
			}
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
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





void MipsOnGenome::genFastas(){
	const VecStr mips = bib::getVecOfMapKeys( mipArms_->mips_);
	const VecStr genomes = bib::getVecOfMapKeys( gMapper_.genomes_	);
	bib::concurrent::LockableQueue<std::string> mipQueue(mips);
	std::mutex logMut;
	Json::Value log;
	log["date"] = bib::getCurrentDateFull();
	auto genFastasFunc = [this, &mipQueue,&genomes,&logMut,&log](){
		std::string mipName;
		while(mipQueue.getVal(mipName)){
			std::stringstream ss;
			bool succes = false;
			auto outOpts = SeqIOOptions::genFastaOut(bib::files::make_path(fastaDir_, mipName + "_withArms"));
			outOpts.out_.overWriteFile_ = true;
			auto trimmedOutOpts = SeqIOOptions::genFastaOut(bib::files::make_path(fastaDir_, mipName));
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
						const std::string genome = bedOpt.first;

						std::vector<GenomicRegion> regions =    bedPtrsToGenomicRegs(getBeds(bedOpt.second->inFilename_.string()));
						std::vector<GenomicRegion> extRegions = bedPtrsToGenomicRegs(getBeds(bib::replaceString(bedOpt.second->inFilename_.string(), ".bed", "-ext.bed")));
						std::vector<GenomicRegion> ligRegions = bedPtrsToGenomicRegs(getBeds(bib::replaceString(bedOpt.second->inFilename_.string(), ".bed", "-lig.bed")));
						if(regions.empty()){
							succes = false;
							ss << "Error in parsing " << bedOpt.second->inFilename_ << ", it was empty\n";
						}else if(regions.size() != extRegions.size() || regions.size() != ligRegions.size()){
							succes = false;
							ss << "Error in parsing " << bedOpt.second->inFilename_ << ", the number of ext regions doesn't equal full regions or lig regions don't equal full regions\n";
						}else{
							for(const auto & regPos : iter::range(regions.size())){
								TwoBit::TwoBitFile twoBitFile(bib::files::make_path(inputParameters_.gMapperPars_.genomeDir_, genome + ".2bit"));
								std::string seq = "";
								twoBitFile[regions[regPos].chrom_]->getSequence(seq, regions[regPos].start_, regions[regPos].end_);
								//consider leaving lower case
								bib::strToUpper(seq);
								if(regions[regPos].reverseSrand_){
									seq = seqUtil::reverseComplement(seq, "DNA");
								}
								MetaDataInName meta(regions[regPos].uid_);
								auto extractionName = genome;
								if(meta.containsMeta("extractionNumber")){
									auto extractionNumber = meta.getMeta("extractionNumber");
									if("0" != extractionNumber){
										extractionName = genome + "." + extractionNumber;
									}
								}
								seqs.emplace_back(seqInfo(extractionName, seq));
								if(extRegions[regPos].overlaps(ligRegions[regPos])){
									ss << "Couldn't trim seq for " << mipName << " in " << genome << " because ligation and ext arms overlap" << "\n";
									ss << "ext: " << extRegions[regPos].chrom_ << ":" << extRegions[regPos].start_ << "-" << extRegions[regPos].end_ << "\n";
									ss << "lig: " << ligRegions[regPos].chrom_ << ":" << ligRegions[regPos].start_ << "-" << ligRegions[regPos].end_ << "\n";
								}else{
									seqInfo trimmedSeq(extractionName, seq);
									readVecTrimmer::trimOffForwardBases(trimmedSeq, extRegions[regPos].getLen());
									readVecTrimmer::trimOffEndBases(trimmedSeq, ligRegions[regPos].getLen());
									trimmedSeqs.emplace_back(trimmedSeq);
								}
							}
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
	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
			threads.emplace_back(std::thread(genFastasFunc));
		}
		bib::concurrent::joinAllJoinableThreads(threads);
	}


	//combine for mip families
	auto mipFamilies = mipArms_->getMipFamilies();

	bib::concurrent::LockableQueue<std::string> mipFamilyQueue(mipFamilies);
	auto combineMipTars = [this,&mipFamilyQueue](){
		std::string mipFamily = "";
		while(mipFamilyQueue.getVal(mipFamily)){
			{
				std::vector<bfs::path> files;
				for(const auto & tar : mipArms_->getMipsForFamily(mipFamily)){
					if(bfs::exists(pathToMipFastaWithoutArms(tar))){
						files.emplace_back(pathToMipFastaWithoutArms(tar));
					}
				}
//				std::cout << mipFamily << std::endl;
//				std::cout << "\t" << bib::conToStr(files, ", ") << std::endl;
				if(!files.empty()){
					OutOptions fastaOut(pathToMipFamilyFastaWithoutArms(mipFamily));
					fastaOut.overWriteFile_ = true;
					concatenateFiles(files, fastaOut);
				}
			}
			{
				std::vector<bfs::path> files;
				for(const auto & tar : mipArms_->getMipsForFamily(mipFamily)){
					if(bfs::exists(pathToMipFasta(tar))){
						files.emplace_back(pathToMipFasta(tar));
					}
				}
				if(!files.empty()){
					OutOptions fastaOut(pathToMipFamilyFasta(mipFamily));
					fastaOut.overWriteFile_ = true;
					concatenateFiles(files, fastaOut);
				}
			}
		}
	};
	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
			threads.emplace_back(std::thread(combineMipTars));
		}
		bib::concurrent::joinAllJoinableThreads(threads);
	}

	OutOptions logOpts(bib::files::make_path(logDir_, "fastaLog-" + bib::getCurrentDate() + ".json"));
	logOpts.outFilename_ = bib::files::findNonexitantFile(logOpts.outFilename_);
	OutputStream outFile(logOpts);
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
	for(const auto & genome : gMapper_.genomes_){
		auto bedFnp = bib::files::make_path(pathToMipBed(tar, genome.first));
		if(bfs::exists(bedFnp)){
			Bed6RecordCore bedCore;
			BioDataFileIO<Bed6RecordCore> bedReader((IoOptions(InOptions(bedFnp))));
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
				Bed6RecordCore bedCore;
				BioDataFileIO<Bed6RecordCore> bedReader((IoOptions(InOptions(bedFnp))));
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
	if (!bib::in(genome, gMapper_.genomes_)) {
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
			Bed6RecordCore bedCore;
			BioDataFileIO<Bed6RecordCore> bedReader((IoOptions(InOptions(bedFnp))));
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
				locs.addRow(genome, mipArms_->mips_.at(tar).regionGroup_, tar,
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

table MipsOnGenome::getMipTarStatsForGenomes(const VecStr & genomes,
		const VecStr & mipTars, bool allRecords) const{
	table ret(getMipTarStatsForGenomeHeader_);
	for (const auto & mipTar : mipTars) {
		auto seqsOpts = SeqIOOptions::genFastaIn(pathToMipFasta(mipTar));
		if(!seqsOpts.inExists()){
			continue;
		}
		seqInfo seq;
		uint32_t totalHapsPossible = 0;
		uint32_t hapNum = 0;
		bool lengthVariation = false;
		std::vector<uint32_t> readLens;
		SeqInput reader(seqsOpts);
		reader.openIn();
		std::unordered_map<std::string, std::shared_ptr<readObject>> readsByAllNames;
		while(reader.readNextRead(seq)){
			readLens.emplace_back(len(seq));
			auto toks = tokenizeString(seq.name_, "-");
			totalHapsPossible+= toks.size();
			++hapNum;
			for(const auto & tok : toks){
				auto outSeq = std::make_shared<readObject>(seq);
				readsByAllNames[tok] = outSeq;
			}
		}
		auto minLen = vectorMinimum(readLens);
		auto maxLen = vectorMaximum(readLens);
		if(maxLen - minLen > 4){
			lengthVariation = true;
		}
		for(const auto & genome : genomes){
			auto bedFnp = bib::files::make_path(pathToMipBed(mipTar, genome));
			if(!bfs::exists(bedFnp)){
				continue;
			}
			//std::cout << mipTar << std::endl;

			Bed6RecordCore bedCore;
			BioDataFileIO<Bed6RecordCore> bedReader((IoOptions(InOptions(bedFnp))));
			bedReader.openIn();
			while(bedReader.readNextRecord(bedCore)){
				MetaDataInName meta(bedCore.name_);
				auto extractionNumber = meta.getMeta("extractionNumber");
				if("0" != extractionNumber && !allRecords){
					break;
				}
				auto genomeName = genome;
				if("0" != extractionNumber){
					genomeName = genome + "." + extractionNumber;
				}
				//bool containsTandems = false;
				auto search = readsByAllNames.find(genomeName);
				std::shared_ptr<readObject> refSeq;
				if(search == readsByAllNames.end()){
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error couldn't find seq for " << genomeName << " for " << mipTar << "\n";
					throw std::runtime_error{ss.str()};
				}else{
					refSeq = search->second;
				}
				double gcContent = 0;
				uint32_t longestHomopolymer = 0;
				if("" != refSeq->seqBase_.name_){
					refSeq->setLetterCount();
					refSeq->counter_.resetAlphabet(true);
					refSeq->counter_.setFractions();
					refSeq->counter_.calcGcContent();
					gcContent = refSeq->counter_.gcContent_;
					refSeq->createCondensedSeq();
					longestHomopolymer = vectorMaximum(refSeq->condensedSeqCount);
//					auto tandems = aligner::findTandemRepeatsInSequence(refSeq->seqBase_.seq_, 2, -2, -7, 20);
//					for(const auto & tandem : tandems){
//						if(tandem.numberOfRepeats_ > 2){
//							containsTandems = true;
//							break;
//						}
//					}
				}
				ret.addRow(mipArms_->mips_[mipTar].regionGroup_,
									mipTar, genome, extractionNumber,
									bedCore.chrom_, bedCore.chromStart_, bedCore.chromEnd_,
									bedCore.strand_, bedCore.length(),
									//containsTandems ? "yes": "no",
									"notChecked",
									lengthVariation ? "yes": "no",
									hapNum, totalHapsPossible,
									static_cast<double>(hapNum)/totalHapsPossible,
									gcContent, longestHomopolymer);
			}
		}
	}
	return ret;
}


const VecStr MipsOnGenome::getMipTarStatsForGenomeHeader_{ "region", "target", "genome", "extractionNumber", "chrom", "start", "end", "strand",
			"length", "containsTandemRepeat", "PossibleLengthVariation", "variantNum",
			"totalVariantsPossible", "variantRatio", "GCContent", "longestHomopolymer" };


table MipsOnGenome::genExtractionNumberTable() const{
	auto genomes = getGenomes();
	table ret{toVecStr("target", genomes)};

	for(const auto & m : getMips()){
		VecStr row{m};
		for(const auto & g : genomes){
			auto bedPath = pathToMipBed(m,g);
			uint32_t count = 0;
			if(bfs::exists(bedPath)){
				auto bRecods = getBeds(bedPath);
				count = bRecods.size();
			}
			row.emplace_back(estd::to_string(count));
		}
		ret.addRow(row);
	}

	return ret;
}

table MipsOnGenome::getMipTarStatsForGenome(const std::string & genome,
		const VecStr & mipTars, bool allRecords) const{
	table ret(getMipTarStatsForGenomeHeader_);
	for (const auto & mipTar : mipTars) {
		auto seqsOpts = SeqIOOptions::genFastaIn(pathToMipFasta(mipTar));
		if(!seqsOpts.inExists()){
			continue;
		}
		seqInfo seq;
		uint32_t totalHapsPossible = 0;
		uint32_t hapNum = 0;
		bool lengthVariation = false;
		std::vector<uint32_t> readLens;
		SeqInput reader(seqsOpts);
		reader.openIn();
		std::unordered_map<std::string, std::shared_ptr<readObject>> readsByAllNames;
		while (reader.readNextRead(seq)) {
			readLens.emplace_back(len(seq));
			auto toks = tokenizeString(seq.name_, "-");
			totalHapsPossible += toks.size();
			++hapNum;
			for (const auto & tok : toks) {
				auto outSeq = std::make_shared<readObject>(seq);
				readsByAllNames[tok] = outSeq;
			}
		}
		auto minLen = vectorMinimum(readLens);
		auto maxLen = vectorMaximum(readLens);
		if (maxLen - minLen > 4) {
			lengthVariation = true;
		}
		auto bedFnp = bib::files::make_path(pathToMipBed(mipTar, genome));
		if (!bfs::exists(bedFnp)) {
			continue;
		}
		//std::cout << mipTar << std::endl;

		Bed6RecordCore bedCore;
		BioDataFileIO<Bed6RecordCore> bedReader((IoOptions(InOptions(bedFnp))));
		bedReader.openIn();
		while (bedReader.readNextRecord(bedCore)) {
			MetaDataInName meta(bedCore.name_);
			auto extractionNumber = meta.getMeta("extractionNumber");
			if ("0" != extractionNumber && !allRecords) {
				break;
			}
			auto genomeName = genome;
			if ("0" != extractionNumber) {
				genomeName = genome + "." + extractionNumber;
			}
			bool containsTandems = false;
			auto search = readsByAllNames.find(genomeName);
			std::shared_ptr<readObject> refSeq;
			if (search == readsByAllNames.end()) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error couldn't find seq for " << genomeName << " for " << mipTar << "\n";
				throw std::runtime_error { ss.str() };
			} else {
				refSeq = search->second;
			}
			double gcContent = 0;
			uint32_t longestHomopolymer = 0;
			if ("" != refSeq->seqBase_.name_) {
				refSeq->setLetterCount();
				refSeq->counter_.resetAlphabet(true);
				refSeq->counter_.setFractions();
				refSeq->counter_.calcGcContent();
				gcContent = refSeq->counter_.gcContent_;
				refSeq->createCondensedSeq();
				longestHomopolymer = vectorMaximum(refSeq->condensedSeqCount);
				auto tandems = aligner::findTandemRepeatsInSequence(
						refSeq->seqBase_.seq_, 2, -2, -7, 20);
				for (const auto & tandem : tandems) {
					if (tandem.numberOfRepeats_ > 2) {
						containsTandems = true;
						break;
					}
				}
			}
			ret.addRow(mipArms_->mips_[mipTar].regionGroup_, mipTar, genome, extractionNumber,
					bedCore.chrom_, bedCore.chromStart_, bedCore.chromEnd_,
					bedCore.strand_, bedCore.length(), containsTandems ? "yes" : "no",
					lengthVariation ? "yes" : "no", hapNum, totalHapsPossible,
					static_cast<double>(hapNum) / totalHapsPossible, gcContent,
					longestHomopolymer);
		}
	}
	return ret;
}


//table MipsOnGenome::getMipTarStatsForGenomes(const VecStr & genomes,
//		const VecStr & mipTars, bool allRecords) const {
//	table ret;
//	for (const auto & genome : genomes) {
//		auto gTable = getMipTarStatsForGenome(genome, mipTars, allRecords);
//		if (ret.empty()) {
//			ret = gTable;
//		} else {
//			ret.rbind(gTable, false);
//		}
//	}
//	return ret;
//}






void MipsOnGenome::genBeds(const comparison & allowableError) {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : gMapper_.genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip{gen.first, m.second.name_});
		}
	}
	bib::concurrent::LockableQueue<GenomeMip> pairsQueue(pairs);
	std::mutex logMut;
	Json::Value log;
	log["date"] = bib::getCurrentDateFull();
	auto genBedsFunc = [this, &pairsQueue,&logMut,&log,&allowableError](){
		GenomeMip pair;
		while(pairsQueue.getVal(pair)){
			std::string outStub = bib::files::make_path(mapDir_,
					pair.genome_ + "_" + pair.mip_).string();
			std::string outCheckExt = outStub + "_ext.sorted.bam";
			std::string outCheckLig = outStub + "_lig.sorted.bam";
			std::stringstream ss;
			bool succes = false;
			if (!bfs::exists(outCheckExt) || !bfs::exists(outCheckLig)) {
				std::lock_guard<std::mutex> lock(logMut);
				ss << "Failed to find both " << outCheckExt  << " and " << outCheckLig << " for " << pair.mip_ << " to "
						<< pair.genome_ << std::endl;
			}else{
				uint32_t insertSizeCutoff = 1000;
				//temporary fix
				if(bib::containsSubString(pair.mip_, "full")){
					insertSizeCutoff = 10000;
				}
				std::vector<std::shared_ptr<AlignmentResults>> alnResultsExt = gatherMapResults(
						outCheckExt, gMapper_.genomes_.at(pair.genome_)->fnpTwoBit_, allowableError);
				std::vector<std::shared_ptr<AlignmentResults>> alnResultsLig = gatherMapResults(
						outCheckLig, gMapper_.genomes_.at(pair.genome_)->fnpTwoBit_, allowableError);

				if(alnResultsExt.empty() || alnResultsLig.empty()){
					if(alnResultsExt.empty()) {
						ss << "Failed to get results from " << outCheckExt << " for " << pair.mip_ << " to "
						<< pair.genome_ << std::endl;
					}
					if(alnResultsLig.empty()) {
						ss << "Failed to get results from " << outCheckLig << " for " << pair.mip_ << " to "
						<< pair.genome_ << std::endl;
					}
				} else {
					OutOptions bedOpts   (bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + ".bed"));
					OutOptions bedExtOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + "-ext.bed"));
					OutOptions bedLigOpts(bib::files::make_path(bedsDir_, pair.genome_ + "_" + pair.mip_ + "-lig.bed"));
					if(bedOpts.outExists() &&
							bib::files::firstFileIsOlder(outCheckExt ,bedOpts.outName()) &&
							bib::files::firstFileIsOlder(outCheckLig ,bedOpts.outName())){
						succes = true;
						ss << bedOpts.outName() << " already up to date" << "\n";
					}else{
						auto extractions = getPossibleGenomeExtracts(alnResultsExt, alnResultsLig, insertSizeCutoff);
						//this sort should put the better matching extraction on top and giving them the lower extraction counts
						bib::sort(extractions, [](const GenomeExtractResult & result1, const GenomeExtractResult & result2){
							return result1.ext_->comp_.distances_.eventBasedIdentity_ + result1.lig_->comp_.distances_.eventBasedIdentity_ >
										 result2.ext_->comp_.distances_.eventBasedIdentity_ + result2.lig_->comp_.distances_.eventBasedIdentity_;
						});
						if(extractions.empty()){
							ss << "Failed to extract any results for " << pair.mip_ << " in " << pair.genome_ << std::endl;
						}else{
							succes = true;
							uint32_t count = 0;
							//full region
							std::ofstream outFile;
							bedOpts.overWriteFile_ = true;
							bedOpts.openFile(outFile);
							//ext
							std::ofstream outExtFile;
							bedExtOpts.overWriteFile_ = true;
							bedExtOpts.openFile(outExtFile);
							//lig
							std::ofstream outLigFile;
							bedLigOpts.overWriteFile_ = true;
							bedLigOpts.openFile(outLigFile);
							for(auto & extraction : extractions){
								MetaDataInName meta;
								meta.addMeta("genome", pair.genome_);
								meta.addMeta("mipTar", pair.mip_);
								meta.addMeta("extractionNumber", count);
								//full region
								extraction.gRegion_->uid_ = meta.createMetaName();
								outFile << extraction.gRegion_->genBedRecordCore().toDelimStr() << std::endl;
								//ext
								extraction.ext_->gRegion_.uid_ = meta.createMetaName();
								outExtFile << extraction.ext_->gRegion_.genBedRecordCore().toDelimStr() << std::endl;
								//lig
								extraction.lig_->gRegion_.uid_ = meta.createMetaName();
								outLigFile << extraction.lig_->gRegion_.genBedRecordCore().toDelimStr() << std::endl;
								++count;
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
	{
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
			threads.emplace_back(std::thread(genBedsFunc));
		}
		bib::concurrent::joinAllJoinableThreads(threads);
	}
	{
		auto genomeNames = getGenomes();
		bib::concurrent::LockableQueue<std::string> genomeQueue(genomeNames);

		auto combineBedsPerGenomes = [this,&genomeQueue](){
			std::string genome = "";
			while(genomeQueue.getVal(genome)){
				std::vector<std::shared_ptr<Bed6RecordCore>> regions;
				OutOptions outOpts(bib::files::make_path(bedsPerGenomeDir_, genome + ".bed"));
				OutputStream outOut(outOpts);
				for(const auto & mip : mipArms_->mips_){
					auto mipBedFnp = bib::files::make_path(bedsDir_, bib::pasteAsStr(genome, "_", mip.first, ".bed"));
					if(bfs::exists(mipBedFnp)){
						addOtherVec(regions, getBeds(mipBedFnp));
					}
				}
				if("" != gMapper_.genomes_.at(genome)->gffFnp_ &&
						bfs::exists(gMapper_.genomes_.at(genome)->gffFnp_)){
					intersectBedLocsWtihGffRecordsPars intersectPars = gMapper_.pars_.gffIntersectPars_;
					intersectPars.gffFnp_ = gMapper_.genomes_.at(genome)->gffFnp_;
					intersectBedLocsWtihGffRecords(regions, intersectPars);
				}
				for(const auto & b : regions){
					outOut << b->toDelimStrWithExtra() << std::endl;
				}
			}
		};
		std::vector<std::thread> threads;
		for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
			threads.emplace_back(std::thread(combineBedsPerGenomes));
		}
		bib::concurrent::joinAllJoinableThreads(threads);

	}

	std::ofstream outFile;
	OutOptions logOpts(bib::files::make_path(logDir_, "bedLog-" + bib::getCurrentDate() + ".json"));
	logOpts.outFilename_ = bib::files::findNonexitantFile(logOpts.outFilename_);
	openTextFile(outFile, logOpts);
	outFile << log << std::endl;

}

std::vector<MipsOnGenome::GenomeMip> MipsOnGenome::genGenomeMipPairs() const {
	std::vector<GenomeMip> pairs;
	for (const auto & gen : gMapper_.genomes_) {
		for (const auto & m : mipArms_->mips_) {
			pairs.emplace_back(GenomeMip { gen.first, m.second.name_ });
		}
	}
	return pairs;
}






bfs::path MipsOnGenome::pathToMipBed(const std::string & mipName, const std::string & genome)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(!bib::in(genome, gMapper_.genomes_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no genome found matching " << genome << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(gMapper_.genomes_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(bedsDir_, genome + "_" + mipName + ".bed");
}


bfs::path MipsOnGenome::pathToAllInfoPrimaryGenome() const {
	return bib::files::make_path(tablesDir_,
			"allTarInfo_" + inputParameters_.gMapperPars_.primaryGenome_ + ".tab.txt");
}

bfs::path MipsOnGenome::pathToAllInfoForGenome(
		const std::string & genome) const {
	return bib::files::make_path(tablesDir_, "allTarInfo_" + genome + ".tab.txt");
}

bfs::path MipsOnGenome::pathToAllInfoAllGenomes() const {
	return bib::files::make_path(tablesDir_,
			"allTarInfo_allGenomes.tab.txt");
}

bfs::path MipsOnGenome::pathToExtractionCounts() const{
	return bib::files::make_path(tablesDir_,
			"extractionCountsTable.tab.txt");
}


bfs::path MipsOnGenome::pathToMipFasta(const std::string & mipName)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaDir_, mipName + "_withArms.fasta");
}

bfs::path MipsOnGenome::pathToMipFastaWithoutArms(const std::string & mipName)const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaDir_, mipName + ".fasta");
}

bfs::path MipsOnGenome::pathToMipFamilyFasta(const std::string & mipFamily)const{
	if(!bib::in(mipFamily, mipArms_->mipFamilies_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipFamily << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaByFamilyDir_, mipFamily + "_withArms.fasta");
}

bfs::path MipsOnGenome::pathToMipFamilyFastaWithoutArms(const std::string & mipFamily)const{
	if(!bib::in(mipFamily, mipArms_->mipFamilies_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipFamily << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(fastaByFamilyDir_, mipFamily + ".fasta");
}



bfs::path MipsOnGenome::pathToMipExtArmFasta(const std::string & mipName) const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(armsDir_, mipName + "_ext-arm.fasta");
}

bfs::path MipsOnGenome::pathToMipLigArmFasta(const std::string & mipName) const{
	if(!bib::in(mipName, mipArms_->mips_)){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ": Error, no mip found matching " << mipName << "\n";
		ss << "Options are: " << bib::conToStr(bib::getVecOfMapKeys(mipArms_->mips_), ", ") << "\n";
		throw std::runtime_error{ss.str()};
	}
	return bib::files::make_path(armsDir_, mipName + "_lig-arm.fasta");
}
//,

VecStr MipsOnGenome::getMips() const {
	auto ret = bib::getVecOfMapKeys(mipArms_->mips_);
	MipNameSorter::sort(ret);
	return ret;
}

VecStr MipsOnGenome::getGenomes() const{
	VecStr genomes = bib::getVecOfMapKeys(gMapper_.genomes_);
	bib::sort(genomes);
	return genomes;
}




void MipsOnGenome::genTables() const{
	auto genomeNames = getGenomes();
	bib::concurrent::LockableQueue<std::string> genomeQueue(genomeNames);
	std::mutex allInfoMut;
	OutOptions allTabOpts (pathToAllInfoAllGenomes());
	allTabOpts.overWriteFile_ = true;
	OutputStream allTabOut(allTabOpts);

	allTabOut << bib::conToStr(getMipTarStatsForGenomeHeader_, "\t") << std::endl;
	auto getGenomeInfo = [this,&genomeQueue,&allInfoMut,&allTabOut](){
		std::string genome = "";
		while(genomeQueue.getVal(genome)){
			auto allTarInfo = getMipTarStatsForGenome(genome, getMips());
			auto tabOpts = TableIOOpts::genTabFileOut(pathToAllInfoForGenome(genome));
			tabOpts.out_.overWriteFile_ = true;
			allTarInfo.outPutContents(tabOpts);
			{
				std::lock_guard<std::mutex> lock(allInfoMut);
				allTarInfo.hasHeader_ = false;
				allTarInfo.outPutContents(allTabOut, "\t");
			}
		}
	};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < inputParameters_.gMapperPars_.numThreads_; ++t){
		threads.emplace_back(std::thread(getGenomeInfo));
	}
	bib::concurrent::joinAllJoinableThreads(threads);

	OutOptions genomeOpts(bib::files::make_path(tablesDir_, "genomes.txt"));
	genomeOpts.overWriteFile_ = true;
	OutOptions targetsOpts(bib::files::make_path(tablesDir_, "mipTargets.txt"));
	targetsOpts.overWriteFile_ = true;

	OutputStream genomeOut(genomeOpts);
	OutputStream targetsOut(targetsOpts);

	printVector(getGenomes(), "\n", genomeOut);
	printVector(getMips(), "\n", targetsOut);

	auto extractionCounts = genExtractionNumberTable();

	auto extractionCountsOpts = TableIOOpts::genTabFileOut(pathToExtractionCounts());
	extractionCountsOpts.out_.overWriteFile_ = true;
	extractionCounts.outPutContents(extractionCountsOpts);


}







}  // namespace bibseq
