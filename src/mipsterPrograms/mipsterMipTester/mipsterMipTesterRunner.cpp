
//  mipsterMipTesterRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
//  Copyright (c) 2015 Nick Hathaway. All rights reserved.
//

    
#include "mipsterMipTesterRunner.hpp"
#include <elucidator/BamToolsUtils.h>
#include <TwoBit.h>


namespace bibseq {

mipsterMipTesterRunner::mipsterMipTesterRunner() :
		bib::progutils::programRunner(
				{
					addFunc("callMircosateliteSizes", callMircosateliteSizes, false)
				},
				"mipsterMipTester") {
}//







int mipsterMipTesterRunner::callMircosateliteSizes(
		const bib::progutils::CmdArgs & inputCommands) {
	mipCorePars pars;
	bfs::path genomeFnp = "";
	bfs::path bedFile = "";
	mipsterMipTesterSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	pars.processDefaults(setUp);
	setUp.setOption(genomeFnp, "--genome", "Genome to compare against", true);
	setUp.setOption(bedFile, "--bedFile", "Bed file to extract info for this region", true);
	setUp.setOption(pars.numThreads, "--numThreads", "Number of threads to utilize");
	setUp.processDirectoryOutputName("mipVarCalling_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	BioCmdsUtils cmdRunner(setUp.pars_.verbose_);
	cmdRunner.RunBowtie2Index(genomeFnp);
	cmdRunner.RunFaToTwoBit(genomeFnp);

	auto genomePrefix = genomeFnp;
	genomePrefix.replace_extension(""	);

	SetUpMaster mipMaster(pars.masterDir);
	mipMaster.setMipArmFnp(pars.mipArmsFileName);
	mipMaster.setMipsSampsNamesFnp(pars.mipsSamplesFile);
	mipMaster.loadMipsSampsInfo(pars.allowableErrors);
	mipMaster.setMetaData(pars.sampleMetaFnp);
	auto warnings = mipMaster.checkDirStruct();
	if(!warnings.empty()){
		std::stringstream ss;
		ss << "Error in directory structure, make sure you are in the correct analysis directory" << std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error{ss.str()};
	}

	auto regions = gatherRegions(bedFile.string(), "", setUp.pars_.verbose_);

	//check for overlaps
	for (const auto & pos : iter::range(regions.size())) {
		for (const auto & subPos : iter::range<uint32_t>(0, pos)) {
			if (regions[pos].chrom_ == regions[subPos].chrom_) {
				if ((regions[pos].start_ > regions[subPos].start_
						&& regions[pos].start_ < regions[subPos].end_)
						|| (regions[pos].end_ > regions[subPos].start_
								&& regions[pos].end_ < regions[subPos].end_)) {
					std::stringstream ss;
					ss << __PRETTY_FUNCTION__ << ", error " << regions[pos].uid_
							<< " overlaps with " << regions[subPos].uid_ << "\n";
					ss << "(regions[pos].start_ > regions[subPos].start_&& regions[pos].start_ < regions[subPos].end_)" << std::endl;
					ss << bib::colorBool((regions[pos].start_ > regions[subPos].start_
						&& regions[pos].start_ < regions[subPos].end_)) << std::endl;
					ss << "(regions[pos].end_ > regions[subPos].start_ && regions[pos].end_ < regions[subPos].end_)" << std::endl;
					ss << bib::colorBool((regions[pos].end_ > regions[subPos].start_
								&& regions[pos].end_ < regions[subPos].end_)) << std::endl;
					ss << "regions[pos].end_:" << regions[pos].end_ << std::endl;
					ss << "regions[subPos].start_:" << regions[subPos].start_ << std::endl;
					ss << "regions[subPos].end_:" << regions[subPos].end_ << std::endl;
					throw std::runtime_error { ss.str() };
				}
			}
		}
	}

	auto finalResultPairs = mipMaster.getPairsWithPopClustered(pars.numThreads);

	std::unordered_map<std::string, std::vector<MipFamSamp>> bySample;

	for(const auto & pair : finalResultPairs){
		bySample[pair.samp_].emplace_back(pair);
	}
	std::unordered_map<std::string, bib::sys::RunOutput> outputs;
	std::mutex outputsMut;
	bib::concurrent::LockableQueue<std::string> sampQueue(getVectorOfMapKeys(bySample));
	auto alignSamp =
			[&sampQueue,&bySample,&setUp,&outputsMut,&outputs,&genomeFnp,&mipMaster,&cmdRunner]() {
				std::string samp = "";
				while(sampQueue.getVal(samp)) {
					auto sampDir = bib::files::makeDir(setUp.pars_.directoryName_, bib::files::MkdirPar(samp));
					SeqIOOptions outOpts = SeqIOOptions::genFastqOut(bib::files::make_path(sampDir, samp));
					//outOpts.out_.append_ = true;
					SeqOutput writer(outOpts);
					writer.openOut();
					seqInfo seq;
					for(const auto & pair : bySample.at(samp)) {
						SeqIOOptions opts = SeqIOOptions::genFastqIn(mipMaster.pathPopClusFinalHaplo(pair));
						SeqInput reader(opts);
						reader.openIn();
						while(reader.readNextRead(seq)) {
							MetaDataInName meta(seq.name_);
							auto mipTar = mipMaster.mips_->getMipsForFamily(meta.getMeta("mipFam")).front();
							//seq.prepend(bib::strToLowerRet(mipMaster.mips_->mips_[mipTar].extentionArm_));
							//seq.append(bib::strToLowerRet(mipMaster.mips_->mips_[mipTar].ligationArm_));
							seq.prepend(bib::strToLowerRet(mipMaster.mips_->mips_[mipTar].extentionArm_));
							seq.append(bib::strToLowerRet(mipMaster.mips_->mips_[mipTar].ligationArm_));
							writer.write(seq);
						}
					}
					writer.closeOut();
					auto opts = SeqIOOptions::genFastaIn(outOpts.out_.outName());
					opts.out_.outFilename_ = bib::files::make_path(sampDir, samp + ".sorted.bam");
					opts.out_.outExtention_ = ".sorted.bam";
					auto runOutput = cmdRunner.bowtie2Align(opts,
							genomeFnp);
					{
						std::lock_guard<std::mutex> lock(outputsMut);
						outputs.emplace(samp, runOutput);
					}
				}
			};

	std::vector<std::thread> threads;
	for(uint32_t t = 0; t < pars.numThreads; ++t){
		threads.emplace_back(alignSamp);
	}

	for(auto & t : threads){
		t.join();
	}

	std::unordered_map<std::string, std::vector<GenomicRegion>> regionsByChrom;
	for(const auto & reg : regions){
		regionsByChrom[reg.chrom_].emplace_back(reg);
	}

	//TwoBit::TwoBitFile gReader(genomePrefix.string() + ".2bit");
	aligner alignerObj(600, setUp.pars_.gapInfo_, setUp.pars_.scoring_, false);
	std::unordered_map<std::string, uint32_t> missingRegionsCounts;
	uint32_t unmapped = 0;
	for (const auto & samp : bySample) {
		auto extractOpts = SeqIOOptions::genFastqOut(bib::files::make_path(setUp.pars_.directoryName_, samp.first,
						samp.first + "_extracted"));
		SeqOutput writer(extractOpts);
		writer.openOut();
		if(setUp.pars_.verbose_){
			std::cout << "Sample: " << samp.first << std::endl;
		}
		BamTools::BamAlignment bAln;
		BamTools::BamReader bReader;
		bReader.Open(
				bib::files::make_path(setUp.pars_.directoryName_, samp.first,
						samp.first + ".sorted.bam").string());
		checkBamOpenThrow(bReader);
		auto refData = bReader.GetReferenceData();
		while(bReader.GetNextAlignment(bAln)){
			if(setUp.pars_.debug_){
				std::cout << "\t" << bAln.Name << std::endl;
			}
			if(!bAln.IsMapped()){
				++unmapped;
				if(setUp.pars_.debug_){
					std::cout << "\t" << bAln.Name << " didn't map, next"<< std::endl;
				}
				continue;
			}
			MetaDataInName meta(bAln.Name);
			auto chromName = refData[bAln.RefID].RefName;
			std::vector<GenomicRegion> regionsContained;
			if(bib::in(chromName, regionsByChrom)){
				for(const auto & reg : regionsByChrom.at(chromName)){
					if(reg.start_ >= bAln.Position && reg.end_ <= bAln.GetEndPosition()){
						regionsContained.emplace_back(reg);
					}
					/*
					if((reg.start_ >= bAln.Position && reg.start_ <= bAln.GetEndPosition() ) ||
							(reg.end_ >= bAln.Position && reg.end_ <= bAln.GetEndPosition() ) ){
						regionsContained.emplace_back(reg);
					}*/
				}
			}
			if (regionsContained.empty()) {
				++missingRegionsCounts[meta.getMeta("mipFam")];
				//std::cout << meta.getMeta("mipFam") << " contains no regions" << std::endl;
			}else{
				auto alnInfos = bamAlnToAlnInfoLocal(bAln);
				auto query = seqInfo(bAln.Name, bAln.QueryBases, bAln.Qualities,
						SangerQualOffset);
				//get the sequence for the reference at this alignment location
				std::string refSeq = std::string(alnInfos.begin()->second.localASize_, 'X');
				auto ref = seqInfo("ref", refSeq);
				//convert bam alignment to bibseq::seqInfo object
				alignerObj.alignObjectA_ = baseReadObject(ref);
				alignerObj.alignObjectB_ = baseReadObject(query);
				//set the alignment info and and rearrange the sequences so they can be profiled with gaps
				alignerObj.parts_.lHolder_ = alnInfos.begin()->second;
				alignerObj.rearrangeObjsLocal(alignerObj.alignObjectA_,
						alignerObj.alignObjectB_);
				for(const auto & region : regionsContained){
					auto start = region.start_ - bAln.Position;
					auto stop = region.end_ - bAln.Position;
					auto metaCopy = meta;
					metaCopy.addMeta("region", region.uid_, true);
					auto extractStart = alignerObj.getAlignPosForSeqAPos(start);
					auto extractStop = alignerObj.getAlignPosForSeqAPos(stop);
					auto extract = alignerObj.alignObjectB_.seqBase_.getSubRead(extractStart, extractStop - extractStart);
					extract.removeGaps();
					metaCopy.resetMetaInName(extract.name_);
					writer.write(extract);
				}
			}
			if (regionsContained.size() > 1) {
				if(setUp.pars_.debug_){
					std::cout << meta.getMeta("mipFam")
							<< " contains more than 1 region has " << regionsContained.size();
					std::cout << '\t';
					for (const auto & reg : regionsContained) {
						std::cout << reg.uid_ << ", ";
					}
					std::cout << std::endl;
				}
			}
			if(setUp.pars_.debug_){
				std::cout << "\tdone: " << bAln.Name << std::endl;
			}
		}
		if(setUp.pars_.verbose_){
			std::cout << "Done Sample: " << samp.first << std::endl;
		}
	}

	table outTab(missingRegionsCounts, {"family", "count"});
	outTab.sortTable("family", false);
	if(setUp.pars_.debug_){
		outTab.outPutContents(std::cout, "\t");
		std::cout << "unmapped: " << unmapped << std::endl;
	}
	table finalResults(concatVecs(VecStr{"sample", "region","popUID", "length", "count", "fraction"}, getVectorOfMapKeys(mipMaster.meta_->groupData_)));
	auto samplesKeys = getVectorOfMapKeys(bySample);
	bib::sort(samplesKeys);
	struct NameBarNum{
		std::string name_ = "";
		uint32_t barNum_ = 0;
	};
	for (const auto & sampKey : samplesKeys) {
		auto extractInOpts = SeqIOOptions::genFastqIn(bib::files::make_path(setUp.pars_.directoryName_, sampKey,
				sampKey + "_extracted.fastq"));
		SeqInput reader(extractInOpts);
		auto allSeqs = reader.readAllReadsPtrs<seqInfo>();
		auto splitByRegions = splitSeqsByMetaField(allSeqs, "region");
		auto regionKeys = getVectorOfMapKeys(splitByRegions);
		MipNameSorter::sortByRegion(regionKeys);
		for(const auto & regionKey : regionKeys){
			const auto & region = splitByRegions.at(regionKey);
			auto splitByFam = splitSeqsByMetaField(region, "mipFam");
			if(splitByFam.size() > 0){
				uint32_t bestCount = 0;
				std::string bestFamily = "";
				for(const auto & fam : splitByFam){
					uint32_t totalBarcodeCnt = 0;
					for(const auto & seq : fam.second){
						totalBarcodeCnt += getBarcodeCntFromMipName(seq->name_);
					}
					if(totalBarcodeCnt > bestCount){
						bestCount = totalBarcodeCnt;
						bestFamily = fam.first;
					}
				}
				double totalBarcodeCnt = 0;
				for(const auto & seq : splitByFam.at(bestFamily)){
					totalBarcodeCnt += getBarcodeCntFromMipName(seq->name_);
				}
				std::unordered_map<uint32_t, NameBarNum> byLength;
				for(const auto & seq : splitByFam.at(bestFamily)){
					auto barNum = getBarcodeCntFromMipName(seq->name_);
					MetaDataInName meta(seq->name_);
					if(byLength.end() == byLength.find(len(*seq))){
						byLength[len(*seq)] = NameBarNum{meta.getMeta("h_popUID"), barNum};
					}else{
						byLength[len(*seq)].barNum_ += barNum;
					}
				}
				for(const auto & microLen : byLength){
					auto info = toVecStr(sampKey, regionKey, microLen.second.name_, microLen.first, microLen.second.barNum_, microLen.second.barNum_/totalBarcodeCnt);
					for(const auto & gData : mipMaster.meta_->groupData_){
						info.emplace_back(gData.second->getGroupForSample(sampKey));
					}
					finalResults.addRow(info);
				}
			}else{
				std::stringstream ss;
				ss << __FILE__ << " " << __LINE__ << " " << __PRETTY_FUNCTION__ << " error, shouldn't be happening splitByFam " << " for " << sampKey << " is empty" << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}
	finalResults.outPutContents(TableIOOpts::genTabFileOut(bib::files::make_path(setUp.pars_.directoryName_, "regionLengthInfos.tab.txt")));
	return 0;

}



                    
} // namespace bibseq
