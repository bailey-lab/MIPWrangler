
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
					addFunc("testingVariationCalling", testingVariationCalling, false)
				},
				"mipsterMipTester") {
}//



bib::sys::RunOutput bowtie2Align(const bfs::path & input,
		const bfs::path & genomePrefix, const bfs::path & output) {
	std::stringstream templateCmd;
	bfs::path outputFnp = bib::appendAsNeededRet(output.string(), ".sorted.bam");
	templateCmd << "bowtie2 -U " << input << " -x " << genomePrefix << " "
			<< "| samtools view - -b " << "| samtools sort - -o " << outputFnp << " "
			<< "&& samtools index " << outputFnp;
	auto ret = bib::sys::run( { templateCmd.str() });
	BioCmdsUtils::checkRunOutThrow(ret, __PRETTY_FUNCTION__);
	return ret;
}



int mipsterMipTesterRunner::testingVariationCalling(
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
			[&sampQueue,&bySample,&setUp,&outputsMut,&outputs,&genomePrefix,&mipMaster]() {
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
							writer.write(seq);
						}
					}
					writer.closeOut();
					auto runOutput = bowtie2Align(outOpts.out_.outName(),
							genomePrefix,
							bib::files::make_path(sampDir, samp + ".sorted.bam"));
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

	std::unordered_map<std::string, uint32_t> missingRegionsCounts;

	for (const auto & samp : bySample) {
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

					if(reg.start_ + 30 >= bAln.Position && reg.end_ -30 <= bAln.GetEndPosition()){
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
			}
			if (regionsContained.size() > 1) {
				std::cout << meta.getMeta("mipFam")
						<< " contains more than 1 region has " << regionsContained.size();
				std::cout << '\t';
				for (const auto & reg : regionsContained) {
					std::cout << reg.uid_ << ", ";
				}
				std::cout << std::endl;
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
	outTab.outPutContents(std::cout, "\t");

	return 0;

}



                    
} // namespace bibseq
