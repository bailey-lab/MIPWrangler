/*
 * mipsterUtilsRunner_ExtractTargetsFromGenomes.cpp
 *
 *  Created on: Aug 29, 2018
 *      Author: nick
 */



#include "mipsterUtilsRunner.hpp"
#include <SeekDeep/utils.h>


namespace bibseq {


int mipsterUtilsRunner::ExtractTargetsFromGenomes(
		const bib::progutils::CmdArgs & inputCommands) {

	bfs::path mipArmsFnp;
	extractBetweenSeqsPars pars;
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mipArmsFnp, "--mipArmsFilename",
				"Name of the mip arms file", true);
	pars.setUpCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

	bib::files::makeDir(pars.outputDirPars);
	auto topDir = pars.outputDirPars.dirName_;
	pars.outputDirPars.dirName_ = bib::files::make_path(pars.outputDirPars.dirName_, "extractions");

	OutOptions outOpts(bib::files::make_path(topDir, "mip_targets_as_primers.tab.txt"));
	MipCollection mips(mipArmsFnp, 6);
	{
		//write out mip arms
		OutputStream out(outOpts);
		out << "target\tforward\treverse" << std::endl;
		auto tarNames = mips.getMipTars();
		for(const auto & mipName : tarNames){
			out << mipName
					<< "\t" << mips.mips_[mipName].extentionArm_
					<< "\t" << mips.mips_[mipName].ligationArm5_to_3prime_ << std::endl;
		}
	}
	//set written out mip arms as the targets
	pars.primersFile = outOpts.outName();

	PrimersAndMids ids(pars.primersFile);
	if(0 == ids.getTargets().size() ){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error in reading in target primers file " << pars.primersFile << "\n";
		ss << "Make sure there is a line that starts with target in file" << "\n";
		throw std::runtime_error{ss.str()};
	}
	ids.initPrimerDeterminator();

	//extract sequences
	extractBetweenSeqs(ids, pars);
	setUp.startARunLog(topDir.string());


	auto forMIPWranglerDir = bib::files::make_path(topDir, "forMIPWrangler");

	bib::files::makeDir(bib::files::MkdirPar{forMIPWranglerDir});
	auto refSeqsDir = bib::files::make_path(forMIPWranglerDir, "refSeqs");
	bib::files::makeDir(bib::files::MkdirPar{refSeqsDir});
	OutOptions lenCutOffsOpts(bib::files::make_path(forMIPWranglerDir, "lenCutOffs.txt"));
	OutputStream lenCutOffsOut(lenCutOffsOpts);
	lenCutOffsOut << "target\tminlen\tmaxlen" << "\n";

	OutOptions overlapStatusOpts(bib::files::make_path(forMIPWranglerDir, "overlapStatuses.txt"));
	OutputStream overlapStatusOut(overlapStatusOpts);
	overlapStatusOut << "target\tstatus" << "\n";

	auto familyNames = mips.getMipFamilies();

	for(const auto & familyName : familyNames){
		auto tarNames = mips.getMipsForFamily(familyName);
		std::vector<uint32_t> readLengths;
		auto primersRemovedFinalFnp = bib::files::make_path(refSeqsDir, familyName + ".fasta");
		auto finalFamilySeqOpts = SeqIOOptions::genFastaOut(primersRemovedFinalFnp);
		SeqOutput finalWriter(finalFamilySeqOpts);
		for(const auto & tar : tarNames){
			auto primersRemovedFnp =      bib::files::make_path(pars.outputDirPars.dirName_, tar, tar + "_primersRemoved.fasta");
			auto extractedSeqsFnp =       bib::files::make_path(pars.outputDirPars.dirName_, tar, tar + ".fasta");
			if(bfs::exists(primersRemovedFnp)){
				{
					SeqInput reader(SeqIOOptions::genFastaIn(primersRemovedFnp));
					reader.openIn();
					seqInfo seq;
					while(reader.readNextRead(seq)){
						finalWriter.openWrite(seq);
					}
				}
				{
					SeqInput reader(SeqIOOptions::genFastaIn(extractedSeqsFnp));
					reader.openIn();
					seqInfo seq;
					while(reader.readNextRead(seq)){
						readLengths.emplace_back(len(seq));
					}
				}
			}else{
				std::cerr << "Warning, no sequences extracted for " << tar << std::endl;
			}
		}
		if(!readLengths.empty()){
			auto minlen = vectorMinimum(readLengths);
			auto maxlen = vectorMaximum(readLengths);
			lenCutOffsOut << familyName
					<< "\t" << (minlen > pars.lenCutOffSizeExpand ? minlen - pars.lenCutOffSizeExpand : 0)
					<< "\t" << maxlen + pars.lenCutOffSizeExpand << std::endl;
			uint32_t finalSize = maxlen + pars.barcodeSize;
			if(finalSize > pars.pairedEndLength && finalSize < (2* pars.pairedEndLength - 10)){
				overlapStatusOut << familyName
						<< "\t" << "R1EndsInR2" << std::endl;
			} else if(finalSize < pars.pairedEndLength){
				overlapStatusOut << familyName
						<< "\t" << "R1BeginsInR2" << std::endl;
			} else {
				overlapStatusOut << familyName
						<< "\t" << "NoOverLap" << std::endl;
			}
		}
	}
	return 0;
}

}  //namespace bibseq

