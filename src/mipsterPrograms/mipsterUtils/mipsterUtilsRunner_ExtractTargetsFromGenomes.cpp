/*
 * mipsterUtilsRunner_ExtractTargetsFromGenomes.cpp
 *
 *  Created on: Aug 29, 2018
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
#include "mipsterUtilsRunner.hpp"
#include <SeekDeep/utils.h>


namespace njhseq {


int mipsterUtilsRunner::ExtractTargetsFromGenomes(
		const njh::progutils::CmdArgs & inputCommands) {

	bfs::path mipArmsFnp;
	extractBetweenSeqsPars pars;
	mipsterUtilsSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(mipArmsFnp, "--mipArmsFilename",
				"Name of the mip arms file", true);
	pars.setUpCoreOptions(setUp);
	setUp.finishSetUp(std::cout);

	njh::files::makeDir(pars.outputDirPars);
	auto topDir = pars.outputDirPars.dirName_;
	pars.outputDirPars.dirName_ = njh::files::make_path(pars.outputDirPars.dirName_, "extractions");

	OutOptions outOpts(njh::files::make_path(topDir, "mip_targets_as_primers.tab.txt"));
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


	auto forMIPWranglerDir = njh::files::make_path(topDir, "forMIPWrangler");

	njh::files::makeDir(njh::files::MkdirPar{forMIPWranglerDir});
	auto refSeqsDir = njh::files::make_path(forMIPWranglerDir, "refSeqs");
	njh::files::makeDir(njh::files::MkdirPar{refSeqsDir});
	OutOptions lenCutOffsOpts(njh::files::make_path(forMIPWranglerDir, "lenCutOffs.txt"));
	OutputStream lenCutOffsOut(lenCutOffsOpts);
	lenCutOffsOut << "target\tminlen\tmaxlen" << "\n";

	OutOptions overlapStatusOpts(njh::files::make_path(forMIPWranglerDir, "overlapStatuses.txt"));
	OutputStream overlapStatusOut(overlapStatusOpts);
	overlapStatusOut << "target\tstatus" << "\n";

	auto familyNames = mips.getMipFamilies();

	for(const auto & familyName : familyNames){
		auto tarNames = mips.getMipsForFamily(familyName);
		std::vector<uint32_t> readLengths;
		auto primersRemovedFinalFnp = njh::files::make_path(refSeqsDir, familyName + ".fasta");
		auto finalFamilySeqOpts = SeqIOOptions::genFastaOut(primersRemovedFinalFnp);
		SeqOutput finalWriter(finalFamilySeqOpts);
		for(const auto & tar : tarNames){
			auto primersRemovedFnp =      njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + "_primersRemoved.fasta");
			auto extractedSeqsFnp =       njh::files::make_path(pars.outputDirPars.dirName_, tar, tar + ".fasta");
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
                    << "\t" << (minlen > pars.minLenCutOffSizeExpand ? minlen - pars.minLenCutOffSizeExpand : 0)
                    << "\t" << maxlen + pars.maxLenCutOffSizeExpand << std::endl;
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

}  //namespace njhseq

