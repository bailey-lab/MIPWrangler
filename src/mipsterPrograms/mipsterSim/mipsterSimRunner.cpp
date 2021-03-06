
//  mipsterSimRunner.cpp
//
//  Created by Nick Hathaway on 2015/08/24.
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

    
#include "mipsterSimRunner.hpp"

#include <njhseq/BamToolsUtils.h>
#include <njhseq/objects/BioDataObject.h>
//#include <elucidator/simulation.h>
//#include <elucidator/seqToolsUtils.h>


namespace njhseq {

mipsterSimRunner::mipsterSimRunner()
    : njh::progutils::ProgramRunner({
		addFunc("extractMipCaptureSequence", extractMipCaptureSequence, false),
		addFunc("simMips", simMips, false),
		addFunc("createArmFastaFiles", createArmFastaFiles, false),
		addFunc("searchForArms", searchForArms, false),
		addFunc("mipSimSetup", mipSimSetup, false),
		addFunc("testMipExtract", testMipExtract, true)
	}, "MIPWranglerSim") {

}//

int mipsterSimRunner::testMipExtract(const njh::progutils::CmdArgs & inputCommands) {
	std::string genomeDirectory = "";
	std::string mipFile = "";
	std::string genomeNames = "";
	std::string regionNames = "";
	bool includeTrim = false;
	uint32_t numThreads = 1;
	seqSetUp setUp(inputCommands);
	setUp.setOption(numThreads, "--numThreads", "Number of threads to use");
	setUp.setOption(includeTrim, "--includeTrim", "Include a fasta with arms trimmed off targets");
	setUp.setOption(genomeDirectory, "--genomeDirectory", "Genome Directory", true);
	setUp.setOption(genomeNames, "--genomes", "Genomes");
	setUp.setOption(regionNames, "--regionNames", "Region Names");
	setUp.setOption(mipFile, "--mipArmsFilename", "Mip Arms File", true);
	setUp.processDirectoryOutputName("mipSimSetup_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	MipCollection mCol(mipFile, 6);
	VecStr regions;
	if("" != regionNames){
		regions = njh::tokenizeString(regionNames, ",");
	}else{
		regions = mCol.getMipRegions();
	}
	for(const auto & reg : regions){
		auto mipTars = mCol.getMipTarsForRegion(reg);
		for(const auto & mipTar : mipTars){
			auto extOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, mipTar + "_ext-arm").string());
			auto ligOpts = SeqIOOptions::genFastaOut(njh::files::make_path(setUp.pars_.directoryName_, mipTar + "_lig-arm").string());
			seqInfo extArm("[mipTar=" + mCol.mips_[mipTar].name_ + ";mipFam=" + mCol.mips_[mipTar].familyName_ +";]", mCol.mips_[mipTar].extentionArm_);
			seqInfo ligArm("[mipTar=" + mCol.mips_[mipTar].name_ + ";mipFam=" + mCol.mips_[mipTar].familyName_ +";]", seqUtil::reverseComplement(mCol.mips_[mipTar].ligationArm_, "DNA"));
			SeqOutput::write(std::vector<seqInfo>{extArm}, extOpts);
			SeqOutput::write(std::vector<seqInfo>{ligArm}, ligOpts);
		}
	}

	VecStr genomes;
	if(genomeNames != ""){
		 genomes = njh::tokenizeString(genomeNames, ",");
	}else{
		auto files = njh::files::filesInFolder(genomeDirectory);
		for(const auto & f : files){
			if(!bfs::is_directory(f)){
				if(njh::endsWith(f.string(), ".fasta")){
					genomes.emplace_back(f.filename().replace_extension("").string());
				}
			}
		}
	}

	VecStr cmds;
	for(const auto & region : regions){
		auto mipTars = mCol.getMipTarsForRegion(region);
		for(const auto & mipTar : mipTars){
			for(const auto & genome : genomes){
				std::string outStub = njh::files::make_path(setUp.pars_.directoryName_, genome + "_" + mipTar).string();
				std::stringstream bowtie2Cmd;
				//bowtie2Cmd << " bowtie2  -p 1 -D 20 -R 3 -N 1 -L 18 -i S,1,0.5 --gbar 1 -k 500 --end-to-end "
				bowtie2Cmd << " bowtie2  -p 1 -D 20 -R 3 -N 1 -L 18 -i S,1,0.5 --gbar 1 --end-to-end "
						<< "-x " << njh::files::make_path(genomeDirectory, genome)
						<< " -f -1 " << njh::files::join(setUp.pars_.directoryName_,mipTar)  << "_ext-arm.fasta"
						<< " -2 " << njh::files::join(setUp.pars_.directoryName_,mipTar) << "_lig-arm.fasta"
						<< " -S " << outStub + ".sam";
				std::stringstream samtoolsCmds;
				samtoolsCmds
						<< "samtools view -Sb " << outStub << ".sam | "
						<< "samtools sort - -o " << outStub << ".sorted.bam && "
						<< "samtools index " << outStub << ".sorted.bam";
				std::stringstream extractMipBedSeqCmd;
				extractMipBedSeqCmd << "MIPWrangler extractMipCaptureSequence -bam "
						<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << mipTar << ".sorted.bam -out "
						<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << mipTar << "-mips.bed";
				std::stringstream getFastaFromBedCmd;
				getFastaFromBedCmd << "elucidator getFastaWithBed "
						<< "--twoBit " << genomeDirectory
						<< genome << ".2bit "
						<< "--file "<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << mipTar << "-mips.bed "
						<< "--outFile "<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << mipTar << "-mips.fasta";
				VecStr currentCmd{bowtie2Cmd.str(), samtoolsCmds.str(), extractMipBedSeqCmd.str(), getFastaFromBedCmd.str()};
				if(includeTrim){
					std::stringstream trimFastaCmd;
					trimFastaCmd << "sequenceTools trimBetweenSeqs  "
							<< "--fasta " << njh::files::join(setUp.pars_.directoryName_,genome) << "_" << mipTar << "-mips.fasta "
							<< "--forwardSeq " << mCol.mips_.at(mipTar).extentionArm_ << " "
							<< "--backSeq " << mCol.mips_.at(mipTar).ligationArm_ << " "
							<< "--hqMismatches " << mCol.allowableArmError_ << " "
							<< "--out "<< njh::files::join(setUp.pars_.directoryName_, "trimmed_" + genome) << "_" << mipTar << "-mips.fasta ";
					currentCmd.emplace_back(trimFastaCmd.str());
				}
				cmds.emplace_back(njh::conToStr(currentCmd, " && "));
			}
		}
	}
	auto runOuts = njh::sys::runCmdsThreaded(cmds, numThreads, setUp.pars_.verbose_, setUp.pars_.debug_);
	Json::Value logJson = njh::json::toJson(runOuts);
	std::ofstream logFile;
	openTextFile(logFile, OutOptions(njh::files::make_path( setUp.pars_.directoryName_,"cmdLogs.json")));
	logFile << logJson;
	return 0;
}








int mipsterSimRunner::mipSimSetup(const njh::progutils::CmdArgs & inputCommands) {
	std::string genomeDirectory = "";
	std::string mipArmsDirectory = "";
	std::string outDir = "";
	std::string mipFile = "";
	std::string genomeNames = "";
	std::string regionNames = "";
	bool includeTrim = false;
	seqSetUp setUp(inputCommands);
	setUp.setOption(includeTrim, "--includeTrim", "Include a fasta with arms trimmed off targets");
	setUp.setOption(genomeDirectory, "--genomeDirectory", "Genome Directory", true);
	setUp.setOption(mipArmsDirectory, "--mipArmsDirectory", "Mip Arms Directory", true);
	setUp.setOption(genomeNames, "--genomes", "Genomes", true);
	setUp.setOption(regionNames, "--regionNames", "Region Names", true);
	setUp.setOption(mipFile, "--mipArmsFilename", "Mip Arms File", true);
	setUp.processDirectoryOutputName("mipSimSetup_TODAY", true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);

	MipCollection mCol(mipFile, 6);
	auto regions = njh::tokenizeString(regionNames, ",");
	auto genomes = njh::tokenizeString(genomeNames, ",");
	Json::Value logJson;
	for(const auto & region : regions){
		/**@todo add saftey check for region */
		for(const auto & genome : genomes){
			std::stringstream searchForArmsCmd;
			searchForArmsCmd << "MIPWrangler searchForArms "
					<< "--genomeDirectory " << genomeDirectory << " "
					<< "--genome " << genome << " "
					<< "--mipGene " << region << " "
					<< "--outDir " << setUp.pars_.directoryName_ << " "
					<< "--mipArmsDirectory " << mipArmsDirectory;
			std::stringstream extractMipBedSeqCmd;
			extractMipBedSeqCmd << "MIPWrangler extractMipCaptureSequence -bam "
					<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << region << ".sorted.bam -out "
					<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << region << "-mips.bed";
			std::stringstream getFastaFromBedCmd;
			getFastaFromBedCmd << "elucidator getFastaWithBed "
					<< "--twoBit " << genomeDirectory
					<< genome << ".2bit "
					<< "--file "<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << region << "-mips.bed "
					<< "--outFile "<< njh::files::join(setUp.pars_.directoryName_,genome) << "_" << region << "-mips.fasta";
			auto searchForArmsCmdOutput = njh::sys::run(VecStr{searchForArmsCmd.str()});
			auto extractMipBedSeqCmdOutput = njh::sys::run(VecStr{extractMipBedSeqCmd.str()});
			auto getFastaFromBedCmdOutput = njh::sys::run(VecStr{getFastaFromBedCmd.str()});
			logJson[region][genome]["searchForArmsCmd"] = searchForArmsCmdOutput.toJson();
			logJson[region][genome]["extractMipBedSeqCmdOutput"] = extractMipBedSeqCmdOutput.toJson();
			logJson[region][genome]["getFastaFromBedCmdOutput"] = getFastaFromBedCmdOutput.toJson();
			if(includeTrim){
				std::stringstream outFastaName;
				outFastaName << njh::files::join(setUp.pars_.directoryName_,genome) << "_" << region << "-mips.fasta";
				auto opts = SeqIOOptions::genFastaIn(outFastaName.str(), false);
				SeqInput reader(opts);
				if(!reader.isFirstEmpty()){
					auto seqs = reader.readAllReads<readObject>();
					for(auto & seq : seqs){
						seq.processNameForMeta();
						auto mipName = seq.getMeta("mipTar");
						readVecTrimmer::trimOffEndBases(seq,mCol.mips_.at(mipName).ligationArm_.size());
						readVecTrimmer::trimOffForwardBases(seq,mCol.mips_.at(mipName).extentionArm_.size());
					}
					std::stringstream outOutFastaName;
					outOutFastaName << njh::files::join(setUp.pars_.directoryName_, "trimmed_" + genome) << "_" << region << "-mips.fasta";
					auto outOpts = SeqIOOptions::genFastaOut(outOutFastaName.str());
					SeqOutput::write(seqs, outOpts);
				}
			}
		}
	}
	std::ofstream logFile;
	openTextFile(logFile,OutOptions(njh::files::make_path(setUp.pars_.directoryName_ ,"cmdLogs.json")));
	logFile << logJson;
	return 0;
}

int mipsterSimRunner::createArmFastaFiles(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string mipFile = "";
	std::string outDir = "";
	bool overWrite = false;
	setUp.setOption(mipFile, "--mipArmsFilename", "Mip Arms File", true);
	setUp.setOption(overWrite, "--overWrite", "Over Write Previous Files");
	setUp.setOption(outDir, "--outDir", "Output Directory, will be created if needed", true);
	setUp.finishSetUp(std::cout);
	MipCollection mCol(mipFile, 6);
	njh::files::makeDirP(njh::files::MkdirPar(outDir));

	auto regions = mCol.getMipRegions();
	for(const auto & region : regions){
		auto mips = mCol.getMipTarsForRegion(region);
		MultiSeqIO out;
		auto extOpts = SeqIOOptions::genFastaOut(njh::files::make_path(outDir, region + "_ext-arm").string());
		extOpts.out_.overWriteFile_ = overWrite;
		auto ligOpts = SeqIOOptions::genFastaOut(njh::files::make_path(outDir, region + "_lig-arm").string());
		ligOpts.out_.overWriteFile_ = overWrite;
		out.addReader("ext", extOpts);
		out.addReader("lig", ligOpts);
		for(const auto & m : mips){
			seqInfo extArm("[mipTar=" + mCol.mips_[m].name_ + ";mipFam=" + mCol.mips_[m].familyName_ +";]", mCol.mips_[m].extentionArm_);
			seqInfo ligArm("[mipTar=" + mCol.mips_[m].name_ + ";mipFam=" + mCol.mips_[m].familyName_ +";]", seqUtil::reverseComplement(mCol.mips_[m].ligationArm_, "DNA"));
			out.openWrite("ext", extArm);
			out.openWrite("lig", ligArm);
		}
	}

	return 0;
}

int mipsterSimRunner::searchForArms(const njh::progutils::CmdArgs & inputCommands) {
	std::string genomeDirectory = "";
	std::string mipArmsDirectory = "";
	std::string mipGene = "";
	std::string genome = "";
	std::string outDir = "";
	OutOptions logOpts("", ".json");

	seqSetUp setUp(inputCommands);
	setUp.setOption(genomeDirectory, "--genomeDirectory", "Genome Directory", true);
	setUp.setOption(mipArmsDirectory, "--mipArmsDirectory", "Mip Arms Directory", true);
	setUp.setOption(mipGene, "--mipGene", "Mip Gene", true);
	setUp.setOption(genome, "--genome", "Genome", true);
	setUp.setOption(outDir, "--outDir", "outDir", true);
	setUp.processWritingOptions(logOpts);
	setUp.finishSetUp(std::cout);

	std::string outStub = njh::files::make_path(outDir, genome + "_" + mipGene).string();
	if("" == logOpts.outFilename_){
		logOpts.outFilename_ = outStub + "Log";
	}
	std::stringstream bowtie2Cmd;
	bowtie2Cmd << " bowtie2  -p 60 -D 20 -R 3 -N 1 -L 18 -i S,1,0.5 --gbar 1 -k 500 --end-to-end "
			<< "-x " << njh::files::make_path(genomeDirectory, genome)
			<< " -f -1 " << mipGene << "_ext-arm.fasta"
			<< " -2 " << mipGene << "_lig-arm.fasta"
			<< " -S " << outStub + ".sam";
	std::stringstream samtoolsCmds;
	samtoolsCmds
			<< "samtools view -Sb " << outStub << ".sam | "
			<< "samtools sort - -o " << outStub << ".sorted.bam && "
			<< "samtools index " << outStub << ".sorted.bam";

	auto bowtie2Run = njh::sys::run(VecStr{bowtie2Cmd.str()});
	auto samtoolsRun = njh::sys::run(VecStr{samtoolsCmds.str()});
	Json::Value log;
	auto & main = log["mainCmd"];
	main["wd"] = setUp.commands_.workingDir_;
	main["cmd"] = setUp.commands_.commandLine_;
	log["bowtie2Cmd"] = bowtie2Run.toJson();
	log["samtoolsCmd"] = samtoolsRun.toJson();
	std::ofstream logFile;
	openTextFile(logFile, logOpts);
	logFile << log;
	return 0;
}









int mipsterSimRunner::extractMipCaptureSequence(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.processDefaultReader({"-bam"}, true);
	setUp.pars_.ioOptions_.out_.outExtention_ = ".bed";
	setUp.pars_.ioOptions_.out_.outFileFormat_ = "bed";
	setUp.finishSetUp(std::cout);

	std::vector<MipMapResult> results = getMipMapResults(setUp.pars_.ioOptions_.firstName_);


	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.ioOptions_.out_);
	for (auto & result : results) {

		if (result.isConcordant() && result.isMapped()) {
			outFile << result.region_.genBedRecordCore().toDelimStr() << std::endl;
		} else {
			if (!result.isMapped()) {
				if (!result.extAln_.IsMapped()) {
					std::cerr << "Result for " << result.mipName_ << " mapped to "
							<< result.genomeName_ << " extension arm didn't map "
							<< std::endl;
				}
				if (!result.ligAln_.IsMapped()) {
					std::cerr << "Result for " << result.mipName_ << " mapped to "
							<< result.genomeName_ << " ligation arm didn't map" << std::endl;
				}
			} else if (!result.isConcordant()) {
				std::cerr << "Result for " << result.mipName_ << " mapped to "
						<< result.genomeName_ << " had discordant arms" << std::endl;
			}
		}

	}
	return 0;
}


struct LibAdundInfo{
	struct GenomeAbund{
		GenomeAbund(const std::string & genome, double relAbund): genome_(genome), relAbund_(relAbund){}
		std::string genome_;
		double relAbund_;
		uint32_t templateAmount_ = 0;
	};
	LibAdundInfo(const std::string & libName): libName_(libName){}
	std::string libName_;
	std::vector<GenomeAbund> genomesRelAbundance_;

	void addGenome(const std::string & genome, double relAbund){
		if(relAbund > 0){
			if(njh::has(genomesRelAbundance_, genome, [](const GenomeAbund & conVal, const std::string & addVal){ return conVal.genome_ == addVal;})){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": trying to add a genome that has already been added, " << genome << "\n";
				throw std::runtime_error{ss.str()};
			}
			genomesRelAbundance_.emplace_back(GenomeAbund{genome, relAbund});
		}else{
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": trying to add a library with a relative abundance of 0\n";
			throw std::runtime_error{ss.str()};
		}
	}

	void setTemplateAmount(const uint32_t templateAmount) {
		for (auto & gAbund : genomesRelAbundance_) {
			gAbund.templateAmount_ = std::round(gAbund.relAbund_ * templateAmount);
		}
	}

	void resetAbundances(){
		if(0 == numOfGenomes()){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, trying to reset abundances when library is empty\n";
			throw std::runtime_error{ss.str()};
		}
		double total = 0;
		for(const auto & gAbund : genomesRelAbundance_){
			total += gAbund.relAbund_;
		}
		for (auto & gAbund : genomesRelAbundance_) {
			gAbund.relAbund_ /= total;
		}
	}

	uint32_t numOfGenomes()const{
		return len(genomesRelAbundance_);
	}

	VecStr genomeNamesForStartingTemplate() const {
		VecStr ret;
		for (const auto & gAbund : genomesRelAbundance_) {
			addOtherVec(ret, VecStr(gAbund.templateAmount_, gAbund.genome_));
		}
		return ret;
	}

};

std::map<std::string, LibAdundInfo> processAbundanceLibaries(const std::string & abundanceFile, const VecStr & availableGenomes){
	std::map<std::string, LibAdundInfo> ret;
	table abundTab(abundanceFile, "\t", true);
	VecStr genomes = abundTab.getColumn("Genomes");
	for(const auto & genome : genomes){
		if(!njh::in(genome, availableGenomes)){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error, trying to set up a library with an unavailable genome, " << genome << "\n";
			ss << "Options are: " << njh::conToStr(availableGenomes, ",") << "\n";
			throw std::runtime_error{ss.str()};
		}
	}
	for(const auto & colPos : iter::range(abundTab.nCol())){
		if("Genomes" != abundTab.columnNames_[colPos]){
			std::string currentLibName = abundTab.columnNames_[colPos];
			if(njh::in(currentLibName, ret)){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, trying to add library with a name anther library that has already been added, " << currentLibName << "\n";
				throw std::runtime_error{ss.str()};
			}
			VecStr currentLibAbunds = abundTab.getColumn(colPos);
			LibAdundInfo lib(currentLibName);
			for(const auto & rowPos : iter::range(currentLibAbunds.size())){
				double abund = njh::lexical_cast<double>(currentLibAbunds[rowPos]);
				if(abund > 0){
					lib.addGenome(genomes[rowPos], abund);
				}
			}
			lib.resetAbundances();
			ret.emplace(currentLibName, lib);
		}
	}
	return ret;
}

void checkLibraryAbundancesThrow(const std::map<std::string, LibAdundInfo> & libAbunds, uint32_t startingTemplate){
	for(const auto & lib : libAbunds){
		for(const auto & gAbund : lib.second.genomesRelAbundance_){
			uint32_t templateAmount = std::round(gAbund.relAbund_ * startingTemplate);
			if(0 == templateAmount){
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ": Error, abundance, "
						<< gAbund.relAbund_ << ", for " << gAbund.genome_
						<< " isn't enough with a starting template of "
						<< startingTemplate << " to simulate any template" << "\n";
				ss << "Minimum starting template for this abundance would be " << std::ceil(1.0/gAbund.relAbund_) << "\n";
				throw std::runtime_error{ss.str()};
			}
		}
	}
}

int mipsterSimRunner::simMips(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
	setUp.finishSetUp(std::cout);

	std::cout <<"Feature coming" << std::endl;
	return 0;
}

//namespace sim{
//
//
//void simMipLib(const LibAdundInfo & libInfo,
//		const MipCollection & mCol,
//		const VecStr & regions,
//		const std::unordered_map<std::string, std::unordered_map<std::string, seqInfo>> & seqs,
//				const std::string & workingDir,
//				double captureEfficiency,
//				uint64_t intErrorRate,
//				uint32_t finalReadAmount,
//				uint32_t pcrRounds,
//				uint32_t initialPcrRounds,
//				uint32_t numThreads,
//				bool verbose) {
//	if(initialPcrRounds >= pcrRounds){
//		std::stringstream ss;
//		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
//		ss << "initialPcrRounds should be less than pcrRounds" << std::endl;
//		ss << "pcrRounds:" << pcrRounds << std::endl;
//		ss << "initialPcrRounds:" << initialPcrRounds << std::endl;
//		throw std::runtime_error{ss.str()};
//	}
//
//	std::ofstream libOutFile;
//	openTextFile(libOutFile, OutOptions(njh::files::make_path(workingDir, libInfo.libName_ + ".fasta")));
//	std::mutex seqFileLock;
//
//
//	if(verbose){
//		std::cout << "Simulating library " << libInfo.libName_ << std::endl;
//	}
//	njh::randomGenerator gen;
//	std::unordered_map<std::string,std::unordered_map<std::string, uint64_t>> allSeqCounts;
//	std::unordered_map<std::string,std::string> barcodedSeqs;
//	std::unordered_map<std::string,std::pair<uint64_t,uint64_t>> templateNonMutated;
//	auto genomesToSample = libInfo.genomeNamesForStartingTemplate();
//	uint32_t captureAmount = std::round(genomesToSample.size() * captureEfficiency);
//	if (0 == captureAmount) {
//		std::stringstream ss;
//		ss << "Error in " << __PRETTY_FUNCTION__ << std::endl;
//		ss << "With capture efficiency of" << captureEfficiency
//				<< " and starting template amount: " << genomesToSample.size()
//				<< " the capture amount would be zero" << std::endl;
//		throw std::runtime_error { ss.str() };
//	}
//	for(const auto & region : regions){
//
//
////		for(const auto & gAbund : libInfo.genomesRelAbundance_){
////			SeqIOOptions inOpts = SeqIOOptions::genFastaIn(njh::files::make_path(mipArmsDir, gAbund.genome_ + "_" + region + "-mips.fasta").string());
////			seqInfo seq;
////			SeqInput reader(inOpts);
////			reader.openIn();
////			while(reader.readNextRead(seq)){
////				std::unordered_map<std::string, std::string> meta;
////				seq.processNameForMeta(meta);
////				seqs[meta["genome"]][meta["mipTar"]] = seq;
////			}
////		}
//		auto mips = mCol.getMipTarsForRegion(region);
//		for(const auto & mipName : mips){
//			auto mip = mCol.mips_.at(mipName);
//			auto capturedGenomes = gen.unifRandSelectionVec(genomesToSample, captureAmount, false);
//			uint32_t count = 0;
//			if(verbose){
//				std::cout << "\tSimulating " << mipName << std::endl;
//			}
//			njh::ProgressBar pBar(capturedGenomes.size());
//			for(const auto & genome : capturedGenomes){
//				if(verbose){
//					pBar.outputProgAdd(std::cout, 1, true);
//				}
//				std::string extBarcode = simulation::evenRandStr(mip.extBarcodeLen_, std::vector<char>{'A', 'C', 'G', 'T'}, gen);
//				std::string ligBarcode = simulation::evenRandStr(mip.ligBarcodeLen_, std::vector<char>{'A', 'C', 'G', 'T'}, gen);
//				std::string seq = extBarcode + seqs.at(genome).at(mipName).seq_ + ligBarcode;
//				std::unordered_map<std::string, uint64_t> seqCounts;
//				std::mutex seqMapLock;
//				std::stringstream nameStream;
//				nameStream << "["
//						<< "genome=" << genome << ";"
//						<< "mipTar=" << mipName << ";"
//						<< "mipFam=" << mip.familyName_ << ";"
//						<< "region=" << region << ";"
//						<< "extBarcode=" << extBarcode << ";"
//						<< "ligBarcode=" << ligBarcode << ";"
//						<< "libName=" << libInfo.libName_ << ";"
//						<< "simCount=" << count << ";"
//						<< "]";
//				std::string name = nameStream.str();
//				auto finalAmount = runPcr(intErrorRate, numThreads, initialPcrRounds, seq,
//						1, name, seqCounts, seqMapLock,
//						false);
//				uint64_t finalPerfectAmount = 1 * std::pow(2, initialPcrRounds);
//				templateNonMutated[name] = {finalAmount, finalPerfectAmount};
//				barcodedSeqs[name] = seq;
//				allSeqCounts[name] = seqCounts;
//				++count;
//			}
//		}
//	}
//
//	auto sampleNumber = sampleReadsWithoutReplacementFinishPCR(barcodedSeqs,
//			allSeqCounts, finalReadAmount, libOutFile, seqFileLock,
//			pcrRounds - initialPcrRounds, intErrorRate, numThreads, verbose);
//	if(verbose){
//		/*
//		{
//
//			std::cout << "PCR amounts: " << std::endl;
//			uint64_t nonMutated = 0;
//			uint64_t mutated = 0;
//			for (const auto & read : reads) {
//				std::cout << read->name_ << std::endl;
//				nonMutated += templateNonMutated[read->name_].first;
//				mutated += templateNonMutated[read->name_].second
//						- templateNonMutated[read->name_].first;
//				std::cout << "\t"
//						<< getPercentageString(templateNonMutated[read->name_].first,
//								templateNonMutated[read->name_].second) << std::endl;
//			}
//			std::cout << "Total: " << nonMutated + mutated << std::endl;
//			std::cout << "\t" << getPercentageString(nonMutated, nonMutated + mutated)
//					<< std::endl;
//			std::cout << "\t" << getPercentageString(mutated, nonMutated + mutated)
//					<< std::endl;
//
//		}
//		std::cout << "Sampling Amounts:" << std::endl;
//		uint64_t nonMutated = 0;
//		uint64_t mutated = 0;
//		for(const auto & read : reads){
//			std::cout << read->name_ << std::endl;
//			auto total = sampleNumber[read->name_].first + sampleNumber[read->name_].second;
//			nonMutated+= sampleNumber[read->name_].first;
//			mutated+= sampleNumber[read->name_].second;
//			std::cout << "\tSampled     : " << getPercentageString(total, finalReadAmountCounts[read->name_])<< std::endl;
//			std::cout << "\tNon-Mutated : " << getPercentageString(sampleNumber[read->name_].first, total) << std::endl;
//			std::cout << "\tMutated     : " << getPercentageString(sampleNumber[read->name_].second, total) << std::endl;
//		}
//		std::cout << "Total\t       : " << nonMutated + mutated << std::endl;
//		std::cout << "\tNon-Mutated : " << getPercentageString(nonMutated, nonMutated + mutated) << std::endl;
//		std::cout << "\tMutated     : " << getPercentageString(mutated, nonMutated + mutated) << std::endl;
//		 */
//	}
//	libOutFile.close();
//}
//
//} // namesapce sim
//
//
//int mipsterSimRunner::simMips(const njh::progutils::CmdArgs & inputCommands) {
//  seqSetUp setUp(inputCommands);
//
//
//  uint32_t startingTemplate = 3000;
//	double captureEfficiency = 0.10;
//	uint32_t finalReadAmount = 5000;
//	uint32_t pcrRounds = 20;
//	uint32_t initialPcrRounds = 10;
//	long double errorRate = 3.5e-06;
//	uint32_t numThreads = 2;
//
//	bfs::path genomeDir = "";
//	std::string genomesName = "";
//	std::string regionsName = "";
//	std::string mgvDirectory = "";
//	std::string abundanceFile = "";
//
//	bool sim454 = false;
//
//	uint32_t pairedEndReadLength = 250;
//	std::string mipFile = "";
//	std::string outDir = "";
//
//
//	setUp.setOption(captureEfficiency, "--captureEfficiency", "Efficiency of capture of the starting template, percent of starting template that actually gets captured");
//	setUp.setOption(genomeDir, "--genomeDir", "Genome Directory containing the genome files", true);
//	setUp.setOption(mgvDirectory, "--mgvDirectory", "Directory with Reference Files to simulate off of", true);
//	setUp.setOption(mipFile, "--mipArmsFilename", "Mip Arms File", true);
//	setUp.setOption(genomesName, "--genomes", "Genomes To Simulate", true);
//	setUp.setOption(regionsName, "--regions", "Regions To Simulate", true);
//	setUp.setOption(abundanceFile, "--abundanceFile", "Abundance File, first column Genome names, each additonal column is genome abundance", true);
//	setUp.setOption(pairedEndReadLength, "--pairedEndReadLength", "Paired End Read Length For illumina simulation");
//	setUp.setOption(sim454, "--sim454", "Simulate 454 as well as illumina");
//	setUp.setOption(startingTemplate, "--startingTemplate", "Staring Template Amount");
//	setUp.setOption(pcrRounds, "--pcrRounds", "Number of PCR rounds");
//	setUp.setOption(initialPcrRounds, "--initialPcrRounds", "Number of Initial rounds of PCR before sampling");
//	setUp.setOption(errorRate, "--errorRate", "Polymerase Error Rate");
//	setUp.setOption(numThreads, "--numThreads", "Number of Threads to Use");
//	setUp.setOption(finalReadAmount, "--finalReadAmount", "Final Read Amount to Sample Per Simulation Library");
//	setUp.processDirectoryOutputName("simMips_TODAY", true);
//	setUp.processVerbose();
//	setUp.processWritingOptions();
//	setUp.finishSetUp(std::cout);
//	setUp.startARunLog(setUp.pars_.directoryName_);
//	setUp.writeParametersFile(setUp.pars_.directoryName_ + "parameters.tab.txt", false, true);
//	MipCollection mCol(mipFile, 6);
//	njh::randomGenerator gen;
//	uint64_t intErrorRate = errorRate * std::numeric_limits<uint64_t>::max();
//
//	std::stringstream ss;
//	bool fail = false;
//	ss << __PRETTY_FUNCTION__ << ", error the following files need to exist" << "\n";
//	auto checkDir = [&fail,&ss](const bfs::path & fnp){
//		if(!bfs::exists(fnp)){
//			ss << fnp << "\n";
//			fail = true;
//		}
//	};
//	checkDir(abundanceFile);
//	checkDir(mgvDirectory);
//	checkDir(njh::files::make_path(mgvDirectory, "beds"));
//	checkDir(genomeDir);
//	if(fail){
//		throw std::runtime_error{ss.str()};
//	}
//
//
//	VecStr genomes = njh::tokenizeString(genomesName, ",");
//	VecStr regions = njh::tokenizeString(regionsName, ",");
//
//	VecStr libNames;
//	VecStr libDirNames;
//	auto libraryAbundances = processAbundanceLibaries(abundanceFile, genomes);
//	checkLibraryAbundancesThrow(libraryAbundances, startingTemplate);
//	//set starting template amounts
//	for (auto & lib : libraryAbundances) {
//		lib.second.setTemplateAmount(startingTemplate);
//	}
//	bfs::path pcrSimsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("pcrSims"));
//	bfs::path sequenceSimsDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("sequenceSims"));
//	bfs::path idFilesDir = njh::files::makeDir(setUp.pars_.directoryName_, njh::files::MkdirPar("id_files"));
//	bfs::copy_file(mipFile, njh::files::make_path(idFilesDir, bfs::path(mipFile).filename().string()));
//	MipsSamplesNames names(mCol.getMipFamsForRegions(regions), getVectorOfMapKeys(libraryAbundances));
//	std::ofstream sampleNamesFile;
//	openTextFile(sampleNamesFile, OutOptions(njh::files::make_path(idFilesDir, "allMipsSamplesNames.tab.txt")));
//	names.write(sampleNamesFile);
//	std::ofstream logFile;
//	Json::Value jsonLog;
//	openTextFile(logFile, OutOptions(bfs::path(setUp.pars_.directoryName_ + "simProgramLogs.json")));
//
//	std::unordered_map<std::string, std::unordered_map<std::string, seqInfo>> seqs;
//	MultiGenomeMapper gMapper(genomeDir, genomes.front());
//	gMapper.init();
//	bool failedREgionExtraction = false;
//	std::stringstream extractionMessage;
//	extractionMessage << __PRETTY_FUNCTION__ << ", found the following errors while extracting regions " << "\n";
//	//check for all bed files first
//	for(const auto & tar : names.mips_){
//		for(const auto & genome : genomes){
//			auto bedFnp = njh::files::make_path(mgvDirectory, "beds", genome + "_" + tar + ".bed");
//			if(!bfs::exists(bedFnp)){
//				extractionMessage << "No bed file for " << genome << " and " << tar << ", should found at " << bedFnp<< "\n";
//				failedREgionExtraction = true;
//			}
//		}
//	}
//	if(failedREgionExtraction){
//		throw std::runtime_error{extractionMessage.str()};
//	}
//	//check for one region only
//	for(const auto & tar : names.mips_){
//		for(const auto & genome : genomes){
//			auto bedFnp = njh::files::make_path(mgvDirectory, "beds", genome + "_" + tar + ".bed");
//			auto beds = getBed3s(bedFnp);
//			if(beds.size() != 1){
//				extractionMessage << "Was expecting only one region in " << bedFnp << " but found: " << beds.size() << "\n";
//				failedREgionExtraction = true;
//			}
//		}
//	}
//	if(failedREgionExtraction){
//		throw std::runtime_error{extractionMessage.str()};
//	}
//	for(const auto & tar : names.mips_){
//		for(const auto & genome : genomes){
//			auto bedFnp = njh::files::make_path(mgvDirectory, "beds", genome + "_" + tar + ".bed");
//			auto beds = getBeds(bedFnp);
//			auto targetRegion = GenomicRegion(*beds.front());
//			seqs[genome][tar] = gMapper.extractGenomeRegion(genome, targetRegion);
//		}
//	}
//	for (const auto & lib : libraryAbundances) {
//
//		sim::simMipLib(lib.second, mCol, regions, seqs,
//				pcrSimsDir.string(), captureEfficiency, intErrorRate,
//				finalReadAmount, pcrRounds, initialPcrRounds, numThreads,
//				setUp.pars_.verbose_);
//
//		//sim 454
//		if (sim454) {
//			std::string simCmd454 = "454sim -d $(echo $(dirname $(which 454sim))/gen) "
//					+ njh::files::make_path(pcrSimsDir, lib.second.libName_).string() + ".fasta -o " + njh::files::make_path(sequenceSimsDir, lib.second.libName_).string() + ".sff";
//			auto simOutPut454 = njh::sys::run(VecStr{simCmd454});
//			jsonLog[lib.first]["454SimLog"] = simOutPut454.toJson();
//		}
//		//sim illumina
//		uint32_t illuminaAttempts = 10; //art fails for no reason sometimes
//		std::string simCmdIllumina = "art_illumina -amp -p -na -i "
//				+ njh::files::make_path(pcrSimsDir, lib.second.libName_).string()+ ".fasta -l " + estd::to_string(pairedEndReadLength) + " -f 1 -o "
//				+ njh::files::make_path(sequenceSimsDir , lib.second.libName_).string() + "_R";
//		auto simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
//		jsonLog[lib.first]["454SimLog"] = simOutPutIllumina.toJson();
//		uint32_t numberOfAttempts = 1;
//		jsonLog[lib.first]["illumina_sim_log"]["attempt-" + estd::to_string(numberOfAttempts)] = simOutPutIllumina.toJson();
//		while(!simOutPutIllumina.success_ && numberOfAttempts <= illuminaAttempts){
//			++numberOfAttempts;
//			simOutPutIllumina = njh::sys::run(VecStr{simCmdIllumina});
//			jsonLog[lib.first]["illumina_sim_log"]["attempt-" + estd::to_string(numberOfAttempts)] = simOutPutIllumina.toJson();
//		}
//		//rename with proper fastq endings
//		bfs::rename(njh::files::make_path(sequenceSimsDir,  lib.second.libName_ ).string()+ "_R1.fq", njh::files::make_path(sequenceSimsDir , lib.second.libName_).string() + "_R1.fastq");
//		bfs::rename(njh::files::make_path(sequenceSimsDir, lib.second.libName_).string() + "_R2.fq", njh::files::make_path(sequenceSimsDir , lib.second.libName_ ).string()+ "_R2.fastq");
//		//gzip files
//		std::stringstream gzipCmd;
//		gzipCmd << "gzip " << njh::files::make_path(sequenceSimsDir , lib.second.libName_).string() + "_R1.fastq"
//				<< "&& gzip " << njh::files::make_path(sequenceSimsDir,  lib.second.libName_ ).string()+ "_R2.fastq";
//		auto gzipRunOuput = njh::sys::run({gzipCmd.str()});
//		jsonLog[lib.first]["illumina_sim_log"]["gzip"] = gzipRunOuput.toJson();
//	}
//	logFile << jsonLog << std::endl;
//
//	std::cout << "Done" << std::endl;
//	setUp.logRunTime(std::cout);
//	return 0;
//}



} // namespace njhseq
