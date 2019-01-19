/*
 * mipsterUtilsRunner_benchmarkingForControlMixtures.cpp
 *
 *  Created on: Jan 18, 2019
 *      Author: nicholashathaway
 */

#include "mipsterUtilsRunner.hpp"
#include <SeekDeep/utils.h>


namespace njhseq {

class ControlMixSetUp {
public:
	ControlMixSetUp(const std::string & name,
			const std::unordered_map<std::string, double> & relativeAbundances) :
			name_(name), rawRelativeAbundances_(relativeAbundances) {
		if (relativeAbundances.empty()) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error "
					<< " relativeAbundances can't be empty" << "\n";
			throw std::runtime_error { ss.str() };
		}

		double total = 0;
		for (const auto & relAbund : rawRelativeAbundances_) {

			if (relAbund.second <= 0) {
				std::stringstream ss;
				ss << __PRETTY_FUNCTION__ << ", error "
						<< " abundances need to be greater than zero, " << relAbund.first
						<< ": " << relAbund.second << "\n";
				throw std::runtime_error { ss.str() };
			}
			total += relAbund.second;
		}
		for (const auto & relAbund : rawRelativeAbundances_) {
			relativeAbundances_[relAbund.first] = relAbund.second/total;
		}

	}


	std::string name_;
	std::unordered_map<std::string, double> rawRelativeAbundances_;
	std::unordered_map<std::string, double> relativeAbundances_;

	MetaDataInName meta_;

	VecStr getStrains() const {
		return njh::getVecOfMapKeys(rawRelativeAbundances_);
	}

};


int mipsterUtilsRunner::benchmarkingForControlMixtures(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path extractionDir = "";
	bfs::path masterAnalysisDir = "";
	bfs::path sampleToMixture = "";
	bfs::path mixtureSetUp = "";

	mipsterUtilsSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(extractionDir,     "--extractionDir",     "Extraction Directory, extractions for MIP arms created by MIPWrangler setUpViewMipsOnGenome", true);
	setUp.setOption(masterAnalysisDir, "--masterAnalysisDir", "Master Analysis Directory, finished analysis of MIPWrangler", true);
	setUp.setOption(sampleToMixture,   "--sampleToMixture",   "Sample To Mixture, 2 columns 1)sample, 2)MixName", true);
	setUp.setOption(mixtureSetUp,      "--mixtureSetUp",      "Mixture Set Up, 3 columns 1)MixName, 2)strain, 3)relative_abundance", true);
	std::string defaultDirName = njh::rstripRet(masterAnalysisDir.string(), '/') + "_benchmarkingForControlMixtures_TODAY";
	setUp.processDirectoryOutputName(defaultDirName, true);
	setUp.finishSetUp(std::cout);
	setUp.startARunLog(setUp.pars_.directoryName_);
	SetUpMaster mipMaster(masterAnalysisDir);
	mipMaster.loadMipsSampsInfo(6);
	mipMaster.checkDirStructThrow(__PRETTY_FUNCTION__);


	//read in mixture setup
	table mixtureSetupTab(mixtureSetUp, "\t", true);
	mixtureSetupTab.checkForColumnsThrow(VecStr{"MixName", "strain", "relative_abundance"}, __PRETTY_FUNCTION__);

	if(mixtureSetupTab.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << mixtureSetUp << " is empty " << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, ControlMixSetUp> mixSetups;
	std::unordered_map<std::string, std::unordered_map<std::string, double>> mixInfos;
	for(const auto & row : mixtureSetupTab){
		std::string mixname = row[mixtureSetupTab.getColPos("MixName")];
		std::string strain = row[mixtureSetupTab.getColPos("strain")];
		double relative_abundance  = njh::StrToNumConverter::stoToNum<double>(row[mixtureSetupTab.getColPos("relative_abundance")]);
		if(njh::in(strain, mixInfos[mixname])){
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ", error " << " already have " << strain << " for mixture " << mixname<< "\n";
			throw std::runtime_error{ss.str()};
		}
		mixInfos[mixname][strain] = relative_abundance;
	}
	for(const auto & mixInfo : mixInfos){
		mixSetups.emplace(mixInfo.first, ControlMixSetUp(mixInfo.first, mixInfo.second));
	}

	//read in samples to mixture
	table sampleToMixtureTab(sampleToMixture, "\t", true);
	sampleToMixtureTab.checkForColumnsThrow(VecStr{"sample", "MixName"}, __PRETTY_FUNCTION__);
	if(sampleToMixtureTab.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << sampleToMixture << " is empty " << "\n";
		throw std::runtime_error{ss.str()};
	}
	std::unordered_map<std::string, std::string> samplesToMix;
	for(const auto & row : sampleToMixtureTab){
		samplesToMix[row[sampleToMixtureTab.getColPos("sample")]] = row[sampleToMixtureTab.getColPos("MixName")];
	}
	VecStr missingMixs;
	for(const auto & sampToMix : samplesToMix){
		if(!njh::in(sampToMix.second, mixSetups)){
			missingMixs.emplace_back(sampToMix.second);
		}
	}
	if(!missingMixs.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " error missing the following mixture information from " << mixtureSetUp << "\n";
		throw std::runtime_error{ss.str()};
	}

	//read in extraction table and determine which mips can bed used per mixture
	std::set<std::string> allStrains;
	for(const auto & mixSetup : mixSetups){
		auto strains = mixSetup.second.getStrains();
		allStrains.insert(strains.begin(), strains.end());
	}

	bfs::path extractionTableFnp = njh::files::make_path(extractionDir, "tables", "extractionCountsTable.tab.txt");
	table extractionTable(extractionTableFnp, "\t", true);
	VecStr missingStrains;
	for(const auto & strain : allStrains){
		if(!njh::in(strain, extractionTable.columnNames_)){
			missingStrains.emplace_back(strain);
		}
	}
	if(!missingStrains.empty()){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error " << " the following strains are missing " << njh::conToStrEndSpecial(missingStrains, ", ", " and ") << " from " << extractionTableFnp<< "\n";
		throw std::runtime_error{ss.str()};
	}
	bfs::copy_file(sampleToMixture, njh::files::make_path(setUp.pars_.directoryName_, "sampleToMixture.tab.txt"));
	bfs::copy_file(mixtureSetUp, njh::files::make_path(setUp.pars_.directoryName_, "mixtureSetUp.tab.txt"));

	std::unordered_map<std::string, VecStr> targetsForMix;
	{
		OutputStream unusedTargetsOut(njh::files::make_path(setUp.pars_.directoryName_, "extractionCountsForUnusedTargets.tab.txt"));
		unusedTargetsOut << njh::conToStr(extractionTable.columnNames_, "\t") << std::endl;
		for(const auto & row : extractionTable){
			auto tarName = row[extractionTable.getColPos("target")];
			for(const auto & mixSetup : mixSetups){
				bool failed = false;
				for(const auto & strain : mixSetup.second.relativeAbundances_){
					if("1" != row[extractionTable.getColPos(strain.first)]){
						failed = true;
						break;
					}
				}
				if(failed){
					unusedTargetsOut << njh::conToStr(row, "\t") << std::endl;
				}else{
					targetsForMix[mixSetup.first].emplace_back(tarName);
				}
			}
		}
	}


	std::unordered_map<std::string, std::set<std::string>> familiesForMix;
	for(const auto & tarsForMix : targetsForMix){
		for(const auto & tar : tarsForMix.second){
			familiesForMix[tarsForMix.first].emplace(mipMaster.mips_->getFamilyForTarget(tar));
		}
	}
	std::set<std::string> allMipFams;
	for(const auto & famForMix : familiesForMix){
		allMipFams.insert(famForMix.second.begin(), famForMix.second.end());
	}

	{
		OutputStream familiesForMixOut(njh::files::make_path(setUp.pars_.directoryName_, "familiesForMix.json"));
		familiesForMixOut << njh::json::toJson(familiesForMix) << std::endl;
	}
	{
		OutputStream targetsForMixOut(njh::files::make_path(setUp.pars_.directoryName_, "targetsForMix.json"));
		targetsForMixOut << njh::json::toJson(targetsForMix) << std::endl;
	}
	std::unordered_map<std::string, std::unordered_map<std::string, uint32_t>> barcodeTotals;
	OutputStream haplotypesClassified(njh::files::make_path(setUp.pars_.directoryName_, "classifiedHaplotypes.tab.txt"));
	haplotypesClassified << "sample\tmix\tmipFamily\tseqName\treadCnt\tbarcodeCnt\tbarcodeFrac\tmatchExpcted\texpectedRef\texpectedFrac" << std::endl;
	OutputStream performanceOut(njh::files::make_path(setUp.pars_.directoryName_, "performancePerTarget.tab.txt"));
	performanceOut << "sample\tmix\tmipFamily\ttotalBarcodes\trecoveredHaps\tfalseHaps\ttotalHaps\ttotalExpectedHaps\thapRecovery\tfalseHapsRate\tRMSE" << std::endl;
	for(const auto & sample : samplesToMix){
		for(const auto & mipFamily : familiesForMix[sample.second]){
			bfs::path resultsFnp = mipMaster.pathPopClusFinalHaplo(MipFamSamp(mipFamily, sample.first));
			if(bfs::exists(resultsFnp)){
				bfs::path expectedSeqeuncesFnp = njh::files::make_path(extractionDir, "fastas", "byFamily", mipFamily + ".fasta");
				SeqInput expReader(SeqIOOptions::genFastaIn(expectedSeqeuncesFnp));
				auto allExpSeqs = expReader.readAllReadsPtrs<seqInfo>();
				std::vector<std::shared_ptr<seqInfo>> expSeqs;
				std::unordered_map<std::string, std::string> nameKey;
				std::unordered_map<std::string, double> expectedFracs;
				for(const auto & seq : allExpSeqs){
					auto toks = tokenizeString(seq->name_, ",");
					std::vector<std::string> expecteds;
					double expectedFrac = 0;
					for(const auto & tok : toks){
						if(njh::in(tok, mixSetups.at(sample.second).relativeAbundances_)){
							expectedFrac+= mixSetups.at(sample.second).relativeAbundances_[tok];
							expecteds.emplace_back(tok);
						}
					}
					if(!expecteds.empty()){
						njh::sort(expecteds);
						expSeqs.push_back(seq);
						expectedFracs[seq->name_] = expectedFrac;
						nameKey[seq->name_] = njh::conToStr(expecteds, ",");
					}
				}
				auto resSeqs = SeqInput::getSeqVec<seqInfo>(SeqIOOptions::genFastqIn(resultsFnp));
				double barCodeTotal = 0;
				for (const auto & resSeq : resSeqs) {
					MetaDataInName seqMeta(resSeq.name_);
					barCodeTotal += seqMeta.getMeta<double>("barcodeCnt");
				}
				uint32_t recoveredHaps = 0;
				uint32_t falseHaps = 0;
				double sumOfSquares = 0;

				for(const auto & resSeq : resSeqs){
					MetaDataInName seqMeta(resSeq.name_);
					std::string matchingRef = "";
					double expectedFrac = 0;
					for(const auto & expSeq : expSeqs){
						if(resSeq.seq_ == expSeq->seq_){
							matchingRef = nameKey[expSeq->name_];
							expectedFrac = expectedFracs[expSeq->name_];
							break;
						}
					}
					if("" != matchingRef){
						++recoveredHaps;
						sumOfSquares += std::pow((seqMeta.getMeta<double>("barcodeCnt")/barCodeTotal) - expectedFrac, 2.0);
					}else{
						++falseHaps;
					}
					haplotypesClassified << sample.first
							<< "\t" << sample.second
							<< "\t" << mipFamily
							<< "\t" << resSeq.name_
							<< "\t" << seqMeta.getMeta("readCnt")
							<< "\t" << seqMeta.getMeta("barcodeCnt")
							<< "\t" << seqMeta.getMeta<double>("barcodeCnt")/barCodeTotal
							<< "\t" << ("" == matchingRef ? "FALSE" : "TRUE")
							<< "\t" << matchingRef
							<< "\t" << expectedFrac
							<< std::endl;
				}
				performanceOut << sample.first
						<< "\t" << sample.second
						<< "\t" << mipFamily
						<< "\t" << barCodeTotal
						<< "\t" << recoveredHaps
						<< "\t" << falseHaps
						<< '\t' << recoveredHaps + falseHaps
						<< "\t" << expSeqs.size()
						<< "\t" << static_cast<double>(recoveredHaps)/expSeqs.size()
						<< "\t" << static_cast<double>(falseHaps)/(recoveredHaps + falseHaps)
						<< "\t" << (0 == sumOfSquares ? 0 : std::sqrt(sumOfSquares)) << std::endl;

				barcodeTotals[sample.first][mipFamily] = barCodeTotal;
			}else{
				barcodeTotals[sample.first][mipFamily] = 0;
			}
		}
	}

	OutputStream barcodeTotalsOut(njh::files::make_path(setUp.pars_.directoryName_, "barcodeTotals.tab.txt"));
	auto samples = getVectorOfMapKeys(samplesToMix);
	njh::sort(samples);
	std::unordered_map<std::string, std::vector<uint32_t>> barcodeCounts;
	barcodeTotalsOut << "Mip\t" << njh::conToStr(samples, "\t") << std::endl;
	for(const auto & mipFam : allMipFams){
		barcodeTotalsOut << mipFam;
		for(const auto & samp : samples){
			barcodeTotalsOut << '\t' << barcodeTotals[samp][mipFam];
			if(barcodeTotals[samp][mipFam] > 0){
				barcodeCounts[samp].emplace_back(barcodeTotals[samp][mipFam]);
			}
		}
		barcodeTotalsOut << std::endl;
	}

	OutputStream sampleStatsOut(njh::files::make_path(setUp.pars_.directoryName_, "sampBarocdeStats.tab.txt"));
	sampleStatsOut << "Sample\tMixName\ttargetsWithBarcodes\ttotalTargets\ttotalBarcodes\tmedianBarcodeCnt" << std::endl;
	for(const auto & samp : samples){
		sampleStatsOut << samp
				<< '\t' << samplesToMix[samp]
				<< '\t' << barcodeCounts[samp].size()
				<< '\t' << familiesForMix[samplesToMix[samp]].size()
				<< "\t" << vectorSum(barcodeCounts[samp])
				<< "\t" << vectorMedianRef(barcodeCounts[samp]) << std::endl;
	}

	return 0;

}


} //namespace njhseq
