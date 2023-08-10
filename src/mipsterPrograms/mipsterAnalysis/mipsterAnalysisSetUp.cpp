/*
 * mipsterAnalysisSetUp.cpp
 *
 *  Created on: Jan 2, 2015
 *      Author: nickhathaway
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
#include "mipsterAnalysisSetUp.hpp"

namespace njhseq {

void mipsterAnalysisSetUp::processDefaults(mipCorePars & pars){
	processVerbose();
	processDebug();
	pars.processDefaults(*this);
}

void mipsterAnalysisSetUp::processMultipleDefaults(mipCorePars & pars){
	processDefaults(pars);
	setOption(pars.numThreads, "--numThreads", "Number of threads to use");
	if(0 == pars.numThreads){
		failed_ = true;
		addWarning("Number of threads can't be zero");
	}
	setOption(pars.logFilename, "--logFile", "Name of a file to log information about the run", pars.logFileRequired);
	setOption(pars.overWriteLog, "--overWriteLog", "Overwrite the log file");
}


void mipsterAnalysisSetUp::setUpExtractFromRawMultiple(extractFromRawParsMultiple & pars){

	pars.infoFilesRequired = true;
	pars.logFileRequired = false;
	pars.logFilename = "extractFromRawLog";
	processMultipleDefaults(pars);
	setOption(pars.dir, "--dir", "Name of raw data directory", true);

	pars.qFilPars_.checkingQFrac_ = true;
	setOption(pars.qFilPars_.qualWindow_, "--qualWindow",
				"Sliding Quality Window, format is WindowSize,WindowStep,Threshold");
	seqUtil::processQualityWindowString(pars.qFilPars_.qualWindow_, pars.qFilPars_.qualityWindowLength_,
			pars.qFilPars_.qualityWindowStep_, pars.qFilPars_.qualityWindowThres_);
	setOption(pars.qFilPars_.qualCheck_, "--qualCheck", "Qual Check Level");
	setOption(pars.qFilPars_.qualCheckCutOff_, "--qualCheckCutOff",
			"Cut Off for fraction of bases above qual check of "
					+ estd::to_string(pars.qFilPars_.qualCheck_));
	if (commands_.hasFlagCaseInsenNoDash("-qualWindow")) {
		pars.qFilPars_.checkingQWindow = true;
		pars.qFilPars_.checkingQFrac_ = false;
	}
	pars.qFilPars_.trimAtQual_ = setOption(pars.qFilPars_.trimAtQualCutOff_, "--trimAtQual",
			"Trim Reads at first occurrence of quality score");
	setOption(pars.allowableErrors, "--allowableErrors",
			"Number of Errors To Allow in Arm");
	setOption(pars.smallFragmentLength, "--smallFragmentLength",
			"Length to consider a read to be small fragment and shouldn't be processed");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms including barcode to be from the front of the read");
	setOption(pars.minCaptureLength, "--minCaptureLength", "Minimum Capture Length cut off, captures below this will be thrown out, can be set for each mip individually in mip table");
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars.fileOpenLimit_, "--fileOpenLimit", "Number of file allowed to open by one process");
	setOption(pars.refDir, "--refDir", "Directory with possible reference sequence to rename popUIDs to");

	setOption(pars.processPairPars_.r1Trim_, "--r1Trim", "Number of Bases to trim off at the end of the r1 read");
	setOption(pars.processPairPars_.r2Trim_, "--r2Trim", "Number of Bases to trim off at the end of the r2 read");

  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");
  setOption(pars.writeOutInitialExtractedPairs, "--writeOutInitialExtractedPairs", "Write Out Initial Extracted Pairs");


	pars.fileOpenLimit_ = (pars.fileOpenLimit_ - pars.numThreads) /pars.numThreads;
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();


	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpRunGzExtractStitch(
		runGzExtractStitchPars & pars) {

	pars.infoFilesRequired = true;
	pars.logFileRequired = false;
	pars.logFilename = "gzStitchLog";
	processMultipleDefaults(pars);
	setOption(pars.trim, "--trim",
			"Whether or not to trim with sickle before stitching");
	setOption(pars.minOverlap, "--minOverlap",
			"The minimum overlap to be used with the stitcher");
	setOption(pars.maxOverlap, "--maxOverlap",
			"The minimum overlap to be used with the stitcher");
	setOption(pars.dir, "--dir", "Name of raw data directory", true);
	setOption(pars.usePear, "--usePear",
			"use pear to stitch, default is to use flash");
	setOption(pars.usePanda, "--usePanda",
			"use pandaseq default to stitch, default is to use flash");
	setOption(pars.mismatchDensity, "--mismatchDensity",
			"The percentage of mismatches allowed in the overlap");
	finishSetUp(std::cout);
}


void mipsterAnalysisSetUp::setUpMipIllumArmExtractionPaired(mipIllumArmExtractionPars & pars){
	processDefaults(pars);
	pars.qFilPars_.checkingQFrac_ = true;
	processQualityFiltering();
	if (commands_.hasFlagCaseInsenNoDash("-qualWindow")) {
		pars.qFilPars_.checkingQWindow = true;
		pars.qFilPars_.checkingQFrac_ = false;
	}
	setOption(pars.qFilPars_.qualWindow_, "--qualWindow",
				"Sliding Quality Window, format is WindowSize,WindowStep,Threshold");
	seqUtil::processQualityWindowString(pars.qFilPars_.qualWindow_, pars.qFilPars_.qualityWindowLength_,
			pars.qFilPars_.qualityWindowStep_, pars.qFilPars_.qualityWindowThres_);
	setOption(pars.qFilPars_.qualCheck_, "--qualCheck", "Qual Check Level");
	setOption(pars.qFilPars_.qualCheckCutOff_, "--qualCheckCutOff",
			"Cut Off for fraction of bases above qual check of "
					+ estd::to_string(pars.qFilPars_.qualCheck_));

	pars.qFilPars_.trimAtQual_ = setOption(pars.qFilPars_.trimAtQualCutOff_, "--trimAtQual",
			"Trim Reads at first occurrence of quality score");
	setOption(pars.sampleName, "--sample", "Name of the sample/sample dir to cluster", true);
	setOption(pars.allowableErrors, "--allowableErrors",
			"Number of Errors To Allow in Arm");
	setOption(pars.smallFragmentLength, "--smallFragmentLength",
			"Length to consider a read to be small fragment and shouldn't be processed");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms to be from the front of the read");
	setOption(pars.minCaptureLength, "--minCaptureLength", "Minimum Capture Length cut off, captures below this will be thrown out, can be set for each mip individually in mip table");
	processReadInNames(VecStr{"--fastq1", "--fastq2"});
	pars_.ioOptions_.revComplMate_ = true;
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipIllumArmExtractionMultiplePaired(
		mipIllumArmExtractionParsMultiple & pars){
	pars.logFilename = "mipArmExtractLog";
	pars.logFileRequired = false;
	processMultipleDefaults(pars);
	pars.qFilPars_.checkingQFrac_ = true;
	processQualityFiltering();
	if (commands_.hasFlagCaseInsenNoDash("-qualWindow")) {
		pars.qFilPars_.checkingQWindow = true;
		pars.qFilPars_.checkingQFrac_ = false;
	}
	pars.qFilPars_.trimAtQual_ = setOption(pars.qFilPars_.trimAtQualCutOff_, "--trimAtQual",
			"Trim Reads at first occurrence of quality score");
	setOption(pars.allowableErrors, "--allowableErrors",
			"Number of Errors To Allow in Arm");
	setOption(pars.smallFragmentLength, "--smallFragmentLength",
			"Length to consider a read to be small fragment and shouldn't be processed");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms including barcode to be from the front of the read");
	setOption(pars.minCaptureLength, "--minCaptureLength", "Minimum Capture Length cut off, captures below this will be thrown out, can be set for each mip individually in mip table");
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars.fileOpenLimit_, "--fileOpenLimit", "Number of file allowed to open by one process");
	pars.fileOpenLimit_ = (pars.fileOpenLimit_ - pars.numThreads) /pars.numThreads;
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipIllumArmExtraction(
		mipIllumArmExtractionPars & pars) {
	processDefaults(pars);
	pars_.qFilPars_.checkingQFrac_ = true;
	processQualityFiltering();
	if (commands_.hasFlagCaseInsenNoDash("-qualWindow")) {
		pars_.qFilPars_.checkingQWindow = true;
		pars_.qFilPars_.checkingQFrac_ = false;
	}
	setOption(pars.sampleName, "--sample", "Name of the sample/sample dir to cluster", true);
	setOption(pars.allowableErrors, "--allowableErrors",
			"Number of Errors To Allow in Arm");
	setOption(pars.smallFragmentLength, "--smallFragmentLength",
			"Length to consider a read to be small fragment and shouldn't be processed");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms to be from the front of the read");
	setOption(pars.minCaptureLength, "--minCaptureLength", "Minimum Capture Length cut off, captures below this will be thrown out, can be set for each mip individually in mip table");
	processReadInNames(VecStr{"--fastq"});
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipIllumArmExtractionMultiple(
		mipIllumArmExtractionParsMultiple & pars) {
	pars.logFilename = "mipArmExtractLog";
	pars.logFileRequired = false;
	processMultipleDefaults(pars);
	pars_.qFilPars_.checkingQFrac_ = true;
	processQualityFiltering();
	if (commands_.hasFlagCaseInsenNoDash("-qualWindow")) {
		pars_.qFilPars_.checkingQWindow = true;
		pars_.qFilPars_.checkingQFrac_ = false;
	}
	setOption(pars.allowableErrors, "--allowableErrors",
			"Number of Errors To Allow in Arm");
	setOption(pars.smallFragmentLength, "--smallFragmentLength",
			"Length to consider a read to be small fragment and shouldn't be processed");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms including barcode to be from the front of the read");
	setOption(pars.minCaptureLength, "--minCaptureLength", "Minimum Capture Length cut off, captures below this will be thrown out, can be set for each mip individually in mip table");
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars.fileOpenLimit_, "--fileOpenLimit", "Number of file allowed to open by one process");
	pars.fileOpenLimit_ = (pars.fileOpenLimit_ - pars.numThreads) /pars.numThreads;
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();
	finishSetUp(std::cout);
}
void mipsterAnalysisSetUp::setUpMipBarcodeCorrection(mipBarcodeCorrectionPars & pars){
	processDefaults(pars);
	setOption(pars.sampleName, "--sample", "Name of the sample/sample dir to cluster", true);
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	setOption(pars.allowableErrors, "--allowableErrors", "Number of errors to allow in arms,"
			" should be the same as was set for extraction");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms to be from the front of the read, should be the same as was set for extraction");
	setOption(pars.barcodeIdentity, "--barcodeIdentity", "The amount of identity between reads to allow for same barcode clustering");
	setOption(pars.qualRep, "--qualCalculation", "How to calculate the per base quality scores");
	processAlignerDefualts();

  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");
	setOption(pars.doNotDownSample_, "--doNotDownSample", "do Not Down Sample");
	setOption(pars.downSampleAmount_, "--downSampleAmount", "downsample read Amount");
	std::random_device rd;
	pars.seed_ = rd();
	setOption(pars.seed_, "--downSampleSeed", "down sample seed");

	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipBarcodeCorrectionMultiple(mipBarcodeCorrectionParsMultiple & pars){
	pars.logFilename = "mipBarCor";
	pars.logFileRequired = false;
	processMultipleDefaults(pars);
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	setOption(pars.allowableErrors, "--allowableErrors", "Number of errors to allow in arms,"
			" should be the same as was set for extraction");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms to be from the front of the read, should be the same as was set for extraction");
	setOption(pars.barcodeIdentity, "--barcodeIdentity", "The amount of identity between reads to allow for same barcode clustering");
	setOption(pars.qualRep, "--qualCalculation", "How to calculate the per base quality scores");
	processAlignerDefualts();

  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");
	setOption(pars.doNotDownSample_, "--doNotDownSample", "do Not Down Sample");
	setOption(pars.downSampleAmount_, "--downSampleAmount", "downsample read Amount");
	std::random_device rd;
	pars.seed_ = rd();
	setOption(pars.seed_, "--downSampleSeed", "down sample seed");


	finishSetUp(std::cout);
	//getOneMipPopSeqsPostHandler
}

void mipsterAnalysisSetUp::setUpMipClustering(mipClusteringPars & pars) {
	processDefaults(pars);
	setOption(pars.sampleName, "--sample", "Name of the sample/sample dir to cluster", true);
	setOption(pars.parameterFile, "--par", "Name of the parameters file for the clustering");
	setOption(pars.useReadLen, "--useReadLen", "Use Read Length large differences to skip alignment comparison");
	setOption(pars.writeOutClusters, "--writeOutClusters", "Write out the sequences that make up the final clusters");
	setOption(pars.cacheAlignments,  "--cacheAlignments", "Cache alignments so if the analysis is re-ran it will be much quicker");
	setOption(pars.qualRep, "--qualCalculation", "How to calculate the per base quality scores");
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "5,1";
	pars_.gapLeft_ = "5,1";
	pars_.qualThres_ = "25,20";
	pars_.qScorePars_.primaryQual_ = 25;
	pars_.qScorePars_.secondaryQual_ = 20;
	pars_.qScorePars_.qualThresWindow_ = 0;
	processAlignerDefualts();
	if(!failed_){
		if ("" == pars.parameterFile) {
			pars.iterMap = CollapseIterations::genIlluminaDefaultPars(100);
		} else {
			pars.iterMap = processIteratorMap(pars.parameterFile);
		}
	}
  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");


	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipClusteringMultiple(mipClusteringParsMultiple & pars) {
	pars.logFilename = "mipClustering";
	pars.logFileRequired = false;
	processMultipleDefaults(pars);
	setOption(pars.parameterFile, "--par", "Name of the parameters file for the clustering");
	setOption(pars.useReadLen, "--useReadLen", "Use Read Length large differences to skip alignment comparison");
	setOption(pars.qualRep, "--qualCalculation", "How to calculate the per base quality scores");
	setOption(pars.writeOutClusters, "--writeOutClusters", "Write out the sequences that make up the final clusters");
	setOption(pars.cacheAlignments,  "--cacheAlignments", "Cache alignments so if the analysis is re-ran it will be much quicker");
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "5,1";
	pars_.gapLeft_ = "5,1";
	pars_.qualThres_ = "25,20";
	pars_.qScorePars_.primaryQual_ = 25;
	pars_.qScorePars_.secondaryQual_ = 20;
	pars_.qScorePars_.qualThresWindow_ = 0;
	processAlignerDefualts();
	if(!failed_){
		if ("" == pars.parameterFile) {
			pars.iterMap = CollapseIterations::genIlluminaDefaultPars(100);
		} else {
			pars.iterMap = processIteratorMap(pars.parameterFile);
		}
	}

  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");


	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipPopulationClustering(
		mipPopulationClusteringPars & pars) {
	processDefaults(pars);
  setOption(pars.parameters, "--par", "ParametersFileName");
  setOption(pars.mipName, "--mipName", "Mip Name", true);
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars_.colOpts_.skipOpts_.skipOnLetterCounterDifference_, "-skip",
			"skipOnLetterCounterDifference");
	setOption(pars_.colOpts_.skipOpts_.fractionDifferenceCutOff_, "-skipCutOff",
			"fractionDifferenceCutOff");
	setOption(pars.previousPopDir, "--refDir", "Directory with possible reference sequence to rename popUIDs to");
	auto refPopFile = njh::files::make_path(pars.previousPopDir, pars.mipName + ".fasta");
	if(bfs::exists(refPopFile)){
		pars.previousPopFilename = refPopFile.string();
		pars.refIoOptions = SeqIOOptions::genFastaIn(refPopFile);
	}

	processComparison(pars.previousPopErrors, "previousPop");
  setOption(pars.runsRequired, "-runsRequired", "runsRequired");
  setOption(pars.clusteringCutOffs.clusterSizeCutOff, "--minimumBarcodeCutOff,-cutoff", "The Minimum molecular Barcode count Cut Off to be inlcuded in final analysis");
  setOption(pars.fracCutoff, "-fraccutoff", "PopulationClusteringFractionCutoff");
  pars.differentPar = setOption(pars.parametersPopulation, "--popPar",
      "ParametersForPopulationCollapse");
  setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
  setOption(pars.keepChimeras, "-keepChimeras", "KeepChimeras");
  setOption(pars_.chiOpts_.parentFreqs_, "-parfreqs", "ParentFrequence_multiplier_cutoff");


  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");


  pars_.ioOptions_.lowerCaseBases_ = "upper";
  pars_.directoryName_ = "analysis";
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "5,1";
	pars_.gapLeft_ = "5,1";
  processGap();
  processQualThres();
  processScoringPars();
  processRefFilename();
  // get the qualities
  if (!failed_) {
		if ("" == pars.parameters) {
			pars.iteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(300);
		} else {
			pars.iteratorMap = processIteratorMap(pars.parameters);
		}
    if(pars.differentPar){
    	pars.popIteratorMap = processIteratorMap(pars.parametersPopulation);
    }else{
    	pars.popIteratorMap = pars.iteratorMap;
    }
    if(pars_.verbose_){
      std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
      std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
      std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
      std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
    }
  }
  // read in the paramteres from the parameters file
  finishSetUp(std::cout);

}

void mipsterAnalysisSetUp::setUpMipPopulationClusteringMultiple(
		mipPopulationClusteringParsMultiple & pars) {
	pars.logFilename = "mipPopClus";
	pars.logFileRequired = false;
	processMultipleDefaults(pars);
  setOption(pars.parameters, "-par", "ParametersFileName");
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars_.colOpts_.skipOpts_.skipOnLetterCounterDifference_, "-skip",
			"skipOnLetterCounterDifference");
	setOption(pars_.colOpts_.skipOpts_.fractionDifferenceCutOff_, "-skipCutOff",
			"fractionDifferenceCutOff");

	setOption(pars.previousPopDir, "-refDir", "Directory with possible reference sequence to rename popUIDs to");

	processComparison(pars.previousPopErrors, "previousPop");
  setOption(pars.runsRequired, "-runsRequired", "runsRequired");
  setOption(pars.clusteringCutOffs.clusterSizeCutOff, "--minimumBarcodeCutOff,-cutoff", "The Minimum molecular Barcode count Cut Off to be inlcuded in final analysis");


  setOption(pars.fracCutoff, "-fraccutoff", "PopulationClusteringFractionCutoff");
  pars.differentPar = setOption(pars.parametersPopulation, "-poppar",
      "ParametersForPopulationCollapse");
  setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
  setOption(pars.keepChimeras, "-keepChimeras", "KeepChimeras");
  setOption(pars_.chiOpts_.parentFreqs_, "-parfreqs", "ParentFrequence_multiplier_cutoff");

  setOption(pars.keepIntermediateFiles, "--keepIntermediateFiles", "Keep Intermediate Files");



  pars_.ioOptions_.lowerCaseBases_ = "upper";
  pars_.directoryName_ = "analysis";
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "5,1";
	pars_.gapLeft_ = "5,1";
  processGap();
  processQualThres();
  processScoringPars();
  processRefFilename();
  // get the qualities
  if (!failed_) {
		if ("" == pars.parameters) {
			pars.iteratorMap = CollapseIterations::genStrictNoErrorsDefaultPars(300);
		} else {
			pars.iteratorMap = processIteratorMap(pars.parameters);
		}
    if(pars.differentPar){
    	pars.popIteratorMap = processIteratorMap(pars.parametersPopulation);
    }else{
    	pars.popIteratorMap = pars.iteratorMap;
    }
    if(pars_.verbose_){
      std::cout << "p: " << pars_.qScorePars_.primaryQual_ << std::endl;
      std::cout << "s: " << pars_.qScorePars_.secondaryQual_ << std::endl;
      std::cout << "go: " << pars_.gapInfo_.gapOpen_ << std::endl;
      std::cout << "ge: " << pars_.gapInfo_.gapExtend_ << std::endl;
    }
  }
  // read in the paramteres from the parameters file
  finishSetUp(std::cout);
}


void mipsterAnalysisSetUp::setUpMipCorrectForContamWithSameBarcodes(mipCorrectForContamWithSameBarcodesPars & pars){
	pars.processDefaults(*this);
	setOption(pars.readCutOff, "--readCutOff",
			"Remove reads that share another barcode with another sample with read coverage of this or below");
	setOption(pars.mipName, "--mipName",
			"Name of mip to correct for contamination", true);
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipCorrectForContamWithSameBarcodesMultiple(mipCorrectForContamWithSameBarcodesParsMultiple & pars){
	processMultipleDefaults(pars);
	setOption(pars.readCutOff, "--readCutOff",
			"Remove reads that share another barcode with another sample with read coverage of this or below");
	finishSetUp(std::cout);
}


} /* namespace njhseq */
