/*
 * mipsterAnalysisSetUp.cpp
 *
 *  Created on: Jan 2, 2015
 *      Author: nickhathaway
 */

#include "mipsterAnalysisSetUp.hpp"

namespace bibseq {

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
	setOption(pars.logFilename, "--logFile", "Name of a file to log information about the run",true);
	setOption(pars.overWriteLog, "--overWriteLog", "Overwrite the log file");
}


void mipsterAnalysisSetUp::setUpRunGzExtractStitch(
		runGzExtractStitchPars & pars) {
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
	setOption(pars.minLen, "--minLen", "Minimum Length of Read to Be Extracted");
	processDefaultReader(true);
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipIllumArmExtractionMultiple(
		mipIllumArmExtractionParsMultiple & pars) {
	pars.logFilename = "mipArmExtractLog";
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
	setOption(pars.minLen, "--minLen", "Minimum Length of Read to Be Extracted");
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
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipBarcodeCorrectionMultiple(mipBarcodeCorrectionParsMultiple & pars){
	pars.logFilename = "mipBarCor";
	processMultipleDefaults(pars);
	pars_.gap_ = "5,1";
	pars_.gapRight_ = "0,0";
	pars_.gapLeft_ = "0,0";
	setOption(pars.allowableErrors, "--allowableErrors", "Number of errors to allow in arms,"
			" should be the same as was set for extraction");
	setOption(pars.wiggleRoom, "--wiggleRoom",
			"Amount of bases to allow the arms to be from the front of the read, should be the same as was set for extraction");
	setOption(pars.barcodeIdentity, "--barcodeIdentity", "The amount of identity between reads to allow for same barcode clustering");
	processAlignerDefualts();
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipClustering(mipClusteringPars & pars) {
	processDefaults(pars);
	setOption(pars.sampleName, "--sample", "Name of the sample/sample dir to cluster", true);
	setOption(pars.parameterFile, "--par", "Name of the parameters file for the clustering");
	setOption(pars.useReadLen, "--useReadLen", "Use Read Length large differences to skip alignment comparison");
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
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipClusteringMultiple(mipClusteringParsMultiple & pars) {
	pars.logFilename = "mipClustering";
	processMultipleDefaults(pars);
	setOption(pars.parameterFile, "--par", "Name of the parameters file for the clustering");
	setOption(pars.useReadLen, "--useReadLen", "Use Read Length large differences to skip alignment comparison");
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
	finishSetUp(std::cout);
}

void mipsterAnalysisSetUp::setUpMipPopulationClustering(
		mipPopulationClusteringPars & pars) {
	processDefaults(pars);
  setOption(pars.parameters, "-par", "ParametersFileName");
  setOption(pars.mipName, "--mipName", "Mip Name", true);
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars_.colOpts_.skipOpts_.skipOnLetterCounterDifference_, "-skip",
			"skipOnLetterCounterDifference");
	setOption(pars_.colOpts_.skipOpts_.fractionDifferenceCutOff_, "-skipCutOff",
			"fractionDifferenceCutOff");
	setOption(pars.previousPopFilename, "-previousPop", "previousPopFilename");
	processComparison(pars.previousPopErrors, "previousPop");
  setOption(pars.runsRequired, "-runsRequired", "runsRequired");
  setOption(pars.cutOff, "-cutoff", "FractionCutoff");


  setOption(pars.fracCutoff, "-fraccutoff", "PopulationClusteringFractionCutoff");
  pars.differentPar = setOption(pars.parametersPopulation, "-poppar",
      "ParametersForPopulationCollapse");
  setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
  setOption(pars.keepChimeras, "-keepChimeras", "KeepChimeras");
  setOption(pars_.chiOpts_.parentFreqs_, "-parfreqs", "ParentFrequence_multiplier_cutoff");

	setOption(pars.groupingsFile, "--groupingsFile",
			"A file to sort samples into different groups and add meta data to sample output");

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
	processMultipleDefaults(pars);
  setOption(pars.parameters, "-par", "ParametersFileName");
	setOption(pars.seqFileSuffix, "--seqFileSuffix", "The ending of the sequence append to sample name");
	setOption(pars_.colOpts_.skipOpts_.skipOnLetterCounterDifference_, "-skip",
			"skipOnLetterCounterDifference");
	setOption(pars_.colOpts_.skipOpts_.fractionDifferenceCutOff_, "-skipCutOff",
			"fractionDifferenceCutOff");
	setOption(pars.previousPopFilename, "-previousPop", "previousPopFilename");
	processComparison(pars.previousPopErrors, "previousPop");
  setOption(pars.runsRequired, "-runsRequired", "runsRequired");
  setOption(pars.cutOff, "-cutoff", "FractionCutoff");


  setOption(pars.fracCutoff, "-fraccutoff", "PopulationClusteringFractionCutoff");
  pars.differentPar = setOption(pars.parametersPopulation, "-poppar",
      "ParametersForPopulationCollapse");
  setOption(pars_.chiOpts_.checkChimeras_, "-markChimeras", "MarkChimeras");
  setOption(pars.keepChimeras, "-keepChimeras", "KeepChimeras");
  setOption(pars_.chiOpts_.parentFreqs_, "-parfreqs", "ParentFrequence_multiplier_cutoff");


	setOption(pars.groupingsFile, "--groupingsFile",
			"A file to sort samples into different groups and add meta data to sample output");


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


} /* namespace bibseq */