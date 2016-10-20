#pragma once
/*
 * parameters.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


#include "mipster/common.h"

namespace bibseq {

struct mipCorePars{

	std::string masterDir = "";
	std::string mipsSamplesFile = "";
	uint32_t allowableErrors = 0;
	uint32_t wiggleRoom = 0;

	std::string seqFileSuffix = ".extendedFrags.fastq";

	std::string mipArmsFileName = "";
	uint32_t numThreads = 1;
	std::string logFilename = "";
	bool overWriteLog = false;
	bool overWriteDirs = false;
	void processDefaults(seqSetUp & setUp);

	void addCorePathsToConfig(Json::Value & config);

};

struct runGzExtractStitchPars : mipCorePars{
	std::string dir = "";
	uint32_t minOverlap = 10;
	uint32_t maxOverlap = 150;
	double mismatchDensity = 0.25;
	bool usePear = false;
	bool usePanda = false;
	bool trim = false;
};


struct mipIllumArmExtractionPars : mipCorePars{
	std::string sampleName = "";
	uint32_t minLen = 150;
  uint32_t smallFragmentLength = 50;


#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	uint32_t fileOpenLimit_ = 200; /**< The maximum number of files to be kept open at one time*/
#else
	uint32_t fileOpenLimit_ = 1000; /**< The maximum number of files to be kept open at one time */
#endif
};

struct mipIllumArmExtractionParsMultiple : public mipIllumArmExtractionPars {


	mipIllumArmExtractionPars createForSample(const std::string & newSampleName)const{
		mipIllumArmExtractionPars ret;
		ret.minLen = minLen;
		ret.smallFragmentLength = smallFragmentLength;
		ret.fileOpenLimit_ = fileOpenLimit_;
		ret.sampleName = newSampleName;

		ret.wiggleRoom = wiggleRoom;
		ret.allowableErrors = allowableErrors;
		ret.masterDir = masterDir;
		ret.mipsSamplesFile = mipsSamplesFile;
		ret.mipArmsFileName = mipArmsFileName;
		ret.numThreads = numThreads;
		ret.logFilename =logFilename;
		ret.overWriteLog = overWriteLog;
		ret.overWriteDirs = overWriteDirs;
		return ret;
	}
};



struct mipBarcodeCorrectionPars : mipCorePars{
	uint32_t clusterCutOff = 5;
  bool useReadLen = false;
  uint32_t readlenDiff = 15;
	double barcodeIdentity = 0.98;
	std::string sampleName = "";

	bool writeExtra = false;

};

struct mipBarcodeCorrectionParsMultiple : public mipBarcodeCorrectionPars  {
	//for when running a sample
	mipBarcodeCorrectionPars createForSample(const std::string & newSampleName)const{
		mipBarcodeCorrectionPars ret;
		ret.clusterCutOff = clusterCutOff;
		ret.useReadLen = useReadLen;
		ret.readlenDiff = readlenDiff;
		ret.barcodeIdentity = barcodeIdentity;
		ret.writeExtra = writeExtra;
		ret.sampleName = newSampleName;

		ret.allowableErrors = allowableErrors;
		ret.wiggleRoom = wiggleRoom;

		ret.masterDir = masterDir;
		ret.mipsSamplesFile = mipsSamplesFile;
		ret.mipArmsFileName = mipArmsFileName;
		ret.numThreads = numThreads;
		ret.logFilename =logFilename;
		ret.overWriteLog = overWriteLog;
		ret.overWriteDirs = overWriteDirs;
		return ret;
	}
};

struct mipClusteringPars : mipCorePars{
	std::string parameterFile = "";
  bool useReadLen = false;
  uint32_t readlenDiff = 15;
  CollapseIterations iterMap;
	std::string sampleName = "";

};

struct mipClusteringParsMultiple : mipClusteringPars{
	mipClusteringPars createForSample(const std::string & newSampleName)const{
		mipClusteringPars ret;
		ret.useReadLen = useReadLen;
		ret.readlenDiff = readlenDiff;
		ret.parameterFile = parameterFile;
		ret.iterMap = iterMap;
		ret.sampleName = newSampleName;

		ret.allowableErrors = allowableErrors;
		ret.masterDir = masterDir;
		ret.mipsSamplesFile = mipsSamplesFile;
		ret.mipArmsFileName = mipArmsFileName;
		ret.numThreads = numThreads;
		ret.logFilename =logFilename;
		ret.overWriteLog = overWriteLog;
		ret.overWriteDirs = overWriteDirs;
		return ret;
	}
};

struct mipPopulationClusteringPars : mipCorePars{

	std::string parameters = "";
	uint32_t cutOff = 1;
	double fracCutoff = 0.005;
	uint32_t runsRequired = 1;
	bool keepChimeras = false;
	std::string parametersPopulation = "";
	bool differentPar = false;
	CollapseIterations popIteratorMap;
	CollapseIterations iteratorMap;
	std::string previousPopFilename = "";
	comparison previousPopErrors;
	std::string seqFileSuffix = "_clustered.fastq";
	std::string mipName = "";
	std::string groupingsFile = "";

};

struct mipPopulationClusteringParsMultiple : public mipPopulationClusteringPars{
	mipPopulationClusteringPars createForMip(const std::string & newMipName)const{
		mipPopulationClusteringPars ret;
		ret.parameters = parameters;
		ret.cutOff = cutOff;
		ret.fracCutoff = fracCutoff;
		ret.runsRequired = runsRequired;
		ret.keepChimeras = keepChimeras;
		ret.parametersPopulation = parametersPopulation;
		ret.differentPar = differentPar;
		ret.popIteratorMap = popIteratorMap;
		ret.iteratorMap = iteratorMap;
		ret.previousPopFilename = previousPopFilename;
		ret.previousPopErrors = previousPopErrors;
		ret.seqFileSuffix = seqFileSuffix;
		ret.groupingsFile = groupingsFile;

		ret.mipName = newMipName;

		ret.allowableErrors = allowableErrors;
		ret.masterDir = masterDir;
		ret.mipsSamplesFile = mipsSamplesFile;
		ret.mipArmsFileName = mipArmsFileName;
		ret.numThreads = numThreads;
		ret.logFilename = logFilename;
		ret.overWriteLog = overWriteLog;
		ret.overWriteDirs = overWriteDirs;
		return ret;
	}

};


struct mipCorrectForContamWithSameBarcodesPars : public mipCorePars{

	std::string mipName = "";
	uint32_t readCutOff = 5;

};

struct mipCorrectForContamWithSameBarcodesParsMultiple : public mipCorrectForContamWithSameBarcodesPars{
	mipCorrectForContamWithSameBarcodesPars createForMip(const std::string & newMipName)const{
		mipCorrectForContamWithSameBarcodesPars ret;
		ret.seqFileSuffix = seqFileSuffix;
		ret.wiggleRoom = wiggleRoom;
		ret.mipName = newMipName;
		ret.readCutOff = readCutOff;
		ret.allowableErrors = allowableErrors;
		ret.masterDir = masterDir;
		ret.mipsSamplesFile = mipsSamplesFile;
		ret.mipArmsFileName = mipArmsFileName;
		ret.numThreads = numThreads;
		ret.logFilename = logFilename;
		ret.overWriteLog = overWriteLog;
		ret.overWriteDirs = overWriteDirs;
		return ret;
	}

};




}  // namespace bibseq

