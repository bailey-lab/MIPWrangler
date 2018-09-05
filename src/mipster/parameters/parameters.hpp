#pragma once
/*
 * parameters.hpp
 *
 *  Created on: Feb 8, 2016
 *      Author: nick
 */


#include "mipster/common.h"
#include <SeekDeep/objects/PairedReadProcessor.hpp>

namespace bibseq {

struct mipCorePars{

	bfs::path masterDir = "";
	bfs::path mipsSamplesFile = "";
	bfs::path sampleMetaFnp = "";
	uint32_t allowableErrors = 6;
	uint32_t wiggleRoom = 0;

	std::string seqFileSuffix = ".extendedFrags.fastq";

	bfs::path mipArmsFileName = "";
	uint32_t numThreads = 1;
	bfs::path logFilename = "";
	bool overWriteLog = false;
	bool overWriteDirs = false;
	bool infoFilesRequired = false;
	bool logFileRequired = true;
	bool verbose_ = false;

	void processDefaults(seqSetUp & setUp);

	void addCorePathsToConfig(Json::Value & config);

	void copyCore(const  mipCorePars & otherPars);

};





struct runGzExtractStitchPars : mipCorePars{
	bfs::path dir = "";
	uint32_t minOverlap = 10;
	uint32_t maxOverlap = 150;
	double mismatchDensity = 0.25;
	bool usePear = false;
	bool usePanda = false;
	bool trim = false;
};


struct mipIllumArmExtractionPars : mipCorePars{
	mipIllumArmExtractionPars(){
	  //quality filtering
	  qFilPars_.checkingQWindow = false;
	  qFilPars_.qualWindow_ = "50,5,20";
	  qFilPars_.qualityWindowLength_ = 50;
	  qFilPars_.qualityWindowStep_ = 5;
	  qFilPars_.qualityWindowThres_ = 20;

	  qFilPars_.checkingQFrac_ = false;
	  qFilPars_.qualCheck_ = 30;
	  qFilPars_.qualCheckCutOff_ = 0.75;
	}
	std::string sampleName = "";
	uint32_t minLen = 150;
  uint32_t smallFragmentLength = 50;
  QualFilteringPars qFilPars_;
	bool cacheAlignments = false;

	PairedReadProcessor::ProcessParams processPairPars_;



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
		ret.qFilPars_ = qFilPars_;
		ret.cacheAlignments = cacheAlignments;
		ret.processPairPars_ = processPairPars_;
		ret.copyCore(*this);
		return ret;
	}
};


struct extractFromRawPars : mipIllumArmExtractionPars {

	bfs::path dir = "";

};

struct extractFromRawParsMultiple : public extractFromRawPars {

	extractFromRawPars createForSample(const std::string & newSampleName)const{
		extractFromRawPars ret;
		ret.dir = dir;
		ret.minLen = minLen;
		ret.smallFragmentLength = smallFragmentLength;
		ret.fileOpenLimit_ = fileOpenLimit_;
		ret.sampleName = newSampleName;
		ret.qFilPars_ = qFilPars_;
		ret.cacheAlignments = cacheAlignments;
		ret.processPairPars_ = processPairPars_;
		ret.copyCore(*this);
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

	std::string qualRep = "median";

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
		ret.qualRep = qualRep;
		ret.sampleName = newSampleName;

		ret.copyCore(*this);
		return ret;
	}
};

struct mipClusteringPars : mipCorePars{
	bfs::path parameterFile = "";
  bool useReadLen = false;
  uint32_t readlenDiff = 15;
  CollapseIterations iterMap;
  std::string qualRep = "median";
	std::string sampleName = "";
	bool writeOutClusters = false;
	bool cacheAlignments = false;

};

struct mipClusteringParsMultiple : mipClusteringPars{
	mipClusteringPars createForSample(const std::string & newSampleName)const{
		mipClusteringPars ret;
		ret.useReadLen = useReadLen;
		ret.readlenDiff = readlenDiff;
		ret.parameterFile = parameterFile;
		ret.iterMap = iterMap;
		ret.qualRep = qualRep;
		ret.writeOutClusters = writeOutClusters;
		ret.cacheAlignments = cacheAlignments;
		ret.sampleName = newSampleName;

		ret.copyCore(*this);
		return ret;
	}
};

struct mipPopulationClusteringPars : mipCorePars{

	bfs::path parameters = "";
	collapse::SampleCollapseCollection::PreFilteringCutOffs clusteringCutOffs;

	double fracCutoff = 0.005;
	uint32_t runsRequired = 1;
	bool keepChimeras = false;
	std::string parametersPopulation = "";
	bool differentPar = false;
	CollapseIterations popIteratorMap;
	CollapseIterations iteratorMap;
	bfs::path previousPopFilename = "";
	bfs::path previousPopDir = "";
	SeqIOOptions refIoOptions;
	comparison previousPopErrors;
	std::string seqFileSuffix = "_clustered.fastq";
	std::string mipName = "";

};

struct mipPopulationClusteringParsMultiple : public mipPopulationClusteringPars{
	mipPopulationClusteringPars createForMip(const std::string & newMipName)const{
		mipPopulationClusteringPars ret;
		ret.parameters = parameters;
		ret.clusteringCutOffs = clusteringCutOffs;
		ret.fracCutoff = fracCutoff;
		ret.runsRequired = runsRequired;
		ret.keepChimeras = keepChimeras;
		ret.parametersPopulation = parametersPopulation;
		ret.differentPar = differentPar;
		ret.popIteratorMap = popIteratorMap;
		ret.iteratorMap = iteratorMap;
		ret.previousPopErrors = previousPopErrors;

		ret.previousPopDir = previousPopDir;


		ret.mipName = newMipName;
		auto refPopFile = bib::files::make_path(ret.previousPopDir, ret.mipName + ".fasta");
		if(bfs::exists(refPopFile)){
			ret.previousPopFilename = refPopFile.string();
			ret.refIoOptions = SeqIOOptions::genFastaIn(refPopFile);
		}
		ret.copyCore(*this);
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
		ret.mipName = newMipName;
		ret.readCutOff = readCutOff;
		ret.copyCore(*this);
		return ret;
	}

};




}  // namespace bibseq

