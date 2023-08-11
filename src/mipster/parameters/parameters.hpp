#pragma once
/*
 * parameters.hpp
 *
 *  Created on: Feb 8, 2016
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
#include "mipster/common.h"
#include <SeekDeep/objects/IlluminaUtils/PairedReadProcessor.hpp>

namespace njhseq {

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
	bool develop_ = false;

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
	  processPairPars_.r1Trim_ = 1;
	  processPairPars_.r2Trim_ = 1;
	}
	std::string sampleName;
	uint32_t minCaptureLength = 100;
  uint32_t smallFragmentLength = 50;
  QualFilteringPars qFilPars_;
	bool cacheAlignments = false;
	bool keepIntermediateFiles = false;
	bool writeOutInitialExtractedPairs = false;

	bfs::path refDir = "";

	PairedReadProcessor::ProcessParams processPairPars_;

	uint32_t seqOutCacheLimit_ = 50000;

#if defined( __APPLE__ ) || defined( __APPLE_CC__ ) || defined( macintosh ) || defined( __MACH__ )
	uint32_t fileOpenLimit_ = 200; /**< The maximum number of files to be kept open at one time*/
#else
	uint32_t fileOpenLimit_ = 1000; /**< The maximum number of files to be kept open at one time */
#endif
};

struct mipIllumArmExtractionParsMultiple : public mipIllumArmExtractionPars {


	mipIllumArmExtractionPars createForSample(const std::string & newSampleName)const{
		mipIllumArmExtractionPars ret;
		ret.minCaptureLength = minCaptureLength;
		ret.smallFragmentLength = smallFragmentLength;
		ret.fileOpenLimit_ = fileOpenLimit_;
		ret.sampleName = newSampleName;
		ret.qFilPars_ = qFilPars_;
		ret.cacheAlignments = cacheAlignments;
		ret.keepIntermediateFiles = keepIntermediateFiles;
		ret.writeOutInitialExtractedPairs = writeOutInitialExtractedPairs;
		ret.processPairPars_ = processPairPars_;
		ret.refDir = refDir;
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
		ret.minCaptureLength = minCaptureLength;
		ret.smallFragmentLength = smallFragmentLength;
		ret.fileOpenLimit_ = fileOpenLimit_;
		ret.sampleName = newSampleName;
		ret.qFilPars_ = qFilPars_;
		ret.cacheAlignments = cacheAlignments;
		ret.processPairPars_ = processPairPars_;
		ret.keepIntermediateFiles = keepIntermediateFiles;
		ret.writeOutInitialExtractedPairs = writeOutInitialExtractedPairs;
		ret.refDir = refDir;

		ret.seqOutCacheLimit_ = seqOutCacheLimit_;

		ret.copyCore(*this);
		return ret;
	}
};


struct mipBarcodeCorrectionPars : mipCorePars{
	uint32_t clusterCutOff = 5;
  bool useReadLen = false;
  uint32_t readlenDiff = 15;
	double barcodeIdentity = 0.98;
	std::string sampleName;

	bool writeExtra = false;

	std::string qualRep = "median";

	bool keepIntermediateFiles = false;
	bool cacheAlignments = false;


	bool doNotDownSample_{false};
	uint32_t downSampleAmount_{20000};
	uint32_t seed_{0};
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
		ret.keepIntermediateFiles = keepIntermediateFiles;
		ret.cacheAlignments = cacheAlignments;

		ret.doNotDownSample_ = doNotDownSample_;
		ret.downSampleAmount_ = downSampleAmount_;
		ret.seed_ = seed_;

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

	bool keepIntermediateFiles = false;


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

		ret.keepIntermediateFiles = keepIntermediateFiles;

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

	bool keepIntermediateFiles = false;
	bool cacheAlignments = false;

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

		ret.keepIntermediateFiles = keepIntermediateFiles;
		ret.cacheAlignments = cacheAlignments;

		ret.mipName = newMipName;
		auto refPopFile = njh::files::make_path(ret.previousPopDir, ret.mipName + ".fasta");
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




}  // namespace njhseq

