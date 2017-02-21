/*
 * parameters.cpp
 *
 *  Created on: Feb 12, 2016
 *      Author: nick
 */

#include "parameters.hpp"
namespace bibseq {

void mipCorePars::processDefaults(seqSetUp & setUp) {
	setUp.setOption(masterDir, "--masterDir", "Name of main analysis directory",
			true);
	//set mip id fnp
	mipArmsFileName = bib::files::make_path(masterDir.string(), "resources",
					"mip_arm_id.tab.txt");
	setUp.setOption(mipArmsFileName, "--mipArmsFilename",
				"Name of the mip arms file", infoFilesRequired);
	//set input mip and sample names
	mipsSamplesFile = bib::files::make_path(masterDir, "resources",
					"allMipsSamplesNames.tab.txt");
	setUp.setOption(mipsSamplesFile, "--mipSampleFile",
				"Mip sample filename, two columns, one column is mips, other is samples", infoFilesRequired);
	//set possible meta data
	auto possibleMeta = bib::files::make_path(masterDir, "resources",
					"samplesMeta.tab.txt");
	if(bfs::exists(possibleMeta)){
		sampleMetaFnp = possibleMeta;
	}
	setUp.setOption(sampleMetaFnp,"--samplesMeta", "A file with meta data for samples");
	setUp.setOption(overWriteDirs, "--overWriteDirs",
			"Over write the barcode correction sample dir if it already exists");
}

void mipCorePars::addCorePathsToConfig(Json::Value & config){
	config["masterDir"] = bib::json::toJson(masterDir);
	config["mipArmsFileName"] = bib::json::toJson(mipArmsFileName);
	config["mipsSamplesFile"] = bib::json::toJson(mipsSamplesFile);
	config["seqFileSuffix"] = bib::json::toJson(seqFileSuffix);
	config["sampleMetaFnp"] = bib::json::toJson(sampleMetaFnp);
}


void mipCorePars::copyCore(const mipCorePars & otherPars) {
	wiggleRoom = otherPars.wiggleRoom;
	allowableErrors = otherPars.allowableErrors;
	masterDir = otherPars.masterDir;
	mipsSamplesFile = otherPars.mipsSamplesFile;
	mipArmsFileName = otherPars.mipArmsFileName;
	numThreads = otherPars.numThreads;
	logFilename = otherPars.logFilename;
	overWriteLog = otherPars.overWriteLog;
	overWriteDirs = otherPars.overWriteDirs;
	sampleMetaFnp = otherPars.sampleMetaFnp;
	seqFileSuffix = otherPars.seqFileSuffix;

	infoFilesRequired = otherPars.infoFilesRequired;
	logFileRequired = otherPars.logFileRequired;

}

}  // namespace bibseq
