/*
 * parameters.cpp
 *
 *  Created on: Feb 12, 2016
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
#include "parameters.hpp"
namespace njhseq {

void mipCorePars::processDefaults(seqSetUp & setUp) {
	setUp.setOption(masterDir, "--masterDir", "Name of main analysis directory",
			true);
	//set mip id fnp
	mipArmsFileName = njh::files::make_path(masterDir.string(), "resources",
					"mip_arm_id.tab.txt");
	setUp.setOption(mipArmsFileName, "--mipArmsFilename",
				"Name of the mip arms file", infoFilesRequired);
	//set input mip and sample names
	mipsSamplesFile = njh::files::make_path(masterDir, "resources",
					"allMipsSamplesNames.tab.txt");
	setUp.setOption(mipsSamplesFile, "--mipSampleFile",
				"Mip sample filename, two columns, one column is mips, other is samples", infoFilesRequired);
	//set possible meta data
	auto possibleMeta = njh::files::make_path(masterDir, "resources",
					"samplesMeta.tab.txt");
	if(bfs::exists(possibleMeta)){
		sampleMetaFnp = possibleMeta;
	}
	setUp.setOption(sampleMetaFnp,"--samplesMeta", "A file with meta data for samples");
	setUp.setOption(overWriteDirs, "--overWriteDirs",
			"Over write the barcode correction sample dir if it already exists");
}

void mipCorePars::addCorePathsToConfig(Json::Value & config){
	config["masterDir"] = njh::json::toJson(masterDir);
	config["mipArmsFileName"] = njh::json::toJson(mipArmsFileName);
	config["mipsSamplesFile"] = njh::json::toJson(mipsSamplesFile);
	config["seqFileSuffix"] = njh::json::toJson(seqFileSuffix);
	config["sampleMetaFnp"] = njh::json::toJson(sampleMetaFnp);
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

	verbose_ = otherPars.verbose_;

}

}  // namespace njhseq
