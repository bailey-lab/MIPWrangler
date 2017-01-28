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
	if (!setUp.setOption(mipArmsFileName, "--mipArmsFilename",
			"Name of the mip arms file", infoFilesRequired)) {
		mipArmsFileName = bib::files::make_path(masterDir.string(), "resources",
				"mip_arm_id.tab.txt");
	}
	if (!setUp.setOption(mipsSamplesFile, "--mipSampleFile",
			"Mip sample filename, two columns, one column is mips, other is samples", infoFilesRequired)) {
		mipsSamplesFile = bib::files::make_path(masterDir, "resources",
				"allMipsSamplesNames.tab.txt");
	}
	setUp.setOption(overWriteDirs, "--overWriteDirs",
			"Over write the barcode correction sample dir if it already exists");
}

void mipCorePars::addCorePathsToConfig(Json::Value & config){
	config["masterDir"] = bib::json::toJson(masterDir);
	config["mipArmsFileName"] = bib::json::toJson(mipArmsFileName);
	config["mipsSamplesFile"] = bib::json::toJson(mipsSamplesFile);
	config["seqFileSuffix"] = bib::json::toJson(seqFileSuffix);
}

}  // namespace bibseq
