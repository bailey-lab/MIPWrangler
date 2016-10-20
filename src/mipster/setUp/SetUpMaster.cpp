/*
 * SetUpMaster.cpp
 *
 *  Created on: Feb 11, 2016
 *      Author: nick
 */


#include "SetUpMaster.hpp"

namespace bibseq {

void SetUpMaster::makeBarcodeCorDirs()const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to create sample directories when when member names_ isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	for (const auto & samp : names_->samples_) {
		SampleDirectoryMaster sampDirMaster(directoryMaster_,
				MipFamSamp("", samp));
		sampDirMaster.ensureBarCorDirectoryExist();
	}
}

void SetUpMaster::makeClusteringDirs()const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to create sample directories when when member names_ isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	for (const auto & samp : names_->samples_) {
		SampleDirectoryMaster sampDirMaster(directoryMaster_, MipFamSamp("", samp));
		sampDirMaster.ensureClusDirectoryExist();
	}
}

SetUpMaster::SetUpMaster(const std::string & masterDir) :
		directoryMaster_(masterDir), mipArmFnp_(
				bib::files::join(directoryMaster_.resourceDir_.string(),
						"mip_arm_id.tab.txt")), mipsSampsNamesFnp_(
				bib::files::join(directoryMaster_.resourceDir_.string(),
						"allMipsSamplesNames.tab.txt")) {
}

void SetUpMaster::setServerName(const std::string & mipServerName){
	mipServerName_ = mipServerName;
}

bfs::path SetUpMaster::getMipSerDir() const{
	return directoryMaster_.getMipSerDir(mipServerName_);
}

void SetUpMaster::setMipArmFnp(const std::string & mipArmFnp) {
	if(!bfs::exists(mipArmFnp)){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " " << mipArmFnp
				<< " doesn't exist" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipArmFnp_ = mipArmFnp;
}

void SetUpMaster::setMipsSampsNamesFnp(const std::string & mipsSampsNamesFnp){
	if(!bfs::exists(mipsSampsNamesFnp)){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " " << mipsSampsNamesFnp
				<< " doesn't exist" << std::endl;
		throw std::runtime_error { ss.str() };
	}
	mipsSampsNamesFnp_ = mipsSampsNamesFnp;
}

void SetUpMaster::createDirStructSkeleton() const {
	if (bfs::exists(directoryMaster_.masterDir_)) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " directory "
				<< directoryMaster_.masterDir_
				<< " already exists, either delete it or choose another directory"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	} else {
		MipAnalysisDirectoryMaster::initMipAnalysisDirectoryStructure(
				directoryMaster_.masterDir_);
		//create readmes
		std::ofstream readme(directoryMaster_.resourceDir_.string() + "README.txt");
		readme << "Need to add three files to this directory" << std::endl;
		readme << "1) mip_arm_id.tab.txt" << std::endl;
		readme << "\t*This file should contain the following columns mip_id,mip_family,extension_arm,ligation_arm,extension_barcode_length,ligation_barcode_length,species,minlen" << std::endl;
	  readme << "\t\t-mip_id: A unique identifier for this mip target" << std::endl;
	  readme << "\t\t-mip_family: A mip family this mip target belongs to, if target does not belong to a family should just be the mip id again" << std::endl;
	  readme << "\t\t-extension_arm: The extension arm sequence, should be as is found in the fastq file" << std::endl;
	  readme << "\t\t-ligation_arm: The ligation arm sequence, should be in the reverse complement of what is in the stitched fastq file" << std::endl;
	  readme << "\t\t-extension_barcode_length: The length of the extension arm barcode for this mip target" << std::endl;
	  readme << "\t\t-ligation_barcode_length: The length of the ligation arm barcode for this mip target" << std::endl;
	  readme << "\t\t-species: The species/genome name for this mip target" << std::endl;
	  readme << "\t\t-minlen: The minimum length expected for this mip target" << std::endl;
		readme << "2) allMipsSamplesNames.tab.txt" << std::endl;
		readme << "\t*This file should contain two columns, mips and samples" << std::endl;
	  readme << "\t\t-mips: All the mip family names for this dataset" << std::endl;
	  readme << "\t\t-samples: All the sample names for this dataset" << std::endl;
		readme << "3) mipInfo.json" << std::endl;
		readme << "\t*This file should contain all important information for each mip target, including genomic location, copies captured,important snps etc" << std::endl;

		std::ofstream readmeMarkDown(directoryMaster_.resourceDir_.string() + "README.md");
		readmeMarkDown << "Need to add three files to this directory, the resource directory  " << std::endl;
		readmeMarkDown << std::endl;
		readmeMarkDown << "+  **mip\\_arm\\_id.tab.txt**  "<< std::endl;
		readmeMarkDown << "This file should contain at least the following columns: "<< std::endl;
		readmeMarkDown << "	+ **mip\\_id**: A unique identifier for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **mip\\_family**: A mip family this mip target belongs to, if target does not belong to a family should just be the mip id again  "<< std::endl;
		readmeMarkDown << "	+ **extension\\_arm**: The extension arm sequence, should be as is found in the fastq file  "<< std::endl;
		readmeMarkDown << "	+ **ligation\\_arm**: The ligation arm sequence, should be in the reverse complement of what is in the stitched fastq file  "<< std::endl;
		readmeMarkDown << "	+ **extension\\_barcode\\_length**: The length of the extension arm barcode for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **ligation\\_barcode\\_length**: The length of the ligation arm barcode for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **species**: The species/genome name for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **minlen**: The minimum length expected for this mip target  "<< std::endl;
		readmeMarkDown << "+ allMipsSamplesNames.tab.txt  "<< std::endl;
		readmeMarkDown << "	This file should contain at least the following two columns, mips and samples  "<< std::endl;
		readmeMarkDown << "	+ **mips**: All the mip family names for this dataset  "<< std::endl;
		readmeMarkDown << "	+ **samples**: All the sample names for this dataset  "<< std::endl;
		readmeMarkDown << "+ mipInfo.json  "<< std::endl;
		readmeMarkDown << "	This file should contain all important information for each mip target, including genomic location, copies captured,important snps etc  "<< std::endl;
		if(!bfs::exists(mipArmFnp_)){
			std::ofstream mipArm(mipArmFnp_);
			mipArm << "REPLACE ME!" << std::endl;
		}
		if(!bfs::exists(mipsSampsNamesFnp_)){
			std::ofstream mipSampNames(mipsSampsNamesFnp_);
			mipSampNames << "REPLACE ME!" << std::endl;
		}
	}
}

void SetUpMaster::createDirStructSkeleton(const std::string & mipSampleFile,
		const std::string & mipArms) {
	if (bfs::exists(directoryMaster_.masterDir_)) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__ << " directory "
				<< directoryMaster_.masterDir_
				<< " already exists, either delete it or choose another directory"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	} else {
		MipAnalysisDirectoryMaster::initMipAnalysisDirectoryStructure(
				directoryMaster_.masterDir_);
		//create readmes
		std::ofstream readme(directoryMaster_.resourceDir_.string() + "README.txt");
		readme << "Need to add three files to this directory" << std::endl;
		readme << "1) mip_arm_id.tab.txt" << std::endl;
		readme << "\t*This file should contain the following columns mip_id,mip_family,extension_arm,ligation_arm,extension_barcode_length,ligation_barcode_length,species,minlen" << std::endl;
	  readme << "\t\t-mip_id: A unique identifier for this mip target" << std::endl;
	  readme << "\t\t-mip_family: A mip family this mip target belongs to, if target does not belong to a family should just be the mip id again" << std::endl;
	  readme << "\t\t-extension_arm: The extension arm sequence, should be as is found in the fastq file" << std::endl;
	  readme << "\t\t-ligation_arm: The ligation arm sequence, should be in the reverse complement of what is in the stitched fastq file" << std::endl;
	  readme << "\t\t-extension_barcode_length: The length of the extension arm barcode for this mip target" << std::endl;
	  readme << "\t\t-ligation_barcode_length: The length of the ligation arm barcode for this mip target" << std::endl;
	  readme << "\t\t-species: The species/genome name for this mip target" << std::endl;
	  readme << "\t\t-minlen: The minimum length expected for this mip target" << std::endl;
		readme << "2) allMipsSamplesNames.tab.txt" << std::endl;
		readme << "\t*This file should contain two columns, mips and samples" << std::endl;
	  readme << "\t\t-mips: All the mip family names for this dataset" << std::endl;
	  readme << "\t\t-samples: All the sample names for this dataset" << std::endl;
		readme << "3) mipInfo.json" << std::endl;
		readme << "\t*This file should contain all important information for each mip target, including genomic location, copies captured,important snps etc" << std::endl;

		std::ofstream readmeMarkDown(directoryMaster_.resourceDir_.string() + "README.md");
		readmeMarkDown << "Need to add three files to this directory, the resource directory  " << std::endl;
		readmeMarkDown << std::endl;
		readmeMarkDown << "+  **mip\\_arm\\_id.tab.txt**  "<< std::endl;
		readmeMarkDown << "This file should contain at least the following columns: "<< std::endl;
		readmeMarkDown << "	+ **mip\\_id**: A unique identifier for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **mip\\_family**: A mip family this mip target belongs to, if target does not belong to a family should just be the mip id again  "<< std::endl;
		readmeMarkDown << "	+ **extension\\_arm**: The extension arm sequence, should be as is found in the fastq file  "<< std::endl;
		readmeMarkDown << "	+ **ligation\\_arm**: The ligation arm sequence, should be in the reverse complement of what is in the stitched fastq file  "<< std::endl;
		readmeMarkDown << "	+ **extension\\_barcode\\_length**: The length of the extension arm barcode for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **ligation\\_barcode\\_length**: The length of the ligation arm barcode for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **species**: The species/genome name for this mip target  "<< std::endl;
		readmeMarkDown << "	+ **minlen**: The minimum length expected for this mip target  "<< std::endl;
		readmeMarkDown << "+ allMipsSamplesNames.tab.txt  "<< std::endl;
		readmeMarkDown << "	This file should contain at least the following two columns, mips and samples  "<< std::endl;
		readmeMarkDown << "	+ **mips**: All the mip family names for this dataset  "<< std::endl;
		readmeMarkDown << "	+ **samples**: All the sample names for this dataset  "<< std::endl;
		readmeMarkDown << "+ mipInfo.json  "<< std::endl;
		readmeMarkDown << "	This file should contain all important information for each mip target, including genomic location, copies captured,important snps etc  "<< std::endl;
		bfs::copy_file(mipArms, mipArmFnp_);
		bfs::copy_file(mipSampleFile, mipsSampsNamesFnp_);
		names_ = std::make_shared<MipsSamplesNames>(mipsSampsNamesFnp_);
		createPopClusMipDirs(1);
		createTopSampleDirs();
	}
}

void SetUpMaster::createTopSampleDirs() const {
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to create sample directories when when member names_ isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (bfs::exists(directoryMaster_.masterDir_)) {
		for (const auto & samp : names_->samples_) {
			bib::files::makeDir(
					bib::files::MkdirPar(
							bib::files::join(directoryMaster_.masterDir_.string(), samp)));
		}
	} else {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to create sample directories when "
				<< directoryMaster_.masterDir_ << " doesn't exist"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
}

void SetUpMaster::createPopClusMipDirs(uint32_t numThreads) const {
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to create mip pop directories when when member names_ isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	if (bfs::exists(directoryMaster_.populationClusteringDir_)) {
		bib::concurrent::LockableQueue<std::string> mipQueue(names_->mips_);
		std::string popDir = directoryMaster_.populationClusteringDir_.string();
		auto mkDirs = [&mipQueue,&popDir]() {
			std::string m = "";
			std::stringstream ss;
			while(mipQueue.getVal(m)) {
				bib::files::makeDir(popDir, bib::files::MkdirPar(m));
			}
		};
		std::vector<std::thread> threads;
		for (uint32_t t = 0; t < numThreads; ++t) {
			threads.emplace_back(std::thread(mkDirs));
		}
		for (auto & t : threads) {
			t.join();
		}
	} else {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to fill the population clustering directory "
				<< directoryMaster_.populationClusteringDir_ << " when it doesn't exist"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
}


void SetUpMaster::loadMipsSampsInfo(uint32_t allowableArmErrors){
	mips_ = std::make_shared<MipCollection>(mipArmFnp_, allowableArmErrors);
	names_ = std::make_shared<MipsSamplesNames>(mipsSampsNamesFnp_);
}

VecStr SetUpMaster::checkDirStruct() const{
	bool passedAllChecks = true;
	VecStr warnings;
	//check for directories
	if (!checkForFileDirExistence(directoryMaster_.masterDir_, warnings)) {
		passedAllChecks = false;
	}
	if (!checkForFileDirExistence(directoryMaster_.logsDir_, warnings)) {
		passedAllChecks = false;
	}
	if (!checkForFileDirExistence(directoryMaster_.populationClusteringDir_, warnings)) {
		passedAllChecks = false;
	}
	if (!checkForFileDirExistence(directoryMaster_.resourceDir_, warnings)) {
		passedAllChecks = false;
	}
	if (!checkForFileDirExistence(directoryMaster_.serResourceDir_, warnings)) {
		passedAllChecks = false;
	}
	//check for files
	bool armFilePass = false;
	bool mipSampFilePass = false;
	if(!checkForFileDirExistence(mipArmFnp_, warnings)){
		passedAllChecks = false;
	}else{
		//check header
		table mipInfo(mipArmFnp_, "whitespace", true);
		VecStr neededColumns = { "mip_id", "extension_arm", "ligation_arm",
				"mip_family", "extension_barcode_length", "ligation_barcode_length" };
		VecStr columnsNotFound;
		for (const auto & col : neededColumns) {
			if (!bib::in(col, mipInfo.columnNames_)) {
				columnsNotFound.emplace_back(col);
			}
		}
		if (!columnsNotFound.empty()) {
			std::stringstream ss;
			ss << "Warning for :" << mipArmFnp_ << std::endl;
			ss << "Need to have " << vectorToString(neededColumns, ",") << std::endl;
			ss << "Did not find " << vectorToString(columnsNotFound, ",") << std::endl;
			ss << "Only have " << vectorToString(mipInfo.columnNames_, ",")
					<< std::endl;
			warnings.emplace_back(ss.str());
		}else{
			armFilePass = true;
		}
	}
	if(!checkForFileDirExistence(mipsSampsNamesFnp_, warnings)){
		passedAllChecks = false;
	}else{
		table mipSampInfo(mipsSampsNamesFnp_, "\t", true);
		VecStr missingCols;
		VecStr neededCols { "mips", "samples" };
		for (const auto & col : neededCols) {
			if (!bib::in(col, mipSampInfo.columnNames_)) {
				missingCols.emplace_back(col);
			}
		}
		if (!missingCols.empty()) {
			std::stringstream ss;
			ss << "Warnings for: "<< mipsSampsNamesFnp_ << std::endl;
			ss << "Missing the following columns, "
					<< bib::conToStr(missingCols, ",") << std::endl;
			ss << "Need the following columns, " << bib::conToStr(neededCols)
					<< std::endl;
			warnings.emplace_back(ss.str());
		}else{
			mipSampFilePass = true;
		}
	}
	if(armFilePass && mipSampFilePass){
		std::stringstream ss;
		MipsSamplesNames names(mipsSampsNamesFnp_);
		MipCollection mips(mipArmFnp_, 0);
		bool fail = false;
		for(const auto & m : names.mips_){
			if (!mips.hasMipFamily(m)) {
				fail = true;
				ss << "Error, have " << m << " in "
						<< directoryMaster_.resourceDir_.string()
								+ "allMipsSamplesNames.tab.txt"
						<< " but have no description for it in "
						<< mipArmFnp_
						<< std::endl;
			}
		}
		if(fail){
			ss << "Possible Families are: " << vectorToString(mips.mipFamilies_, ", ") << std::endl;;
			warnings.emplace_back(ss.str());
		}
	}
	return warnings;
}

bool SetUpMaster::checkForExtractedMipFamForSamp(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	//check to see if family has any mips with extracted reads;
	bool familyHasExtractedReads = false;
	auto mipsForFam = mips_->getMipsForFamily(mipSampName.mipFam_);
	for (const auto & mipName : mipsForFam) {
		if (bfs::exists(
				bib::files::join(
						VecStr { sampDirMaster.extractDir_.string(), mipName, mipName + ".fastq" }))) {
			familyHasExtractedReads = true;
			break;
		}
	}
	return familyHasExtractedReads;
}

bool SetUpMaster::checkForBarCorMipFamForSamp(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	return bfs::exists(bib::files::join(VecStr { sampDirMaster.barCorDir_.string(),
		mipSampName.mipFam_, mipSampName.mipFam_ + "_all.fastq" }));
}

bool SetUpMaster::checkForExtractedMipFamForSampThrow(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	sampDirMaster.checkForExtractDirectoryThrow();
	sampDirMaster.checkForBarCorDirectoryThrow();
	//check to see if family has any mips with extracted reads;
	bool familyHasExtractedReads = false;
	auto mipsForFam = mips_->getMipsForFamily(mipSampName.mipFam_);
	for (const auto & mipName : mipsForFam) {
		if (bfs::exists(
				bib::files::join(
						VecStr { sampDirMaster.extractDir_.string(), mipName, mipName + ".fastq" }))) {
			familyHasExtractedReads = true;
			break;
		}
	}
	return familyHasExtractedReads;
}

bool SetUpMaster::checkForBarCorMipFamForSampThrow(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	sampDirMaster.checkForExtractDirectoryThrow();
	sampDirMaster.checkForBarCorDirectoryThrow();
	return bfs::exists(bib::files::join(VecStr { sampDirMaster.barCorDir_.string(),
		mipSampName.mipFam_, mipSampName.mipFam_ + "_all.fastq" }));
}




std::vector<MipFamSamp> SetUpMaster::getPairsWithClustered(uint32_t numThreads)const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	if(!mips_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member mips_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	bib::concurrent::LockableQueue<std::string> sampQueue(names_->samples_);
	std::vector<MipFamSamp> pairs;

	std::mutex pairsMut;
	auto checkSample = [&pairsMut,&pairs](bib::concurrent::LockableQueue<std::string>& sampQueue,
			const SetUpMaster & mipMaster){
		std::string samp = "";
		std::vector<MipFamSamp> currentPairs;
		while(sampQueue.getVal(samp)){
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", samp));
			for(const auto & mipFam : mipMaster.names_->mips_){
					if (bfs::exists(bib::files::join(VecStr { sampDirMaster.clusDir_.string(),
						mipFam, mipFam + "_clustered.fastq" }))) {
						currentPairs.emplace_back(MipFamSamp(mipFam, samp));
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentPairs);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t threadNum = 0; threadNum < numThreads; ++threadNum){
		threads.emplace_back(std::thread(checkSample, std::ref(sampQueue), std::cref(*this)));
	}
	for(auto & t : threads){
		t.join();
	}
	return pairs;
}
std::vector<MipFamSamp> SetUpMaster::getMipFamsWithPopClustered(uint32_t numThreads) const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	bib::concurrent::LockableQueue<std::string> mipQueue(names_->mips_);
	std::vector<MipFamSamp> pairs;

	std::mutex pairsMut;
	auto checkSample = [&pairsMut,&pairs](bib::concurrent::LockableQueue<std::string>& mipQueue,
			const SetUpMaster & mipMaster){
		std::string mipFam = "";
		std::vector<MipFamSamp> currentPairs;
		while(mipQueue.getVal(mipFam)){
			auto popClusPath = mipMaster.pathMipPopClusHaplo(MipFamSamp(mipFam, ""));
			if(bfs::exists(popClusPath)){
				currentPairs.emplace_back(MipFamSamp(mipFam, ""));
			}
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentPairs);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t threadNum = 0; threadNum < numThreads; ++threadNum){
		threads.emplace_back(std::thread(checkSample, std::ref(mipQueue), std::cref(*this)));
	}
	for(auto & t : threads){
		t.join();
	}
	return pairs;
}


std::vector<MipFamSamp> SetUpMaster::getPairsWithBarCor(uint32_t numThreads)const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	bib::concurrent::LockableQueue<std::string> sampQueue(names_->samples_);
	std::vector<MipFamSamp> pairs;
	std::mutex pairsMut;

	auto checkSample = [&pairsMut,&pairs](bib::concurrent::LockableQueue<std::string>& sampQueue,
			const SetUpMaster & mipMaster){
		std::string samp = "";
		std::vector<MipFamSamp> currentPairs;
		while(sampQueue.getVal(samp)){
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", samp));
			for(const auto & mipFam : mipMaster.names_->mips_){
					if (bfs::exists(bib::files::join(VecStr { sampDirMaster.barCorDir_.string(),
						mipFam, mipFam + "_all.fastq" }))) {
						currentPairs.emplace_back(MipFamSamp(mipFam, samp));
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentPairs);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t threadNum = 0; threadNum < numThreads; ++threadNum){
		threads.emplace_back(std::thread(checkSample, std::ref(sampQueue), std::cref(*this)));
	}
	for(auto & t : threads){
		t.join();
	}
	return pairs;
}

void SetUpMaster::setRawDataSuffix(const std::string & rawDataSuffix){
	rawDataSuffix_ = rawDataSuffix;
}

std::vector<MipFamSamp> SetUpMaster::getPairsWithExtracted(uint32_t numThreads) const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	bib::concurrent::LockableQueue<std::string> sampQueue(names_->samples_);
	std::vector<MipFamSamp> pairs;
	std::mutex pairsMut;

	auto checkSample = [&pairsMut,&pairs](bib::concurrent::LockableQueue<std::string>& sampQueue,
			const SetUpMaster & mipMaster){
		std::string samp = "";
		std::vector<MipFamSamp> currentPairs;
		while(sampQueue.getVal(samp)){
			SampleDirectoryMaster sampDirMaster(mipMaster.directoryMaster_, MipFamSamp("", samp));
			for(const auto & mipFam : mipMaster.names_->mips_){
				auto mipsForFam = mipMaster.mips_->getMipsForFamily(mipFam);
				for (const auto & mipName : mipsForFam) {
					if (bfs::exists(
							bib::files::join(
									VecStr { sampDirMaster.extractDir_.string(), mipName, mipName + ".fastq" }))) {
						currentPairs.emplace_back(MipFamSamp(mipFam, samp));
						break;
					}
				}
			}
		}
		{
			std::lock_guard<std::mutex> lock(pairsMut);
			addOtherVec(pairs, currentPairs);
		}
	};
	std::vector<std::thread> threads;
	for(uint32_t threadNum = 0; threadNum < numThreads; ++threadNum){
		threads.emplace_back(std::thread(checkSample, std::ref(sampQueue), std::cref(*this)));
	}
	for(auto & t : threads){
		t.join();
	}
	return pairs;
}

std::vector<MipFamSamp> SetUpMaster::getSamplesWithRawData(uint32_t numThreads)const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::vector<MipFamSamp> ret;
	for(const auto & samp : names_->samples_){
		MipFamSamp mipSamp("", samp);
		if(checkForRawDataForSamp(mipSamp)){
			ret.emplace_back(mipSamp);
		}
	}
	return ret;
}


bool SetUpMaster::checkForClusteredMipFamForSamp(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	return bfs::exists(bib::files::join(VecStr { sampDirMaster.clusDir_.string(),
		mipSampName.mipFam_, mipSampName.mipFam_ + "_clustered.fastq" }));
}

bool SetUpMaster::checkForRawDataForSamp(const MipFamSamp & mipSampName) const{
	SampleDirectoryMaster sampDirMaster(directoryMaster_, mipSampName);
	return bfs::exists(bib::files::join(VecStr { sampDirMaster.masterSampleDir_.string(),
		mipSampName.samp_ + rawDataSuffix_}));
}

VecStr SetUpMaster::getAllMipTargets() const {
	if (!names_) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	if (!mips_) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member mips_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}

	VecStr ret;
	for (const auto & mip : names_->mips_) {
		addOtherVec(ret, mips_->getMipsForFamily(mip));
	}
	return ret;
}

VecStr SetUpMaster::getAllMipFamilies() const{
	if (!names_) {
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error { ss.str() };
	}
	return names_->mips_;
}

VecStr SetUpMaster::getMipGroupings() const {
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	std::set<std::string> groupings;
	for(const auto & mip : names_->mips_){
		groupings.emplace(getGroupForMipFam(mip));
	}
	return VecStr{groupings.begin(), groupings.end()};
}

VecStr SetUpMaster::getMipFamiliesForMipGroup(const std::string & groupName) const{
	if(!names_){
		std::stringstream ss;
		ss << "Error in " << __PRETTY_FUNCTION__
				<< " attempting to use member names_ when it isn't a valid pointer"
				<< std::endl;
		throw std::runtime_error{ss.str()};
	}
	VecStr mipFamilies;
	auto prefix = groupName + "_";
	for(const auto & mipFam : names_->mips_){
		if(bib::beginsWith(mipFam, prefix)){
			mipFamilies.emplace_back(mipFam);
		}
	}
	return mipFamilies;
}

std::string SetUpMaster::getGroupForMipFam(const std::string & mipFamName)const{
	return mipFamName.substr(0,mipFamName.find("_"));
}

bfs::path SetUpMaster::pathMipPopClusSampInfo(
		const MipFamSamp & mipSampName) const {
	return bib::files::make_path(directoryMaster_.populationClusteringDir_,
			mipSampName.mipFam_, "analysis", "selectedClustersInfo.tab.txt");
}

bfs::path SetUpMaster::pathSampPopClusSampInfo(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(getMipSerDir(),"popClusInfo",
			mipSampName.samp_ + ".tab.txt");
}

bfs::path SetUpMaster::pathMipExtractInfo(const std::string & mipTar) const{
	return bib::files::make_path(getMipSerDir(),"extractionInfo",
			mipTar + ".tab.txt");
}

bfs::path SetUpMaster::pathMipPopClusPopInfo(
		const MipFamSamp & mipSampName) const {
	return bib::files::make_path(directoryMaster_.populationClusteringDir_,
			mipSampName.mipFam_, "analysis", "population",
			"populationCluster.tab.txt");
}

bfs::path SetUpMaster::pathMipPopClusHaplo(
		const MipFamSamp & mipSampName) const {
	/*return bib::files::make_path(directoryMaster_.populationClusteringDir_,
			mipSampName.mipFam_, "analysis", "population",
			mipSampName.mipFam_ + ".fastq");*/
	return bib::files::make_path(directoryMaster_.populationClusteringDir_,
				mipSampName.mipFam_, "analysis", "population",
				"PopSeqs.fastq");
}

bfs::path SetUpMaster::pathSampleExtractInfo(const MipFamSamp & mipSampName) const {
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
			mipSampName.samp_ + "_mipExtraction", "extractInfoByTarget.txt");
}

bfs::path SetUpMaster::pathSampleRawData(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
				mipSampName.samp_ + rawDataSuffix_);
}

bfs::path SetUpMaster::pathPopClusFinalHaplo(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.populationClusteringDir_,
				mipSampName.mipFam_, "analysis", "samplesOutput",mipSampName.samp_, "final",
				mipSampName.samp_ + ".fastq");
}
/*
bfs::path SetUpMaster::pathPopClusOriginalHaplo(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.populationClusteringDir_,
				mipSampName.mipFam_, "analysis", "originals",
				mipSampName.samp_ + ".fastq");
}*/

bfs::path SetUpMaster::pathMipSampBarCorHap(
		const MipFamSamp & mipSampName) const {
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
			mipSampName.samp_ + "_mipBarcodeCorrection", mipSampName.mipFam_,
			mipSampName.mipFam_ + "_all.fastq");
}
bfs::path SetUpMaster::pathMipSampBarCorBars(
		const MipFamSamp & mipSampName, const std::string & mipTarName) const {
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
			mipSampName.samp_ + "_mipBarcodeCorrection", mipSampName.mipFam_,
			mipTarName + "_barcodes.tab.txt");
}

//

bfs::path SetUpMaster::pathMipSampClusDir(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
					mipSampName.samp_ + "_mipClustering");
}

bfs::path SetUpMaster::pathMipSampBarCorDir(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_,
				mipSampName.samp_ + "_mipBarcodeCorrection");
}

bfs::path SetUpMaster::pathSampDir(const MipFamSamp & mipSampName) const{
	return bib::files::make_path(directoryMaster_.masterDir_, mipSampName.samp_);
}

bfs::path SetUpMaster::pathMipSampExtractDir(const MipTarFamSamp & mipTarSampName) const{
	return bib::files::make_path(directoryMaster_.masterDir_, mipTarSampName.mipFamSamp_.samp_, mipTarSampName.mipTar_);
}

void SetUpMaster::prepareMipAnalysisServer(uint32_t numThreads) const{
	auto warnings = checkDirStruct();
	if (!warnings.empty()) {
		std::stringstream ss;
		ss
				<< "Error in directory structure, make sure you are in the correct analysis directory"
				<< std::endl;
		ss << "Following warnings;" << std::endl;
		ss << bib::conToStr(warnings, "\n") << std::endl;
		throw std::runtime_error { ss.str() };
	}
	auto samplesExtracted = getSamplesWithRawData(numThreads);
	bib::files::makeDirP(
			bib::files::MkdirPar(
					getMipSerDir().string() + "extractionInfo"));
	bib::files::makeDirP(
			bib::files::MkdirPar(getMipSerDir().string() + "popClusInfo"));
	//extraction
	{
		auto extractionMasterTab = gatherExtractStats(
				samplesExtracted,numThreads);
		if (!extractionMasterTab.content_.empty()) {
			TableIOOpts allExtractInfoOpts(InOptions(), "\t",
					OutOptions(
							getMipSerDir().string()
									+ "extractionInfo/allExtractInfo.tab.txt"), "\t", true);
			allExtractInfoOpts.out_.overWriteFile_ = true;
			extractionMasterTab.outPutContents(allExtractInfoOpts);
			std::vector<bfs::path> extractInfoFilepaths;
			for (const auto & samp : samplesExtracted) {
				extractInfoFilepaths.emplace_back(
						pathSampleExtractInfo(samp));
			}

			MasterTableStaticCache allExtractInfo(
					TableIOOpts(InOptions(), "\t", OutOptions(), "\t", true),
					extractInfoFilepaths);
			auto byTarget = allExtractInfo.get().splitTableOnColumn("mipTarget");
			auto allTargets = mips_->getMipsForFamily(names_->mips_);
			for (auto & tar : byTarget) {
				if(bib::in(tar.first, allTargets)){
					tar.second.trimElementsAtFirstOccurenceOf("(");
					TableIOOpts tarOpts(
							OutOptions(
									getMipSerDir().string() + "extractionInfo/"
											+ tar.first + ".tab.txt"), "\t", true);
					tarOpts.out_.overWriteFile_ = true;
					tar.second.outPutContents(tarOpts);
				}
			}
		}
	}
	//pop clustering
	{
		auto mipsPopClustered = getMipFamsWithPopClustered(
				numThreads);
		std::vector<bfs::path> popClusInfofilepaths;
		for (const auto & mipFam : mipsPopClustered) {
			popClusInfofilepaths.emplace_back(
					pathMipPopClusSampInfo(mipFam));
		}
		TableIOOpts allPopInfoOpts(InOptions(), "\t",
				OutOptions(
						getMipSerDir().string() + "popClusInfo/allInfo.tab.txt"),
				"\t", true);
		allPopInfoOpts.out_.overWriteFile_ = true;
		MasterTableStaticCache allPopInfo(allPopInfoOpts, popClusInfofilepaths);
		auto bySample = allPopInfo.get().splitTableOnColumn("s_Sample");
		for (auto & samp : bySample) {
			if(bib::in(samp.first, names_->samples_)){
				TableIOOpts sampOpts(
						OutOptions(
								getMipSerDir().string() + "popClusInfo/" + samp.first
										+ ".tab.txt"), "\t", true);
				sampOpts.out_.overWriteFile_ = true;
				samp.second.sortTable("p_geneName", "p_targetName", "c_clusterID", false);
				samp.second.outPutContents(sampOpts);
			}
		}
		allPopInfo.writeTab();
	}
}


table SetUpMaster::gatherExtractStats(const std::vector<MipFamSamp> & samplesExtracted, uint32_t numThreads)const {

	TableIOOpts allExtractInfoOpts(InOptions(), "\t",
			OutOptions(
					getMipSerDir().string()
							+ "extractionInfo/allExtractInfo.tab.txt"), "\t", true);
	allExtractInfoOpts.out_.overWriteFile_ = true;
	bool needsUpdate = false;
	if (!bfs::exists(allExtractInfoOpts.out_.outFilename_)) {
		needsUpdate = true;
	} else {
		auto extractStatTime = bib::files::last_write_time(
				allExtractInfoOpts.out_.outFilename_);
		for (const auto & samp : samplesExtracted) {
			if (bib::files::last_write_time(pathSampleExtractInfo(samp))
					> extractStatTime) {
				needsUpdate = true;
				break;
			}
		}
	}
	table extractionMasterTab(VecStr { "Sample", "raw", "assembled", "discarded",
			"unassembled", "matchingExtArm", "UnmatachedExtArm",
			"readsFailing_LigationArm", "readsFailing_Minlen", "readsFailing_Quality",
			"smallFragment" });
	if (needsUpdate) {
		bib::concurrent::LockableQueue<MipFamSamp> samplesExtractQueue(
				samplesExtracted);
		std::mutex extractionMasterTabMut;
		auto gatherSampExtractInfo =
				[&extractionMasterTabMut,&extractionMasterTab](bib::concurrent::LockableQueue<MipFamSamp> & samplesExtractQueue,
						const SetUpMaster & mipMaster) {
					MipFamSamp samp("", "");
					while(samplesExtractQueue.getVal(samp)) {
						auto filePath = mipMaster.pathSampleRawData(samp);
						if (bfs::exists(filePath)) {
							auto assembled = countSeqs(
									SeqIOOptions(filePath.string(),
											SeqIOOptions::getInFormat(
													bib::files::getExtension(filePath.string())), false), false);
							uint32_t discarded = countSeqs(
									SeqIOOptions(bib::files::make_path(mipMaster.directoryMaster_.masterDir_,samp.samp_,samp.samp_ + ".discarded.fastq").string(),
											SeqIOOptions::inFormats::FASTQ, false), false);
							uint32_t unassembled = countSeqs(
									SeqIOOptions(bib::files::make_path(mipMaster.directoryMaster_.masterDir_,samp.samp_,samp.samp_ + ".notCombined_1.fastq").string(),
											SeqIOOptions::inFormats::FASTQ, false), false);
							uint32_t raw = assembled + discarded + unassembled;
							table extractTab(mipMaster.pathSampleExtractInfo(samp).string(), "\t", true);
							VecStr rowsNeeded {"unmatched", "smallFragment", "total"};
							extractTab = extractTab.extractByComp("mipTarget",
									[&rowsNeeded](const std::string & row) {
										return bib::in(row, rowsNeeded);
									});
							extractTab.trimElementsAtFirstOccurenceOf("(");
							uint32_t matchingExtArm = 0;
							uint32_t unmatchedExtArm = 0;
							uint32_t readsFailingLigationArm = 0;
							uint32_t raedsFailingMinLen = 0;
							uint32_t readsFailingQuality = 0;
							uint32_t readsSmallFragment = 0;
							for(const auto & row : extractTab.content_) {
								if("unmatched" == row[extractTab.getColPos("mipTarget")]) {
									unmatchedExtArm = bib::lexical_cast<uint32_t>(row[extractTab.getColPos("readNumber")]);
								} else if("smallFragment" == row[extractTab.getColPos("mipTarget")]) {
									readsSmallFragment = bib::lexical_cast<uint32_t>(row[extractTab.getColPos("readNumber")]);
								} else if("total" == row[extractTab.getColPos("mipTarget")]) {
									matchingExtArm = bib::lexical_cast<uint32_t>(row[extractTab.getColPos("goodReads")]);
									readsFailingLigationArm = bib::lexical_cast<uint32_t>(row[extractTab.getColPos("failedLigationArm")]);
									raedsFailingMinLen = bib::lexical_cast<uint32_t>(row[6]);
									readsFailingQuality = bib::lexical_cast<uint32_t>(row[7]);
								}
							}
							{
								std::lock_guard<std::mutex> lock(extractionMasterTabMut);
								extractionMasterTab.content_.emplace_back(toVecStr(samp.samp_,raw, assembled, discarded,
												unassembled, matchingExtArm, unmatchedExtArm, readsFailingLigationArm,
												raedsFailingMinLen,readsFailingQuality, readsSmallFragment));
							}
						}
					}
				};
		std::vector<std::thread> threads;
		for (uint32_t tNum = 0; tNum < numThreads; ++tNum) {
			threads.emplace_back(
					std::thread(gatherSampExtractInfo, std::ref(samplesExtractQueue),
							std::cref(*this)));
		}
		for (auto & t : threads) {
			t.join();
		}
		extractionMasterTab.sortTable("Sample", false);
	}
	return extractionMasterTab;
}



}  // namespace bibseq
