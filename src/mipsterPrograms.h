#pragma once

// Created on 2014/12/30
// Including headers in programs
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

#include "mipsterPrograms/mipsterAnalysis.h"
#include "mipsterPrograms/mipsterServer.h"
#include "mipsterPrograms/mipsterUtils.h"
#include "mipsterPrograms/mipsterSim.h"
#include "mipsterPrograms/mipsterSetUp.h"
#include "mipsterPrograms/mipsterMipExplorer.h"
#include "mipsterPrograms/mipsterMipTester.h"


namespace njhseq {

class mipsterRunner: public njh::progutils::OneRing {
public:
	mipsterRunner();
};
mipsterRunner::mipsterRunner() :
		njh::progutils::OneRing(
				{ addRing<mipsterAnalysisRunner>(),
					addRing<mipsterServerRunner>(),
					addRing<mipsterUtilsRunner>(),
					addRing<mipsterSetUpRunner>(),
				  addRing<mipsterSimRunner>(),
					addRing<mipsterMipExplorerRunner>(),
					addRing<mipsterMipTesterRunner>()}, { }, "MIPWrangler",
					"1", "0", "1-dev") {
}
} // namespace njhseq

