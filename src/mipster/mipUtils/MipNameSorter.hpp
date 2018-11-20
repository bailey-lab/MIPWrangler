#pragma once
/*
 * MipNameSorter.hpp
 *
 *  Created on: Feb 21, 2017
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

namespace njhseq {

class MipNameSorter{
public:
	static bool compareNames(const std::string & name1,
			const std::string & name2, const std::regex & namePat);

	static bool compareNames(const std::string & name1,
			const std::string & name2,
			const std::regex & namePat,
			const std::regex & secondaryPat);


	static const std::regex mipNamePat;
	static const std::regex regionNamePat;

	template<typename T>
	static void sort(std::vector<T> & seqs);

	template<typename T>
	static void sortByRegion(std::vector<T> & seqs);

	template<typename T>
	static void sort(std::vector<T> & seqs,  std::function<std::string(const T&)> nameFunc);

	static void sort(VecStr & names, const std::regex & namePat);
	static void sort(VecStr & names, const std::regex & namePat, const std::regex & secondaryPat);

	template<typename T>
	static void sort(std::vector<T> & seqs, const std::regex & namePat);
	template<typename T>
	static void sort(std::vector<T> & seqs, const std::regex & namePat, const std::regex & secondaryPat);

	template<typename T>
	static void sort(std::vector<T> & input, const std::regex & namePat, std::function<std::string(const T&)> nameFunc);
	template<typename T>
	static void sort(std::vector<T> & input, const std::regex & namePat, const std::regex & secondaryPat, std::function<std::string(const T&)> nameFunc);

};


template<typename T>
void MipNameSorter::sortByRegion(std::vector<T> & seqs){
	sort(seqs,regionNamePat);
}

template<typename T>
void MipNameSorter::sort(std::vector<T> & seqs){

	sort(seqs, mipNamePat, regionNamePat);
}

template<typename T>
void MipNameSorter::sort(std::vector<T> & seqs,  std::function<std::string(const T&)> nameFunc){
	sort(seqs, mipNamePat, regionNamePat, nameFunc);
}


template<typename T>
void MipNameSorter::sort(std::vector<T> & seqs, const std::regex & namePat){
	if (namePat.mark_count() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error regex pattern must have at least 2 sub patterns capture to sort by"
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	njh::sort(seqs, [&namePat](const T & seq1,
			const T & seq2) {
		return compareNames(getSeqBase(seq1).name_, getSeqBase(seq2).name_, namePat);
	});
}


template<typename T>
void MipNameSorter::sort(std::vector<T> & seqs, const std::regex & namePat,
		const std::regex & secondaryPat){
	if(namePat.mark_count() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error regex pattern, namePat, must have at least 2 sub patterns capture to sort by, found less: " << namePat.mark_count() << "\n";
		throw std::runtime_error{ss.str()};
	}
	if(secondaryPat.mark_count() < 2){
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << ", error regex pattern, secondaryPat, must have at least 2 sub patterns capture to sort by, found less: " << secondaryPat.mark_count() << "\n";
		throw std::runtime_error{ss.str()};
	}

	njh::sort(seqs, [&namePat,&secondaryPat](const T & seq1,
			const T & seq2){
		return compareNames(getSeqBase(seq1).name_, getSeqBase(seq2).name_, namePat, secondaryPat);
	});
}

template<typename T>
void MipNameSorter::sort(std::vector<T> & input, const std::regex & namePat,
		const std::regex & secondaryPat,
		std::function<std::string(const T&)> nameFunc) {
	if (namePat.mark_count() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error regex pattern, namePat, must have at least 2 sub patterns capture to sort by, found less: "
				<< namePat.mark_count() << "\n";
		throw std::runtime_error { ss.str() };
	}
	if (secondaryPat.mark_count() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error regex pattern, secondaryPat, must have at least 2 sub patterns capture to sort by, found less: "
				<< secondaryPat.mark_count() << "\n";
		throw std::runtime_error { ss.str() };
	}
	njh::sort(input, [&namePat,&secondaryPat,&nameFunc](const T & seq1,
			const T & seq2) {
		return compareNames(nameFunc(seq1), nameFunc(seq2), namePat, secondaryPat);
	});
}

template<typename T>
void MipNameSorter::sort(std::vector<T> & input,
		const std::regex & namePat,
		std::function<std::string(const T&)> nameFunc){
	if (namePat.mark_count() < 2) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__
				<< ", error regex pattern must have at least 2 sub patterns capture to sort by"
				<< "\n";
		throw std::runtime_error { ss.str() };
	}
	njh::sort(input, [&namePat,&nameFunc](const T & input1,
			const T & input2) {
		return compareNames(nameFunc(input1), nameFunc(input2), namePat);
	});
}


}  // namespace njhseq




