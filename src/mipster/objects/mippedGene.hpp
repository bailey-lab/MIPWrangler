#pragma once
/*
 * mippedGene.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: nickhathaway
 */


#include "mipster/common.h"
#include <njhseq/objects/seqObjects/readObject.hpp>
#include <njhseq/alignment/aligner.h>
#include <njhseq/objects/BioDataObject/RefSeqGeneRecord.hpp>


namespace njhseq {



class mipRead : public readObject {

public:
	mipRead(const seqInfo & seqBase, bool processed, bool finalAnalysisNameFormat): readObject(seqBase, processed){
		if(finalAnalysisNameFormat){
			auto rPos = seqBase.name_.rfind("_R");
			auto bPos = seqBase.name_.rfind("_B");
			auto fPos = seqBase.name_.rfind("_f");
			readNumber_ = estd::stou(seqBase.name_.substr(rPos + 2, bPos - 2 - rPos));
			seqBase_.cnt_ = estd::stou(seqBase.name_.substr(bPos + 2, fPos - 2 - bPos));
			seqBase_.frac_ = njh::lexical_cast<double>(seqBase.name_.substr(fPos + 2));
		}else{
			auto rPos = seqBase.name_.rfind("_R");
			auto tPos = seqBase.name_.rfind("_t");
			readNumber_ = estd::stou(seqBase.name_.substr(rPos + 2, tPos - 2 - rPos));
		}


		auto firstUnderscore = seqBase.name_.find("_");
		auto firstPeriod = seqBase.name_.find(".");
		mipSampName_ = seqBase.name_.substr(0, firstUnderscore);
		mipTargetName_ = seqBase.name_.substr(firstUnderscore +1, firstPeriod - 1 - firstUnderscore);
		auto targetUnderscore = mipTargetName_.rfind("_");
		targetId_ = estd::stou(mipTargetName_.substr(targetUnderscore + 1));
		auto otherUnderScore = seqBase.name_.find("_", firstPeriod);
		clusterId_ = estd::stou(seqBase.name_.substr(firstPeriod + 1, otherUnderScore - 1 - firstPeriod));
	}

	uint32_t readNumber_;
	double readFraction_ = 0;

	std::string mipSampName_;
	std::string mipTargetName_;
	uint32_t targetId_;
	uint32_t clusterId_;

	void setReadFraction(double totalCnt){
		readFraction_ = readNumber_/totalCnt;
	}

};


class mipOverlap {
public:
	class edge;
	class node {
	public:
		node(const std::shared_ptr<mipRead> & read): value_(read){}
		std::shared_ptr<mipRead> value_;
		std::vector<std::shared_ptr<edge>> parentEdges_;
		std::vector<std::shared_ptr<edge>> childrenEdges_;

		uint32_t groupByParents_ = std::numeric_limits<uint32_t>::max();
		uint32_t groupByChildren_ = std::numeric_limits<uint32_t>::max();


		void printDownwardPaths(std::ostream & out, std::string currentPath)const{
			currentPath.append(value_->seqBase_.name_ + " "  + estd::to_string(value_->seqBase_.frac_));
			if(!childrenEdges_.empty()){
				currentPath.append(" -> ");
				for(const auto & child : childrenEdges_){
					auto cn = child->child_.lock();
					cn->printDownwardPaths(out, currentPath);
				}
			}else{
				out << currentPath << "\n";
			}
		}

		void printPossibleHap(std::ostream & out, std::string currentPath, std::string pathDelimiter = " -> ")const{
			currentPath.append("L" + estd::to_string(value_->targetId_+ 1) + ".A"  + estd::to_string(value_->clusterId_ + 1));
			if(!childrenEdges_.empty()){
				currentPath.append(pathDelimiter);
				for(const auto & child : childrenEdges_){
					auto cn = child->child_.lock();
					cn->printPossibleHap(out, currentPath, pathDelimiter);
				}
			}else{
				out << currentPath << "\n";
			}
		}

		VecStr getChildrenNames()const {
			VecStr ret;
			for(const auto & c : childrenEdges_){
				auto child = c->child_.lock();
				ret.emplace_back(child->value_->seqBase_.name_);
			}
			return ret;
		}
		VecStr getParentNames()const {
			VecStr ret;
			for(const auto & p : parentEdges_){
				auto par = p->par_.lock();
				ret.emplace_back(par->value_->seqBase_.name_);
			}
			return ret;
		}


	};

	class edge {
	public:
		/*edge(const std::shared_ptr<node> & node1, const std::shared_ptr<node> & node2){
		 	 con_[node1->value_->seqBase_.name_] = node2;
		 	 con_[node2->value_->seqBase_.name_] = node1;
		 }
		 std::unordered_map<std::string, std::weak_ptr<node>> con_;
		 */
		edge(const std::shared_ptr<node> & parent,
				const std::shared_ptr<node> & child) :
				par_(parent), child_(child) {
		}
		std::weak_ptr<node> par_;
		std::weak_ptr<node> child_;

	};
	std::vector<std::shared_ptr<node>> nodes_;
	std::unordered_map<std::string, std::shared_ptr<node>> nameToNode_;
	/**@b Connect two nodes, node 1 should be the parent of node 2
	 *
	 * @param nodeName1
	 * @param nodeName2
	 */
	void connectNodes(const std::string & parentNode, const std::string & childNode){
		auto con = std::make_shared<edge>(nameToNode_[parentNode], nameToNode_[childNode]);
		nameToNode_[parentNode]->childrenEdges_.emplace_back(con);
		nameToNode_[childNode]->parentEdges_.emplace_back(con);
	}

	void addNode(const std::shared_ptr<mipRead> & read){
		auto n = std::make_shared<node>(read);
		nameToNode_[read->seqBase_.name_] = n;
		nodes_.push_back(n);
	}





};


class mipTargetReads {
public:

	mipTargetReads(const std::string & targetName, bool reverseStrand):targetName_(targetName), reverseStrand_(reverseStrand){

	}

	std::vector<std::shared_ptr<mipRead>> reads_;
	std::string targetName_;
	uint32_t start_ = 0;
	uint32_t stop_ = 0;
	std::vector<uint32_t> allStarts_;
	std::vector<uint32_t> allStops_;
	bool reverseStrand_;


	void addRead(const std::shared_ptr<mipRead> & read, uint32_t start, uint32_t stop){
		reads_.emplace_back(read);
		addStartStop(start, stop);
	}

	void addStartStop(uint32_t start, uint32_t stop){
		allStarts_.emplace_back(start);
		allStops_.emplace_back(stop);
	}

	void setGeneralStartStop(){
		start_ = std::round(vectorMedianCopy(allStarts_));
		stop_ = std::round(vectorMedianCopy(allStops_));
	}

	void setFractionsOfReads(){
		double totalReadCnt = 0;
		double totalBarcodeCnt = 0;
		for(const auto & read : reads_){
			totalReadCnt += read->readNumber_;
			totalBarcodeCnt += read->seqBase_.cnt_;
		}
		for(auto & read : reads_){
			read->setFractionByCount(totalBarcodeCnt);
			read->setReadFraction(totalReadCnt);
		}
	}

	/**
	 *
	 * @param otherReads
	 * @param minOverlapRequired
	 * @return
	 *
	 */
	bool overLap(const mipTargetReads & otherReads, uint32_t minOverlapRequired = 1)const{
		if(otherReads.start_ > start_ && otherReads.start_ < stop_){
			uint32_t overlap = stop_ - otherReads.start_;
			if(overlap >=minOverlapRequired){
				return true;
			}else{
				return false;
			}
		}else if(start_ > otherReads.start_ && start_ < otherReads.stop_){
			uint32_t overlap = otherReads.stop_ - start_;
			if(overlap >=minOverlapRequired){
				return true;
			}else{
				return false;
			}
		}else {
			return false;
		}
	}

};


class mippedGene {
public:

	mippedGene(const seqInfo & seqBase,
			const std::unordered_map<std::string, mipTargetReads> & mReads);

	seqInfo genomicDna_;

	std::shared_ptr<RefSeqGeneRecord> refGeneRecord_ ;//= nullptr;

	std::unordered_map<std::string, mipTargetReads> processedReads_;

	mipOverlap overlapGraph_;

	seqInfo protein_;
	bool proteinCalculated = false;

	void setUpGraph(aligner & alignerObj, comparison comp);

	void printEdges(std::ostream & out)const ;

	void printAllPaths(std::ostream & out)const ;

	void printLociInfo(std::ostream & out)const ;

	void printPossibleHaps(std::ostream & out)const ;

	void setGroups();

	void printGroupInfo()const;

	std::vector<readObject> getAlignedTargets(aligner & alignerObj) ;

	std::unordered_map<std::string, table> findSnps(aligner & alignerObj);
	std::unordered_map<std::string, table> findProteinSnps(aligner & alignerObj);

	std::map<uint32_t, mismatch> getProteinSnps(aligner & alignerObj, const std::map<uint32_t, mismatch> & mismatches);

	void setProtein();

	seqInfo getMutatedProtein(const std::map<uint32_t, mismatch> & mismatches)const;

};


std::unordered_map<std::string, mipTargetReads> processMipReads(
		const std::vector<readObject> & reads, const readObject & genomicDna,
		aligner & alignerObj, bool processed, bool finalAnalysisNameFormat);

} /* namespace njhseq */
