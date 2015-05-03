//    Copyright 2015 Christina Teflioudi
// 
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
// 
//        http://www.apache.org/licenses/LICENSE-2.0
// 
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

/*
 * Bucket.h
 *
 *  Created on: Nov 19, 2013
 *      Author: chteflio
 */

#ifndef PROBEBUCKET_H_
#define PROBEBUCKET_H_


#include <boost/unordered_set.hpp>
using boost::unordered_set;

namespace ta{

class Retriever;
typedef boost::shared_ptr<Retriever> retriever_ptr;


enum Index_Type{
	PLAIN = 0,
	SL = 1,
	INT_SL = 2,
	TREE = 3,
        AP = 4,
        BLSH = 5,
        LSH = 6

};


class ProbeBucket{


public:
	void* ptrIndexes[NUM_INDEXES];
	std::pair<double, double> normL2, invNormL2; // min and max length information
	double bucketScanThreshold, runtime, t_b; // all theta_b < t_b do LENGTH

	col_type colNum, numLists;
	row_type startPos, endPos, rowNum;
	retriever_ptr ptrRetriever;
        
        row_type activeQueries; // active queries for this bucket. If multiple threads this will just be an estimation



	rg::Timer bucketTimer;

	inline ProbeBucket(): numLists(1), t_b(1), runtime(0),  activeQueries(0){ ///////////////////
		ptrIndexes[PLAIN] = 0;
		ptrIndexes[SL] = 0;
		ptrIndexes[INT_SL] = 0;
		ptrIndexes[TREE] = 0;
                ptrIndexes[AP] = 0;
                ptrIndexes[BLSH] = 0;
                ptrIndexes[LSH] = 0;

	};

	inline ~ProbeBucket(){
		if(ptrIndexes[SL] != 0){
			delete static_cast<QueueElementLists*>(ptrIndexes[SL]);
		}
		if(ptrIndexes[INT_SL] != 0){
			delete static_cast<IntLists*>(ptrIndexes[INT_SL]);
		}
// 		if(ptrIndexes[TREE] != 0){
// 			delete static_cast<TreeIndex*>(ptrIndexes[TREE]);
// 		}
//                 if(ptrIndexes[AP] != 0){
// 			delete static_cast<L2apIndex*>(ptrIndexes[AP]);
// 		}
//                 if(ptrIndexes[BLSH] != 0){
// 			delete static_cast<BlshIndex*>(ptrIndexes[BLSH]);
// 		}
//                 if(ptrIndexes[LSH] != 0){
//                     delete static_cast<LshIndex*>(ptrIndexes[LSH]);
//                 }
	}

	inline  void init(const VectorMatrix& matrix, row_type startInd, row_type endInd, LEMPArg& args){
		startPos = startInd;
		endPos = endInd;
		rowNum = endPos - startPos;
		normL2.second = matrix.getVectorLength(startPos);
		normL2.first = matrix.getVectorLength(endPos-1);
                
                invNormL2.first = 1 / normL2.first;
                invNormL2.second = 1 / normL2.second;
                
		colNum = matrix.colNum;
	}

	inline void setAfterTuning(col_type lists, double thres){
		numLists = lists;
		t_b = thres;
	}

	inline bool hasIndex(Index_Type type){
		if (ptrIndexes[type] == 0)
			return false;
		else
			return true;
	}

	inline void* getIndex(Index_Type type){
		return ptrIndexes[type];
	}



};





}


#endif /* PROBEBUCKET_H_ */
