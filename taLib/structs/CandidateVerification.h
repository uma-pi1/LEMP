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
 * CandidateVerification.h
 *
 *  Created on: May 12, 2014
 *      Author: chteflio
 */

#ifndef CANDIDATEVERIFICATION_H_
#define CANDIDATEVERIFICATION_H_


#include <iterator>

#include "RetrievalArguments.h" 


namespace ta {

    inline void verifyCandidates_lengthTest(const double * query, row_type numCandidatesToVerify, RetrievalArguments* arg) {
        std::pair<bool, double> p;

        for (row_type i = 0; i < numCandidatesToVerify; ++i) {
            row_type row = arg->candidatesToVerify[i];
            p = arg->probeMatrix->passesThreshold(row, query, arg->theta);

            if (p.first) {
                arg->results.emplace_back(p.second, arg->queryId, arg->probeMatrix->getId(row));
            }
        }
        arg->comparisons += numCandidatesToVerify;

    }

    inline void verifyCandidates_noLengthTest(const double * query, row_type numCandidatesToVerify,  RetrievalArguments* arg) {       

        for (row_type i = 0; i < numCandidatesToVerify; ++i) {
            row_type row = arg->candidatesToVerify[i];
            double ip = arg->probeMatrix->innerProduct(row, query);

            if (ip >= arg->theta) {               
                 arg->results.emplace_back(ip, arg->queryId, arg->probeMatrix->getId(row));
            }
        }
        arg->comparisons += numCandidatesToVerify;

    }

    inline void verifyCandidatesTopK_noLengthTest(const double * query,row_type numCandidatesToVerify, RetrievalArguments* arg) {

        double minScore = arg->heap.front().data;

        for (row_type i = 0; i < numCandidatesToVerify; ++i) {
            row_type row = arg->candidatesToVerify[i];
            double ip = arg->probeMatrix->innerProduct(row, query);

            if (ip > minScore) {
                std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                arg->heap.pop_back();
                arg->heap.emplace_back(ip,  arg->probeMatrix->getId(row));
                std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                minScore = arg->heap.front().data;
            }
        }
        arg->comparisons += numCandidatesToVerify;

    }

    inline void verifyCandidatesTopK_lengthTest(const double * query, row_type numCandidatesToVerify,  RetrievalArguments* arg) {


        double minScore = arg->heap.front().data;

        for (row_type i = 0; i < numCandidatesToVerify; ++i) {
            row_type row = arg->candidatesToVerify[i];

            if (arg->probeMatrix->getVectorLength(row) <= minScore)
                continue;

            double ip = arg->probeMatrix->innerProduct(row, query);
            arg->comparisons++;

            if (ip > minScore) {
                std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                arg->heap.pop_back();
                arg->heap.emplace_back(ip, arg->probeMatrix->getId(row));
                std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                minScore = arg->heap.front().data;
            }
        }
    }


    // These are used by TA
    // examines if the item should be included in the result and does the corresponding housekeeping

    inline void verifyCandidate(row_type posMatrix, const double* query, RetrievalArguments* arg) {
        arg->comparisons++;

        std::pair<bool, double> p;
        p = arg->probeMatrix->passesThreshold(posMatrix, query, arg->theta);

        if (p.first) {
             arg->results.emplace_back(p.second, arg->queryId, arg->probeMatrix->getId(posMatrix));
        }
    }

    inline void verifyCandidateTopk(row_type posMatrix, const double* query, RetrievalArguments* arg) {
        std::pair<bool, double> p;
        arg->comparisons++;       
        p = arg->probeMatrix->passesThreshold(posMatrix, query, arg->heap.front().data);

        if (p.first) {
            // remove min element from the heap
            pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>()); // Yes! I need to use greater to get a min heap!
            arg->heap.pop_back();
            // and push new element inside the heap
            arg->heap.emplace_back(p.second,  arg->probeMatrix->getId(posMatrix));
            push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());

        }

    }



}


#endif /* CANDIDATEVERIFICATION_H_ */
