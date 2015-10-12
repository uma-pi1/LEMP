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
 * QueryBucket_withTuning2.h
 *
 *  Created on: Mar 25, 2014
 *      Author: chteflio
 */

#ifndef QUERYBATCH_H_
#define QUERYBATCH_H_

namespace ta {

    double matrixToMatrixTime2 = 0;


    class QueryBucket_withTuning {
        col_type* queues;


    public:

        LshIndex * lshIndex;
        row_type startPos, endPos, rowNum;
        row_type inactiveCounter;
        bool initializedQueues;
        std::vector<bool> inactiveQueries;
        std::pair<double, double> normL2; // first: min second: max
        

        inline void preprocess(const VectorMatrix& userMatrix, col_type maxLists);

        inline col_type* getQueue(row_type user, col_type maxLists) const {
            return &queues[user * maxLists];
        }
  

        inline void createLshIndex(const VectorMatrix& matrix) {
            lshIndex = new LshIndex();
            lshIndex->initializeLists(matrix, false, startPos, endPos);
        }


        inline QueryBucket_withTuning() : initializedQueues(false), queues(nullptr),
        inactiveCounter(0),  lshIndex(nullptr){
        };

        inline ~QueryBucket_withTuning() {
            if (queues != nullptr)
                delete[] queues;

            if (lshIndex != nullptr)
                delete lshIndex;
        }

        inline void init(const VectorMatrix& matrix, row_type startInd, row_type endInd, LEMPArg& args) {
            startPos = startInd;
            endPos = endInd;
            rowNum = endPos - startPos;
            normL2.second = matrix.getVectorLength(startPos);
            normL2.first = matrix.getVectorLength(endPos - 1);

            if (args.k > 0) {
                inactiveQueries.resize(rowNum);
            }
        }

    };

    inline void QueryBucket_withTuning::preprocess(const VectorMatrix& userMatrix, col_type maxLists) {

        queues = new col_type[rowNum * maxLists];

        std::vector<QueueElement> tmp;
        tmp.resize(maxLists);

        for (row_type j = 0; j < rowNum; ++j) {
            const double* query = userMatrix.getMatrixRowPtr(j + startPos);

            for (col_type i = 0; i < maxLists; ++i) {
                double value = fabs(query[i]);

                tmp[i] = QueueElement(value, i);
            }
            std::make_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());

            for (col_type i = maxLists; i < userMatrix.colNum; ++i) {
                double value = fabs(query[i]);


                if (value > tmp.front().data) {
                    std::pop_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                    tmp.pop_back();
                    tmp.emplace_back(value,i);
                    std::push_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                }
            }

            std::sort(tmp.begin(), tmp.end(), std::greater<QueueElement>());
            for (col_type i = 0; i < maxLists; ++i) {
                queues[j * maxLists + i] = tmp[i].id;
            }
        }

    }

}

#endif /* QUERYBATCH_H_ */
