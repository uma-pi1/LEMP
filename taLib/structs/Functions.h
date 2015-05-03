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
 * Functions.h
 *
 *  Created on: Oct 17, 2013
 *      Author: chteflio
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_



namespace ta {

    inline bool compareResultsForUser(std::vector<QueueElement> & check, std::vector<QueueElement> & right, row_type user) {
        std::sort(check.begin(), check.end(), std::less<QueueElement>());
        std::sort(right.begin(), right.end(), std::less<QueueElement>());


        bool allOk = true;

        for (int i = 0; i < right.size(); i++) {
            bool thisOk = false;

            //            std::cout <<"i: "<<i<<" --> "<< r1[user][i].i << " " << r2[user][i].i << " --> " << ( r1[user][i].i ==  r2[user][i].i) << std::endl;

            for (int j = 0; j < check.size(); j++) {

                if (right[i].data == check[j].data || right[i].id == check[j].id) {
                    thisOk = true;
                    break;
                }
            }
            if (!thisOk) {
                allOk = false;
                break;
            }

        }

        if (!allOk) {
            std::cout << "NOT EQUAL at user: " << user << std::endl;
            for (int i = 0; i < check.size(); i++) {
                std::cout << std::setprecision(15) << check[i] << " " << std::setprecision(15) << right[i] << std::endl;
            }
            exit(1);
        }

        return allOk;

    }

//    inline bool compareResults(std::vector<std::vector<QueueElement> >& r1, std::vector<std::vector<QueueElement> >& r2) {
//        bool userOk = true;
//
//        for (int i = 0; i < r1.size(); i++) {
//            userOk = compareResultsForUser(r1, r2, i);
//        }
//
//        if (userOk) {
//            std::cout << "equal results" << std::endl;
//        } else {
//            std::cout << "NOT equal results" << std::endl;
//        }
//
//        return userOk;
//    }

    inline bool compareResults(std::vector<QueueElement> & r1, std::vector<QueueElement> & r2, int k) {
        bool userOk = true;
        int numUsers = r1.size() / k;

        std::vector<QueueElement> tmp1, tmp2;

        for (int i = 0; i < numUsers; i++) {

            tmp1.assign(r1.begin() + i*k, r1.begin() + (i + 1) * k);
            tmp2.assign(r2.begin() + i*k, r2.begin() + (i + 1) * k);

            userOk = compareResultsForUser(tmp1, tmp2, i);
        }

        if (userOk) {
            std::cout << "equal results" << std::endl;
        } else {
            std::cout << "NOT equal results" << std::endl;
        }

        return userOk;
    }

    inline bool compareResultsForUserRecall(std::vector<QueueElement>& check, std::vector<QueueElement>& right, row_type user) {
        std::sort(check.begin(), check.end(), std::less<QueueElement>());
        std::sort(right.begin(), right.end(), std::less<QueueElement>());
        comp_type nonSame = 0;

        for (int i = 0; i < right.size(); i++) {
            bool thisOk = false;

            for (int j = 0; j < check.size(); j++) {

                if (right[i].id == check[j].id || right[i].data == check[j].data) {
                    thisOk = true;
                    break;
                }
            }
            if (!thisOk) {
                nonSame++;
            }

        }

        return nonSame;

    }

    inline bool compareResultsForRecall(std::vector<QueueElement> & r1, std::vector<QueueElement> & r2, int k) {
        comp_type nonSame = 0;    
        int numUsers = r1.size() / k;

        std::vector<QueueElement> tmp1, tmp2;

        for (int i = 0; i < numUsers; i++) {

            tmp1.assign(r1.begin() + i*k, r1.begin() + (i + 1) * k);
            tmp2.assign(r2.begin() + i*k, r2.begin() + (i + 1) * k);

           
            nonSame += compareResultsForUserRecall(tmp1, tmp2, i);
        }
        
        
        
//        for (int i = 0; i < r1.size(); i++) {
//            nonSame += compareResultsForUserRecall(r1, r2, i);
//        }

        std::cout << "nonSame: " << nonSame << std::endl;

    }

    template<typename T>
    inline bool calculateIntervals(const double* query, const col_type* listsQueue, T& invLists, std::vector<IntervalElement>& intervals,
            double localTheta, col_type lists) {


        bool shouldScan = true;
        std::pair<row_type, row_type> necessaryIndices;

        for (col_type i = 0; i < lists; i++) {

            getBounds(query[listsQueue[i]], localTheta, invLists.getList(listsQueue[i]), necessaryIndices);

            intervals[i].col = listsQueue[i];
            intervals[i].start = necessaryIndices.first;
            intervals[i].end = necessaryIndices.second;

            if (intervals[i].end <= intervals[i].start) {
                shouldScan = false;
                break;
            }
        }

        std::sort(intervals.begin(), intervals.begin() + lists);

        return shouldScan;
    }

    inline void preprocess(const VectorMatrix* matrix, int id, col_type numLists, std::vector<QueueElement>& queue) {
        queue.resize(numLists);
        double* query = matrix->getMatrixRowPtr(id);

        for (col_type i = 0; i < numLists; i++) {
            double value = fabs(query[i]);
            queue[i] = QueueElement(value, i);
        }
        std::make_heap(queue.begin(), queue.end(), std::greater<QueueElement>());

        for (col_type i = numLists; i < matrix->colNum; i++) {
            double value = fabs(query[i]);

            if (value > queue.front().data) {
                std::pop_heap(queue.begin(), queue.end(), std::greater<QueueElement>());
                queue.pop_back();
                queue.push_back(QueueElement(value, i));
                std::push_heap(queue.begin(), queue.end(), std::greater<QueueElement>());
            }
        }

        std::sort(queue.begin(), queue.end(), std::greater<QueueElement>());

    }



}


#endif /* FUNCTIONS_H_ */
