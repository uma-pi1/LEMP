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
 * TAState.h
 *
 *  Created on: Sep 18, 2012
 *      Author: chteflio
 */

#ifndef TASTATE_H_
#define TASTATE_H_

#include "Lists.h"


namespace ta {

    class TAState {        
        
    public:
        const double* query;
        QueueElementLists* invLists;
        boost::dynamic_bitset<> exploredItems; 
        col_type colNum;
        row_type rowNum;
        std::vector<long> fringePos;
        bool allSeen; // true if you should stop scanning


        TAState(col_type colNum) : colNum(colNum), invLists(nullptr) {
            fringePos.resize(colNum);
        }

        ///////////////// virtual methods ////////////////////////

        inline virtual void initForNewQuery(const double* q) {
            query = q;
            exploredItems.reset();
            allSeen = false;
        }

        /*
         * given that the step happened on col update the structures
         * returns the value of the next non-explored entry in this column
         */
        inline virtual void updateState(col_type col) {
            if (!allSeen) {
                row_type rowId;

                if (query[col] < 0) {//scan downwards //////////////////
                    // keep moving downwards
                    do {
                        fringePos[col]++;
                        if (fringePos[col] == (col + 1) * invLists->size) {
                            allSeen = true;
                            break;
                        }

                        rowId = invLists->getElement(fringePos[col])->id;
                    } while (exploredItems[rowId]);
                } else {//scan upwards
                    // keep moving upwards
                    do {
                        fringePos[col]--;
                        if (fringePos[col] < col * invLists->size) {
                            allSeen = true;
                            break;
                        }
                        rowId = invLists->getElement(fringePos[col])->id;
                    } while (exploredItems[rowId]);
                }
            }
        }

        inline virtual col_type chooseStep() {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        ///////////////// concrete methods ////////////////////////

        inline void initializeForNewBucket(QueueElementLists* m1) {// first initialization
            invLists = m1;
            rowNum = m1->getRowNum();
            exploredItems.resize(rowNum);
        }

        inline double initializeThreshold() {
            double stopThreshold = 0;

            for (col_type i = 0; i < colNum; ++i) {
                stopThreshold += invLists->getElement(fringePos[i])->data * query[i];
            }

            return stopThreshold;
        }

        inline bool isThresholdUnterTheta(double& stopThreshold, double localTheta, col_type stepOnCol, double oldValue, bool forCosine) {
            double piNew = invLists->getElement(fringePos[stepOnCol])->data;

            stopThreshold += query[stepOnCol] * (piNew - oldValue);

            if (stopThreshold < localTheta) {
                return true;
            }

            return false;
        }

        inline void maintainExploredLists(const VectorMatrix& matrix, row_type countSeen, col_type stepOnCol) {
            if (countSeen == rowNum) {
                allSeen = true;
            } else {
                row_type rowId = invLists->getElement(fringePos[stepOnCol])->id;
                exploredItems[rowId] = true;
            }
        }
    };

    class TAStateMAX : public TAState {
    public:

        /*
         * this is for the scheduler
         */
        maxHeap piqi;

        TAStateMAX(col_type colNum) : TAState(colNum) {
        };

        inline void initForNewQuery(const double* q) {// initialize to specific rows
            TAState::initForNewQuery(q);
            piqi.clear();

            for (col_type i = 0; i < colNum; ++i) {
                if (query[i] < 0) { //scan downwards///////////////////////
                    fringePos[i] = i * invLists->size + 0;
                } else {
                    fringePos[i] = i * invLists->size + rowNum - 1;
                }
                if (query[i] != 0)
                    piqi.add(QueueElement(invLists->getElement(fringePos[i])->data * query[i], i));
            }
        }

        /*
         * given that the step happened on col update the structures
         * returns the value of the next non-explored entry in this column
         */
        inline void updateState(col_type col) {
            TAState::updateState(col);
            if (!allSeen) {// update the front item                                
                piqi.updateRoot(invLists->getElement(fringePos[col])->data * query[col]);
            }
        }

        inline col_type chooseStep() {
            col_type stepOnCol = piqi.getMaxId();
            return stepOnCol;
        }

    };

    class TAStateRR : public TAState {
    public:

        /*
         * this is for the scheduler
         */
        std::vector<col_type> activeCols;
        col_type colIndx;

        TAStateRR(col_type colNum) : TAState(colNum), colIndx(0) {
        };

        inline void initForNewQuery(const double* q) {// initialize to specific rows
            TAState::initForNewQuery(q);
            colIndx = 0;            
            activeCols.clear();

            for (col_type i = 0; i < colNum; ++i) {
                if (query[i] < 0) { //scan downwards///////////////////////
                    fringePos[i] = i * invLists->size + 0;
                } else {
                    fringePos[i] = i * invLists->size + rowNum - 1;
                }
                if (query[i] != 0) {
                    activeCols.push_back(i);
                }
            }
        }

        /*
         * given that the step happened on col update the structures
         * returns the value of the next non-explored entry in this column
         */
        inline void updateState(col_type col) {
            TAState::updateState(col);
        }

        inline col_type chooseStep() {
            if (colIndx == activeCols.size())
                colIndx = 0;

            col_type stepOnCol = activeCols[colIndx];
            colIndx++;
            return stepOnCol;
        }

    };
}



#endif /* TASTATE_H_ */
