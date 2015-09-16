// Copyright 2015 Christina Teflioudi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
/* 
 * File:   TANRAState.h
 * Author: chteflio
 *
 * Created on August 24, 2015, 9:10 AM
 */

#ifndef TANRASTATE_H
#define	TANRASTATE_H

namespace ta {

    struct TANRAState {
        // fields for the 1-sided version
        const double* query;
        QueueElementLists* invLists;

        unordered_map<row_type, std::pair<double, std::vector<col_type> > > worstScores;

        col_type colNum;
        row_type rowNum;
        std::vector<long> fringePos;
        col_type inactiveCols;

        /*
         * this is for the scheduler
         */
        maxHeap piqi;

        TANRAState(col_type colNum) : colNum(colNum), invLists(0) {
            fringePos.resize(colNum);
        }

        inline void initializeForNewBucket(QueueElementLists* m1) {// first initialization
            invLists = m1;
            rowNum = m1->getRowNum();
        }

        inline void initForNewQuery(const double* q) {// initialize to specific rows

            inactiveCols = 0;
            worstScores.erase(worstScores.begin(), worstScores.end());

            query = q;
            piqi.clear();
            double value;


            for (col_type i = 0; i < colNum; i++) {
                if (query[i] < 0) { //scan downwards///////////////////////
                    fringePos[i] = i * invLists->size + 0;

                } else {
                    fringePos[i] = i * invLists->size + rowNum - 1;
                }
                value = invLists->getElement(fringePos[i])->data * query[i];
                piqi.add(QueueElement(value, i));
            }
        }

        /*
         * given that the step happened on col update the structures
         * returns the value of the next non-explored entry in this column
         */

        inline void updateScores(row_type key, col_type stepOnCol) {
            worstScores[key].first += invLists->getElement(fringePos[stepOnCol])->data * query[stepOnCol];
            worstScores[key].second.push_back(stepOnCol);
        }

        inline void getCandidates(std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify, double localTheta, double stopThreshold) {
            unordered_map<row_type, std::pair<double, std::vector<col_type> > >::iterator it;

            for (it = worstScores.begin(); it != worstScores.end(); it++) {
                row_type key = it->first;
                double ws = it->second.first;

                if (ws >= localTheta) {
                    candidatesToVerify[numCandidatesToVerify] = key;
                    numCandidatesToVerify++;
                    continue;
                }
                ws += stopThreshold;

                for (int i = 0; i < it->second.second.size(); i++) {
                    col_type col = it->second.second[i];
                    ws -= invLists->getElement(fringePos[col])->data * query[col];
                }

                if (ws >= localTheta) {
                    candidatesToVerify[numCandidatesToVerify] = key;
                    numCandidatesToVerify++;
                }
            }

        }

        inline void updateState(col_type col) {// given that the step happened on col update the structures

            row_type rowId;
            bool inactivate = false;

            if (query[col] < 0) {//scan downwards //////////////////

                fringePos[col]++;
                if (fringePos[col] == (col + 1) * invLists->size) {
                    inactivate = true;
                } else {
                    rowId = invLists->getElement(fringePos[col])->id;
                }

            } else {//scan upwards

                // keep moving upwards
                fringePos[col]--;
                if (fringePos[col] < col * invLists->size) {
                    inactivate = true;
                } else {
                    rowId = invLists->getElement(fringePos[col])->id;
                }
            }

            if (inactivate) {
                piqi.updateRoot((-1) * std::numeric_limits<double>::max());
                inactiveCols++;
            } else {
                piqi.updateRoot(invLists->getElement(fringePos[col])->data * query[col]);
            }
        }

   

        ///////////////// methods for threshold update ////////////////////////

        inline double initializeThreshold() {
            double stopThreshold = 0;
            for (col_type i = 0; i < colNum; i++) {
                stopThreshold += invLists->getElement(fringePos[i])->data * query[i];
            }
            return stopThreshold;
        }

        inline void updateThreshold(double& stopThreshold, col_type stepOnCol, double oldValue) {
            double queryPart = query[stepOnCol];
            double piNew = invLists->getElement(fringePos[stepOnCol])->data;
            stopThreshold += queryPart * (piNew - oldValue);
        }

        inline bool isThresholdUnterTheta(double& stopThreshold, double localTheta, col_type stepOnCol, double oldValue, bool forCosine) {
            if (colNum <= inactiveCols) {
                return true;
            }

            updateThreshold(stopThreshold, stepOnCol, oldValue);

            if (stopThreshold < localTheta) {
                return true;
            }

            if (forCosine) {
                if (stopThreshold > 1 && 1 < localTheta) {
                    return true;
                }
            }

            return false;
        }

        ///////////////// methods for choosing next step ////////////////////////

        inline col_type chooseStepGreatestPiQi() {
            col_type stepOnCol = piqi.getMaxId();
            return stepOnCol;
        }

        inline col_type chooseStep() {
            col_type stepOnCol = chooseStepGreatestPiQi();
            return stepOnCol;
        }
        ////////////////////////////////////////


    };
}

#endif	/* TANRASTATE_H */

