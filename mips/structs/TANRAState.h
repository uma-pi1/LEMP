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

namespace mips {

    class TANRAState {
    protected:
        const double* query;
        QueueElementLists* invLists;

        std::unordered_map<row_type, std::pair<double, std::vector<col_type> > > worstScores;

        col_type colNum;
        row_type rowNum;
        std::vector<long> fringePos;
        col_type inactiveCols;

    public:

        TANRAState(col_type colNum) : colNum(colNum), invLists(0) {
            fringePos.resize(colNum);
        }

        virtual ~TANRAState() {
        }

        inline long getDepthInFringeInCol(col_type col) {
            return fringePos[col];
        }

        inline void initializeForNewBucket(QueueElementLists* m1) {// first initialization
            invLists = m1;
            rowNum = m1->getRowNum();
        }

        inline virtual void initForNewQuery(const double* q) {// initialize to specific rows
            inactiveCols = 0;
            query = q;
        }

        /*
         * given that the step happened on col update the structures
         * returns the value of the next non-explored entry in this column
         */

        inline void updateScores(row_type key, col_type stepOnCol) {
            worstScores[key].first += invLists->getElement(fringePos[stepOnCol])->data * query[stepOnCol];
            worstScores[key].second.push_back(stepOnCol);
        }

        inline void getCandidates(row_type* candidatesToVerify, row_type& numCandidatesToVerify, double localTheta, double stopThreshold) {
            std::unordered_map<row_type, std::pair<double, std::vector<col_type> > >::iterator it;

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
            worstScores.erase(worstScores.begin(), worstScores.end());
        }

        inline virtual void updateState(col_type col) {// given that the step happened on col update the structures

            if (query[col] < 0) {//scan downwards //////////////////

                fringePos[col]++;
                if (fringePos[col] >= (col + 1) * invLists->getRowNum()) {
                    inactiveCols++;
                }

            } else {//scan upwards

                // keep moving upwards
                fringePos[col]--;
                if (fringePos[col] < col * invLists->getRowNum()) {
                    inactiveCols++;
                }
            }

        }



        ///////////////// methods for threshold update ////////////////////////

        inline double getInitialThreshold() {
            double stopThreshold = 0;
            for (col_type i = 0; i < colNum; i++) {
                stopThreshold += invLists->getElement(fringePos[i])->data * query[i];
            }
            return stopThreshold;
        }


        inline bool isThresholdUnterTheta(double& stopThreshold, double localTheta, col_type stepOnCol, double oldValue, bool forCosine) {

            if (colNum <= inactiveCols) {
                return true;
            }
            double piNew = invLists->getElement(fringePos[stepOnCol])->data;

            if (piNew == oldValue)
                return false;

            //updateThreshold
            stopThreshold += query[stepOnCol] * (piNew - oldValue);


            if (stopThreshold < localTheta) {
                return true;
            }

            return false;
        }

        inline virtual col_type chooseStep() {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }
   


    };

    class TANRAStateMax : public TANRAState {
        maxHeap piqi;



    public:

        TANRAStateMax(col_type colNum) : TANRAState(colNum) {
        }

        virtual ~TANRAStateMax() {
        }

        inline void initForNewQuery(const double* q) {// initialize to specific rows

            TANRAState::initForNewQuery(q);

            double value;
            
            

            for (col_type i = 0; i < colNum; i++) {
                if (query[i] < 0) { //scan downwards///////////////////////
                    fringePos[i] = i * invLists->getRowNum() + 0;

                } else {
                    fringePos[i] = i * invLists->getRowNum() + rowNum - 1;
                }
                value = invLists->getElement(fringePos[i])->data * query[i];
                piqi.add(QueueElement(value, i));
            }
        }

        inline virtual void updateState(col_type col) {// given that the step happened on col update the structures

            row_type rowId;
            bool inactivate = false;

            if (query[col] < 0) {//scan downwards //////////////////

                fringePos[col]++;
                if (fringePos[col] == (col + 1) * invLists->getRowNum()) {
                    inactivate = true;
                } else {
                    rowId = invLists->getElement(fringePos[col])->id;
                }

            } else {//scan upwards

                // keep moving upwards
                fringePos[col]--;
                if (fringePos[col] < col * invLists->getRowNum()) {
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


        ///////////////// methods for choosing next step ////////////////////////

        inline col_type chooseStep() {
            col_type stepOnCol = piqi.getMaxId();
            return stepOnCol;
        }
        ////////////////////////////////////////


    };

    class TANRAStateRR : public TANRAState {
        std::vector<col_type> activeCols;
        col_type colIndx;

    public:

        TANRAStateRR(col_type colNum) : TANRAState(colNum) {
        }

        virtual ~TANRAStateRR() {
        }

        inline virtual void initForNewQuery(const double* q) {// initialize to specific rows

            TANRAState::initForNewQuery(q);
            activeCols.clear();
            colIndx = 0;

            for (col_type i = 0; i < colNum; i++) {
                if (query[i] < 0) { //scan downwards///////////////////////
                    fringePos[i] = i * invLists->getRowNum() + 0;
                } else {
                    fringePos[i] = i * invLists->getRowNum() + rowNum - 1;
                }
                if (query[i] != 0) {
                    activeCols.push_back(i);
                }
            }
        }

        using TANRAState::updateState;

        ///////////////// methods for choosing next step ////////////////////////

        inline virtual col_type chooseStep() {
                      
            if (colIndx == activeCols.size())
                colIndx = 0;

            col_type stepOnCol = activeCols[colIndx];
            colIndx++;
            return stepOnCol;
        }
        ////////////////////////////////////////


    };
}

#endif	/* TANRASTATE_H */

