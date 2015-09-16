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
        std::vector<col_type> activeCols;
        col_type colIndx;


        TANRAState(col_type colNum) : colNum(colNum), invLists(0){
            fringePos.resize(colNum);
        }

        inline void initializeForNewBucket(QueueElementLists* m1) {// first initialization
            invLists = m1;
            rowNum = m1->getRowNum();	    
        }

        inline void initForNewQuery(const double* q) {// initialize to specific rows

            inactiveCols = 0;
            query = q;
	    activeCols.clear();


            for (col_type i = 0; i < colNum; i++) {
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

     

        inline void updateScores(row_type key, col_type stepOnCol) {   	    
	    worstScores[key].first += invLists->getElement(fringePos[stepOnCol])->data * query[stepOnCol];
            worstScores[key].second.push_back(stepOnCol);
        }

  

        inline void getCandidates(std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify, double localTheta, double stopThreshold) {
            unordered_map<row_type, std::pair<double, std::vector<col_type> >  >::iterator it;

            for (it = worstScores.begin(); it != worstScores.end(); it++) {
                row_type key = it->first;
                double ws = it->second.first;

                if (ws >= localTheta) {
                    candidatesToVerify[numCandidatesToVerify] = key;
                    numCandidatesToVerify++;
                    continue;
                }
                ws += stopThreshold;
                
                for (int i=0; i< it->second.second.size() ; i++) {
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
        
        


        inline void updateState(col_type col) {// given that the step happened on col update the structures
	    double newPiqi;

            if (query[col] < 0) {//scan downwards //////////////////

                fringePos[col]++;
                if (fringePos[col] >= (col + 1) * invLists->size) {
		    inactiveCols++;
                } 

            } else {//scan upwards

                // keep moving upwards
                fringePos[col]--;
                if (fringePos[col] < col * invLists->size) {		  
		    inactiveCols++;	
                } 
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


        inline bool isThresholdUnterTheta(double& stopThreshold, double localTheta, col_type stepOnCol, double oldValue, bool forCosine) {
            if (colNum <= inactiveCols) {
                return true;
            }
	    double piNew = invLists->getElement(fringePos[stepOnCol])->data;
	    
	    if(piNew == oldValue)
	      return false;

	    //updateThreshold
	    stopThreshold += query[stepOnCol] * (piNew - oldValue);
         

            if (stopThreshold < localTheta) {
                return true;
            }

//             if (forCosine) {
//                 if (stopThreshold > 1 && 1 < localTheta) {
//                     return true;
//                 }
//             }

            return false;
        }

        ///////////////// methods for choosing next step ////////////////////////


        inline col_type chooseStep() {
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

