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
 * File:   tanra.h
 * Author: chteflio
 *
 * Created on August 24, 2015, 9:15 AM
 */

#ifndef TANRA_H
#define	TANRA_H

namespace ta {

    class tanraRetriever : public Retriever {
    public:

        tanraRetriever() = default;

        ~tanraRetriever() = default;

        //lists is dummy

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            col_type stepOnCol = 0;
            row_type posMatrix;
            double stopThreshold;
            double oldValue = 0;
            double localTheta;

            unordered_map<row_type, double>::iterator it;

            if (arg->forCosine) {
                localTheta = probeBucket.bucketScanThreshold / query[-1];
            } else {
                localTheta = arg->theta;
            }

#ifdef TIME_IT
            arg->t.start();
#endif   
            //initialize scheduler and state          
            arg->tanraState->initForNewQuery(query);
            stopThreshold = arg->tanraState->getInitialThreshold();
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif
	  

	    
            //check if you need to stop. Otherwise continue;
            while (stopThreshold >= localTheta ) {
#ifdef TIME_IT
                arg->t.start();
#endif
                //choose next step
                stepOnCol = arg->tanraState->chooseStep();
                //pick up the item
                oldValue = invLists->getElement(arg->tanraState->getDepthInFringeInCol(stepOnCol))->data;
                posMatrix = invLists->getElement(arg->tanraState->getDepthInFringeInCol(stepOnCol))->id;
		row_type key = posMatrix + probeBucket.startPos;    
	
              
		

                arg->tanraState->updateScores(key, stepOnCol);
                arg->tanraState->updateState(stepOnCol);
	#ifdef TIME_IT	
                arg->t.stop();
                arg->scanTime += arg->t.elapsedTime().nanos();


                arg->t.start();
		#endif
		arg->tanraState->isThresholdUnterTheta(stopThreshold, localTheta, stepOnCol, oldValue, arg->forCosine);
#ifdef TIME_IT                 
                arg->t.stop();
                arg->boundsTime += arg->t.elapsedTime().nanos();
#endif

            }
            #ifdef TIME_IT
            arg->t.start();
	    #endif

            row_type numCandidatesToVerify = 0;
            arg->tanraState->getCandidates(arg->candidatesToVerify, numCandidatesToVerify, localTheta, stopThreshold);

 #ifdef TIME_IT
            arg->t.stop();
            arg->filterTime += arg->t.elapsedTime().nanos();
            arg->t.start();
    #endif
            verifyCandidates_noLengthTest(query, numCandidatesToVerify, arg);
          #ifdef TIME_IT

            arg->t.stop();
            arg->ipTime += arg->t.elapsedTime().nanos();
	    #endif

        }

        inline virtual void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            
             row_type sampleSize = probeBucket.xValues->size();

            if (sampleSize > 0) {
                sampleTimes.resize(sampleSize);
                QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

                retrArg[0].tanraState->initializeForNewBucket(invLists);

                for (row_type i = 0; i < sampleSize; i++) {

                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    retrArg[0].tunerTimer.start();
                    run(query, probeBucket, &retrArg[0]);
                    retrArg[0].tunerTimer.stop();
                    sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                }

            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }


        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket,  std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            arg->tanraState->initializeForNewBucket(invLists);

            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }

                QueryBatch& queryBatch = arg->queryBatches[q];

                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;
                    arg->queryId = arg->queryMatrix->getId(i);
                    run(query, probeBucket, arg);
                }

            }

        }

    };


}

#endif	/* TANRA_H */

