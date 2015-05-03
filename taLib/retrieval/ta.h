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
 * ta.h
 *
 *  Created on: Jun 25, 2014
 *      Author: chteflio
 */

#ifndef TA_H_
#define TA_H_


namespace ta {

    class taRetriever : public Retriever {
    public:

        //double bestTimeLTA;
        row_type t_b_indx;
        std::vector<double> high;
        std::vector<double> low;

        taRetriever() : Retriever(LEMP_TA) {
        };

        ~taRetriever() {
        };

        //lists is dummy

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            col_type stepOnCol = 0;
            row_type posMatrix, row;
            bool stopFlag = false;
            double stopThreshold;
            row_type prevRow;
            double oldValue, newValue;
            bool showFlag = false;
            double localTheta;

            if (arg->forCosine) {
                localTheta = probeBucket.bucketScanThreshold / query[-1];
            } else {
                localTheta = arg->theta;
            }


            //initialize scheduler and state
            arg->state->initForNewQuery(query);
            //stepOnCol = state.setQuery();
            stopThreshold = arg->state->initializeThreshold();



            row_type countSeen = 0;
            //check if you need to stop. Otherwise continue;
            while (!arg->state->allSeen) {

                //choose next step

                stepOnCol = arg->state->chooseStep();
                //pick up the item
                //                row = arg->state->fringeDepth[stepOnCol];
                //                oldValue = invLists->getValue(row, stepOnCol);
                //                posMatrix = invLists->getRowPointer(row, stepOnCol);
                oldValue = invLists->getElement(arg->state->fringePos[stepOnCol])->data;
                posMatrix = invLists->getElement(arg->state->fringePos[stepOnCol])->id;

                //if I haven't explored the item already, do RA
                if (!arg->state->exploredItems[posMatrix]) {

                    verifyCandidate(posMatrix + probeBucket.startPos, query, arg);
                    // maintain the explored lists
                    countSeen++;
                    arg->state->maintainExploredLists(*(arg->probeMatrix), countSeen, stepOnCol);
                }

                // update the state
                arg->state->updateState(stepOnCol);


                //check if you need to stop. Otherwise continue;
                if (!arg->state->allSeen) {
                    stopFlag = arg->state->isThresholdUnterTheta(stopThreshold, localTheta, stepOnCol, oldValue, arg->forCosine);
                }// otherwise you stop in  any case.


                if (stopFlag || arg->state->allSeen) {
                    break;
                }

            }
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
#ifdef TIME_IT
            //std::cout<<" TA-based retrieval running "<<std::endl;
#endif
            // this can go out of th loop. I run incr for all

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            arg->state->initializeForNewBucket(invLists);


            for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
                const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                    break;
                arg->queryId = arg->queryMatrix->getId(i);
                run(query, probeBucket, arg);
            }

        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            col_type stepOnCol = 0;
            row_type posMatrix, row;
            bool stopFlag = false;
            double stopThreshold;
            row_type prevRow;
            double oldValue, newValue;
            bool showFlag = false;
            bool initializeThresholder = false;

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));


            double x1 = 1;
            double x2 = 1;

            if (arg->forCosine) {
                x1 = probeBucket.invNormL2.second;
                x2 = probeBucket.invNormL2.first;
            }



            //initialize scheduler and state
            arg->state->initForNewQuery(query);
            //stepOnCol = state.setQuery();
            //std::cout<<"done set query"<<std::endl;


            row_type countSeen = 0;
            row_type i = 0;
            //check if you need to stop. Otherwise continue;
            while (!arg->state->allSeen) {


                //choose next step

                stepOnCol = arg->state->chooseStep();
                //pick up the item
                //                row = arg->state->fringeDepth[stepOnCol];
                //                posMatrix = invLists->getRowPointer(row, stepOnCol);
                //                oldValue = invLists->getValue(row, stepOnCol);

                oldValue = invLists->getElement(arg->state->fringePos[stepOnCol])->data;
                posMatrix = invLists->getElement(arg->state->fringePos[stepOnCol])->id;


                //if I haven't explored the item already, do RA
                if (!arg->state->exploredItems[posMatrix]) {
                    verifyCandidateTopk(posMatrix + probeBucket.startPos, query, arg);
                    // maintain the explored lists
                    countSeen++;
                    arg->state->maintainExploredLists(*(arg->probeMatrix), countSeen, stepOnCol);
                }


                if (!initializeThresholder) {
                    stopThreshold = arg->state->initializeThreshold();
                    initializeThresholder = true;
                }

                // update the state
                arg->state->updateState(stepOnCol);


                /* Optimization for skipping threshold computation for top-k */
                if (arg->heap.size() < arg->k)
                    continue;


                //check if you need to stop. Otherwise continue;
                if (!arg->state->allSeen) {
                    double localTheta = arg->heap.front().data * (arg->heap.front().data > 0 ? x1 : x2);
                    stopFlag = arg->state->isThresholdUnterTheta(stopThreshold, localTheta, stepOnCol, oldValue, arg->forCosine);
                }


                if (stopFlag || arg->state->allSeen) {
                    break;
                }
                i++;
            }
        }

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            if (!invLists->initialized) {
#ifdef TIME_IT
                arg->t.start();
#endif
                invLists->initializeLists(*(arg->probeMatrix), probeBucket.startPos, probeBucket.endPos);
#ifdef TIME_IT
                arg->t.stop();
                arg->initializeListsTime += arg->t.elapsedTime().nanos();
#endif
            }

            arg->state->initializeForNewBucket(invLists);

            row_type user = queryBatch.startPos;
            int start = queryBatch.startPos * arg->k;
            int end = queryBatch.endPos * arg->k;
            for (row_type i = start; i < end; i += arg->k) {
                if (queryBatch.inactiveQueries[user - queryBatch.startPos]) {
                    user++;
                    continue;
                }


                const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                double minScore = arg->topkResults[i].data;

                if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                    queryBatch.inactiveQueries[user - queryBatch.startPos] = true;
                    queryBatch.inactiveCounter++;
                    user++;
                    continue;
                }

                arg->moveTopkToHeap(i);

                arg->queryId = arg->queryMatrix->getId(user);
                runTopK(query, probeBucket, arg);

                arg->writeHeapToTopk(user);

                user++;
            }



        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            if (xValues->size() > 0) {
                sampleTimes.resize(xValues->size());
                QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

                retrArg[0].state->initializeForNewBucket(invLists);


                for (row_type i = 0; i < xValues->size(); i++) {

                    int t = xValues->at(i).i;
                    int ind = xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    retrArg[0].tunerTimer.start();
                    run(query, probeBucket, &retrArg[0]);
                    retrArg[0].tunerTimer.stop();
                    sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                }

            }


        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type b = retrArg[0].bucketInd;

            if (xValues->size() > 0) {
                sampleTimes.resize(xValues->size());
                QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

                // todo shouldn't I initialize the lists? no! already initialized from Incremental

                retrArg[0].state->initializeForNewBucket(invLists);

                for (row_type i = 0; i < xValues->size(); i++) {

                    int t = xValues->at(i).i;
                    int ind = xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    std::vector<QueueElement>& prevResults = retrArg[0].globalData[b - 1][t][ind].results;

                    retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                    std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                    retrArg[0].tunerTimer.start();
                    runTopK(query, probeBucket, &retrArg[0]);
                    retrArg[0].tunerTimer.stop();
                    sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                }


            }

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));


            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                    continue;

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];




                if (!invLists->initialized) {
#ifdef TIME_IT
                    arg->t.start();
#endif
                    invLists->initializeLists(*(arg->probeMatrix), probeBucket.startPos, probeBucket.endPos);
#ifdef TIME_IT
                    arg->t.stop();
                    arg->initializeListsTime += arg->t.elapsedTime().nanos();
#endif
                }

                arg->state->initializeForNewBucket(invLists);

                row_type user = queryBatch.startPos;
                int start = queryBatch.startPos * arg->k;
                int end = queryBatch.endPos * arg->k;
                for (row_type i = start; i < end; i += arg->k) {
                    if (queryBatch.inactiveQueries[user - queryBatch.startPos]) {
                        user++;
                        continue;
                    }


                    const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                    double minScore = arg->topkResults[i].data;

                    if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                        queryBatch.inactiveQueries[user - queryBatch.startPos] = true;
                        queryBatch.inactiveCounter++;
                        user++;
                        continue;
                    }

                    arg->moveTopkToHeap(i);

                    arg->queryId = arg->queryMatrix->getId(user);
                    runTopK(query, probeBucket, arg);

                    arg->writeHeapToTopk(user);

                    user++;
                }

            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            arg->state->initializeForNewBucket(invLists);

            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];

#ifdef TIME_IT
                //std::cout<<" TA-based retrieval running "<<std::endl;
#endif
     
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
#endif /* TA_H_ */
