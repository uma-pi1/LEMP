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

        taRetriever() = default;
        ~taRetriever() = default;

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            col_type stepOnCol = 0;
            row_type posMatrix;
            double stopThreshold;
            double oldValue;
            double localTheta;

            if (arg->forCosine) {
                localTheta = probeBucket.bucketScanThreshold / query[-1];
            } else {
                localTheta = arg->theta;
            }
#ifdef TIME_IT
            arg->t.start();
#endif
            //initialize scheduler and state
            arg->state->initForNewQuery(query);
            stopThreshold = arg->state->getInitialThreshold();
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            row_type countSeen = 0;
            //check if you need to stop. Otherwise continue;
            while (stopThreshold >= localTheta) {

                //choose next step
                stepOnCol = arg->state->chooseStep();
                //pick up the item
                
                oldValue = invLists->getElement(arg->state->getDepthInFringeInCol(stepOnCol))->data;
                posMatrix = invLists->getElement(arg->state->getDepthInFringeInCol(stepOnCol))->id;
#ifdef TIME_IT
                arg->t.start();
#endif
                //if I haven't explored the item already, do RA
                if (!arg->state->isExplored(posMatrix)) {
                    verifyCandidate(posMatrix + probeBucket.startPos, query, arg);
                    countSeen++;
                    arg->state->maintainExploredLists(*(arg->probeMatrix), countSeen, stepOnCol);
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->ipTime += arg->t.elapsedTime().nanos();

                arg->t.start();
#endif
                // update the state
                arg->state->updateState(stepOnCol);
#ifdef TIME_IT
                arg->t.stop();
                arg->scanTime += arg->t.elapsedTime().nanos();
#endif
                if (arg->state->areAllSeen()) {
                    break;
                }

#ifdef TIME_IT
                arg->t.start();
#endif
                arg->state->isThresholdUnterTheta(stopThreshold, localTheta, stepOnCol, oldValue, arg->forCosine);
#ifdef TIME_IT
                arg->t.stop();
                arg->boundsTime += arg->t.elapsedTime().nanos();
#endif
            }


        }

        inline virtual void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

            col_type stepOnCol = 0;
            row_type posMatrix;
            double oldValue;

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));


            double x1 = 1;
            double x2 = 1;

            if (arg->forCosine) {
                x1 = probeBucket.invNormL2.second;
                x2 = probeBucket.invNormL2.first;
            }
            double localTheta = arg->heap.front().data * (arg->heap.front().data > 0 ? x1 : x2);
#ifdef TIME_IT
            arg->t.start();
#endif
            arg->state->initForNewQuery(query);
            double stopThreshold = arg->state->getInitialThreshold();
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            row_type countSeen = 0;
            while (stopThreshold > localTheta) {

                //choose next step
                stepOnCol = arg->state->chooseStep();
                //pick up the item
                oldValue = invLists->getElement(arg->state->getDepthInFringeInCol(stepOnCol))->data;
                posMatrix = invLists->getElement(arg->state->getDepthInFringeInCol(stepOnCol))->id;
#ifdef TIME_IT
                arg->t.start();
#endif
                //if I haven't explored the item already, do RA
                if (!arg->state->isExplored(posMatrix)) {
                    verifyCandidateTopk(posMatrix + probeBucket.startPos, query, arg);
                    countSeen++;
                    arg->state->maintainExploredLists(*(arg->probeMatrix), countSeen, stepOnCol);
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->ipTime += arg->t.elapsedTime().nanos();

                arg->t.start();
#endif
                // update the state
                arg->state->updateState(stepOnCol);
#ifdef TIME_IT
                arg->t.stop();
                arg->scanTime += arg->t.elapsedTime().nanos();
#endif
                if (arg->state->areAllSeen())
                    break;
#ifdef TIME_IT
                arg->t.start();
#endif
                localTheta = arg->heap.front().data * (arg->heap.front().data > 0 ? x1 : x2);
                arg->state->isThresholdUnterTheta(stopThreshold, localTheta, stepOnCol, oldValue, arg->forCosine);
#ifdef TIME_IT
                arg->t.stop();
                arg->boundsTime += arg->t.elapsedTime().nanos();
#endif
            }
        }

        inline virtual void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type sampleSize = probeBucket.xValues->size();
            
            if (sampleSize > 0) {
                sampleTimes.resize(sampleSize);
                QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

                retrArg[0].state->initializeForNewBucket(invLists);

                for (row_type i = 0; i < sampleSize; ++i) {

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

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type sampleSize = (probeBucket.xValues!=nullptr ? probeBucket.xValues->size() : 0);

            if (sampleSize > 0) {
                sampleTimes.resize(sampleSize);
                QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

                // todo shouldn't I initialize the lists? no! already initialized from Incremental

                retrArg[0].state->initializeForNewBucket(invLists);

                for (row_type i = 0; i < sampleSize; ++i) {

                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    const std::vector<QueueElement>& prevResults = prevBucket.sampleThetas->at(t).at(ind).results;

                    retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                    std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                    retrArg[0].tunerTimer.start();
                    runTopK(query, probeBucket, &retrArg[0]);
                    retrArg[0].tunerTimer.stop();
                    sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                }
            }else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.isWorkDone())
                    continue;

                if (!invLists->isInitialized()) {
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
                    if (queryBatch.isQueryInactive(user)) {
                        user++;
                        continue;
                    }


                    const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                    double minScore = arg->topkResults[i].data;

                    if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                        queryBatch.inactivateQuery(user);
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

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            arg->state->initializeForNewBucket(invLists);

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }


                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;
                    arg->queryId = arg->queryMatrix->getId(i);
                    run(query, probeBucket, arg);
                }

            }

        }
        
               inline virtual void cleanupAfterTuning() {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

    };


}
#endif /* TA_H_ */
