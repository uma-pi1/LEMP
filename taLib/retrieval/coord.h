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
 * coord.h
 *
 *  Created on: Jun 25, 2014
 *      Author: chteflio
 */

#ifndef COORD_H_
#define COORD_H_

namespace ta {

    class CoordRetriever : public Retriever {
    public:
        std::unique_ptr<ListTuneData> dataForTuning;

        CoordRetriever() = default;
        ~CoordRetriever() = default;

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

            IntLists* invLists = static_cast<IntLists*> (probeBucket.getIndex(INT_SL));
            row_type numItemsToVerify = 0;
            double localTheta = probeBucket.bucketScanThreshold / query[-1];

#ifdef TIME_IT
            arg->t.start();
#endif
            bool shouldScan = invLists->calculateIntervals(query, arg->listsQueue, arg->intervals, localTheta, arg->numLists);
#ifdef TIME_IT
            arg->t.stop();
            arg->boundsTime += arg->t.elapsedTime().nanos();
#endif



            if (shouldScan) {
                //initialize         
#ifdef TIME_IT
                arg->t.start();
#endif
                row_type* entry = invLists->getElement(arg->intervals[0].start);
                row_type length = arg->intervals[0].end - arg->intervals[0].start;

                for (row_type i = 0; i < length; ++i) {
                    arg->cp_array[entry[i]] = 1;
                }

                // add Candidates
                for (col_type j = 1; j < arg->numLists; ++j) {

                    entry = invLists->getElement(arg->intervals[j].start);
                    length = arg->intervals[j].end - arg->intervals[j].start;

                    for (row_type i = 0; i < length; ++i) {
                        arg->cp_array[entry[i]]++;
                    }
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->scanTime += arg->t.elapsedTime().nanos();
                arg->t.start();
#endif
                // final scan
                entry = invLists->getElement(arg->intervals[0].start);
                length = arg->intervals[0].end - arg->intervals[0].start;

                for (row_type i = 0; i < length; ++i) {
                    if (arg->cp_array[entry[i]] == arg->numLists) {
                        arg->candidatesToVerify[numItemsToVerify] = entry[i] + probeBucket.startPos;
                        numItemsToVerify++;
                    }
                }

#ifdef TIME_IT
                arg->t.stop();
                arg->filterTime += arg->t.elapsedTime().nanos();
                arg->t.start();
#endif
                verifyCandidates_lengthTest(query, numItemsToVerify, arg);
#ifdef TIME_IT
                arg->t.stop();
                arg->ipTime += arg->t.elapsedTime().nanos();
#endif
            }
        }

        inline virtual void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

#ifdef TIME_IT
            arg->t.start();
#endif
            if (!queryBatch.hasInitializedQueues()) { //preprocess
                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            arg->numLists = probeBucket.numLists;
            for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                    break;

                col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);
                arg->setQueues(localQueue);
                arg->queryId = arg->queryMatrix->getId(i);

                run(query, probeBucket, arg);
            }

        }

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const{

            row_type numItemsToVerify = 0;
            IntLists* invLists = static_cast<IntLists*> (probeBucket.getIndex(INT_SL));
            double localTheta = arg->heap.front().data * (arg->heap.front().data > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);

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
#ifdef TIME_IT
            arg->t.start();
#endif

            bool shouldScan = invLists->calculateIntervals(query, arg->listsQueue, arg->intervals, localTheta, arg->numLists);

#ifdef TIME_IT
            arg->t.stop();
            arg->boundsTime += arg->t.elapsedTime().nanos();
#endif
            if (shouldScan) {
#ifdef TIME_IT
                arg->t.start();
#endif
                //initialize
                row_type* entry = invLists->getElement(arg->intervals[0].start);
                row_type length = arg->intervals[0].end - arg->intervals[0].start;

                for (row_type i = 0; i < length; ++i) {
                    arg->cp_array[entry[i]] = 1;
                }

                // add Candidates
                for (col_type j = 1; j < arg->numLists; ++j) {

                    entry = invLists->getElement(arg->intervals[j].start);
                    length = arg->intervals[j].end - arg->intervals[j].start;

                    for (row_type i = 0; i < length; ++i) {
                        arg->cp_array[entry[i]]++;
                    }
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->scanTime += arg->t.elapsedTime().nanos();
                arg->t.start();
#endif
                // final scan
                entry = invLists->getElement(arg->intervals[0].start);
                length = arg->intervals[0].end - arg->intervals[0].start;

                for (row_type i = 0; i < length; ++i) {
                    if (arg->cp_array[entry[i]] == arg->numLists) {
                        arg->candidatesToVerify[numItemsToVerify] = entry[i] + probeBucket.startPos;
                        numItemsToVerify++;
                    }
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->filterTime += arg->t.elapsedTime().nanos();
                arg->t.start();
#endif
                verifyCandidatesTopK_lengthTest(query, numItemsToVerify, arg);
#ifdef TIME_IT
                arg->t.stop();
                arg->ipTime += arg->t.elapsedTime().nanos();
#endif
            }

        }

        inline virtual void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const{
#ifdef TIME_IT
            arg->t.start();
#endif
            if (!queryBatch.hasInitializedQueues()) { //preprocess
                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            arg->numLists = probeBucket.numLists;

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
                col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                arg->setQueues(localQueue);

                runTopK(query, probeBucket, arg);

                arg->writeHeapToTopk(user);
                user++;
            }


        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            if (probeBucket.xValues->size() > 0) {
                dataForTuning = std::unique_ptr<ListTuneData>(new ListTuneData());
                dataForTuning->tune(probeBucket, prevBucket, retrArg, this);
            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            
                dataForTuning = std::unique_ptr<ListTuneData>(new ListTuneData());
                dataForTuning->tuneTopk(probeBucket, prevBucket, retrArg, this);
            

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            arg->numLists = probeBucket.numLists;

            for (auto& queryBatch : arg->queryBatches) {
                if (queryBatch.isWorkDone())
                    continue;
#ifdef TIME_IT
                arg->t.start();
#endif
                if (!queryBatch.hasInitializedQueues()) { //preprocess
                    queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif                

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
                    col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                    arg->setQueues(localQueue);

                    runTopK(query, probeBucket, arg);

                    arg->writeHeapToTopk(user);
                    user++;
                }




            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const{
            arg->numLists = probeBucket.numLists;

            for (auto& queryBatch : arg->queryBatches) {
                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }
#ifdef TIME_IT
                arg->t.start();
#endif
                if (!queryBatch.hasInitializedQueues()) { //preprocess
                    queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif                
                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;

                    col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);
                    arg->setQueues(localQueue);
                    arg->queryId = arg->queryMatrix->getId(i);

                    run(query, probeBucket, arg);
                }
            }

        }
       inline virtual void cleanupAfterTuning() {
           dataForTuning.reset(nullptr);
        }


    };





}
#endif /* COORD_H_ */
