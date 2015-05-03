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
 * File:   incr2.h
 * Author: chteflio
 *
 * Created on December 14, 2014, 11:59 AM
 */

#ifndef INCR2_H
#define	INCR2_H

#include "ListsTuneData.h"


namespace ta {

    class IncrRetriever : public Retriever {
    public:

        ListTuneData dataForTuning;

        IncrRetriever() : Retriever(LEMP_I) {
        };

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));

            double qi;
            row_type numCandidatesToVerify = 0;


            double seenQi2 = 1;

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

                if (arg->numLists == 1) {// ||arg->intervals[0].end-arg->intervals[0].start < 10 || arg->intervals[0].end-arg->intervals[0].start < 0.05 * (probeBucket.endPos-probeBucket.startPos)){ // scan only one list

                    scanOnly1List(query, probeBucket, arg, invLists, false);

                } else {

#ifdef TIME_IT
                    arg->t.start();
#endif
                    //initialize                   
                    qi = query[arg->intervals[0].col];
                    seenQi2 -= qi * qi;

                    QueueElement* entry = invLists->getElement(arg->intervals[0].start);
                    row_type length = arg->intervals[0].end - arg->intervals[0].start;

                    for (row_type i = 0; i < length; i++) {
                        arg->ext_cp_array[entry[i].id].addFirst(qi, entry[i].data);
                    }

                    // add Candidates
                    for (col_type j = 1; j < arg->numLists; j++) {
                        qi = query[arg->intervals[j].col];

                        if (qi == 0)
                            continue;

                        seenQi2 -= qi * qi;

                        entry = invLists->getElement(arg->intervals[j].start);
                        length = arg->intervals[j].end - arg->intervals[j].start;

                        for (row_type i = 0; i < length; i++) {
                            arg->ext_cp_array[entry[i].id].add(qi, entry[i].data);
                        }

                    }


#ifdef TIME_IT
                    arg->t.stop();
                    arg->scanTime += arg->t.elapsedTime().nanos();
                    arg->t.start();
#endif
                    // run first scan again to find items to verify

                    entry = invLists->getElement(arg->intervals[0].start);
                    length = arg->intervals[0].end - arg->intervals[0].start;

                    for (row_type i = 0; i < length; i++) {
                        row_type row = entry[i].id;
                        double len = query[-1] * arg->probeMatrix->lengthInfo[row + probeBucket.startPos].data; //.getVectorLength(row+probeBucket.startPos);
                        
                        if (!arg->ext_cp_array[row].prune(len, arg->theta, seenQi2) ) {
                            arg->candidatesToVerify[numCandidatesToVerify] = row + probeBucket.startPos;
                            numCandidatesToVerify++;
                        }
                        
                        
                        
                    }

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
            }
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
#ifdef TIME_IT
            //std::cout<<" INCR-based retrieval running "<<std::endl;
            arg->t.start();
#endif
            // this can go out of the loop. I run incr for all
            if (!queryBatch.initializedQueues) { //preprocess              
                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);

                //                arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);                
                queryBatch.initializedQueues = true;
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            arg->numLists = probeBucket.numLists;


            for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
                const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                    break;

                col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);

                //                col_type* localQueue = arg->queryMatrix->getQueue(i);           
                arg->setQueues(localQueue);
                arg->queryId = arg->queryMatrix->getId(i);

                run(query, probeBucket, arg);

            }
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {



            double localTheta = arg->heap.front().data * (arg->heap.front().data > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);

            double qi;
            row_type numCandidatesToVerify = 0;

            double seenQi2 = 1;
            col_type validLists = arg->numLists;
            double prevSeenQi2 = 1;

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




#ifdef TIME_IT
            arg->t.start();
#endif
            bool shouldScan = invLists->calculateIntervals(query, arg->listsQueue, arg->intervals, localTheta, validLists);
#ifdef TIME_IT
            arg->t.stop();
            arg->boundsTime += arg->t.elapsedTime().nanos();
#endif

            if (shouldScan) {


                if (validLists == 1) {// || arg->intervals[0].end-arg->intervals[0].start < 10 || arg->intervals[0].end-arg->intervals[0].start < 0.05 * (probeBucket.endPos-probeBucket.startPos)){ // scan only one list

                    scanOnly1ListTopk(query, probeBucket, arg, invLists, true);

                } else {
#ifdef TIME_IT
                    arg->t.start();
#endif
                    //initialize

                    qi = query[arg->intervals[0].col];
                    seenQi2 -= qi * qi;

                    QueueElement* entry = invLists->getElement(arg->intervals[0].start);
                    row_type length = arg->intervals[0].end - arg->intervals[0].start;

                    for (row_type i = 0; i < length; i++) {
                        arg->ext_cp_array[entry[i].id].addFirst(qi, entry[i].data);
                    }



                    // add Candidates
                    for (int j = 1; j < validLists; j++) {
                        qi = query[arg->intervals[j].col];

                        if (qi == 0)
                            continue;

                        seenQi2 -= qi * qi;

                        entry = invLists->getElement(arg->intervals[j].start);
                        length = arg->intervals[j].end - arg->intervals[j].start;

                        for (row_type i = 0; i < length; i++) {
                            arg->ext_cp_array[entry[i].id].add(qi, entry[i].data);
                        }

                    }
#ifdef TIME_IT
                    arg->t.stop();
                    arg->scanTime += arg->t.elapsedTime().nanos();
                    arg->t.start();
#endif
                    // run first scan again to find items to verify

                    entry = invLists->getElement(arg->intervals[0].start);
                    length = arg->intervals[0].end - arg->intervals[0].start;

                    for (row_type i = 0; i < length; i++) {
                        row_type row = entry[i].id;
                        double len = arg->probeMatrix->lengthInfo[row + probeBucket.startPos].data; //.getVectorLength(row+probeBucket.startPos);

                        if (!arg->ext_cp_array[row].prune(len, arg->heap.front().data, seenQi2)) {
                            arg->candidatesToVerify[numCandidatesToVerify] = row + probeBucket.startPos;
                            numCandidatesToVerify++;
                        }
                        
                    }

#ifdef TIME_IT
                    arg->t.stop();
                    arg->filterTime += arg->t.elapsedTime().nanos();
                    arg->t.start();
#endif
                    verifyCandidatesTopK_noLengthTest(query, numCandidatesToVerify, arg);
#ifdef TIME_IT
                    arg->t.stop();
                    arg->ipTime += arg->t.elapsedTime().nanos();
#endif
                }

            }
        }

        // this runs practically COORD on 1 list.
        // Use it when the shortest necessary interval is REALLY short

        inline void scanOnly1List(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg,
                QueueElementLists* invLists, bool topk) {
#ifdef TIME_IT
            arg->t.start();
#endif           

            row_type numCandidatesToVerify = 0;

            QueueElement* entry = invLists->getElement(arg->intervals[0].start);
            row_type length = arg->intervals[0].end - arg->intervals[0].start;

            for (row_type i = 0; i < length; i++) {
                arg->candidatesToVerify[numCandidatesToVerify] = entry[i].id + probeBucket.startPos;
                numCandidatesToVerify++;
            }


#ifdef TIME_IT
            arg->t.stop();
            arg->scanTime += arg->t.elapsedTime().nanos();
            arg->t.start();
#endif
            verifyCandidates_lengthTest(query, numCandidatesToVerify, arg);
#ifdef TIME_IT
            arg->t.stop();
            arg->ipTime += arg->t.elapsedTime().nanos();
#endif

        }

        inline void scanOnly1ListTopk(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg,
                QueueElementLists* invLists, bool topk) {
#ifdef TIME_IT
            arg->t.start();
#endif          
            row_type numCandidatesToVerify = 0;

            QueueElement* entry = invLists->getElement(arg->intervals[0].start);
            row_type length = arg->intervals[0].end - arg->intervals[0].start;

            for (row_type i = 0; i < length; i++) {
                arg->candidatesToVerify[numCandidatesToVerify] = entry[i].id + probeBucket.startPos;
                numCandidatesToVerify++;
            }

#ifdef TIME_IT
            arg->t.stop();
            arg->scanTime += arg->t.elapsedTime().nanos();
            arg->t.start();
#endif
            verifyCandidatesTopK_lengthTest(query, numCandidatesToVerify, arg);
#ifdef TIME_IT
            arg->t.stop();
            arg->ipTime += arg->t.elapsedTime().nanos();
#endif

        }

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {

#ifdef TIME_IT
            arg->t.start();
#endif
            // this can go out of the loop. I run incr for all
            if (!queryBatch.initializedQueues) { //preprocess
                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                //                arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);       
                queryBatch.initializedQueues = true;
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

            arg->numLists = probeBucket.numLists;

            ///////////////////////////////////
            row_type user = queryBatch.startPos;
            int start = queryBatch.startPos * arg->k;
            int end = queryBatch.endPos * arg->k;
            for (row_type i = start; i < end; i += arg->k) {

                //                std:cout<<"i: "<<i<<" user: "<<user
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
                col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                //                col_type* localQueue = arg->queryMatrix->getQueue(user);
                arg->setQueues(localQueue);

                runTopK(query, probeBucket, arg);

                arg->writeHeapToTopk(user);

                user++;
            }


        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            if (xValues->size() > 0) {
                dataForTuning.tune(probeBucket, retrArg, this);
            }
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            if (xValues->size() > 0) {
                dataForTuning.tuneTopk(probeBucket, retrArg, this);
            }

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {

            arg->numLists = probeBucket.numLists;


            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                    continue;

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];

#ifdef TIME_IT
                arg->t.start();
#endif
                // this can go out of the loop. I run incr for all
                if (!queryBatch.initializedQueues) { //preprocess
                    queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                    //                arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);       
                    queryBatch.initializedQueues = true;
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif



                ///////////////////////////////////
                row_type user = queryBatch.startPos;
                int start = queryBatch.startPos * arg->k;
                int end = queryBatch.endPos * arg->k;
                for (row_type i = start; i < end; i += arg->k) {

                    //                std:cout<<"i: "<<i<<" user: "<<user
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

                    col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                    //                col_type* localQueue = arg->queryMatrix->getQueue(user);
                    arg->setQueues(localQueue);

                    runTopK(query, probeBucket, arg);

                    arg->writeHeapToTopk(user);

                    user++;
                }

            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {

            arg->numLists = probeBucket.numLists;

            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }


                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];


#ifdef TIME_IT
                //std::cout<<" INCR-based retrieval running "<<std::endl;
                arg->t.start();
#endif
                // this can go out of the loop. I run incr for all
                if (!queryBatch.initializedQueues) { //preprocess              
                    queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);

                    //                arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);                
                    queryBatch.initializedQueues = true;
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif




                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;

                    col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);

                    //                col_type* localQueue = arg->queryMatrix->getQueue(i);           
                    arg->setQueues(localQueue);
                    arg->queryId = arg->queryMatrix->getId(i);

                    run(query, probeBucket, arg);

                }


            }

        }

    };
}

#endif	/* INCR2_H */

