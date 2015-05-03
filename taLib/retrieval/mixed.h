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
 * mixed.h
 *
 *  Created on: Jun 27, 2014
 *      Author: chteflio
 */

#ifndef MIXED_H_
#define MIXED_H_


namespace ta {

    template<class X> // where X can be C or I
    class LX_Retriever : public Retriever {
    public:
        LengthRetriever plainRetriever;
        X otherRetriever;

        inline LengthRetriever* getLengthRetriever() {
            return &plainRetriever;
        }

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {


#ifdef TIME_IT
            //std::cout<<" LI-based retrieval running "<<std::endl;
#endif
            if (probeBucket.t_b == 1 || (probeBucket.t_b * queryBatch.normL2.first > probeBucket.bucketScanThreshold)) {
                plainRetriever.run(queryBatch, probeBucket, arg);
            } else if (probeBucket.t_b * queryBatch.normL2.second <= probeBucket.bucketScanThreshold) {
                otherRetriever.run(queryBatch, probeBucket, arg);
            } else { // do it per query

#ifdef TIME_IT
                arg->t.start();
#endif

                // this can go out of th loop. I run incr for all
                if (!queryBatch.initializedQueues) { //preprocess
                    queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                    //                     arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);                 
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

                    arg->queryId = arg->queryMatrix->getId(i);
                    if (probeBucket.t_b * query[-1] > probeBucket.bucketScanThreshold) {// do length-based
                        plainRetriever.run(query, probeBucket, arg);
                    } else {

                        col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);

                        //                         col_type* localQueue = arg->queryMatrix->getQueue(i);

                        arg->setQueues(localQueue);
                        otherRetriever.run(query, probeBucket, arg);
                    }

                }
            }



        }

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {

#ifdef TIME_IT
            //std::cout<<" LI-based retrieval running "<<std::endl;
#endif
            if (probeBucket.t_b == 1) {
                plainRetriever.runTopK(queryBatch, probeBucket, arg);
            } else { // do it per query
                arg->numLists = probeBucket.numLists;

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

                    if (probeBucket.t_b * probeBucket.normL2.second > minScore) {// do length-based
                        plainRetriever.runTopK(query, probeBucket, arg);

                    } else {
#ifdef TIME_IT
                        arg->t.start();
#endif
                        // this can go out of th loop. I run incr for all
                        if (!queryBatch.initializedQueues) { //preprocess
                            queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                            //                             arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);
                            queryBatch.initializedQueues = true;
                        }
#ifdef TIME_IT
                        arg->t.stop();
                        arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif
                        col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                        //                         col_type* localQueue = arg->queryMatrix->getQueue(user);
                        arg->setQueues(localQueue);
                        otherRetriever.runTopK(query, probeBucket, arg);

                    }


                    arg->writeHeapToTopk(user);

                    user++;
                }
            }
        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            if (xValues->size() > 0) {

                plainRetriever.xValues = xValues;
                otherRetriever.xValues = xValues;

                plainRetriever.tune(probeBucket, retrArg);
                retrArg[0].competitorMethod = &plainRetriever.sampleTimes;
                otherRetriever.tune(probeBucket, retrArg);

                if (plainRetriever.sampleTotalTime < otherRetriever.dataForTuning.bestTimeX) {
                    probeBucket.setAfterTuning(1, 1);
                } else {
                    double value = (otherRetriever.dataForTuning.t_b_indx == 0 ? -1 : xValues->at(otherRetriever.dataForTuning.t_b_indx).result);
                    probeBucket.setAfterTuning(otherRetriever.dataForTuning.bestNumLists + 1, value);
                }
            }
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type b = retrArg[0].bucketInd;

            otherRetriever.xValues = xValues;

            for (row_type i = 0; i < xValues->size(); i++) {
                int t = xValues->at(i).i;
                int ind = xValues->at(i).j;

                plainRetriever.sampleTimes.push_back(retrArg[0].globalData[b][t][ind].lengthTime);
                plainRetriever.sampleTotalTime += retrArg[0].globalData[b][t][ind].lengthTime;
            }
            retrArg[0].competitorMethod = &plainRetriever.sampleTimes;

            otherRetriever.tuneTopk(probeBucket, retrArg);

//            std::cout<<"LENGTH: "<<plainRetriever.sampleTotalTime<<" INCR: "<<otherRetriever.dataForTuning.bestTimeX<<std::endl;
            
            if (plainRetriever.sampleTotalTime < otherRetriever.dataForTuning.bestTimeX) {
                
                probeBucket.setAfterTuning(1, 1);

            } else {
                double value = (otherRetriever.dataForTuning.t_b_indx == 0 ? -1 : xValues->at(otherRetriever.dataForTuning.t_b_indx).result);
                probeBucket.setAfterTuning(otherRetriever.dataForTuning.bestNumLists + 1, value);

                //                std::cout<<(int)otherRetriever.dataForTuning.bestNumLists + 1<<std::endl;

            }
        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {

#ifdef TIME_IT
            //std::cout<<" LI-based retrieval running "<<std::endl;
#endif
            if (probeBucket.t_b == 1) {
                plainRetriever.runTopK(probeBucket, arg);
            } else { // do it per query
                arg->numLists = probeBucket.numLists;

                for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                    if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                        continue;

                    QueryBucket_withTuning& queryBatch = arg->queryBatches[q];

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

                        if (probeBucket.t_b * probeBucket.normL2.second > minScore) {// do length-based
                            plainRetriever.runTopK(query, probeBucket, arg);

                        } else {
#ifdef TIME_IT
                            arg->t.start();
#endif
                            // this can go out of th loop. I run incr for all
                            if (!queryBatch.initializedQueues) { //preprocess
                                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                                //                             arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);
                                queryBatch.initializedQueues = true;
                            }
#ifdef TIME_IT
                            arg->t.stop();
                            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif
                            col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                            //                         col_type* localQueue = arg->queryMatrix->getQueue(user);
                            arg->setQueues(localQueue);
                            otherRetriever.runTopK(query, probeBucket, arg);

                        }


                        arg->writeHeapToTopk(user);

                        user++;
                    }
                }
            }
        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {


            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];

#ifdef TIME_IT
                //std::cout<<" LI-based retrieval running "<<std::endl;
#endif
                if (probeBucket.t_b == 1 || (probeBucket.t_b * queryBatch.normL2.first > probeBucket.bucketScanThreshold)) {
                    plainRetriever.run(queryBatch, probeBucket, arg);
                } else if (probeBucket.t_b * queryBatch.normL2.second <= probeBucket.bucketScanThreshold) {
                    otherRetriever.run(queryBatch, probeBucket, arg);
                } else { // do it per query

#ifdef TIME_IT
                    arg->t.start();
#endif

                    // this can go out of th loop. I run incr for all
                    if (!queryBatch.initializedQueues) { //preprocess
                        queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                        //                     arg->queryMatrix->preprocessQueues(queryBatch.startPos, queryBatch.endPos);                 
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

                        arg->queryId = arg->queryMatrix->getId(i);
                        if (probeBucket.t_b * query[-1] > probeBucket.bucketScanThreshold) {// do length-based
                            plainRetriever.run(query, probeBucket, arg);
                        } else {

                            col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);

                            //                         col_type* localQueue = arg->queryMatrix->getQueue(i);

                            arg->setQueues(localQueue);
                            otherRetriever.run(query, probeBucket, arg);
                        }

                    }
                }
            }

        }

    };



}


#endif /* MIXED_H_ */
