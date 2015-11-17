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

        LX_Retriever() = default;
        ~LX_Retriever() = default;

        inline LengthRetriever* getLengthRetriever() {
            return &plainRetriever;
        }

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
            std::cout << "You should not call run for single query on LEMP_LI" << std::endl;
            exit(1);
        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            if (probeBucket.xValues->size() > 0) {

                plainRetriever.tune(probeBucket, prevBucket, retrArg);
                retrArg[0].competitorMethod = &plainRetriever.sampleTimes;
                otherRetriever.tune(probeBucket, prevBucket, retrArg);

                if (plainRetriever.sampleTotalTime < otherRetriever.dataForTuning->bestTime) {
                    probeBucket.setAfterTuning(1, 1);
                } else {
                    double value = (otherRetriever.dataForTuning->t_b_indx == 0 ? -1 : probeBucket.xValues->at(otherRetriever.dataForTuning->t_b_indx).result);
                    probeBucket.setAfterTuning(otherRetriever.dataForTuning->bestPhi + 1, value);
                }
            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            row_type sampleSize = (probeBucket.xValues != nullptr ? probeBucket.xValues->size() : 0);

            if (sampleSize > 0) {
                plainRetriever.sampleTimes.reserve(sampleSize);

                for (row_type i = 0; i < sampleSize; ++i) {
                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;

                    plainRetriever.sampleTimes.push_back(probeBucket.sampleThetas->at(t)[ind].lengthTime);
                    plainRetriever.sampleTotalTime += probeBucket.sampleThetas->at(t)[ind].lengthTime;
                }
                retrArg[0].competitorMethod = &plainRetriever.sampleTimes;

                otherRetriever.tuneTopk(probeBucket, prevBucket, retrArg);

                if (plainRetriever.sampleTotalTime < otherRetriever.dataForTuning->bestTime) {
                    probeBucket.setAfterTuning(1, 1);
                } else {
                    double value = (otherRetriever.dataForTuning->t_b_indx == 0 ? -1 : probeBucket.xValues->at(otherRetriever.dataForTuning->t_b_indx).result);
                    probeBucket.setAfterTuning(otherRetriever.dataForTuning->bestPhi + 1, value);
                }

            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const {


            if (probeBucket.t_b == 1) {
                plainRetriever.runTopK(probeBucket, arg);
            } else { // do it per query
                arg->numLists = probeBucket.numLists;

                for (auto& queryBatch : arg->queryBatches) {

                    if (queryBatch.isWorkDone())
                        continue;

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
                        double minScoreAppr = minScore;

#ifdef RELATIVE_APPROX
                        if (minScoreAppr >= 0) {
                            minScoreAppr *= (1 + arg->gamma);
                            arg->currGammaAppr = (1 + arg->gamma);
                        } else {
                            minScoreAppr *= (1 - arg->gamma);
                            arg->currGammaAppr = (1 - arg->gamma);
                        }


#else 
#ifdef         ABS_APPROX             
                        minScoreAppr += arg->queryMatrix->gammaEquivalents[user];
                        arg->currGammaAppr = arg->queryMatrix->gammaEquivalents[user];

#endif
#endif


                        if (probeBucket.normL2.second < minScoreAppr) {// skip this bucket and all other buckets
                            queryBatch.inactivateQuery(user);
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
                            if (!queryBatch.hasInitializedQueues()) { //preprocess
                                queryBatch.preprocess(*(arg->queryMatrix), arg->maxLists);
                            }
#ifdef TIME_IT
                            arg->t.stop();
                            arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif
                            col_type* localQueue = queryBatch.getQueue(user - queryBatch.startPos, arg->maxLists);
                            arg->setQueues(localQueue);
                            otherRetriever.runTopK(query, probeBucket, arg);

                        }


                        arg->writeHeapToTopk(user);
                        user++;
                    }
                }
            }
        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const {
            arg->numLists = probeBucket.numLists;

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }

                if (probeBucket.t_b == 1 || (probeBucket.t_b * queryBatch.minLength() > probeBucket.bucketScanThreshold)) {
                    plainRetriever.run(queryBatch, probeBucket, arg);
                } else if (probeBucket.t_b * queryBatch.maxLength() <= probeBucket.bucketScanThreshold) {
                    otherRetriever.run(queryBatch, probeBucket, arg);
                } else { // do it per query

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

                        arg->queryId = arg->queryMatrix->getId(i);
                        if (probeBucket.t_b * query[-1] > probeBucket.bucketScanThreshold) {// do length-based
                            plainRetriever.run(query, probeBucket, arg);
                        } else {

                            col_type* localQueue = queryBatch.getQueue(i - queryBatch.startPos, arg->maxLists);
                            arg->setQueues(localQueue);
                            otherRetriever.run(query, probeBucket, arg);
                        }

                    }
                }
            }

        }

        inline virtual void cleanupAfterTuning() {
            otherRetriever.cleanupAfterTuning();
        }

    };



}


#endif /* MIXED_H_ */
