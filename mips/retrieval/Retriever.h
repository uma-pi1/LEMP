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
 * naive.h
 *
 *  Created on: Jun 25, 2014
 *      Author: chteflio
 */

#ifndef RETRIEVER_H_
#define RETRIEVER_H_


namespace mips {

    /*
     * Responsible for naive and length-based retrieval
     */
    class Retriever {
    public:

        std::vector<double> sampleTimes;
        double sampleTotalTime = 0;

        Retriever() {
        }

        virtual ~Retriever() {
        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        // this is to be used only by the 1st bucket in Row-Top-k. It just initializes the top-k elements

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const {


            for (auto& queryBatch : arg->queryBatches) {

                //                if (queryBatch.isWorkDone())
                //                    continue;

                row_type user = queryBatch.startPos;
                int start = queryBatch.startPos * arg->k;
                int end = queryBatch.endPos * arg->k;


                for (row_type i = start; i < end; i += arg->k) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                    arg->queryId = arg->queryMatrix->getId(user);

                    for (row_type j = probeBucket.startPos; j < probeBucket.endPos; j++) {
                        arg->comparisons++;
                        double ip = arg->probeMatrix->innerProduct(j, query);

                        arg->heap[j] = QueueElement(ip, arg->probeMatrix->getId(j));

                    }
                    std::make_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());

                    if (arg->worstMinScore > arg->heap.front().data) {
                        arg->worstMinScore = arg->heap.front().data;
                    }

                    arg->writeHeapToTopk(user);
                    user++;
                }
            }

        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }


    };

    class LengthRetriever : public Retriever {
    public:

        inline LengthRetriever() {
        }

        inline ~LengthRetriever() {
        }

        /*
         * scans itemMatrix from position start to position end for inner products above args.theta. Method: Naive
         */
        inline void naive(const double *query, row_type start, row_type end, RetrievalArguments* arg) const {

            for (row_type j = start; j < end; ++j) {
                arg->comparisons++;
                double ip = arg->probeMatrix->innerProduct(j, query);

                if (ip >= arg->theta) {
                    arg->results.emplace_back(ip, arg->queryId, arg->probeMatrix->getId(j));
                }
            }
        }

        inline void naiveTopk(const double *query, row_type start, row_type end, RetrievalArguments* arg) const {
            double minScore = arg->heap.front().data;

            for (row_type j = start; j < end; ++j) {
                arg->comparisons++;
                double ip = arg->probeMatrix->innerProduct(j, query);

                if (ip > minScore) {
                    std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    arg->heap.pop_back();
                    arg->heap.emplace_back(ip, arg->probeMatrix->getId(j));
                    std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    minScore = arg->heap.front().data;
                }
            }
        }

        inline void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const {
#ifdef TIME_IT
            arg->t.start();
#endif
            if (query[-1] * probeBucket.normL2.first < arg->theta) { // LENGTH
                for (row_type j = probeBucket.startPos; j < probeBucket.endPos; ++j) {

                    double* item = arg->probeMatrix->getMatrixRowPtr(j);

                    double len = query[-1] * item[-1];

                    if (len < arg->theta) { // stop scanning for this user
                        break;
                    }

                    arg->comparisons++;

                    double ip = len * arg->probeMatrix->cosine(j, query);

                    if (ip >= arg->theta) {
                        arg->results.emplace_back(ip, arg->queryId, arg->probeMatrix->getId(j));
                    }
                }
            } else {// NAIVE
                naive(query, probeBucket.startPos, probeBucket.endPos, arg);

            }
#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const {
#ifdef TIME_IT
            arg->t.start();
#endif

            for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                const double* query = arg->queryMatrix->getMatrixRowPtr(i);
                if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                    break;

                arg->queryId = arg->queryMatrix->getId(i);
                run(query, probeBucket, arg);
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg)const {
#ifdef TIME_IT
            arg->t.start();
#endif

            double minScore = arg->heap.front().data;
            double minScoreAppr = minScore;

#if defined(RELATIVE_APPROX)     
            minScoreAppr *= arg->currEpsilonAppr;
#else 
#if defined(ABS_APPROX)           
            minScoreAppr += arg->currEpsilonAppr;
#endif
#endif

            for (row_type j = probeBucket.startPos; j < probeBucket.endPos; ++j) {

                double* item = arg->probeMatrix->getMatrixRowPtr(j);

                if (item[-1] < minScoreAppr) { // stop scanning for this user
                    break;
                }
                arg->comparisons++;

                double ip = item[-1] * arg->probeMatrix->cosine(j, query);

                if (ip > minScore) {

                    std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    arg->heap.pop_back();
                    arg->heap.emplace_back(ip, arg->probeMatrix->getId(j));
                    std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    minScore = arg->heap.front().data;

                    minScoreAppr = minScore;

#if defined(RELATIVE_APPROX)     
                    minScoreAppr *= arg->currEpsilonAppr;
#else 
#if defined(ABS_APPROX)     
                    minScoreAppr += arg->currEpsilonAppr;
#endif
#endif

                }
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const {


#ifdef TIME_IT
            arg->t.start();
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
                double minScoreAppr = minScore;

#if defined(RELATIVE_APPROX) 
                if (minScoreAppr >= 0) {
                    minScoreAppr *= (1 + arg->epsilon);
                    arg->currEpsilonAppr = (1 + arg->epsilon);
                } else {
//                    minScoreAppr *= (1 - arg->epsilon);
                    arg->currEpsilonAppr = 1;//(1 - arg->epsilon);
                }
#else 
#if defined(ABS_APPROX)
                minScoreAppr += arg->queryMatrix->epsilonEquivalents[user];
                arg->currEpsilonAppr = arg->queryMatrix->epsilonEquivalents[user];
#endif
#endif


                if (probeBucket.normL2.second < minScoreAppr) {// skip this bucket and all other buckets
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

#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif


        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg)const {

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.isWorkDone())
                    continue;
#ifdef TIME_IT
                arg->t.start();
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
                    double minScoreAppr = minScore;

#if defined(RELATIVE_APPROX)
                    if (minScoreAppr >= 0) {
                        minScoreAppr *= (1 + arg->epsilon);
                        arg->currEpsilonAppr = (1 + arg->epsilon);
                    } else {
//                        minScoreAppr *= (1 - arg->epsilon);
                        arg->currEpsilonAppr = 1;//(1 - arg->epsilon);
                    }
#else 
#if defined(ABS_APPROX)     
                    minScoreAppr += arg->queryMatrix->epsilonEquivalents[user];
                    arg->currEpsilonAppr = arg->queryMatrix->epsilonEquivalents[user];
#endif
#endif


                    if (probeBucket.normL2.second < minScoreAppr) {// skip this bucket and all other buckets
                        queryBatch.inactivateQuery(user);
                        user++;

#ifdef         RELATIVE_APPROX                     

                        arg->moveTopkToHeap(i);
                        arg->totalErrorAfterResults += approximateRelError(*arg);
#endif
                        continue;
                    }

                    arg->moveTopkToHeap(i);

                    arg->queryId = arg->queryMatrix->getId(user);
                    runTopK(query, probeBucket, arg);

                    arg->writeHeapToTopk(user);

                    user++;
                }

#ifdef TIME_IT
                arg->t.stop();
                arg->lengthTime += arg->t.elapsedTime().nanos();
#endif

            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }
#ifdef TIME_IT
                arg->t.start();
#endif

                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;

                    arg->queryId = arg->queryMatrix->getId(i);


                    run(query, probeBucket, arg);
                }
#ifdef TIME_IT
                arg->t.stop();
                arg->lengthTime += arg->t.elapsedTime().nanos();
#endif

            }
        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            //case 1: all sample queries from the same queryMatrix
            // just call this function with retrArg[0]
            row_type sampleSize = probeBucket.xValues->size();

            if (sampleSize > 0) {
                sampleTimes.reserve(sampleSize);
                for (row_type i = 0; i < sampleSize; ++i) {

                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    retrArg[t].tunerTimer.start();
                    run(query, probeBucket, &retrArg[t]);
                    retrArg[t].tunerTimer.stop();
                    sampleTimes.emplace_back(retrArg[t].tunerTimer.elapsedTime().nanos());
                    sampleTotalTime += sampleTimes[i];

                }
            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }



        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
        }
    };




}


#endif /* RETRIEVER_H_ */
