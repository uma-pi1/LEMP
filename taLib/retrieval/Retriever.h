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


namespace ta {

    //typedef boost::shared_ptr< std::vector<QueueElement> > xValues_ptr;

    /*
     * Responsible for naive and length-based retrieval
     */
    class Retriever {
    public:

        xValues_ptr xValues; // data: theta_b(q) id: sampleId
        std::vector<double> sampleTimes;
        double sampleTotalTime;
        LEMP_Method name;

        Retriever() : sampleTotalTime(0) {
        };

        Retriever(LEMP_Method name) : sampleTotalTime(0), name(name) {
        };

        virtual ~Retriever() {
        };

        inline void calculateTimeInCutoff(const std::vector<double>& method1, const std::vector<double>& method2, row_type cutoffSample, double& time) {
            if (cutoffSample == 0) {
                time = 0;
                for (row_type i = 0; i < method1.size(); i++) {
                    time += method1[i];
                }
            } else {
                //std::cout<<"I should not be here"<<std::endl;
                time -= method1[cutoffSample - 1];
                time += method2[cutoffSample - 1]; // other can be LENGTH for example
            }
        }

        inline void findCutOffPoint(const std::vector<double>& method1, const std::vector<double>& method2, double& bestTime, row_type& best_t_b_ind) {
            double time;
            for (row_type i = 0; i < method1.size(); i++) {
                calculateTimeInCutoff(method1, method2, i, time);

                if ((i == 0) || bestTime > time) {
                    bestTime = time;
                    best_t_b_ind = i;
                }

            }
        }

        inline void findCutOffPointInv(const std::vector<double>& method1, const std::vector<double>& method2, double& bestTime, row_type& best_t_b_ind) {
            double time = bestTime;
            for (int i = method1.size() - 1; i >= 0; i--) {
                time -= method1[i];
                time += method2[i];
                if (bestTime > time) {
                    bestTime = time;
                    best_t_b_ind = i;
                }
            }
        }

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        // this is to be used only by the 1st bucket in Row-Top-k. It just initializes the top-k elements

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch,
                ProbeBucket& probeBucket, RetrievalArguments* arg) {

            //            QueueElement* point = &(arg->topkResults[queryBatch.startPos * arg->k]);
            //            row_type p = 0;
            //
            //            for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {
            //
            //                const double* query = arg->queryMatrix->getMatrixRowPtr(i);
            //                arg->queryId = arg->queryMatrix->getId(i);
            //
            //
            //                for (row_type j = probeBucket.startPos; j < probeBucket.endPos; j++) {
            //                    arg->comparisons++;
            //                    double ip = arg->probeMatrix->innerProduct(j, query);
            //
            //                    point[p] = QueueElement(ip, arg->probeMatrix->getId(j));
            //                    p++;
            //
            //                    if (arg->worstMinScore > ip) {
            //                        arg->worstMinScore = ip;
            //                    }
            //
            //                }
            //            }

            ////////////////
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

            //////////////
        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {


            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                    continue;

                row_type user = arg->queryBatches[q].startPos;
                int start = arg->queryBatches[q].startPos * arg->k;
                int end = arg->queryBatches[q].endPos * arg->k;
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

        inline void sampling(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            xValues_ptr ptr(new std::vector<MatItem> ());
            xValues = ptr;
            row_type sampleSize;
            // find the queries
            double thresForQ = probeBucket.bucketScanThreshold;
            std::vector<row_type> activeQueries(retrArg.size());

            for (int t = 0; t < retrArg.size(); t++) {
                std::vector<QueueElement>::const_iterator up = std::lower_bound(retrArg[t].queryMatrix->lengthInfo.begin(),
                        retrArg[t].queryMatrix->lengthInfo.end(), QueueElement(thresForQ, 0), std::greater<QueueElement>());
                activeQueries[t] = up - retrArg[t].queryMatrix->lengthInfo.begin();

                probeBucket.activeQueries += activeQueries[t];
            }


            if (probeBucket.activeQueries < LOWER_LIMIT_PER_BUCKET * 3) {
                sampleSize = 0;
            } else {
                sampleSize = 0.02 * probeBucket.activeQueries;

                if (sampleSize > UPPER_LIMIT_PER_BUCKET) {
                    sampleSize = UPPER_LIMIT_PER_BUCKET;
                }

                if (sampleSize < LOWER_LIMIT_PER_BUCKET) {// if very few elements qualify for this bucket, it does not pay off to tune.
                    sampleSize = LOWER_LIMIT_PER_BUCKET;
                }
            }

            rg::Random32& random = retrArg[0].random;


            if (sampleSize > 0) {
                xValues->reserve(sampleSize);
                sampleSize /= retrArg.size();

                for (int t = 0; t < retrArg.size(); t++) {
                    std::vector<row_type> sampleIndx = rg::sample(random, sampleSize, activeQueries[t]);

                    for (row_type i = 0; i < sampleIndx.size(); i++) {
                        double theta_b_q = probeBucket.bucketScanThreshold / retrArg[t].queryMatrix->getVectorLength(sampleIndx[i]);
                        xValues->push_back(MatItem(theta_b_q, t, sampleIndx[i]));
                    }
                }
                std::sort(xValues->begin(), xValues->end(), std::less<MatItem>());
            }
        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
        }

        inline void setup_xValues_topk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            xValues_ptr ptr(new std::vector<MatItem>());
            xValues = ptr;
            row_type b = retrArg[0].bucketInd;

            for (int t = 0; t < retrArg.size(); t++) {

                row_type lastInd = retrArg[0].globalData[b][t].size();
                unordered_map<row_type, GlobalTopkTuneData >::iterator it = retrArg[0].globalData[b][t].begin();

                while (lastInd > 0) {
                    double localTheta = retrArg[0].globalData[b - 1][t][it->first].results.front().data;
                    localTheta *= (localTheta > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);
                    xValues->push_back(MatItem(localTheta, t, it->first));
                    lastInd--;
                    ++it;
                }
            }
            std::sort(xValues->begin(), xValues->end(), std::less<MatItem>());
        }

    };

    class LengthRetriever : public Retriever {

        /*
         * scans itemMatrix from position start to position end for inner products above args.theta. Method length
         */
        inline void length(const double *query, row_type start, row_type end, RetrievalArguments* arg) {

            for (row_type j = start; j < end; j++) {

                double* item = arg->probeMatrix->getMatrixRowPtr(j);

                double len = query[-1] * item[-1];

                if (len < arg->theta) { // stop scanning for this user
                    break;
                }

                arg->comparisons++;

                double ip = len * arg->probeMatrix->cosine(j, query);

                if (ip >= arg->theta) {
                    arg->results.push_back(MatItem(ip, arg->queryId, arg->probeMatrix->getId(j)));
                }
            }
        }





    public:
        double t_b_low, t_b_high;

        inline LengthRetriever() : Retriever(LEMP_L) {
        }

        /*
         * scans itemMatrix from position start to position end for inner products above args.theta. Method: Naive
         */
        inline void naive(const double *query, row_type start, row_type end, RetrievalArguments* arg) {

            for (row_type j = start; j < end; j++) {
                arg->comparisons++;
                double ip = arg->probeMatrix->innerProduct(j, query);

                if (ip >= arg->theta) {
                    arg->results.push_back(MatItem(ip, arg->queryId, arg->probeMatrix->getId(j)));
                }
            }
        }

        inline void naiveTopk(const double *query, row_type start, row_type end, RetrievalArguments* arg) {
            double minScore = arg->heap.front().data;

            for (row_type j = start; j < end; j++) {
                arg->comparisons++;
                double ip = arg->probeMatrix->innerProduct(j, query);

                if (ip > minScore) {
                    std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    arg->heap.pop_back();
                    arg->heap.push_back(QueueElement(ip, arg->probeMatrix->getId(j)));
                    std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    minScore = arg->heap.front().data;
                }
            }
        }

        /*
         * scans itemMatrix from position start to position end for inner products above current minScore. Method length
         * assumes that  results already contains args.k elements
         */
        inline void lengthTopk(const double *query, row_type start, row_type end, RetrievalArguments* arg) {

            double minScore = arg->heap.front().data;
            for (row_type j = start; j < end; j++) {

                double* item = arg->probeMatrix->getMatrixRowPtr(j);

                if (item[-1] < minScore) { // stop scanning for this user
                    break;
                }
                arg->comparisons++;

                double ip = item[-1] * arg->probeMatrix->cosine(j, query);

                if (ip > minScore) {

                    std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    arg->heap.pop_back();
                    arg->heap.push_back(QueueElement(ip, arg->probeMatrix->getId(j)));
                    std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                    minScore = arg->heap.front().data;
                }
            }

        }

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
#ifdef TIME_IT
            //std::cout<<" LENGTH-based retrieval running "<<std::endl;
            arg->t.start();
#endif
            if (query[-1] * probeBucket.normL2.first < arg->theta) { // LENGTH
                length(query, probeBucket.startPos, probeBucket.endPos, arg);
            } else {// NAIVE
                naive(query, probeBucket.startPos, probeBucket.endPos, arg);
            }
#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
#ifdef TIME_IT
            //std::cout<<" LENGTH-based retrieval running "<<std::endl;
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

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
#ifdef TIME_IT
            //std::cout<<" LENGTH-based retrieval running "<<std::endl;
            arg->t.start();
#endif
            lengthTopk(query, probeBucket.startPos, probeBucket.endPos, arg);
#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {


#ifdef TIME_IT
            //std::cout<<" LENGTH-based retrieval running "<<std::endl;
            arg->t.start();
#endif

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

                //                arg->moveTopkToHeap(i);

                // get the minScore
                //                double minScore = arg->heap.front().data;

                if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                    queryBatch.inactiveQueries[user - queryBatch.startPos] = true;
                    queryBatch.inactiveCounter++;
                    user++;
                    continue;
                }

                arg->moveTopkToHeap(i);

                arg->queryId = arg->queryMatrix->getId(user);
                lengthTopk(query, probeBucket.startPos, probeBucket.endPos, arg);

                arg->writeHeapToTopk(user);

                user++;
            }

#ifdef TIME_IT
            arg->t.stop();
            arg->lengthTime += arg->t.elapsedTime().nanos();
#endif


        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {
            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                    continue;

#ifdef TIME_IT
                //std::cout<<" LENGTH-based retrieval running "<<std::endl;
                arg->t.start();
#endif

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
                    lengthTopk(query, probeBucket.startPos, probeBucket.endPos, arg);

                    arg->writeHeapToTopk(user);

                    user++;
                }

#ifdef TIME_IT
                arg->t.stop();
                arg->lengthTime += arg->t.elapsedTime().nanos();
#endif
            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {


            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];

#ifdef TIME_IT
                //std::cout<<" LENGTH-based retrieval running "<<std::endl;
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

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            //case 1: all sample queries from the same queryMatrix
            // just call this function with retrArg[0]
            sampleTimes.resize(xValues->size());
            for (row_type i = 0; i < xValues->size(); i++) {

                int t = xValues->at(i).i;
                int ind = xValues->at(i).j;
                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                retrArg[t].tunerTimer.start();
                run(query, probeBucket, &retrArg[t]);
                retrArg[t].tunerTimer.stop();
                sampleTimes[i] = retrArg[t].tunerTimer.elapsedTime().nanos();
                sampleTotalTime += sampleTimes[i];

            }


        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
        }
   };




}


#endif /* RETRIEVER_H_ */
