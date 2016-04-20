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
 * ListsTuneData.h
 *
 *  Created on: Jul 4, 2014
 *      Author: chteflio
 */

#ifndef LISTSTUNEDATA_H_
#define LISTSTUNEDATA_H_

namespace mips {

    struct ListTuneData {
        bool upperActive, lowerActive;
        std::unordered_map<col_type, double> bestTimeForPhi; // for each phi explored (key) contains the best achieved total time (value) 
        std::unordered_map<col_type, std::vector<double> > timeX; //key: phi value, value: runtimes for each sample query
        row_type bestPhi; // the phi value that managed to achieve the best total runtime 
        double bestTime; // the best total runtime corresponding to the bestPhi. 
        col_type* queues;
        row_type t_b_indx;

        inline ListTuneData() : queues(nullptr), upperActive(true), lowerActive(true) {
        }

        inline ~ListTuneData() {
            if (queues != nullptr) {
                delete [] queues;
            }
        }

        inline col_type* getQueue(row_type user, col_type numLists) const {
            return &queues[user * numLists];
        }


        // I drop and recreate the queues for each phi I try out

        inline void preprocess(const std::vector<RetrievalArguments>& retrArg, const ProbeBucket& probeBucket, col_type numLists) {

            if (queues != nullptr) {
                delete[] queues;
                queues = nullptr;
            }
            queues = new col_type[numLists * probeBucket.xValues->size()];

            std::vector<QueueElement> tmp;
            tmp.resize(numLists);


            for (row_type j = 0; j < probeBucket.xValues->size(); ++j) {
                int t = probeBucket.xValues->at(j).i;
                int ind = probeBucket.xValues->at(j).j;

                double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                for (col_type i = 0; i < numLists; ++i) {
                    double value = fabs(query[i]);
                    tmp[i] = QueueElement(value, i);
                }
                std::make_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());

                for (col_type i = numLists; i < retrArg[t].queryMatrix->colNum; ++i) {
                    double value = fabs(query[i]);
                    if (value > tmp.front().data) {
                        std::pop_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                        tmp.pop_back();
                        tmp.emplace_back(value, i);
                        std::push_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                    }
                }

                std::sort(tmp.begin(), tmp.end(), std::greater<QueueElement>());

                for (col_type i = 0; i < numLists; ++i) {
                    queues[j * numLists + i] = tmp[i].id;
                }
            }
        }

        inline void findCutOffPointForList(col_type list, row_type sampleSize, double bestTimeOfPrevPhi, const std::vector<double>* competitorMethod, bool upper) {
            double time;


            for (row_type i = 0; i < sampleSize; ++i) {
                calculateTimeInCutoff(timeX[list], *competitorMethod, i, time);

                if ((i == 0) || bestTimeForPhi[list] > time) {
                    bestTimeForPhi[list] = time;
                }
                if ((bestTimeOfPrevPhi < 0 && i == 0) || bestTime > 1.1 * time) {// take into consideration caching effects
                    bestTime = time;
                    t_b_indx = i;
                    bestPhi = list;
                }
                if (competitorMethod == NULL) { // we have LEMP_C or LEMP_I, competitor method is in practice LENGTH
                    t_b_indx = 0;
                    break;
                }
            }

            if (list == 1) { // inactivate exploration left and right if you have hit borders
                lowerActive = false;
            } else if (list == NUM_LISTS) {
                upperActive = false;
            }


            if (bestTimeOfPrevPhi > 0) {// if this is not the first phi value to be explored
                if (bestTimeForPhi[list] > 1.1 * bestTimeOfPrevPhi) { // if you are more than 1.1 times worse than the previous list, stop exploring

                    if (upper)
                        upperActive = false;
                    else
                        lowerActive = false;
                }
            }
        }

        template<typename R>
        inline void tuneBucketForList(ProbeBucket& probeBucket, col_type numLists, double otherTime, bool upper,
                std::vector<RetrievalArguments>& retrArg, R& retriever) {

            // I use retrArg[0] for all running
            retrArg[0].setIntervals(numLists);

            probeBucket.numLists = numLists;

            row_type sampleSize = probeBucket.xValues->size();


            timeX[numLists - 1] = std::vector<double>();
            timeX[numLists - 1].reserve(sampleSize);

            retrArg[0].tunerTimer.start();
            preprocess(retrArg, probeBucket, numLists);
            retrArg[0].tunerTimer.stop();
            double preprocessTime = retrArg[0].tunerTimer.elapsedTime().nanos() / (probeBucket.xValues->size() == 0 ? 1 : probeBucket.xValues->size());


            retrArg[0].numLists = probeBucket.numLists;
            for (row_type i = 0; i < sampleSize; ++i) {

                int t = probeBucket.xValues->at(i).i;
                int ind = probeBucket.xValues->at(i).j;
                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);


                retrArg[0].results.clear();
                const col_type * localQueue = getQueue(i, numLists);
                retrArg[0].setQueues(localQueue);

                retrArg[0].tunerTimer.start();
                retriever->run(query, probeBucket, &retrArg[0]);
                retrArg[0].tunerTimer.stop();
                timeX[numLists - 1].emplace_back(retrArg[0].tunerTimer.elapsedTime().nanos() + preprocessTime);
            }

            findCutOffPointForList(numLists - 1, sampleSize, otherTime, retrArg[0].competitorMethod, upper);

        }

        template<typename R>
        void tuneBucketForListTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg,
                col_type numLists, double otherTime, bool upper, R& retriever) {

            // I use retrArg[0] for all running
            retrArg[0].setIntervals(numLists);
            probeBucket.numLists = numLists;

            row_type sampleSize = probeBucket.xValues->size();

            timeX[numLists - 1] = std::vector<double>();
            timeX[numLists - 1].reserve(sampleSize);

            retrArg[0].tunerTimer.start();
            preprocess(retrArg, probeBucket, numLists);
            retrArg[0].tunerTimer.stop();
            double preprocessTime = retrArg[0].tunerTimer.elapsedTime().nanos() / (probeBucket.xValues->size() == 0 ? 1 : probeBucket.xValues->size());

            retrArg[0].numLists = probeBucket.numLists;
            for (row_type i = 0; i < sampleSize; ++i) {

                int t = probeBucket.xValues->at(i).i;
                int ind = probeBucket.xValues->at(i).j;
                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                const std::vector<QueueElement>& prevResults = prevBucket.sampleThetas[t].at(ind).results; //just reading
                std::copy(prevResults.begin(), prevResults.end(), retrArg[0].heap.begin());

                std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                const col_type * localQueue = getQueue(i, numLists);
                retrArg[0].setQueues(localQueue);

#if defined(RELATIVE_APPROX) 
                if (retrArg[0].heap.front().data >= 0) {
                    retrArg[0].currEpsilonAppr = (1 + retrArg[0].epsilon);
                } else {
                    retrArg[0].currEpsilonAppr = 1;
                }
#else 
#if defined(ABS_APPROX)               
                retrArg[0].currEpsilonAppr = retrArg[t].queryMatrix->epsilonEquivalents[ind];
#endif
#endif

                retrArg[0].tunerTimer.start();
                retriever->runTopK(query, probeBucket, &retrArg[0]); ///////////////
                retrArg[0].tunerTimer.stop();
                timeX[numLists - 1].emplace_back(retrArg[0].tunerTimer.elapsedTime().nanos() + preprocessTime);
            }
            findCutOffPointForList(numLists - 1, sampleSize, otherTime, retrArg[0].competitorMethod, upper);
        }

        template<typename R>
        inline void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg, R* retriever) {

            upperActive = true;
            lowerActive = true;


            col_type list = prevBucket.numLists; // start with the best numLists of the previous bucket

            // explore an QueryBatchial phi value
            tuneBucketForList(probeBucket, list, -1, true, retrArg, retriever);

            // Here I keep track of the best total runtime 
            double bestTimeOfPrevPhiOnLeft = bestTimeForPhi[list - 1];
            double bestTimeOfPrevPhiOnRight = bestTimeForPhi[list - 1];

            // explore phi values left and right of the QueryBatchial phi value
            int tries = NUM_LISTS / 2;
            for (col_type i = 1; i <= tries; ++i) {

                if (upperActive) {
                    tuneBucketForList(probeBucket, list + i, bestTimeOfPrevPhiOnLeft, true, retrArg, retriever);
                    bestTimeOfPrevPhiOnLeft = bestTimeForPhi[list + i - 1];

                }
                if (lowerActive) {
                    tuneBucketForList(probeBucket, list - i, bestTimeOfPrevPhiOnRight, false, retrArg, retriever);
                    bestTimeOfPrevPhiOnRight = bestTimeForPhi[list - i - 1];
                }
                if (!upperActive && !lowerActive)
                    break;
            }

            if (queues != nullptr) {
                delete [] queues;
                queues = nullptr;
            }


            double value = (t_b_indx == 0 ? -1 : probeBucket.xValues->at(t_b_indx).result);
            probeBucket.setAfterTuning(bestPhi + 1, value);
            retriever->sampleTotalTime = bestTime;

        }

        template<typename R>
        void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg, R* retriever) {
            upperActive = true;
            lowerActive = true;

            col_type list = prevBucket.numLists; // start with the best numLists of the previous bucket
            tuneBucketForListTopk(probeBucket, prevBucket, retrArg, list, -1, true, retriever);

            double previousForUpper = bestTimeForPhi[list - 1];
            double previousForLower = bestTimeForPhi[list - 1];


            int tries = NUM_LISTS / 2;

            for (col_type i = 1; i <= tries; ++i) {

                if (upperActive) {
                    tuneBucketForListTopk(probeBucket, prevBucket, retrArg, list + i, previousForUpper, true, retriever);
                    previousForUpper = bestTimeForPhi[list + i - 1];
                }
                if (lowerActive) {
                    tuneBucketForListTopk(probeBucket, prevBucket, retrArg, list - i, previousForLower, false, retriever);
                    previousForLower = bestTimeForPhi[list - i - 1];
                }
                if (!upperActive && !lowerActive)
                    break;
            }

            if (queues != nullptr) {
                delete [] queues;
                queues = nullptr;
            }

            double value = (t_b_indx == 0 ? -1 : probeBucket.xValues->at(t_b_indx).result);
            probeBucket.setAfterTuning(bestPhi + 1, value);
            retriever->sampleTotalTime = bestTime;
        }

    };




}



#endif /* LISTSTUNEDATA_H_ */
