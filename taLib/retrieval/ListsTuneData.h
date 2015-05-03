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

namespace ta {

    struct ListTuneData {
        bool active;
        bool upperActive, lowerActive;
        unordered_map<col_type, double> bestTimeForList;
        unordered_map<col_type, std::vector<double> > timeX; //1d:list 2d: query
        row_type bestNumLists;
        double bestTimeX;
        col_type* queues;
        row_type t_b_indx;

        inline ListTuneData() : queues(0), active(true), upperActive(true), lowerActive(true) {

        }

        inline ~ListTuneData() {
            if (!queues) {
                delete [] queues;
            }
        }

        inline col_type* getQueue(row_type user, col_type numLists) const {
            return &queues[user * numLists];
        }
        
       

        inline void preprocess(std::vector<RetrievalArguments>& retrArg, col_type numLists, Retriever* retriever) {


            if (queues != 0) {
                delete[] queues;
                queues = 0;
            }
            queues = new col_type[numLists * retriever->xValues->size()];

            std::vector<QueueElement> tmp;
            tmp.resize(numLists);


            for (row_type j = 0; j < retriever->xValues->size(); j++) {
                int t = retriever->xValues->at(j).i;
                int ind = retriever->xValues->at(j).j;

                double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                for (col_type i = 0; i < numLists; i++) {
                    double value = fabs(query[i]);
                    tmp[i] = QueueElement(value, i);
                }
                std::make_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());

                for (col_type i = numLists; i < retrArg[t].queryMatrix->colNum; i++) {
                    double value = fabs(query[i]);


                    if (value > tmp.front().data) {
                        std::pop_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                        tmp.pop_back();
                        tmp.push_back(QueueElement(value, i));
                        std::push_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                    }
                }

                std::sort(tmp.begin(), tmp.end(), std::greater<QueueElement>());

                for (col_type i = 0; i < numLists; i++) {
                    queues[j * numLists + i] = tmp[i].id;
                }
            }

        }

        inline void findCutOffPointForList(col_type list, double otherTime, std::vector<double>* competitorMethod, bool upper, Retriever* retriever) {
            double time;


            for (row_type i = 0; i < retriever->xValues->size(); i++) {
                retriever->calculateTimeInCutoff(timeX[list], *competitorMethod, i, time);


                if ((i == 0) || bestTimeForList[list] > time) {
                    bestTimeForList[list] = time;
                }

                if ((otherTime < 0 && i == 0) || bestTimeX > 1.1 * time) {// take into consideration caching effects
                    bestTimeX = time;
                    t_b_indx = i;
                    bestNumLists = list;
                }


                if (competitorMethod == NULL) { // we have LEMP_C or LEMP_I
                    t_b_indx = 0;
                    break;
                }

            }

            if (list == 1) {
                lowerActive = false;
            } else if (list == NUM_LISTS) {
                upperActive = false;
            }

            if (otherTime > 0) {
                if (bestTimeForList[list] > 1.1 * otherTime) {
                    if (upper)
                        upperActive = false;
                    else
                        lowerActive = false;
                }
            }
        }

        inline void tuneBucketForList(ProbeBucket& probeBucket, col_type numLists, double otherTime, bool upper,
                std::vector<RetrievalArguments>& retrArg, Retriever* retriever) {

//            std::cout<<"list "<<(int)numLists<<std::endl;

            // I use retrArg[0] for all running
            retrArg[0].setIntervals(numLists);

            probeBucket.numLists = numLists;
            

            timeX[numLists - 1] = std::vector<double>();
            timeX[numLists - 1].resize(retriever->xValues->size());

            retrArg[0].tunerTimer.start();
            preprocess(retrArg, numLists, retriever);
            retrArg[0].tunerTimer.stop();
            double preprocessTime = retrArg[0].tunerTimer.elapsedTime().nanos() / (retriever->xValues->size() == 0 ? 1 : retriever->xValues->size());


            retrArg[0].numLists = probeBucket.numLists;
            for (row_type i = 0; i < retriever->xValues->size(); i++) {

                int t = retriever->xValues->at(i).i;
                int ind = retriever->xValues->at(i).j;
                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);


               retrArg[0].results.clear();
                const col_type * localQueue = getQueue(i, numLists);
                retrArg[0].setQueues(localQueue);

                retrArg[0].tunerTimer.start();
                retriever->run(query, probeBucket, &retrArg[0]);
                retrArg[0].tunerTimer.stop();
                timeX[numLists - 1][i] = retrArg[0].tunerTimer.elapsedTime().nanos() + preprocessTime;
            }

            findCutOffPointForList(numLists - 1, otherTime, retrArg[0].competitorMethod, upper, retriever);

        }

        void tuneBucketForListTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg,
                col_type numLists, double otherTime, bool upper, Retriever* retriever) {

//            std::cout<<"list "<<(int)numLists<<std::endl;
            // I use retrArg[0] for all running
            retrArg[0].setIntervals(numLists);
            probeBucket.numLists = numLists;

            timeX[numLists - 1] = std::vector<double>();
            timeX[numLists - 1].resize(retriever->xValues->size());

            retrArg[0].tunerTimer.start();
            preprocess(retrArg, numLists, retriever);
            retrArg[0].tunerTimer.stop();
            double preprocessTime = retrArg[0].tunerTimer.elapsedTime().nanos() / (retriever->xValues->size() == 0 ? 1 : retriever->xValues->size());
            //std::cout<<"list"<<(int)numLists<<" preprocess: "<<preprocessTime<<std::endl;
            row_type b = retrArg[0].bucketInd;


            retrArg[0].numLists = probeBucket.numLists;
            for (row_type i = 0; i < retriever->xValues->size(); i++) {

                int t = retriever->xValues->at(i).i;
                int ind = retriever->xValues->at(i).j;
                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                std::vector<QueueElement>& prevResults = retrArg[0].globalData[b - 1][t][ind].results; //just reading
      
                retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());
          

                const col_type * localQueue = getQueue(i, numLists);
                retrArg[0].setQueues(localQueue);
                
                

                retrArg[0].tunerTimer.start();
                retriever->runTopK(query, probeBucket, &retrArg[0]);///////////////
                retrArg[0].tunerTimer.stop();
                timeX[numLists - 1][i] = retrArg[0].tunerTimer.elapsedTime().nanos() + preprocessTime;
            }
            findCutOffPointForList(numLists - 1, otherTime, retrArg[0].competitorMethod, upper, retriever);
        }

        inline void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg, Retriever* retriever) {

            active = true;
            upperActive = true;
            lowerActive = true;
            
            col_type list = (retrArg[0].prevList < 0 ? 1 : retrArg[0].prevList);

//            std::cout<<"Start with lists: "<<(int) list<<std::endl;

            tuneBucketForList(probeBucket, list, -1, true, retrArg, retriever);

            double previousForUpper = bestTimeForList[list - 1];
            double previousForLower = bestTimeForList[list - 1];

            int tries = NUM_LISTS / 2;
            for (col_type i = 1; i <= tries; i++) {

                if (upperActive) {
                    tuneBucketForList(probeBucket, list + i, previousForUpper, true, retrArg, retriever);
                    previousForUpper = bestTimeForList[list + i - 1];
                }
                if (lowerActive) {
                    tuneBucketForList(probeBucket, list - i, previousForLower, false, retrArg, retriever);
                    previousForLower = bestTimeForList[list - i - 1];
                }
                if (!upperActive && !lowerActive)
                    break;
            }


            retrArg[0].prevList = bestNumLists + 1; // ready for the next bucket to start from

            if (queues != 0) {
                delete [] queues;
                queues = 0;
            }


            double value = (t_b_indx == 0 ? -1 : retriever->xValues->at(t_b_indx).result);
            probeBucket.setAfterTuning(bestNumLists + 1, value);

        }

        void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg, Retriever* retriever) {


            active = true;
            upperActive = true;
            lowerActive = true;
        

            col_type list = (retrArg[0].prevList < 0 ? 1 : retrArg[0].prevList);

           
            tuneBucketForListTopk(probeBucket, retrArg, list, -1, true, retriever);

            double previousForUpper = bestTimeForList[list - 1];
            double previousForLower = bestTimeForList[list - 1];
//            std::cout<<"list: "<<(int) list<<" time: "<<previousForLower<<std::endl;


            int tries = NUM_LISTS / 2;

            for (col_type i = 1; i <= tries; i++) {

                if (upperActive) {
                    tuneBucketForListTopk(probeBucket, retrArg, list + i, previousForUpper, true, retriever);
                    previousForUpper = bestTimeForList[list + i - 1];
                }
                if (lowerActive) {
                    tuneBucketForListTopk(probeBucket, retrArg, list - i, previousForLower, false, retriever);
                    previousForLower = bestTimeForList[list - i - 1];
                }
                if (!upperActive && !lowerActive)
                    break;
            }

            retrArg[0].prevList = bestNumLists + 1; // ready for the next bucket to start from

            if (queues != 0) {
                delete [] queues;
                queues = 0;
            }

            double value = (t_b_indx == 0 ? -1 : retriever->xValues->at(t_b_indx).result);
            probeBucket.setAfterTuning(bestNumLists + 1, value);
        }

    };




}



#endif /* LISTSTUNEDATA_H_ */
