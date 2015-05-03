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
 * TuneTopk.h
 *
 *  Created on: Jul 7, 2014
 *      Author: chteflio
 */

#ifndef TUNETOPK_H_
#define TUNETOPK_H_
using boost::unordered_set;
using boost::unordered_map;


namespace ta {

    std::pair<row_type, row_type> findSampleForTuningTopk(std::vector<ProbeBucket>& probeBuckets, std::vector<RetrievalArguments>& retrArg) {
        row_type activeBuckets = 0, bucketsForInit = 0;

        // 1: Bucket 2: valid sample points for bucket
        retrArg[0].globalData.resize(probeBuckets.size());
        row_type id;
        std::pair<row_type, row_type> p;


        LengthRetriever plainRetriever;

        row_type sampleSize = 0;
        for (int t = 0; t < retrArg.size(); t++) {
            sampleSize += retrArg[t].queryMatrix->rowNum;
        }

        sampleSize *= 0.02;

        if (sampleSize > UPPER_LIMIT_PER_BUCKET) {
            sampleSize = UPPER_LIMIT_PER_BUCKET;
        }

        sampleSize /= retrArg.size();

        rg::Random32& random = retrArg[0].random;


        retrArg[0].globalData[0].resize(retrArg.size());
        retrArg[0].heap.resize(retrArg[0].k);



        for (int t = 0; t < retrArg.size(); t++) {

            std::vector<row_type> sampleIndx = rg::sample(random, sampleSize, retrArg[t].queryMatrix->rowNum);


            for (int i = 0; i < sampleSize; i++) {
                id = sampleIndx[i];



                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(id);

                retrArg[0].globalData[0][t][id] = GlobalTopkTuneData();

                for (row_type j = probeBuckets[0].startPos; j < probeBuckets[0].endPos; j++) {
                    double ip = retrArg[0].probeMatrix->innerProduct(j, query);
                    retrArg[0].globalData[0][t][id].results.push_back(QueueElement(ip, retrArg[0].probeMatrix->getId(j)));
                }

                // and now make the heap
                std::make_heap(retrArg[0].globalData[0][t][id].results.begin(), retrArg[0].globalData[0][t][id].results.end(), std::greater<QueueElement>()); // non thread safe

                for (row_type b = 1; b < probeBuckets.size(); b++) {

                    std::vector<QueueElement>& prevResults = retrArg[0].globalData[b - 1][t][id].results;

                    if (prevResults.front().data >= probeBuckets[b].normL2.second) { // bucket check
                        break;

                    } else {

                        if (retrArg[0].globalData[b].size() == 0) {
                            retrArg[0].globalData[b].resize(retrArg.size());
                        }

                        retrArg[0].globalData[b][t][id] = GlobalTopkTuneData();

                        retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                        std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                        retrArg[0].tunerTimer.start();
                        plainRetriever.runTopK(query, probeBuckets[b], &retrArg[0]);
                        retrArg[0].tunerTimer.stop();
                        retrArg[0].globalData[b][t][id].lengthTime = retrArg[0].tunerTimer.elapsedTime().nanos();

                        retrArg[0].globalData[b][t][id].results.assign(retrArg[0].heap.begin(), retrArg[0].heap.end());

                        if (b > bucketsForInit)
                            bucketsForInit = b;
                    }
                }
            }
        }



        for (row_type b = 1; b < probeBuckets.size(); b++) {
            row_type counter = 0;
            for (int t = 0; t < retrArg.size(); t++) {
                counter += retrArg[0].globalData[b][t].size();
            }
            if (counter >= LOWER_LIMIT_PER_BUCKET) {
                activeBuckets++;
            } else {
                break;
            }

        }
        activeBuckets++;
        bucketsForInit++;

        //        std::cout<<"activeBuckets: "<<activeBuckets<<" bucketsForInit: "<<bucketsForInit<<std::endl;

        p.first = activeBuckets;
        p.second = bucketsForInit;

        //        std::cout << "size 0: " << retrArg[0].globalData[0][0].size() << " size activeBuckets: " << retrArg[0].globalData[activeBuckets - 1][0].size() << std::endl;


        return p;
    }

   
    
}



#endif /* TUNETOPK_H_ */
