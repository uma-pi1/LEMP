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


        std::pair<row_type, row_type> p;


        LengthRetriever plainRetriever;

        // calculate how large the sample should be
        row_type sampleSize = 0;
        for (auto& arg : retrArg) {
            sampleSize += arg.queryMatrix->rowNum;
        }
        sampleSize *= 0.02;
        if (sampleSize > UPPER_LIMIT_PER_BUCKET) {
            sampleSize = UPPER_LIMIT_PER_BUCKET;
        }
        sampleSize /= retrArg.size();

        rg::Random32& random = retrArg[0].random;


        for (auto& b : probeBuckets)
            b.sampleThetas.resize(retrArg.size());

        retrArg[0].heap.resize(retrArg[0].k);


        // from the query matrix that corresponds to each thread...
        for (int t = 0; t < retrArg.size(); ++t) {

            //... get a sample
            std::vector<row_type> sampleIndx = rg::sample(random, sampleSize, retrArg[t].queryMatrix->rowNum);


            for (auto id : sampleIndx) {

                const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(id);

                // I will need to keep track of the topk results of each query for each bucket. The kth value in the topk list will give me the theta_b(q))
                // that corresponds to this query
                // I keep this info in the retrArg because it will be needed later in the ListsTuneData.h                
                probeBuckets[0].sampleThetas[t][id] = GlobalTopkTuneData();


                for (row_type j = probeBuckets[0].startPos; j < probeBuckets[0].endPos; ++j) {
                    double ip = retrArg[0].probeMatrix->innerProduct(j, query);
                    probeBuckets[0].sampleThetas[t][id].results.emplace_back(ip, retrArg[0].probeMatrix->getId(j));
                }

                // and now make the heap
                std::make_heap(probeBuckets[0].sampleThetas[t][id].results.begin(), probeBuckets[0].sampleThetas[t][id].results.end(), std::greater<QueueElement>()); // non thread safe

                for (row_type b = 1; b < probeBuckets.size(); ++b) {


                    std::vector<QueueElement>& prevResults = probeBuckets[b - 1].sampleThetas[t][id].results;


                    if (prevResults.front().data >= probeBuckets[b].normL2.second) { // bucket check
                        break;

                    } else {// run LENGTH and measure the time

                        if (probeBuckets[b].sampleThetas.size() == 0) {
                            probeBuckets[b].sampleThetas.resize(retrArg.size());
                        }

                        probeBuckets[b].sampleThetas[t][id] = GlobalTopkTuneData();

                        std::copy(prevResults.begin(), prevResults.end(), retrArg[0].heap.begin());
                        std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                        retrArg[0].tunerTimer.start();
                        plainRetriever.runTopK(query, probeBuckets[b], &retrArg[0]);
                        retrArg[0].tunerTimer.stop();
                        probeBuckets[b].sampleThetas[t][id].lengthTime = retrArg[0].tunerTimer.elapsedTime().nanos();

                        probeBuckets[b].sampleThetas[t][id].results.reserve(retrArg[0].k);
                        std::copy_n(retrArg[0].heap.begin(), retrArg[0].k, std::back_inserter(probeBuckets[b].sampleThetas[t][id].results));

                        if (b > bucketsForInit)
                            bucketsForInit = b;
                    }
                }
            }
        }



        for (row_type b = 1; b < probeBuckets.size(); ++b) {
            row_type counter = 0;

            for (int t = 0; t < retrArg.size(); ++t) {
                counter += probeBuckets[b].sampleThetas[t].size();
            }

            if (counter >= LOWER_LIMIT_PER_BUCKET) {
                activeBuckets++;
            } else { // if I do not have enough sample queries for a bucket, it makes no sense to try  to tune
                break;
            }

        }
        activeBuckets++;        
        bucketsForInit++;
        
        
        for (row_type b = activeBuckets; b < probeBuckets.size(); ++b) {
            if (probeBuckets[b].xValues != nullptr)
                probeBuckets[b].xValues->clear();
            probeBuckets[b].sampleThetas.clear();
        }


        p.first = activeBuckets; // I will try to tune t_b, phi for all these buckets
        p.second = bucketsForInit; // I will create indexes for all these buckets (in multi-threaded it is better to create indexes before you enter the retrieval phase)


        return p;
    }



}



#endif /* TUNETOPK_H_ */
