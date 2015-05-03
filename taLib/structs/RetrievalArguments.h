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
 * RetrievalArguments.h
 *
 *  Created on: Jun 26, 2014
 *      Author: chteflio
 */

#ifndef RETRIEVALARGUMENTS_H_
#define RETRIEVALARGUMENTS_H_

#include <taLib/structs/QueryBatch.h>

namespace ta {

    typedef boost::shared_ptr< std::vector<MatItem> > xValues_ptr; // data: localTheta i: thread j: posInMatrix

    struct GlobalTopkTuneData {
        double lengthTime;
        std::vector<QueueElement> results;

        GlobalTopkTuneData() : lengthTime(0) {
        }
    };

    struct RetrievalArguments {
        std::vector<IntervalElement> intervals;
        std::vector<MatItem > results;
        std::vector<QueueElement> topkResults;
        std::vector<QueueElement> heap;


//         std::vector<bool> done; // for LSH
//         std::vector<float> sums; // for LSH
//         std::vector<row_type> countsOfBlockValues; // for LSH

        std::vector<row_type> candidatesToVerify;
        std::vector<row_type> cp_array;
        std::vector<Candidate_incr> ext_cp_array;
        std::vector<QueryBucket_withTuning> queryBatches;

//         std::vector<double> accum, hashlen, hashval, hashwgt; // for L2AP
        std::vector<double>* competitorMethod;
        std::vector< std::vector< unordered_map< row_type, GlobalTopkTuneData > > > globalData; // 1: Bucket 2: thread 3: valid sample points for bucket -->result
        std::vector<row_type> sampleInd;
//         std::vector<QueueElement> thetaForActiveBlocks; // for LSH
//         std::vector<uint8_t> sketches;//for LSH no need to keep the actual sketches for the probe vectors. just keep the buckets with the ids
//         std::vector<int>* clusterIds;

        VectorMatrix* probeMatrix;
        VectorMatrix* queryMatrix;
        TAState* state; //for TA
//         da_rig_t *rig; // BLSH
//         gsl_rng *r; // random number generator BLSH
//         TreeIndex* tree; //for Tree
        const col_type* listsQueue;

        rg::Timer t, tunerTimer;
        rg::Random32 random;
        
        


        int threads, k, numClusters;
        int prevList; // to be used during tuning

        // for INCR or COORD
        double theta, t_b, epsilon;
        double worstMinScore; // for L2AP
        double boundsTime, ipTime, scanTime, preprocessTime, filterTime, initializeListsTime, lengthTime, onlyRunTime;

        comp_type comparisons;

        LEMP_Method method;

        row_type queryId, bucketInd;
        row_type queryPos; //for Tree

        col_type colnum, maxLists, numLists;
        bool forCosine; // to be used for TA

        RetrievalArguments(col_type colnum, VectorMatrix& queryMatrix, VectorMatrix& probeMatrix, LEMP_Method method,
                bool forCosine = true) :
        colnum(colnum), comparisons(0), probeMatrix(&probeMatrix), queryMatrix(&queryMatrix), forCosine(forCosine), method(method),
        boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0), initializeListsTime(0), lengthTime(0), state(0),
        prevList(-1), threads(1), worstMinScore(std::numeric_limits<double>::max()), //startPosOfPrevBucket(0),
        competitorMethod(0), //rig(NULL), r(NULL), 
        onlyRunTime(0) {
            random = rg::Random32(123); // PSEUDO-RANDOM

        }

        inline RetrievalArguments() : comparisons(0), forCosine(true), boundsTime(0), ipTime(0), scanTime(0),
        preprocessTime(0), filterTime(0), initializeListsTime(0), lengthTime(0), state(0),
        prevList(-1), threads(1), colnum(0), queryMatrix(0), probeMatrix(0), competitorMethod(0), onlyRunTime(0) {
            random = rg::Random32(123);
        }

        inline void initializeBasics(
                VectorMatrix& queryMatrix1, VectorMatrix& probeMatrix1,
                LEMP_Method method1, double theta1, int k1,
                int threads1, double epsilon1, int numClusters1, bool forCosine1 = true) {

            queryMatrix = &queryMatrix1;
            probeMatrix = &probeMatrix1;

            colnum = queryMatrix->colNum;
            method = method1;
            forCosine = forCosine1;

            theta = theta1;
            k = k1;
            threads = threads1;
            epsilon = epsilon1;
            
            numClusters = numClusters1;
        }

        RetrievalArguments(RetrievalArguments* other) : boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0),
        initializeListsTime(0), lengthTime(0), state(0), prevList(-1), comparisons(0), probeMatrix(other->probeMatrix),
        queryMatrix(other->queryMatrix), forCosine(other->forCosine), theta(other->theta), k(other->k), colnum(other->colnum),
        maxLists(other->maxLists), method(other->method), threads(other->threads), // startPosOfPrevBucket(0),
        competitorMethod(0), worstMinScore(std::numeric_limits<double>::max()), onlyRunTime(0), epsilon(other->epsilon){
            random = rg::Random32(123); // PSEUDO-RANDOM

            intervals.resize(other->intervals.size());

            candidatesToVerify.resize(other->candidatesToVerify.size());
            cp_array.resize(other->cp_array.size());
            ext_cp_array.resize(other->ext_cp_array.size());

//             accum.resize(other->accum.size(), -1);
//             hashlen.resize(other->hashlen.size());
//             hashval.resize(other->hashval.size(), 0);
//             hashwgt.resize(other->hashwgt.size());
// 
//             done.resize(other->done.size());




            if (other->state != 0) {
                state = new TAState(colnum);
            }

//             if (other->rig != NULL) {
//                 setupSketches();
//             }
// 
//             if (method == LEMP_LSH) {
//                 //copy
//                 std::copy(other->thetaForActiveBlocks.begin() , other->thetaForActiveBlocks.end(), thetaForActiveBlocks.begin());
//                 sketches.resize(other->sketches.size(),0);
//             }

        }

        ~RetrievalArguments() {
            if (state != 0) {
                delete state;
            }

//             if (method == LEMP_BLSH || method == LEMP_LSH) {
//                 rig->r = NULL;
// //                 gsl_rng_free(r);
//                 gk_free((void**) &(rig->i2fCache), &(rig->intGaussians), &rig, LTERM);
// 
//             }


        }

        inline void clear() {
            lengthTime = 0;
            boundsTime = 0;
            scanTime = 0;
            filterTime = 0;
            preprocessTime = 0;
            ipTime = 0;
            onlyRunTime = 0;
            initializeListsTime = 0;
            comparisons = 0;

            results.clear();
        }

        inline void moveTopkToHeap(row_type pos) {
            std::copy(topkResults.begin() + pos, topkResults.begin() + pos + k, heap.begin());
        }

        inline void writeHeapToTopk(row_type queryPos) {
            row_type p = queryPos*k;
            std::copy(heap.begin(), heap.end(), topkResults.begin() + p);
        }

//         inline void setupSketches() {
// 
//             int nHashBits = LSH_SIGNATURES * 8;
//             sums.resize(nHashBits, 0);
//             countsOfBlockValues.resize(256, 0);
//             int i, j, range;
//             float *i2fCache;
// 
//             // set up rig and sketch parameters & temp memory
// 
//             rig = (da_rig_t*) gk_malloc(sizeof (da_rig_t), NULL); //"allocateBayesLSHMemory: rig"
//             r = gsl_rng_alloc(gsl_rng_taus2);
//             gsl_rng_set(r, 123); // I just put something
// 
//             range = rig->range = 16;
//             rig->leftLimit = (0 - range) / 2.0;
//             rig->rightLimit = rig->leftLimit + range;
//             rig->f2iFactor = F2I2F * 1.0 / range;
//             rig->i2fFactor = range * 1.0 / F2I2F;
//             rig->r = r;
//             rig->size = colnum * nHashBits;
//             rig->intGaussians = da_ui16malloc(rig->size, NULL); //"setupSketches: rig->intGaussians"
//             i2fCache = rig->i2fCache = da_fmalloc(F2I2F + 1, NULL); //"setupSketches: rig->i2fCache"
// 
//             // cache conversion values between int16 and float
//             for (i2fCache[0] = rig->leftLimit, i = 1; i <= F2I2F; i++)
//                 i2fCache[i] = i2fCache[i - 1] + rig->i2fFactor;
// 
//             // generate random vectors
//             fillRig(rig);
//         }

//         inline void initThetaForActiveBlocks() {
//             thetaForActiveBlocks.resize(LSH_SIGNATURES);
//             double logNom = log(epsilon);
//             double exponent = 0.125; // 1/8
// 
//             for (int b = 0; b < LSH_SIGNATURES; b++) {
//                 double logDenom = logNom / (b+1);
//                 double denom = exp(logDenom);
//                 double thres8 = 1 - denom;            
//                 double thres = pow(thres8, exponent);
//                 double acosThres = (1 - thres)* PI;
//                 double theta = cos(acosThres);
//                 
//                 thetaForActiveBlocks[b].data = theta;
//                 thetaForActiveBlocks[b].id = b+1;
//                 
// //                std::cout<<"b: "<<b+1<<"-->"<<theta<<std::endl;
// 
//             }
//             std::sort(thetaForActiveBlocks.begin(), thetaForActiveBlocks.end(), std::less<QueueElement>());
//             
//         }
//         inline int findActiveBlocks(double theta){
//             std::vector<QueueElement>::iterator it = std::upper_bound(thetaForActiveBlocks.begin(), thetaForActiveBlocks.end(), QueueElement(theta,0));
//             
//             int pos =   it - thetaForActiveBlocks.begin();   
//             if (pos >= thetaForActiveBlocks.size()){
//                 pos = thetaForActiveBlocks.size()-1;
//                 return thetaForActiveBlocks[pos].id;
//             }
//                 
//             int b =  thetaForActiveBlocks[pos].id +1; 
//             return b;
//         
//         }

        inline void init(row_type maxProbeBucketSize) {

            if (method == LEMP_LI || method == LEMP_I || method == LEMP_NB) {
                ext_cp_array.reserve(maxProbeBucketSize + 1);
                ext_cp_array.resize(maxProbeBucketSize);
            }


            //maxProbeBucketSize = maxBucketSize;
            if (method == LEMP_LI || method == LEMP_I ||
                    method == LEMP_LC || method == LEMP_C ||
                    method == LEMP_AP || method == LEMP_NB ||
                    method == LEMP_BLSH || method == LEMP_LSH) {
                candidatesToVerify.reserve(maxProbeBucketSize + 2);
                candidatesToVerify.resize(maxProbeBucketSize);
            }

//             if (method == LEMP_AP) {
//                 accum.resize(maxProbeBucketSize, -1);
//                 hashlen.resize(colnum);
//                 hashval.resize(colnum, 0);
//                 hashwgt.resize(colnum);
//             }



            if (method == LEMP_LC || method == LEMP_C) {
                cp_array.resize(maxProbeBucketSize);
            }

            if (method == LEMP_TA) {
                state = new TAState(colnum);
            }

//             if (method == LEMP_BLSH || method == LEMP_LSH) {
//                 setupSketches();
//             }

//             if (method == LEMP_LSH) {
//                 done.resize(maxProbeBucketSize);
//                 initThetaForActiveBlocks();
//                 sketches.resize(maxProbeBucketSize * LSH_SIGNATURES, 0);
//             }
        }

        inline void setIntervals(col_type lists) {
            maxLists = lists;
            intervals.reserve(maxLists + 1);
            intervals.resize(maxLists);
            comparisons = 0;
        }

        inline void setQueues(const col_type* queue) {
            listsQueue = queue;
        }

        void printTimes() {
#ifdef TIME_IT
            std::cout << "-------------" << std::endl;

            std::cout << "lengthTime: " << lengthTime / 1E9 << std::endl;
            std::cout << "boundsTime: " << boundsTime / 1E9 << std::endl;
            std::cout << "scanTime: " << scanTime / 1E9 << std::endl;
            std::cout << "filterTime: " << filterTime / 1E9 << std::endl;
            std::cout << "ipTime: " << ipTime / 1E9 << std::endl;
            std::cout << "preprocessTime: " << preprocessTime / 1E9 << std::endl;
            std::cout << "initializeListsTime: " << initializeListsTime / 1E9 << std::endl;
            std::cout << "OnlyRunTime: " << onlyRunTime / 1E9 << std::endl;

            std::cout << "-------------" << std::endl;
#endif
        }

    };



}

#endif /* RETRIEVALARGUMENTS_H_ */
