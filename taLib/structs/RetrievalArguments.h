/* 
 * File:   RetrievalArguments2.h
 * Author: chteflio
 *
 * Created on June 26, 2015, 3:02 PM
 */

#ifndef RETRIEVALARGUMENTS_H
#define	RETRIEVALARGUMENTS_H
#include <taLib/structs/QueryBatch.h>
#include <boost/dynamic_bitset.hpp>

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


        boost::dynamic_bitset<> done; // for LSH
        std::vector<float> sums; // for LSH
        std::vector<row_type> countsOfBlockValues; // for LSH

        std::vector<row_type> candidatesToVerify;
        std::vector<row_type> cp_array;// for coord
        std::vector<Candidate_incr> ext_cp_array;// for icoord
        std::vector<QueryBucket_withTuning> queryBatches;

        std::vector<double> accum, hashlen, hashval, hashwgt; // for L2AP
        
        // for tuning
        std::vector<double>* competitorMethod;
        std::vector< std::vector< unordered_map< row_type, GlobalTopkTuneData > > > globalData; // 1: Bucket 2: thread 3: valid sample points for bucket -->result
        std::vector<row_type> sampleInd;
        std::vector<QueueElement> thetaForActiveBlocks; // for LSH
        
        
        uint8_t* sketches; //for LSH no need to keep the actual sketches for the probe vectors. just keep the buckets with the ids
        

        VectorMatrix* probeMatrix;
        VectorMatrix* queryMatrix;
        TAState* state; //for TA
        TANRAState* tanraState; //for TANRA
        RandomIntGaussians* rig;


        TreeIndex* tree; //for Tree
        const col_type* listsQueue;// for coord or icoord

        rg::Timer t, tunerTimer;
        rg::Random32 random;

        int threads, k;
        int prevBucketBestPhi; // to be used during tuning

        // for INCR or COORD
        double theta, t_b, R, epsilon, currEpsilonAppr;
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
        boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0), initializeListsTime(0), lengthTime(0), state(0),tanraState(0),
        prevBucketBestPhi(-1), threads(1), worstMinScore(std::numeric_limits<double>::max()), 
        competitorMethod(0), onlyRunTime(0), sketches(0), rig(NULL) {
            random = rg::Random32(123); // PSEUDO-RANDOM

        }

        inline RetrievalArguments() : comparisons(0), forCosine(true), boundsTime(0), ipTime(0), scanTime(0),
        preprocessTime(0), filterTime(0), initializeListsTime(0), lengthTime(0), state(0),tanraState(0),
        prevBucketBestPhi(-1), threads(1), colnum(0), queryMatrix(0), probeMatrix(0), competitorMethod(0), onlyRunTime(0), sketches(0), rig(NULL) {
            random = rg::Random32(123);
        }

        inline void initializeBasics(
                VectorMatrix& queryMatrix1, VectorMatrix& probeMatrix1,
                LEMP_Method method1, double theta1, int k1,
                int threads1, double R1, double epsilon1,  bool forCosine1 = true) {

            queryMatrix = &queryMatrix1;
            probeMatrix = &probeMatrix1;

            colnum = queryMatrix->colNum;
            method = method1;
            forCosine = forCosine1;

            theta = theta1;
            k = k1;
            threads = threads1;
            R = R1;
            epsilon = epsilon1;     
	    
#ifdef RELATIVE_APPROX
	    epsilon = epsilon/(1-epsilon);
#endif

        }

        ~RetrievalArguments() {
            if (state != 0) {
                delete state;
            }
            
            if (tanraState != 0) {
                delete tanraState;
            }

            if (method == LEMP_LSH) {
                if (sketches)
                    delete[] sketches;

                if (rig)
                    delete rig;

            }
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


        inline void initThetaForActiveBlocks() {
            thetaForActiveBlocks.resize(LSH_SIGNATURES);
            double logNom = log(1-R);
            double exponent = 1/((double) LSH_CODE_LENGTH);//0.125; // 1/8 LSH_CODE_LENGTH

            for (int b = 0; b < LSH_SIGNATURES; b++) {
                double logDenom = logNom / (b + 1);
                double denom = exp(logDenom);
                double thres8 = 1 - denom;
                double thres = pow(thres8, exponent);
                double acosThres = (1 - thres) * PI;
                double theta = cos(acosThres);

                thetaForActiveBlocks[b].data = theta;
                thetaForActiveBlocks[b].id = b + 1;
            }
            std::sort(thetaForActiveBlocks.begin(), thetaForActiveBlocks.end(), std::less<QueueElement>());

        }

        inline int findActiveBlocks(double theta) {
            std::vector<QueueElement>::iterator it = std::upper_bound(thetaForActiveBlocks.begin(), thetaForActiveBlocks.end(), QueueElement(theta, 0));

            int pos = it - thetaForActiveBlocks.begin();
            if (pos >= thetaForActiveBlocks.size()) {
                pos = thetaForActiveBlocks.size() - 1;
                return thetaForActiveBlocks[pos].id;
            }

            int b = thetaForActiveBlocks[pos].id + 1;

            return b;

        }

        inline void init(row_type maxProbeBucketSize) {

            if (method == LEMP_LI || method == LEMP_I) {
                ext_cp_array.reserve(maxProbeBucketSize + 1);
                ext_cp_array.resize(maxProbeBucketSize);
            }


            //maxProbeBucketSize = maxBucketSize;
            if (method == LEMP_LI || method == LEMP_I ||
                    method == LEMP_LC || method == LEMP_C ||
                    method == LEMP_AP ||
                    method == LEMP_LSH ) {
                candidatesToVerify.reserve(maxProbeBucketSize + 2);
                candidatesToVerify.resize(maxProbeBucketSize);
            }

            if (method == LEMP_AP) {
                accum.resize(maxProbeBucketSize, -1);
                hashlen.resize(colnum);
                hashval.resize(colnum, 0);
                hashwgt.resize(colnum);
            }

            if (method == LEMP_LC || method == LEMP_C) {
                cp_array.resize(maxProbeBucketSize);
            }

            if (method == LEMP_TA) {
                state = new TAState(colnum);
            }
            
            if (method == LEMP_TANRA) {
                tanraState = new TANRAState(colnum);
                candidatesToVerify.resize(maxProbeBucketSize);
            }
            
            

            if (method == LEMP_LSH) {
                //                setupSketches(LSH_SIGNATURES, 8);
                rig = new RandomIntGaussians(colnum, LSH_SIGNATURES * LSH_CODE_LENGTH);
                done.resize(maxProbeBucketSize);
                
		
		if(forCosine){
			initThetaForActiveBlocks();
			long long totalSketchSize = ((long)maxProbeBucketSize) * (LSH_CODE_LENGTH/8) * 	((long)LSH_SIGNATURES);
			sketches = new uint8_t[totalSketchSize]();
		}else{ // forCosine=false here means that I am running the simple-LSH therefore I will try to reduce the cost of LSH signature creation
		  	long long totalSketchSize = (LSH_CODE_LENGTH/8) * ((long)LSH_SIGNATURES);
			sketches = new uint8_t[totalSketchSize]();	
		}

                sums.resize(LSH_SIGNATURES * LSH_CODE_LENGTH, 0);		
                if (LSH_CODE_LENGTH == 8)
                    countsOfBlockValues.resize(256, 0);

            }

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
#endif	/* RETRIEVALARGUMENTS2_H */

