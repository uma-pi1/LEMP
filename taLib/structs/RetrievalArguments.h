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

#include "Candidates.h"

namespace ta {

    typedef boost::shared_ptr< std::vector<MatItem> > xValues_ptr; // data: localTheta i: thread j: posInMatrix

    struct RetrievalArguments {
        std::vector<IntervalElement> intervals;
        std::vector<MatItem > results;
        std::vector<QueueElement> topkResults;
        std::vector<QueueElement> heap;


        boost::dynamic_bitset<> done; // for LSH
        std::vector<float> sums; // for LSH
        std::vector<row_type> countsOfBlockValues; // for LSH

        row_type* candidatesToVerify;
        row_type* cp_array; // for coord
        Candidate_incr* ext_cp_array; // for icoord

        std::vector<QueryBatch> queryBatches;

        std::vector<double> accum, hashval; // for L2AP
        double* hashlen; // for L2AP
        double* hashwgt; // for L2AP

        // for tuning
        std::vector<double>* competitorMethod;


        uint8_t* sketches; //for LSH no need to keep the actual sketches for the probe vectors. just keep the buckets with the ids

        VectorMatrix* probeMatrix;
        VectorMatrix* queryMatrix;
        TAState* state; //for TA
        TANRAState* tanraState; //for TANRA

        TreeIndex* tree; //for Tree
        const col_type* listsQueue; // for ICOORD or COORD

        rg::Timer t, tunerTimer;
        rg::Random32 random;

        int threads, k;

        double theta, t_b, R, epsilon, currEpsilonAppr; // for ICOORD or COORD
        double worstMinScore; // for L2AP
        double boundsTime, ipTime, scanTime, preprocessTime, filterTime, initializeListsTime, lengthTime;

        mutable comp_type comparisons; // pass  RetrievalArguments as const argument in a function, but still get the comparisons modified 

        LEMP_Method method;

        row_type queryId;
        row_type queryPos; //for Tree

        col_type colnum, maxLists, numLists;
        bool forCosine; // to be used for TA
        bool isTARR; // true: TAStateRR  false: TAStateMAX


        // with constructor delegation

        inline RetrievalArguments() : RetrievalArguments(0, nullptr, nullptr, LEMP_L) {
        }

        RetrievalArguments(col_type colnum, VectorMatrix* queryMatrix, VectorMatrix* probeMatrix, LEMP_Method method,
                bool forCosine = true, bool isTARR = false) :
        colnum(colnum), comparisons(0), probeMatrix(probeMatrix), queryMatrix(queryMatrix), forCosine(forCosine), method(method),
        boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0), initializeListsTime(0), lengthTime(0), state(nullptr), tanraState(nullptr),
        threads(1), worstMinScore(std::numeric_limits<double>::max()), hashwgt(nullptr), hashlen(nullptr),
        competitorMethod(nullptr), sketches(nullptr), isTARR(isTARR), cp_array(nullptr), ext_cp_array(nullptr), candidatesToVerify(nullptr) {
            random = rg::Random32(123); // PSEUDO-RANDOM
        }

        inline void initializeBasics(
                VectorMatrix& queryMatrix1, VectorMatrix& probeMatrix1,
                LEMP_Method method1, double theta1, int k1,
                int threads1, double R1, double epsilon1, bool forCosine1 = true, bool isTARR1 = false) {

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
            isTARR = isTARR1;

#ifdef RELATIVE_APPROX
            epsilon = epsilon / (1 - epsilon);
#endif

        }

        ~RetrievalArguments() {
            if (state != nullptr)
                delete state;
            if (tanraState != nullptr)
                delete tanraState;
            if (cp_array != nullptr)
                delete[] cp_array;
            if (ext_cp_array != nullptr)
                delete[] cp_array;
            if (candidatesToVerify != nullptr)
                delete[] candidatesToVerify;
            if (hashlen != nullptr)
                delete[] hashlen;
            if (hashwgt != nullptr)
                delete[] hashwgt;
            if (sketches != nullptr)
                delete[] sketches;
        }

        inline void clear() {
            lengthTime = 0;
            boundsTime = 0;
            scanTime = 0;
            filterTime = 0;
            preprocessTime = 0;
            ipTime = 0;
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

        inline void init(row_type maxProbeBucketSize) {

            if (method == LEMP_LI || method == LEMP_I) {
                ext_cp_array = new Candidate_incr[maxProbeBucketSize];
            }


            //maxProbeBucketSize = maxBucketSize;
            if (method == LEMP_LI || method == LEMP_I ||
                    method == LEMP_LC || method == LEMP_C ||
                    method == LEMP_AP ||
                    method == LEMP_LSH) {
                candidatesToVerify = new row_type[maxProbeBucketSize];
            }

            if (method == LEMP_AP) {
                accum.resize(maxProbeBucketSize, -1);
                hashval.resize(colnum, 0);
                hashlen = new double[colnum];
                hashwgt = new double[colnum];
            }

            if (method == LEMP_LC || method == LEMP_C) {
                cp_array = new row_type[maxProbeBucketSize];
            }

            if (method == LEMP_TA) {
                if (isTARR) {
                    state = new TAStateRR(colnum);
                } else {
                    state = new TAStateMAX(colnum);
                }                
            }

            if (method == LEMP_TANRA) {
                tanraState = new TANRAState(colnum);

                if (isTARR) {
                    tanraState = new TANRAStateRR(colnum);
                } else {
                    tanraState = new TANRAStateMax(colnum);
                }

                candidatesToVerify = new row_type[maxProbeBucketSize];
            }

            if (method == LEMP_LSH) {
                done.resize(maxProbeBucketSize);
                long long totalSketchSize = ((long) maxProbeBucketSize) * (LSH_CODE_LENGTH / 8) * ((long) LSH_SIGNATURES);
                sketches = new uint8_t[totalSketchSize]();
                sums.resize(LSH_SIGNATURES * LSH_CODE_LENGTH, 0);
                if (LSH_CODE_LENGTH == 8)
                    countsOfBlockValues.resize(256, 0);

            }

        }

        inline void allocTopkResults() {
            heap.resize(k);
            topkResults.resize(queryMatrix->rowNum * k);
        }

        inline void setIntervals(col_type lists) {
            maxLists = lists;
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

            std::cout << "-------------" << std::endl;
#endif
        }

    };



}
#endif	/* RETRIEVALARGUMENTS2_H */

