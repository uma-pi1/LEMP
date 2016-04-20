
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
 * algo_with_tuning2.h
 *
 *  Created on: Mar 25, 2014
 *      Author: chteflio
 */

#ifndef LEMP_H_
#define LEMP_H_


namespace mips {


    // I can change between runs: theta, method, queryMatrix
    // I cannot change between runs: probeMatrix, k

    class Lemp : public Mip {
        std::vector<VectorMatrix> queryMatrices;
        std::vector<double> cweights; //for L2AP

        std::vector<ProbeBucket> probeBuckets;
        std::vector<RetrievalArguments> retrArg;

        row_type maxProbeBucketSize;
        LempArguments args;
        row_type activeBuckets;

        inline row_type initProbeBuckets(VectorMatrix& rightMatrix);
        inline void initializeRetrievers();
        inline void initQueryBatches(VectorMatrix& leftMatrix, row_type maxBlockSize, std::vector<RetrievalArguments>& retrArg);
        inline void initListsInBuckets();
        inline void tune(std::vector<RetrievalArguments>& retrArg, row_type allQueries);
        inline void printAlgoName(const VectorMatrix& queryMatrix);

    public:

        inline void setTheta(double theta) {
            args.theta = theta;
        }

        inline void setMethod(LEMP_Method method) {
            args.method = method;
        }

        inline Lemp(InputArguments& in, int cacheSizeinKB, LEMP_Method method, bool isTARR, double R, double epsilon) :
        maxProbeBucketSize(0) {
            args.copyInputArguments(in);
            args.cacheSizeinKB = cacheSizeinKB;
            args.method = method;
            args.isTARR = isTARR;
            args.R = R;
            args.epsilon = epsilon;

            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (!logging.is_open()) {
                std::cout << "[WARNING] No log will be created!" << std::endl;
            } else {
                std::cout << "[INFO] Logging in " << args.logFile << std::endl;
            }

            omp_set_num_threads(args.threads);
            retrArg.resize(args.threads);
        }

        inline ~Lemp() {
            logging.close();
        }

        inline void initialize(VectorMatrix& rightMatrix) {
            std::cout << "[INIT] ProbeMatrix contains " << rightMatrix.rowNum << " vectors with dimensionality " << (0 + rightMatrix.colNum) << std::endl;
            timer.start();
            maxProbeBucketSize = initProbeBuckets(rightMatrix);
            timer.stop();
            dataPreprocessingTimeRight += timer.elapsedTime().nanos();
        }

        inline void runAboveTheta(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            timer.start();
            //            std::vector<RetrievalArguments> retrArg; //one argument for each thread
            initQueryBatches(leftMatrix, maxProbeBucketSize, retrArg);
            initializeRetrievers();
            results.resultsVector.resize(args.threads);

            for (auto& argument : retrArg)
                argument.init(maxProbeBucketSize);

            timer.stop();
            dataPreprocessingTimeLeft += timer.elapsedTime().nanos();

            timer.start();
            initListsInBuckets();
            timer.stop();
            dataPreprocessingTimeRight += timer.elapsedTime().nanos();


            tune(retrArg, leftMatrix.rowNum);

            std::cout << "[RETRIEVAL] Retrieval (theta = " << args.theta << ") starts ..." << std::endl;
            logging << "theta(" << args.theta << ")\t";

            timer.start();
            col_type maxLists = 1;


            switch (args.method) {
                case LEMP_I:
                case LEMP_LI:
                case LEMP_C:
                case LEMP_LC:

                    std::for_each(probeBuckets.begin(), probeBuckets.begin() + activeBuckets, [&maxLists](const ProbeBucket & b) {
                        if (maxLists < b.numLists) maxLists = b.numLists;
                    });

                    for (auto& argument : retrArg)
                        argument.setIntervals(maxLists);


                    break;
            }

            for (auto& argument : retrArg)
                argument.clear();


            comp_type comparisons = 0;


#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();

                for (row_type b = 0; b < activeBuckets; ++b) {
                    probeBuckets[b].ptrRetriever->run(probeBuckets[b], &retrArg[tid]);
                }
                comparisons += retrArg[tid].comparisons;
                results.moveAppend(retrArg[tid].results, tid);
            }


            int totalSize = results.getResultSize();

            timer.stop();
            retrievalTime += timer.elapsedTime().nanos();
            totalComparisons += comparisons;

            std::cout << "[RETRIEVAL] ... and is finished with " << totalSize << " results" << std::endl;
            logging << totalSize << "\t";

            outputStats();
        }

        inline void runTopK(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            // initialize Retrievers
            timer.start();
            initQueryBatches(leftMatrix, maxProbeBucketSize, retrArg);
            initializeRetrievers();
            for (auto& argument : retrArg) {
                argument.init(maxProbeBucketSize);
                if (args.k > 0) {
                    argument.allocTopkResults();
                }
            }

            results.resultsVector.resize(args.threads);
            timer.stop();
            dataPreprocessingTimeLeft += timer.elapsedTime().nanos();

            tune(retrArg, leftMatrix.rowNum);

            std::cout << "[RETRIEVAL] Retrieval (k = " << args.k << ") starts ..." << std::endl;
            logging << "k(" << args.k << ")\t";

            timer.start();
            col_type maxLists = 1;

            switch (args.method) {
                case LEMP_I:
                case LEMP_LI:
                case LEMP_C:
                case LEMP_LC:

                    std::for_each(probeBuckets.begin(), probeBuckets.begin() + activeBuckets, [&maxLists](const ProbeBucket & b) {
                        if (maxLists < b.numLists) maxLists = b.numLists;
                    });

                    for (auto& argument : retrArg)
                        argument.setIntervals(maxLists);

                    break;
            }
            for (auto& argument : retrArg)
                argument.clear();


            comp_type comparisons = 0;
            double worstMinScore = std::numeric_limits<double>::max();
            double totalError = 0;

#pragma omp parallel reduction(+ : comparisons, totalError)
            {
                row_type tid = omp_get_thread_num();


                for (row_type b = 0; b < probeBuckets.size(); ++b) {//

                    probeBuckets[b].ptrRetriever->runTopK(probeBuckets[b], &retrArg[tid]);

                    if (args.method == LEMP_AP || args.method == LEMP_BLSH) { // synchronize the worstMinScore      

                        // obviously you should avoid using these guys in a parallel setting
#pragma omp critical
                        {
                            if (worstMinScore > retrArg[tid].worstMinScore) {
                                worstMinScore = retrArg[tid].worstMinScore;
                            }
                        }
#pragma omp barrier
                        retrArg[tid].worstMinScore = worstMinScore;
#pragma omp barrier
                        worstMinScore = std::numeric_limits<double>::max();

                    }
                }
                retrArg[tid].extendIncompleteResultItems();
                results.moveAppend(retrArg[tid].results, tid);
                comparisons += retrArg[tid].comparisons;
                totalError += retrArg[tid].totalErrorAfterResults;
            }


            timer.stop();
            retrievalTime += timer.elapsedTime().nanos();
            totalComparisons += comparisons;
            std::cout << "[RETRIEVAL] ... and is finished with " << results.getResultSize() << " results" << std::endl;

//             std::cout << "TOTAL ERROR: " << totalError / leftMatrix.rowNum << " countABOVE: " << countABOVE << std::endl;
            logging << totalError / leftMatrix.rowNum << "\t" << results.getResultSize() << "\t";

            outputStats();
        }


    };

    inline row_type Lemp::initProbeBuckets(VectorMatrix& rightMatrix) {
        std::vector<row_type> probeBucketOffsets;

        probeMatrix.init(rightMatrix, true, false); // normalize and sort

        row_type maxBlockSize = computeBlockOffsetsByFactorCacheFittingForItems(probeMatrix.lengthInfo,
                probeMatrix.rowNum, probeBucketOffsets, FACTOR, ITEMS_PER_BLOCK, args.cacheSizeinKB, probeMatrix.colNum, args);

        bucketize(probeBuckets, probeMatrix, probeBucketOffsets, args);

        std::cout << "[INIT] ProbeBuckets = " << probeBucketOffsets.size() << std::endl;
        return maxBlockSize;
    }

    inline void Lemp::initializeRetrievers() {

        activeBuckets = probeBuckets.size();
        row_type b0 = 0;
        maxProbeBucketSize = 0;

        if (args.k == 0) { // Above-theta
            double maxUserLength = 0;


            std::for_each(queryMatrices.begin(), queryMatrices.end(),
                    [&maxUserLength](const VectorMatrix & m) {
                        if (maxUserLength < m.lengthInfo[0].data) maxUserLength = m.lengthInfo[0].data;
                    });



            for (row_type i = 0; i < probeBuckets.size(); ++i) {
                probeBuckets[i].bucketScanThreshold = args.theta * probeBuckets[i].invNormL2.second;

                // find maxProbeBucketLength
                if (maxUserLength > probeBuckets[i].bucketScanThreshold) {
                    if (maxProbeBucketSize < probeBuckets[i].endPos - probeBuckets[i].startPos)
                        maxProbeBucketSize = probeBuckets[i].endPos - probeBuckets[i].startPos;
                } else {
                    activeBuckets = i;
                    break;
                }
            }
        } else { // Row-top-k

            std::for_each(probeBuckets.begin() + 1, probeBuckets.end(), [this](const ProbeBucket & b) {
                if (maxProbeBucketSize < b.rowNum) maxProbeBucketSize = b.rowNum;
            });

            b0 = 1;
            retriever_ptr firstPtr(new Retriever());
            probeBuckets[0].ptrRetriever = firstPtr;

        }

        switch (args.method) {
            case LEMP_LI:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new LX_Retriever<IncrRetriever>());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[SL] == 0) {
                        probeBuckets[b].ptrIndexes[SL] = new QueueElementLists();
                    }
                }
                break;

            case LEMP_I:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new IncrRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[SL] == 0)
                        probeBuckets[b].ptrIndexes[SL] = new QueueElementLists();
                }
                break;

            case LEMP_LC:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new LX_Retriever<CoordRetriever>());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[INT_SL] == 0)
                        probeBuckets[b].ptrIndexes[INT_SL] = new IntLists();
                }
                break;

            case LEMP_C:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new CoordRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[INT_SL] == 0)
                        probeBuckets[b].ptrIndexes[INT_SL] = new IntLists();
                }
                break;

            case LEMP_TA:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {

                    retriever_ptr rPtr(new taRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[SL] == 0)
                        probeBuckets[b].ptrIndexes[SL] = new QueueElementLists();
                }
                break;

            case LEMP_L:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new LengthRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                }
                break;

            case LEMP_TREE:
#pragma omp parallel for schedule(static,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new SingleTree());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[TREE] == 0)
                        probeBuckets[b].ptrIndexes[TREE] = new TreeIndex();
                }
                break;

            case LEMP_AP:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new apRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[AP] == 0)
                        probeBuckets[b].ptrIndexes[AP] = new L2apIndex();
                }
                break;

            case LEMP_LSH:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new LshRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[LSH] == 0)
                        probeBuckets[b].ptrIndexes[LSH] = new LshIndex();
                }
                break;

            case LEMP_BLSH:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < activeBuckets; ++b) {
                    retriever_ptr rPtr(new BlshRetriever());
                    probeBuckets[b].ptrRetriever = rPtr;
                    if (probeBuckets[b].ptrIndexes[BLSH] == 0)
                        probeBuckets[b].ptrIndexes[BLSH] = new BlshIndex();
                }
                break;
        }
    }

    inline void Lemp::initQueryBatches(VectorMatrix& leftMatrix, row_type maxBlockSize, std::vector<RetrievalArguments>& retrArg) {

        std::cout << "[RETRIEVAL] QueryMatrix contains " << leftMatrix.rowNum << " vectors with dimensionality " << (0 + leftMatrix.colNum) << std::endl;


        row_type nCount = 0;
        row_type myNumThreads = args.threads;

        if (leftMatrix.rowNum < args.threads) {
            myNumThreads = leftMatrix.rowNum;
            std::cout << "[WARNING] Query matrix contains too few elements. Suboptimal running with " << myNumThreads << " thread(s)" << std::endl;
        }
        omp_set_num_threads(myNumThreads);
        queryMatrices.resize(myNumThreads);
        retrArg.resize(myNumThreads);

        if (args.k > 0) { // this is a top-k version
            initializeMatrices(leftMatrix, queryMatrices, false, true, args.epsilon); // normalize but don't sort
        } else {
            initializeMatrices(leftMatrix, queryMatrices, true, false); // normalize and sort
        }

#pragma omp parallel reduction(+ : nCount)
        {

            row_type tid = omp_get_thread_num();
            std::vector<row_type> blockOffsets;
            computeBlockOffsetsForUsersFixed(queryMatrices[tid].rowNum, blockOffsets, args.cacheSizeinKB, queryMatrices[tid].colNum, args, maxBlockSize);
            bucketize(retrArg[tid].queryBatches, queryMatrices[tid], blockOffsets, args);
            nCount += retrArg[tid].queryBatches.size();
            retrArg[tid].initializeBasics(queryMatrices[tid], probeMatrix, args.method, args.theta, args.k, myNumThreads, args.R, args.epsilon, args.numTrees, args.search_k, true, args.isTARR);

        }

        std::cout << "[RETRIEVAL] QueryBatches = " << nCount << std::endl;

        switch (args.method) {
            case LEMP_AP:
                calculateAPneededForQuery(queryMatrices, 0, args.k, cweights);
                break;

            case LEMP_BLSH:

                if (((double) LSH_CODE_LENGTH * LSH_SIGNATURES / 32) < 1) {
                    std::cerr << "LSH_CODE_LENGTH * LSH_SIGNATURES should be at least 32. Change their value in  Definitions.h and recompile!" << std::endl;
                    exit(1);
                }
                // no break here
            case LEMP_LSH:

                if (LSH_CODE_LENGTH != 8 && LSH_CODE_LENGTH != 16 && LSH_CODE_LENGTH != 32 && LSH_CODE_LENGTH != 64) {
                    std::cerr << "LSH_CODE_LENGTH should be a 8, 16, 32 or 64. Change its value in  Definitions.h and recompile!" << std::endl;
                    exit(1);
                }

                repetitionsForTheta.init(args.R);
                break;


        }

    }

    inline void Lemp::initListsInBuckets() {

        double maxQueryLength = 0;
        // in the case of topk the tuning part has modified the activeBuckets
        row_type b0 = (args.k == 0 ? 0 : 1);

        std::cout << "[RETRIEVAL] ProbeBuckets (active) = " << activeBuckets << std::endl;

        double worstCaseTheta;

        switch (args.method) {
            case LEMP_LI:
            case LEMP_I:
            case LEMP_TA:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    static_cast<QueueElementLists*> (probeBuckets[b].ptrIndexes[SL])->initializeLists(probeMatrix, probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_LC:
            case LEMP_C:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    static_cast<IntLists*> (probeBuckets[b].ptrIndexes[INT_SL])->initializeLists(probeMatrix, probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_TREE:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    static_cast<TreeIndex*> (probeBuckets[b].ptrIndexes[TREE])->initializeTree(probeMatrix, args.threads, probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_AP:

                std::for_each(queryMatrices.begin(), queryMatrices.end(),
                        [&maxQueryLength](const VectorMatrix & m) {
                            if (maxQueryLength < m.getVectorLength(0)) maxQueryLength = m.getVectorLength(0);
                        });

                worstCaseTheta = args.theta / maxQueryLength;

#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    static_cast<L2apIndex*> (probeBuckets[b].ptrIndexes[AP])->initializeLists(probeMatrix, worstCaseTheta, cweights,
                            probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_LSH:

#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < activeBuckets; ++b) {
                    static_cast<LshIndex*> (probeBuckets[b].ptrIndexes[LSH])->initializeLists(probeMatrix, true, probeBuckets[b].startPos, probeBuckets[b].endPos);

                    //                    row_type signatures = LSH_SIGNATURES; //LSH_SIGNATURES/(b+1);
                    //                    if (retrArg.size() > 1 && signatures > 0) {
                    //                        row_type tid = omp_get_thread_num();
                    //
                    //                        static_cast<LshIndex*> (probeBuckets[b].ptrIndexes[LSH])->checkAndReallocateAll(&probeMatrix, true, probeBuckets[b].startPos, probeBuckets[b].endPos, signatures,
                    //                                retrArg[tid].sums, retrArg[tid].countsOfBlockValues, retrArg[tid].sketches);
                    //                    }
                }

                break;


            case LEMP_BLSH:
                if (args.k == 0) { // only for the 
                    std::for_each(queryMatrices.begin(), queryMatrices.end(),
                            [&maxQueryLength](const VectorMatrix & m) {
                                if (maxQueryLength < m.getVectorLength(0)) maxQueryLength = m.getVectorLength(0);
                            });

                    worstCaseTheta = args.theta / maxQueryLength;
                    // what happens here for tuning?


#pragma omp parallel for schedule(dynamic,1) 
                    for (row_type b = b0; b < activeBuckets; ++b) {
                        static_cast<BlshIndex*> (probeBuckets[b].ptrIndexes[BLSH])->initializeLists(probeMatrix, worstCaseTheta, true, args.R, probeBuckets[b].startPos, probeBuckets[b].endPos);

                    }

                }
                break;

        }
//         std::cout << "Done creating lists" << std::endl;
    }

    inline void Lemp::tune(std::vector<RetrievalArguments>& retrArg, row_type allQueries) {

        if (activeBuckets > 0) {
            switch (args.method) {
                case LEMP_LI:
                case LEMP_I:
                case LEMP_LC:
                case LEMP_C:
                case LEMP_LSH:
                case LEMP_BLSH:

                    if (probeBuckets[0].isTunable(allQueries)) {
                        if (args.k == 0) {

                            timer.start();
                            // first set-up the xValues in each retriever
                            for (row_type b = 0; b < activeBuckets; ++b) {
                                probeBuckets[b].sampling(retrArg);
                            }


                            // then do the actual tuning
                            for (row_type b = 0; b < activeBuckets; ++b) {
                                probeBuckets[b].ptrRetriever->tune(probeBuckets[b], (b == 0 ? probeBuckets[b] : probeBuckets[b - 1]), retrArg);                             
                            }

                            timer.stop();
                            tuningTime += timer.elapsedTime().nanos();

                        } else {
                            timer.start();
                            std::pair<row_type, row_type> p(probeBuckets.size(), probeBuckets.size());

                            p = findSampleForTuningTopk(probeBuckets, retrArg);


                            activeBuckets = p.first;
#pragma omp parallel for schedule(dynamic,1) 
                            for (row_type b = 1; b < activeBuckets; ++b) {
                                probeBuckets[b].setup_xValues_topk(retrArg, probeBuckets[b - 1].sampleThetas);
                            }

                            timer.stop();
                            tuningTime += timer.elapsedTime().nanos();

                            timer.start();
                            activeBuckets = p.second; // this change is mainly for multithreading. It will force us to initialize more lists up-front --> less locking later
                            initListsInBuckets();
                            activeBuckets = p.first;
                            timer.stop();
                            dataPreprocessingTimeRight += timer.elapsedTime().nanos();

                            timer.start();

                            for (row_type b = 1; b < probeBuckets.size(); ++b) {
                                probeBuckets[b].ptrRetriever->tuneTopk(probeBuckets[b], probeBuckets[b - 1], retrArg);
                            }

                            timer.stop();
                            tuningTime += timer.elapsedTime().nanos();
                        }
                    } else {
                        std::cout << "[WARNING] Too few queries (" << allQueries << ") for tuning (at least " << LOWER_LIMIT_PER_BUCKET * 3 << " needed)" << std::endl;
                        std::cout << "[WARNING] Using default (t_b=1, lists=1) or previous tuning values for all probe buckets " << std::endl;
                        std::cout << "[WARNING] You can either reduce LOWER_LIMIT_PER_BUCKET and recompile or use larger batches of queries" << std::endl;
                    }



                    break;
            }
        }


    }

    inline void Lemp::printAlgoName(const VectorMatrix& queryMatrix) {
        switch (args.method) {
            case LEMP_L:
                logging << "LEMP_L" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_L with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_LI:
                logging << "LEMP_LI" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_LI with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_I:
                logging << "LEMP_I" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_I with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_LC:
                logging << "LEMP_LC" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_LC with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_C:
                logging << "LEMP_C" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_C with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_TA:
                if (args.isTARR) {
                    logging << "LEMP_TA (RR)" << "\t" << args.threads << "\t";
                    std::cout << "[ALGORITHM] LEMP_TA (ROUND ROBIN) with " << args.threads << " thread(s)" << std::endl;
                } else {
                    logging << "LEMP_TA (MAX)" << "\t" << args.threads << "\t";
                    std::cout << "[ALGORITHM] LEMP_TA (MAX HEAP) with " << args.threads << " thread(s)" << std::endl;
                }
                break;
            case LEMP_TREE:
                logging << "LEMP_TREE" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_TREE with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_AP:
                logging << "LEMP_AP" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_AP with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_LSH:
                logging << "LEMP_LSH" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_LSH with " << args.threads << " thread(s)" << std::endl;
                break;
            case LEMP_BLSH:
                logging << "LEMP_BLSH" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] LEMP_BLSH with " << args.threads << " thread(s)" << std::endl;
                break;
       
        }

        logging << "P(" << probeMatrix.rowNum << "x" << (0 + probeMatrix.colNum) << ")\t";
        logging << "Q^T(" << queryMatrix.rowNum << "x" << (0 + queryMatrix.colNum) << ")\t";
    }


}


#endif /* ALGO_WITH_TUNING_H_ */
