
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

#ifndef ALGO_WITH_TUNING_H_
#define ALGO_WITH_TUNING_H_


namespace ta {

    class Algo_withTuning {
        // 1dim: number of threads
        // i: leftId, j: rightId 1 dim:users 2 dim: correct matches for user
        std::vector< std::vector<MatItem >* > results; //threads - results

        rg::Timer t;

        std::vector<VectorMatrix> queryMatrices;
        std::vector<double> cweights; //for L2AP
        VectorMatrix probeMatrix;

        std::vector<ProbeBucket> probeBuckets;
        std::vector<std::unique_ptr<Retriever> > retrievers;
        std::vector<RetrievalArguments> retrArg; //one argument for each thread

        row_type maxProbeBucketSize;
        LEMPArg args;
        std::ofstream logging;
        row_type numBucketsToFireQueries; // these buckets should have non-null retrievers and indexes
        row_type numBucketsToTune;
        row_type numBucketsToInitIndexes;


        inline row_type initProbeBuckets(VectorMatrix& rightMatrix);
        inline void initializeRetrievers();
        inline void initQueryBatches(VectorMatrix& leftMatrix, row_type maxBlockSize);
        inline void initListsInBuckets();
        inline void tune();
        inline void printTimes(rg::Timer& tAll);

    public:

        inline Algo_withTuning(LEMPArg& args) : args(args), maxProbeBucketSize(0) {

            VectorMatrix leftMatrix, rightMatrix;

            if (args.querySideLeft) {
                leftMatrix.readFromFile(args.usersFile, args.r, args.m, true);
                rightMatrix.readFromFile(args.itemsFile, args.r, args.n, false);
            } else {
                leftMatrix.readFromFile(args.itemsFile, args.r, args.n, false);
                rightMatrix.readFromFile(args.usersFile, args.r, args.m, true);
            }


            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (logging.is_open()) {

                std::cout << "Writing output to " << args.logFile << std::endl;
            } else {

                std::cout << "Problem with opening log-file. No log-file will be created" << std::endl;
            }

            t.start();
            omp_set_num_threads(args.threads);

            results.resize(args.threads);

            queryMatrices.resize(args.threads);
            retrArg.resize(args.threads);
            row_type maxBlockSize = initProbeBuckets(rightMatrix);
            initQueryBatches(leftMatrix, maxBlockSize);
            t.stop();
            Algo_withTuning::args.dataManipulationTime += t.elapsedTime().nanos();


            std::cout << "Cache Size (KB) per Processor: " << args.cacheSizeinKB << std::endl;

            switch (args.method) {
                case LEMP_L:
                    logging << "LEMP_L";
                    std::cout << "ALGO: LEMP_L" << std::endl;
                    break;
                case LEMP_LI:
                    logging << "LEMP_LI";
                    std::cout << "ALGO: LEMP_LI" << std::endl;
                    break;
                case LEMP_I:
                    logging << "LEMP_I";
                    std::cout << "ALGO: LEMP_I" << std::endl;
                    break;
                case LEMP_LC:
                    logging << "LEMP_LC";
                    std::cout << "ALGO: LEMP_LC" << std::endl;
                    break;
                case LEMP_C:
                    logging << "LEMP_C";
                    std::cout << "ALGO: LEMP_C" << std::endl;
                    break;
                case LEMP_TA:
                    if (args.isTARR) {
                        logging << "LEMP_TA RR";
                        std::cout << "ALGO: LEMP_TA with RR " << std::endl;
                    } else {
                        logging << "LEMP_TA MaxPiQi";
                        std::cout << "ALGO: LEMP_TA with MaxPiQi " << std::endl;
                    }
                    break;
                case LEMP_TREE:
                    logging << "LEMP_TREE";
                    std::cout << "ALGO: LEMP_TREE" << std::endl;
                    break;
                case LEMP_AP:
                    logging << "LEMP_AP";
                    std::cout << "ALGO: LEMP_AP" << std::endl;
                    t.start();
                    calculateAPneededForQuery(queryMatrices, 0, args.k, cweights);
                    t.stop();
                    Algo_withTuning::args.dataManipulationTime += t.elapsedTime().nanos();
                    break;
                case LEMP_LSH:
                    logging << "LEMP_LSH";
                    std::cout << "ALGO: LEMP_LSH" << std::endl;
                    std::cout << "Max LSH signatures: " << LSH_SIGNATURES << std::endl;
                    std::cout << "Recall rate: " << args.R << std::endl;
                    t.start();
                    static_assert(LSH_CODE_LENGTH % 8 == 0, "Require that LSH_CODE_LENGTH is multiple of 8 bits");
                    repetitionsForTheta.init(args.R);
                    coreLshInfo.init(LSH_CODE_LENGTH, LSH_SIGNATURES);
                    t.stop();
                    Algo_withTuning::args.dataManipulationTime += t.elapsedTime().nanos();
                    break;
            }



            std::cout << "Threads used " << args.threads << std::endl;

            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";
        }

        inline ~Algo_withTuning() {
        }

        inline void multiply() {

            t.start();
            initializeRetrievers();

            for (auto& argument : retrArg)
                argument.init(maxProbeBucketSize);

            initListsInBuckets();

            t.stop();
            Algo_withTuning::args.dataManipulationTime += t.elapsedTime().nanos();

            tune();

            std::cout << "dataManipulationTime: " << Algo_withTuning::args.dataManipulationTime / 1E9 << std::endl;
            std::cout << "Multiplication starts! theta: " << args.theta << std::endl;


            t.start();
            col_type maxLists = 1;


            switch (args.method) {
                case LEMP_I:
                case LEMP_LI:
                case LEMP_C:
                case LEMP_LC:

                    std::for_each(probeBuckets.begin(), probeBuckets.begin() + numBucketsToTune, [&maxLists](const ProbeBucket & b) {
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

                for (row_type b = 0; b < numBucketsToFireQueries; ++b) {
                    retrievers[b]->run(probeBuckets[b], &retrArg[tid]);
                }
                comparisons += retrArg[tid].comparisons;

                results[tid] = &retrArg[tid].results;

            }


            int totalSize = getResultSetSize();

            t.stop();

            std::cout << "Time for Retrieval: " << t << std::endl;
            std::cout << "Size of result: " << totalSize << std::endl;
            std::cout << "Comparisons: " << comparisons << std::endl;

            logging << "\t" << args.theta << "\t" << comparisons << "\t" << totalSize << "\t";
            printTimes(t);

            if (args.resultsFile != "") {
                writeResults(results, args.resultsFile);
            }


            logging.close();
        }

        inline void runTopkPerUser() {

            // initialize Retrievers
            t.start();
            initializeRetrievers();
            for (auto& argument : retrArg) {
                argument.init(maxProbeBucketSize);
                if (args.k > 0) {
                    argument.allocTopkResults();
                }
            }

            t.stop();
            Algo_withTuning::args.dataManipulationTime += t.elapsedTime().nanos();

            tune();

            std::cout << "dataManipulationTime: " << Algo_withTuning::args.dataManipulationTime / 1E9 << std::endl;
            std::cout << "Multiplication starts! k: " << args.k << std::endl;



            t.start();
            col_type maxLists = 1;

            switch (args.method) {
                case LEMP_I:
                case LEMP_LI:
                case LEMP_C:
                case LEMP_LC:

                    std::for_each(probeBuckets.begin(), probeBuckets.begin() + numBucketsToTune, [&maxLists](const ProbeBucket & b) {
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


#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();


                for (row_type b = 0; b < numBucketsToFireQueries; ++b) {

                    retrievers[b]->runTopK(probeBuckets[b], &retrArg[tid]);

                    if (args.method == LEMP_AP) { // synchronize the worstMinScore      

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

                localToGlobalIds(retrArg[tid].topkResults, args.k, retrArg[tid].results, queryMatrices[tid]);
                results[tid] = &retrArg[tid].results;
                comparisons += retrArg[tid].comparisons;

            }


            t.stop();


            std::cout << "Time for Retrieval: " << t << std::endl;
            std::cout << "Size of result: " << getResultSetSize() << std::endl;
            std::cout << "Comparisons: " << comparisons << std::endl;

            logging << "\t" << args.k << "\t" << comparisons << "\t" << getResultSetSize() << "\t";
            printTimes(t);

            if (args.resultsFile != "") {
                writeResults(results, args.resultsFile);
            }

            logging.close();
        }

        // this is mainly for debugging

        std::vector < QueueElement> & getResultsTopk() {
            return retrArg[0].topkResults;
        }

        std::vector< std::vector<MatItem >* >& getResults() {
            return results;
        }

        inline int getResultSetSize() {
            int count = 0;
            for (auto& element : results)
                count += element->size();

            return count;
        }

    };

    inline row_type Algo_withTuning::initProbeBuckets(VectorMatrix& rightMatrix) {
        std::vector<row_type> probeBucketOffsets;

        probeMatrix.init(rightMatrix, true, false); // normalize and sort

        row_type maxBlockSize = computeBlockOffsetsByFactorCacheFittingForItems(probeMatrix.lengthInfo,
                probeMatrix.rowNum, probeBucketOffsets, FACTOR, ITEMS_PER_BLOCK, args.cacheSizeinKB, probeMatrix.colNum, args);

        bucketize(probeBuckets, probeMatrix, probeBucketOffsets, args);

        std::cout << "Probe Buckets: " << probeBucketOffsets.size() << std::endl;

        logging << probeBucketOffsets.size() << "\t";

        return maxBlockSize;
    }

    inline void Algo_withTuning::initializeRetrievers() {

        numBucketsToFireQueries = probeBuckets.size();
        row_type b0 = 0;

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
                    numBucketsToFireQueries = i;
                    break;
                }
            }
            retrievers.reserve(numBucketsToFireQueries);
        } else { // Row-top-k

            std::for_each(probeBuckets.begin() + 1, probeBuckets.end(), [this](const ProbeBucket & b) {
                if (maxProbeBucketSize < b.rowNum) maxProbeBucketSize = b.rowNum;
            });

            b0 = 1;
            retrievers.reserve(numBucketsToFireQueries);
            retrievers.push_back(std::unique_ptr<Retriever >(new Retriever()));

        }

        numBucketsToInitIndexes = numBucketsToFireQueries;



        rg::Timer t;
        t.start();

        switch (args.method) {
            case LEMP_LI:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<LX_Retriever<IncrRetriever> >(new LX_Retriever<IncrRetriever> ()));
                }
                break;

            case LEMP_I:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<IncrRetriever >(new IncrRetriever()));
                }
                break;

            case LEMP_LC:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<LX_Retriever<CoordRetriever> >(new LX_Retriever<CoordRetriever> ()));
                }
                break;

            case LEMP_C:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<CoordRetriever >(new CoordRetriever()));
                }
                break;

            case LEMP_TA:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<taRetriever >(new taRetriever()));
                }
                break;

            case LEMP_L:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<LengthRetriever >(new LengthRetriever()));
                }
                break;

            case LEMP_TREE:
#pragma omp parallel for schedule(static,1) 
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<SingleTree >(new SingleTree()));
                }
                break;

            case LEMP_AP:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<apRetriever >(new apRetriever()));
                }
                break;

            case LEMP_LSH:
#pragma omp parallel for schedule(static,1)
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    retrievers.push_back(std::unique_ptr<LshRetriever >(new LshRetriever()));
                }
                break;
        }

        t.stop();

        std::cout << "Retriever: " << t << std::endl;
    }

    inline void Algo_withTuning::initQueryBatches(VectorMatrix& leftMatrix, row_type maxBlockSize) {
        row_type nCount = 0;

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
            retrArg[tid].initializeBasics(queryMatrices[tid], probeMatrix, args.method, args.theta, args.k, args.threads, args.R, args.epsilon, true, args.isTARR);

        }
        logging << nCount << "\t";
        std::cout << "Query Batches: " << nCount << std::endl;

    }

    inline void Algo_withTuning::initListsInBuckets() {

        double maxQueryLength = 0;
        row_type b0 = (args.k == 0 ? 0 : 1);
        std::cout << "Active Buckets used for index initialization: " << numBucketsToInitIndexes << std::endl;
        double worstCaseTheta;

        rg::Timer t1;
        t1.start();
        switch (args.method) {
            case LEMP_LI:
            case LEMP_I:
            case LEMP_TA:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    probeBuckets[b].ptrIndexes[SL] = new QueueElementLists();

                    if (b < numBucketsToInitIndexes)
                        static_cast<QueueElementLists*> (probeBuckets[b].ptrIndexes[SL])->initializeLists(probeMatrix, probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_LC:
            case LEMP_C:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    probeBuckets[b].ptrIndexes[INT_SL] = new IntLists();

                    if (b < numBucketsToInitIndexes)
                        static_cast<IntLists*> (probeBuckets[b].ptrIndexes[INT_SL])->initializeLists(probeMatrix, probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_TREE:
#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    probeBuckets[b].ptrIndexes[TREE] = new TreeIndex();
                    if (b < numBucketsToInitIndexes)
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
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    probeBuckets[b].ptrIndexes[AP] = new L2apIndex();
                    if (b < numBucketsToInitIndexes)
                        static_cast<L2apIndex*> (probeBuckets[b].ptrIndexes[AP])->initializeLists(probeMatrix, worstCaseTheta, cweights,
                            probeBuckets[b].startPos, probeBuckets[b].endPos);
                }
                break;

            case LEMP_LSH:


#pragma omp parallel for schedule(dynamic,1) 
                for (row_type b = b0; b < numBucketsToFireQueries; ++b) {
                    probeBuckets[b].ptrIndexes[LSH] = new LshIndex();
                    if (b < numBucketsToInitIndexes) {
                        static_cast<LshIndex*> (probeBuckets[b].ptrIndexes[LSH])->initializeLists(probeMatrix, true, probeBuckets[b].startPos, probeBuckets[b].endPos);
                        //                    row_type signatures = LSH_SIGNATURES; //LSH_SIGNATURES/(b+1);
                        //                    if (retrArg.size() > 1 && signatures > 0) {
                        //                        row_type tid = omp_get_thread_num();
                        //
                        //                        static_cast<LshIndex*> (probeBuckets[b].ptrIndexes[LSH])->checkAndReallocateAll(&probeMatrix, true, probeBuckets[b].startPos, probeBuckets[b].endPos, signatures,
                        //                                retrArg[tid].sums, retrArg[tid].countsOfBlockValues, retrArg[tid].sketches);
                        //                    }

                    }



                }

                break;

        }
        t1.stop();
        std::cout << "Time for initializing lists: " << t1 << std::endl;


    }

    inline void Algo_withTuning::tune() {

        if (numBucketsToFireQueries > 0) {
            switch (args.method) {
                case LEMP_LI:
                case LEMP_I:
                case LEMP_LC:
                case LEMP_C:
                case LEMP_LSH:

                    if (args.k == 0) {

                        numBucketsToTune = numBucketsToFireQueries;

                        t.start();
                        // first set-up the xValues in each retriever
                        for (row_type b = 0; b < numBucketsToTune; ++b) {
                            probeBuckets[b].sampling(retrArg);
                        }

                        // then do the actual tuning
                        for (row_type b = 0; b < numBucketsToTune; ++b) {
                            retrievers[b]->tune(probeBuckets[b], (b == 0 ? probeBuckets[b] : probeBuckets[b - 1]), retrArg);
                        }

                        // then release memory not needed any more
                        for (row_type b = 0; b < numBucketsToTune; ++b) {
                            retrievers[b]->cleanupAfterTuning();
                            probeBuckets[b].cleanupAfterTuning();
                        }

                        t.stop();
                        args.tuningTime += t.elapsedTime().nanos();

                    } else {
                        t.start();
                        std::pair<row_type, row_type> p = findSampleForTuningTopk(probeBuckets, retrArg);


                        numBucketsToTune = p.first;
#pragma omp parallel for schedule(dynamic,1) 
                        for (row_type b = 1; b < numBucketsToTune; ++b) {
                            probeBuckets[b].setup_xValues_topk(retrArg, probeBuckets[b - 1].sampleThetas);
                        }

                        t.stop();
                        args.tuningTime += t.elapsedTime().nanos();

                        t.start();
                        numBucketsToInitIndexes = p.second; // this change is mainly for multithreading. It will force us to initialize more lists up-front --> less locking later
                        initListsInBuckets();
                        t.stop();
                        args.dataManipulationTime += t.elapsedTime().nanos();

                        t.start();

                        for (row_type b = 1; b < probeBuckets.size(); ++b) {

                            if (b < numBucketsToTune) {
                                retrievers[b]->tuneTopk(probeBuckets[b], probeBuckets[b - 1], retrArg);
                            } else {
                                probeBuckets[b].setAfterTuning(probeBuckets[b - 1].numLists, probeBuckets[b - 1].t_b);
                            }

                        }


                        // then release memory not needed any more
                        for (row_type b = 1; b < numBucketsToTune; ++b) {
                            retrievers[b]->cleanupAfterTuning();
                            probeBuckets[b].cleanupAfterTuning();
                        }

                        t.stop();
                        args.tuningTime += t.elapsedTime().nanos();
                    }

                    break;
            }
        }


    }

    inline void Algo_withTuning::printTimes(rg::Timer& tAll) {


        args.printTimes();
        retrArg[0].printTimes();

        std::cout << "Total Time: " << (args.dataManipulationTime / 1E9) + tAll.elapsedTime().seconds()+(args.tuningTime / 1E9) << std::endl;
        args.dataManipulationTime -= matrixToMatrixTime; // being nice to L2AP

        double retrTime = tAll.elapsedTime().nanos();
        double timeAll = args.dataManipulationTime + args.tuningTime + retrTime;

        logging << ((args.dataManipulationTime) / 1E9) << "\t" << (args.tuningTime / 1E9) << "\t" << retrTime / 1E9 << "\t" << timeAll / 1E9 << "\n";
        std::cout << ((args.dataManipulationTime) / 1E9) << "\t" << (args.tuningTime / 1E9) << "\t" << retrTime / 1E9 << "\t" << timeAll / 1E9 << std::endl;


    }

}


#endif /* ALGO_WITH_TUNING_H_ */
