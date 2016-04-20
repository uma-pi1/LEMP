// Copyright 2015 Christina Teflioudi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/* 
 * File:   TANRA_all.h
 * Author: chteflio
 *
 * Created on August 24, 2015, 9:53 AM
 */

#ifndef TANRA_ALL_H
#define	TANRA_ALL_H

namespace mips {

    class TaNra : public Mip {
        ProbeBucket probeBucket, probeBucketK;

        LempArguments args;
        std::vector<RetrievalArguments> retrArg;

        inline void printAlgoName(const VectorMatrix& leftMatrix) {
            if (args.isTARR) {
                logging << "TA_NRA (RR)" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] TA_NRA (ROUND ROBIN) with " << args.threads << " thread(s)" << std::endl;
            } else {
                logging << "TA_NRA (MAX)" << "\t" << args.threads << "\t";
                std::cout << "[ALGORITHM] TA_NRA (MAX HEAP) with " << args.threads << " thread(s)" << std::endl;
            }
            logging << "P(" << probeMatrix.rowNum << "x" << (0 + probeMatrix.colNum) << ")\t";
            logging << "Q^T(" << leftMatrix.rowNum << "x" << (0 + leftMatrix.colNum) << ")\t";
        }
        
        
        inline void initializeInternal(std::vector<VectorMatrix>& queryMatrices, VectorMatrix& leftMatrix) {

            std::cout << "[RETRIEVAL] QueryMatrix contains " << leftMatrix.rowNum << " vectors with dimensionality " << (0 + leftMatrix.colNum) << std::endl;
            row_type myNumThreads = args.threads;

            if (leftMatrix.rowNum < args.threads) {
                myNumThreads = leftMatrix.rowNum;
                std::cout << "[WARNING] Query matrix contains too few elements. Suboptimal running with " << myNumThreads << " thread(s)" << std::endl;
            }
            omp_set_num_threads(myNumThreads);
            queryMatrices.resize(myNumThreads);
            splitMatrices(leftMatrix, queryMatrices);



            for (row_type i = 0; i < myNumThreads; i++) {
                retrArg[i].initializeBasics(queryMatrices[i], probeMatrix, LEMP_TANRA, args.theta, args.k, args.threads, args.R, args.epsilon, 0, 0, false, args.isTARR);
                retrArg[i].init(probeMatrix.rowNum);
                retrArg[i].clear();
            }
        }
    public:

        inline void setTheta(double theta) {
            args.theta = theta;
        }

        inline TaNra(InputArguments& input, bool isTARR) {
            args.copyInputArguments(input);
            args.isTARR = isTARR;

            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);


            if (!logging.is_open()) {
                std::cout << "[WARNING] No log will be created!" << std::endl;
            } else {
                std::cout << "[INFO] Logging in " << args.logFile << std::endl;
            }

            omp_set_num_threads(args.threads);
            retrArg.resize(args.threads);

        };


        inline void initialize(VectorMatrix& rightMatrix) {
            std::cout << "[INIT] ProbeMatrix contains " << rightMatrix.rowNum << " vectors with dimensionality " << (0 + rightMatrix.colNum) << std::endl;


            probeMatrix = rightMatrix;

            timer.start();
            probeBucket.init(probeMatrix, 0, probeMatrix.rowNum, args); // initialize
            probeBucket.bucketScanThreshold = 0;
            probeBucket.normL2.second = std::numeric_limits<double>::max();

            retriever_ptr rPtr(new tanraRetriever());
            probeBucket.ptrRetriever = rPtr;
            if (probeBucket.ptrIndexes[SL] == 0)
                probeBucket.ptrIndexes[SL] = new QueueElementLists();

            static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL])->initializeLists(probeMatrix, 0, probeMatrix.rowNum);

            timer.stop();
            dataPreprocessingTimeRight += timer.elapsedTime().nanos();
        }

        inline ~TaNra() {
            logging.close();
        }

        inline void runAboveTheta(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            std::vector<VectorMatrix> queryMatrices;


            initializeInternal(queryMatrices, leftMatrix);
            results.resultsVector.resize(args.threads);

            std::vector<row_type> blockOffsets;
            blockOffsets.push_back(0);
            std::cout << "[RETRIEVAL] Retrieval (theta = " << args.theta << ") starts ..." << std::endl;
            logging << "theta(" << args.theta << ")\t";

            timer.start();

            comp_type comparisons = 0;
#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();
                bucketize(retrArg[tid].queryBatches, queryMatrices[tid], blockOffsets, args);

                probeBucket.ptrRetriever->run(probeBucket, &retrArg[tid]);

                results.moveAppend(retrArg[tid].results, tid);
                comparisons += retrArg[tid].comparisons;

            }

            timer.stop();
            retrievalTime += timer.elapsedTime().nanos();
            totalComparisons += comparisons;
            std::cout << "[RETRIEVAL] ... and is finished with " << results.getResultSize() << " results" << std::endl;
            logging << results.getResultSize() << "\t";

            outputStats();
        }

    };

}

#endif	/* TANRA_ALL_H */
