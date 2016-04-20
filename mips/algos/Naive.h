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
 * naive.h
 *
 *  Created on: Oct 10, 2013
 *      Author: chteflio
 */

#ifndef NAIVEALGO_H_
#define NAIVEALGO_H_



namespace mips {

    class Naive : public Mip {
        InputArguments args;
        std::vector<RetrievalArguments> retrArg;

        inline void printAlgoName(const VectorMatrix& leftMatrix) {
            logging << "NAIVE" << "\t" << args.threads << "\t";
            std::cout << "[ALGORITHM] NAIVE with " << args.threads << " thread(s)" << std::endl;

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

            if (retrArg.size() < args.threads)
                retrArg.resize(myNumThreads);

            for (row_type i = 0; i < retrArg.size(); i++) {
                retrArg[i].initializeBasics(queryMatrices[i], probeMatrix, LEMP_L, args.theta, args.k, myNumThreads, 1, 0, 0,0, false, false);
                retrArg[i].clear();
            }
        }

    public:

        inline Naive(InputArguments& args) : args(args) {
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (!logging.is_open()) {
                std::cout << "[WARNING] No log will be created!" << std::endl;
            } else {
                std::cout << "[INFO] Logging in " << args.logFile << std::endl;
            }

            omp_set_num_threads(args.threads);
        }

        inline ~Naive() {
            logging.close();
        }

        inline void setTheta(double theta) {
            args.theta = theta;
        }

        inline void initialize(VectorMatrix& rightMatrix) {
            std::cout << "[INIT] ProbeMatrix contains " << rightMatrix.rowNum << " vectors with dimensionality " << (0 + rightMatrix.colNum) << std::endl;
            probeMatrix = rightMatrix;
        }


        void runAboveTheta(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            std::vector<VectorMatrix> queryMatrices;
            initializeInternal(queryMatrices, leftMatrix);
            results.resultsVector.resize(args.threads);

            std::cout << "[RETRIEVAL] Retrieval (theta = " << args.theta << ") starts ..." << std::endl;
            logging << "theta(" << args.theta << ")\t";

            timer.start();

            comp_type comparisons = 0;
#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();

                LengthRetriever plainRetriever;
                for (row_type i = 0; i < queryMatrices[tid].rowNum; ++i) {
                    double* query = queryMatrices[tid].getMatrixRowPtr(i);
                    retrArg[tid].queryId = i;
                    plainRetriever.naive(query, 0, probeMatrix.rowNum, &retrArg[tid]);
                }

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

        void runTopK(VectorMatrix& leftMatrix, Results& results) {

            printAlgoName(leftMatrix);

            std::vector<VectorMatrix> queryMatrices;
            initializeInternal(queryMatrices, leftMatrix);
            results.resultsVector.resize(args.threads);

            std::cout << "[RETRIEVAL] Retrieval (k = " << args.k << ") starts ..." << std::endl;
            logging << "k(" << args.k << ")\t";

            timer.start();

            comp_type comparisons = 0;
#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();

                retrArg[tid].allocTopkResults();

                LengthRetriever plainRetriever;
                for (row_type i = 0; i < queryMatrices[tid].rowNum; ++i) {
                    double* query = queryMatrices[tid].getMatrixRowPtr(i);
                    retrArg[tid].queryId = i;

                    for (row_type j = 0; j < args.k; ++j) {
                        double ip = queryMatrices[tid].innerProduct(i, probeMatrix.getMatrixRowPtr(j));
                        retrArg[tid].comparisons++;
                        retrArg[tid].heap[j] = QueueElement(ip, j);
                    }

                    std::make_heap(retrArg[tid].heap.begin(), retrArg[tid].heap.end(), std::greater<QueueElement>()); //make the heap;
                    plainRetriever.naiveTopk(query, retrArg[tid].k, probeMatrix.rowNum, &retrArg[tid]);
                    retrArg[tid].writeHeapToTopk(i);
                }


                retrArg[tid].extendIncompleteResultItems();
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



#endif /* NAIVE_H_ */
