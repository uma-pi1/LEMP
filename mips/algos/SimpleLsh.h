/* 
 * File:   LSH_all.h
 * Author: chteflio
 *
 * Created on March 26, 2015, 9:16 AM
 */

#ifndef LSH_ALL_H
#define	LSH_ALL_H

namespace mips {

    class SimpleLsh : public Mip {
        ProbeBucket probeBucket, probeBucketK;
        std::vector<RetrievalArguments> retrArg;
        LempArguments args;

        inline void transformQueryMatrix(VectorMatrix& leftMatrix, VectorMatrix& queryMatrix) {
            // transform queryMatrix (transform ||q|| to 1 and q = [0;q])

            queryMatrix.rowNum = leftMatrix.rowNum;
            queryMatrix.colNum = leftMatrix.colNum + 1;
            queryMatrix.initializeBasics(queryMatrix.colNum, queryMatrix.rowNum, false);

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < queryMatrix.rowNum; i++) {
                double* dQuery = leftMatrix.getMatrixRowPtr(i);
                double* dTmp = queryMatrix.getMatrixRowPtr(i);
                // division is expensive! multiply with inverse instead
                double invLen = 1 / calculateLength(leftMatrix.getMatrixRowPtr(i), leftMatrix.colNum);
                // all but last coordinates
                scaleAndCopy(dTmp + 1, dQuery, invLen, leftMatrix.colNum);
                // last coordinate
                dTmp[0] = 0;
                queryMatrix.setLengthInData(i, 1); // ||q|| = 1
            }
        }

        inline void transformProbeMatrix(VectorMatrix& rightMatrix) {
            //  transform probeMatrix (transform ||p|| to less than 1 and p = [sqrt(1- ||p|| * ||p||);p])

            // we need the longest vector from probeMatrix
            double maxLen = 0;
            for (row_type i = 0; i < rightMatrix.rowNum; i++) {
                double len = calculateLength(rightMatrix.getMatrixRowPtr(i), rightMatrix.colNum);
                if (len > maxLen) {
                    maxLen = len;
                }
            }
            double invMaxLen = 1 / maxLen;

            probeMatrix.rowNum = rightMatrix.rowNum;
            probeMatrix.colNum = rightMatrix.colNum + 1;
            probeMatrix.initializeBasics(probeMatrix.colNum, probeMatrix.rowNum, false);

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < probeMatrix.rowNum; i++) {
                double* dProbe = rightMatrix.getMatrixRowPtr(i);
                double* dTmp = probeMatrix.getMatrixRowPtr(i);
                // all but last coordinates
                scaleAndCopy(dTmp + 1, dProbe, invMaxLen, rightMatrix.colNum); // multiply with inverse
                // last coordinate
                double len = calculateLength(dTmp, rightMatrix.colNum); // use new values but without the new coordinate          
                dTmp[0] = ((1 - len * len) < 0) ? 0 : sqrt(1 - len * len);
                probeMatrix.setLengthInData(i, 1); // set to 1 by the transformation
            }
        }

        inline void printAlgoName(const VectorMatrix& leftMatrix) {
            logging << "SIMPLE_LSH" << "\t" << args.threads << "\t";
            std::cout << "[ALGORITHM] SIMPLE_LSH with " << args.threads << " thread(s)" << std::endl;

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

            timer.start();
            if (!isTransformed) {
                std::cout << "[RETRIEVAL] QueryMatrix will be transformed" << std::endl;

                VectorMatrix queryMatrix;
                SimpleLsh::transformQueryMatrix(leftMatrix, queryMatrix);
                splitMatrices(queryMatrix, queryMatrices);
            } else {
                splitMatrices(leftMatrix, queryMatrices);
            }
            timer.stop();
            dataPreprocessingTimeLeft += timer.elapsedTime().nanos();


            for (row_type i = 0; i < myNumThreads; i++) {
                retrArg[i].initializeBasics(queryMatrices[i], probeMatrix, LEMP_LSH, args.theta, args.k, myNumThreads, 1, 0, 0, 0, false, false);
                retrArg[i].init(probeMatrix.rowNum);
                retrArg[i].clear();
            }
        }

        inline void processIndexesTopk(double * query, row_type queryId,
                LshIndex* index, LshIndex* queryIndex, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            row_type numCandidatesToVerify = 0;

            index->lshBins->getCandidates(queryIndex->cosSketches->sketches, queryId, arg->candidatesToVerify, numCandidatesToVerify,
                    arg->done, LSH_SIGNATURES, probeBucket.startPos);
            verifyCandidatesTopK_noLengthTest(query, numCandidatesToVerify, arg);
        }


    public:

        bool isTransformed;

        inline SimpleLsh(InputArguments& input, bool isTransformed) : isTransformed(isTransformed) {
            args.copyInputArguments(input);

            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (!logging.is_open()) {
                std::cout << "[WARNING] No log will be created!" << std::endl;
            } else {
                std::cout << "[INFO] Logging in " << args.logFile << std::endl;
            }

            omp_set_num_threads(args.threads);
            retrArg.resize(args.threads);
        }

        inline ~SimpleLsh() {
            logging.close();
        }

        void initialize(VectorMatrix& rightMatrix) {
            std::cout << "[INIT] ProbeMatrix contains " << rightMatrix.rowNum << " vectors with dimensionality " << (0 + rightMatrix.colNum) << std::endl;
            logging << "P(" << rightMatrix.rowNum << "x" << (0 + rightMatrix.colNum) << ")\t";

            if (!isTransformed) {
                std::cout << "[INIT] ProbeMatrix will be transformed" << std::endl;
                timer.start();
                transformProbeMatrix(rightMatrix);
                timer.stop();
                dataPreprocessingTimeRight += timer.elapsedTime().nanos();
            } else {
                probeMatrix = rightMatrix;
            }


            probeBucketK.init(probeMatrix, 0, args.k, args);
            probeBucket.init(probeMatrix, args.k, probeMatrix.rowNum, args); // initialize

            if (probeBucket.ptrIndexes[LSH] == 0)
                probeBucket.ptrIndexes[LSH] = new LshIndex;

            static_cast<LshIndex*> (probeBucket.ptrIndexes[LSH])->initializeLists(probeMatrix, true, args.k, probeMatrix.rowNum);
        }

        inline void runTopK(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            std::vector<VectorMatrix> queryMatrices;
            initializeInternal(queryMatrices, leftMatrix);
            results.resultsVector.resize(args.threads);

            LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

            timer.start();

            if (LSH_SIGNATURES > index->initializedSketchesForIndex) {
                index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, LSH_SIGNATURES,
                        retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, false);
            }

            timer.stop();
            dataPreprocessingTimeRight += timer.elapsedTime().nanos();

            std::cout << "[RETRIEVAL] Retrieval (k = " << args.k << ") starts ..." << std::endl;
            logging << "k(" << args.k << ")\t";

            timer.start();

            comp_type comparisons = 0;
#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();

                LshIndex queryIndex; // separate for each thread
                queryIndex.initializeLists(queryMatrices[tid], false, 0, queryMatrices[tid].rowNum);
                queryIndex.checkAndReallocateAll(retrArg[tid].queryMatrix, false, 0, queryMatrices[tid].rowNum, LSH_SIGNATURES,
                        retrArg[tid].sums, retrArg[tid].countsOfBlockValues, retrArg[tid].sketches, false);

                retrArg[tid].allocTopkResults();

                for (row_type i = 0; i < queryMatrices[tid].rowNum; i++) {
                    double* query = queryMatrices[tid].getMatrixRowPtr(i);
                    retrArg[tid].queryId = i;

                    for (row_type j = 0; j < args.k; j++) {
                        double ip = queryMatrices[tid].innerProduct(i, probeMatrix.getMatrixRowPtr(j));
                        retrArg[tid].comparisons++;
                        retrArg[tid].heap[j] = QueueElement(ip, j);
                    }

                    std::make_heap(retrArg[tid].heap.begin(), retrArg[tid].heap.end(), std::greater<QueueElement>()); //make the heap;
                    processIndexesTopk(query, i, index, &queryIndex, probeBucket, &retrArg[tid]);
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

#endif	/* LSH_ALL_H */
