/* 
 * File:   LSH_all.h
 * Author: chteflio
 *
 * Created on March 26, 2015, 9:16 AM
 */

#ifndef LSH_ALL_H
#define	LSH_ALL_H


namespace ta {

    /* this files are for running the SimpleLSH method AFTER the transformation of the input vectors. I.e. the datasets read need to be the transformed ones.
     */

    class SimpleLSH {
    public:
        rg::Timer t;
        VectorMatrix queryMatrix, probeMatrix;
        ProbeBucket probeBucket, probeBucketK;

        std::vector< MatItem >* thetaResults; // for a specific query holds the itemIDs + the score
        std::vector<QueueElement> * topkResults;


        LEMPArg args;
        RetrievalArguments* retrArg;
        std::ofstream logging;

        double dataManipulationTime;

        inline SimpleLSH(LEMPArg& args) : args(args), dataManipulationTime(0) {
            
             if (args.querySideLeft) {
                queryMatrix.readFromFile(args.usersFile, args.r, args.m, true);
                probeMatrix.readFromFile(args.itemsFile, args.r, args.n, false);
            } else {
                queryMatrix.readFromFile(args.itemsFile, args.r, args.n, false);
                probeMatrix.readFromFile(args.usersFile, args.r, args.m, true);
            }

            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (logging.is_open()) {

                std::cout << "Writing output to " << args.logFile << std::endl;
            } else {

                std::cout << "Problem with opening output file" << std::endl;
            }

            std::cout << "ALGO: SimpleLSH" << std::endl;
            std::cout << "Threads: " << args.threads << std::endl;
            std::cout << "LSH_SIGNATURES: " << LSH_SIGNATURES << std::endl;
            std::cout << "LSH_CODE_LENGTH: " << LSH_CODE_LENGTH << std::endl;

            retrArg = new RetrievalArguments(probeMatrix.colNum, &queryMatrix, &probeMatrix, LEMP_LSH); 
            retrArg->k = args.k;
            retrArg->R = args.R;

            retrArg->init(probeMatrix.rowNum);

            t.start();

            probeBucketK.init(probeMatrix, 0, args.k, args);
            probeBucket.init(probeMatrix, args.k, probeMatrix.rowNum, args); // initialize


            if (probeBucket.ptrIndexes[LSH] == 0)
                probeBucket.ptrIndexes[LSH] = new LshIndex;

            static_cast<LshIndex*> (probeBucket.ptrIndexes[LSH])->initializeLists(probeMatrix, true, args.k, probeMatrix.rowNum);

            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();

            logging << "SimpleLSH";

            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";

        }


        inline void processIndexesTopk(double * query, row_type queryId,
                LshIndex* index, LshIndex* queryIndex, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            row_type numCandidatesToVerify = 0;

            index->lshBins->getCandidates(queryIndex->cosSketches->sketches, queryId, arg->candidatesToVerify, numCandidatesToVerify,
                    arg->done, LSH_SIGNATURES, probeBucket.startPos);
            verifyCandidatesTopK_noLengthTest(query, numCandidatesToVerify, arg);
        }



        inline void runTopkPerUser() {
            retrArg->comparisons = 0;

            LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

            t.start();

            if (LSH_SIGNATURES > index->initializedSketchesForIndex)
                index->checkAndReallocateAll(retrArg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, LSH_SIGNATURES,
                    retrArg->sums, retrArg->countsOfBlockValues, retrArg->sketches, retrArg->rig);
  		
            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();   
            
            std::cout << "dataManipulationTime: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Multiplication starts! k: " << retrArg->k << std::endl;
            
            t.start();
            LshIndex queryIndex;
            queryIndex.initializeLists(queryMatrix, false, 0, queryMatrix.rowNum);
            queryIndex.checkAndReallocateAll(retrArg->queryMatrix, false, 0, queryMatrix.rowNum, LSH_SIGNATURES,
                    retrArg->sums, retrArg->countsOfBlockValues, retrArg->sketches, retrArg->rig);

            retrArg->topkResults.resize(queryMatrix.rowNum * args.k);
            retrArg->heap.resize(args.k);

            for (row_type i = 0; i < queryMatrix.rowNum; i++) {

                double* query = queryMatrix.getMatrixRowPtr(i);
                retrArg->queryId = i;

                for (row_type j = 0; j < args.k; j++) {
                    double ip = queryMatrix.innerProduct(i, probeMatrix.getMatrixRowPtr(j));
                    retrArg->comparisons++;
                    retrArg->heap[j] = QueueElement(ip, j);
                }

                std::make_heap(retrArg->heap.begin(), retrArg->heap.end(), std::greater<QueueElement>()); //make the heap;

                processIndexesTopk(query, i, index, &queryIndex, probeBucket, retrArg);
            
                retrArg->writeHeapToTopk(i);

            }

            topkResults = &(retrArg->topkResults);
            t.stop();

            std::cout << "Time for retrieval: " << t.elapsedTime().seconds() << std::endl;
            std::cout << "Comparisons: " << retrArg->comparisons << std::endl;
            std::cout << "Size of result: " << getResultSetSize() << std::endl;
            std::cout << "Preprocessing time: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Total time: " << (dataManipulationTime / 1E9) + t.elapsedTime().seconds() << std::endl;

            std::vector<MatItem> results;
            localToGlobalIds(*topkResults, args.k, results, queryMatrix);
            std::vector< std::vector<MatItem >* > resultsForWriting;
            resultsForWriting.push_back(&results);

            if (args.resultsFile != "") {
                writeResults(resultsForWriting, args.resultsFile);
            }


            logging << "\t" << args.k << "\t" << LSH_SIGNATURES << "\t" << retrArg->comparisons << "\t" << getResultSetSize() << "\t";
            printTimes(t);


            logging.close();
        }

        inline void printTimes(rg::Timer& tAll) {

            logging << (dataManipulationTime / 1E9) << "\t" << tAll.elapsedTime().seconds() << "\t" << tAll.elapsedTime().seconds()+(dataManipulationTime / 1E9) << "\n";
        }

        inline int getResultSetSize() {
            if (topkResults->size() > 0) {
                return topkResults->size();
            } else {
                return thetaResults->size();
            }
        }

    };

}


#endif	/* LSH_ALL_H */

