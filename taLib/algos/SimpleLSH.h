/* 
 * File:   LSH_all.h
 * Author: chteflio
 *
 * Created on March 26, 2015, 9:16 AM
 */

#ifndef LSH_ALL_H
#define	LSH_ALL_H




namespace ta {

    class SimpleLSH {
    public:
        rg::Timer t;
        std::shared_ptr<VectorMatrix> queryMatrix, probeMatrix;
        ProbeBucket probeBucket, probeBucketK;
   

        std::vector< MatItem >* thetaResults; // for a specific query holds the itemIDs + the score
        std::vector<QueueElement> * topkResults;

        LEMPArg args;
        RetrievalArguments* retrArg;
        std::ofstream logging;

        double dataManipulationTime;

        inline SimpleLSH(LEMPArg& args) : args(args), dataManipulationTime(0) {
            queryMatrix = std::make_shared<VectorMatrix>();
            probeMatrix = std::make_shared<VectorMatrix>();

            if (args.querySideLeft) {
                queryMatrix->readFromFile(args.usersFile, args.r, args.m, true);
                probeMatrix->readFromFile(args.itemsFile, args.r, args.n, false);
            } else {
                queryMatrix->readFromFile(args.itemsFile, args.r, args.n, false);
                probeMatrix->readFromFile(args.usersFile, args.r, args.m, true);
            }

            if (!args.isTransformed) {
                t.start();
                SimpleLSH::transform();
                t.stop();
                std::cout << "transform time: " << t << std::endl;
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

            retrArg = new RetrievalArguments(probeMatrix->colNum, queryMatrix.get(), probeMatrix.get(), LEMP_LSH);
            retrArg->k = args.k;
            retrArg->R = args.R;

            retrArg->init(probeMatrix->rowNum);

            t.start();

            probeBucketK.init(*probeMatrix, 0, args.k, args);
            probeBucket.init(*probeMatrix, args.k, probeMatrix->rowNum, args); // initialize

            if (probeBucket.ptrIndexes[LSH] == 0)
                probeBucket.ptrIndexes[LSH] = new LshIndex;

            static_cast<LshIndex*> (probeBucket.ptrIndexes[LSH])->initializeLists(*probeMatrix, true, args.k, probeMatrix->rowNum);

            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();

            logging << "SimpleLSH";

            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";
        }

        inline void transform() {
            // first transform queryMatrix (transform ||q|| to 1 and q = [0;q])
            std::cout << "Transform queryMatrix" << std::endl;

            std::shared_ptr<VectorMatrix> tmpQueryMatrix(new VectorMatrix());

            tmpQueryMatrix->rowNum = queryMatrix->rowNum;
            tmpQueryMatrix->colNum = queryMatrix->colNum + 1;
            tmpQueryMatrix->initializeBasics(tmpQueryMatrix->colNum, tmpQueryMatrix->rowNum, false);
            for (row_type i = 0; i < tmpQueryMatrix->rowNum; i++) {
                double* dQuery = queryMatrix->getMatrixRowPtr(i);
                double* dTmp = tmpQueryMatrix->getMatrixRowPtr(i);
                // division is expensive! multiply with inverse instead
                double invLen = 1 / calculateLength(queryMatrix->getMatrixRowPtr(i), queryMatrix->colNum);
                // all but last coordinates
                scaleAndCopy(dTmp + 1, dQuery, invLen, queryMatrix->colNum);
                // last coordinate
                dTmp[0] = 0;
                tmpQueryMatrix->setLengthInData(i, 1); // ||q|| = 1
            }

            // then transform probeMatrix (transform ||p|| to less than 1 and p = [sqrt(1- ||p|| * ||p||);p])
            std::cout << "Transform probeMatrix" << std::endl;

            // we need the longest vector from probeMatrix
            double maxLen = 0;
            for (row_type i = 0; i < probeMatrix->rowNum; i++) {
                double len = calculateLength(probeMatrix->getMatrixRowPtr(i), probeMatrix->colNum);
                if (len > maxLen) {
                    maxLen = len;
                }
            }
            double invMaxLen = 1 / maxLen;
            std::shared_ptr<VectorMatrix> tmpProbeMatrix(new VectorMatrix());

            tmpProbeMatrix->rowNum = probeMatrix->rowNum;
            tmpProbeMatrix->colNum = probeMatrix->colNum + 1;
            tmpProbeMatrix->initializeBasics(tmpProbeMatrix->colNum, tmpProbeMatrix->rowNum, false);
            for (row_type i = 0; i < tmpProbeMatrix->rowNum; i++) {
                double* dProbe = probeMatrix->getMatrixRowPtr(i);
                double* dTmp = tmpProbeMatrix->getMatrixRowPtr(i);
                // all but last coordinates
                scaleAndCopy(dTmp + 1, dProbe, invMaxLen, probeMatrix->colNum); // multiply with inverse
                // last coordinate
                double len = calculateLength(dTmp, probeMatrix->colNum); // use new values but without the new coordinate          
                dTmp[0] = ((1 - len * len) < 0) ? 0 : sqrt(1 - len * len);
                tmpProbeMatrix->setLengthInData(i, 1); // set to 1 by the transformation
            }
            // change pointers
            queryMatrix = tmpQueryMatrix;
            probeMatrix = tmpProbeMatrix;
        }

        inline void processIndexesTopk(double * query, row_type queryId,
                LshIndex* index, LshIndex* queryIndex, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            row_type numCandidatesToVerify = 0;

            index->getCandidates(queryIndex->getSketch(), queryId, arg->candidatesToVerify, numCandidatesToVerify,
                    arg->done, LSH_SIGNATURES, probeBucket.startPos);
            verifyCandidatesTopK_noLengthTest(query, numCandidatesToVerify, arg);
        }

        inline void runTopkPerUser() {

            retrArg->comparisons = 0;

            LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

            t.start();

            if (LSH_SIGNATURES > index->initializedSketchesForIndex)
                index->checkAndReallocateAll(retrArg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, LSH_SIGNATURES,
                    retrArg->sums, retrArg->countsOfBlockValues, retrArg->sketches);

            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();

            std::cout << "dataManipulationTime: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Multiplication starts! k: " << retrArg->k << std::endl;

            t.start();
            LshIndex queryIndex;
            queryIndex.initializeLists(*queryMatrix, false, 0, queryMatrix->rowNum);
            queryIndex.checkAndReallocateAll(retrArg->queryMatrix, false, 0, queryMatrix->rowNum, LSH_SIGNATURES,
                    retrArg->sums, retrArg->countsOfBlockValues, retrArg->sketches);

            retrArg->topkResults.resize(queryMatrix->rowNum * args.k);
            retrArg->heap.resize(args.k);

            for (row_type i = 0; i < queryMatrix->rowNum; i++) {

                double* query = queryMatrix->getMatrixRowPtr(i);
                retrArg->queryId = i;


                for (row_type j = 0; j < args.k; j++) {
                    double ip = queryMatrix->innerProduct(i, probeMatrix->getMatrixRowPtr(j));
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
            localToGlobalIds(*topkResults, args.k, results, *queryMatrix);
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

