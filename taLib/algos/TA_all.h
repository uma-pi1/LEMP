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
 * TA_all.h
 *
 *  Created on: Mar 10, 2014
 *      Author: chteflio
 */

#ifndef TA_ALL_H_
#define TA_ALL_H_



namespace ta {

    class TA_all {
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

        inline TA_all(LEMPArg& args) : args(args), dataManipulationTime(0) {

            if (args.querySideLeft) {
                queryMatrix.readFromFile(args.usersFile, true); // do not normalize and do not sort
                probeMatrix.readFromFile(args.itemsFile, false); // do not normalize and do not sort
            } else {
                queryMatrix.readFromFile(args.itemsFile, false); // do not normalize and do not sort
                probeMatrix.readFromFile(args.usersFile, true); // do not normalize and do not sort
            }

            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (logging.is_open()) {
                std::cout << "Writing output to " << args.logFile << std::endl;
            } else {
                std::cout << "Problem with opening log-file. No log-file will be created" << std::endl;
            }
            if (args.isTARR){
                std::cout << "ALGO: TA_all with RR" << std::endl;
            }else{
                std::cout << "ALGO: TA_all with MaxPiQi" << std::endl;
            }
            
            std::cout << "Threads: " << args.threads << std::endl;

            retrArg = new RetrievalArguments(probeMatrix.colNum, queryMatrix, probeMatrix, LEMP_TA, false, args.isTARR);
            retrArg->k = args.k;
            retrArg->theta = args.theta;

            retrArg->init(0);

            t.start();
            if (args.k == 0) {
                probeBucket.init(probeMatrix, 0, probeMatrix.rowNum, args); // initialize
                probeBucket.bucketScanThreshold = args.theta / probeBucket.normL2.second;

                retriever_ptr rPtr(new taRetriever());
                probeBucket.ptrRetriever = rPtr;
                if (probeBucket.ptrIndexes[SL] == 0)
                    probeBucket.ptrIndexes[SL] = new QueueElementLists();

                static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL])->initializeLists(probeMatrix, 0, probeMatrix.rowNum);

            } else {
                probeBucketK.init(probeMatrix, 0, args.k, args);
                probeBucket.init(probeMatrix, args.k, probeMatrix.rowNum, args); // initialize

                retriever_ptr firstPtr(new Retriever());
                probeBucketK.ptrRetriever = firstPtr;

                retriever_ptr rPtr(new taRetriever());
                probeBucket.ptrRetriever = rPtr;
                if (probeBucket.ptrIndexes[SL] == 0)
                    probeBucket.ptrIndexes[SL] = new QueueElementLists();

                static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL])->initializeLists(probeMatrix, args.k, probeMatrix.rowNum);
            }
            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();

            logging << "TA_all";

            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";


        };

        inline ~TA_all() {
        }

        inline void multiply() {
            args.comparisons = 0;
            std::cout << "Multiplication starts! theta: " << retrArg->theta << std::endl;

            t.start();
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL]);
            retrArg->state->initializeForNewBucket(invLists);

            for (row_type i = 0; i < queryMatrix.rowNum; i++) {
                const double* query = queryMatrix.getMatrixRowPtr(i);
                retrArg->queryId = i;
                probeBucket.ptrRetriever->run(query, probeBucket, retrArg);
            }
            thetaResults = &(retrArg->results);            
            t.stop();

            std::cout << "Time for retrieval: " << t << std::endl;
            std::cout << "Comparisons: " << retrArg->comparisons << std::endl;
            std::cout << "Size of result: " << getResultSetSize() << std::endl;
            std::cout << "Preprocessing time: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Total time: " << (dataManipulationTime / 1E9) + t.elapsedTime().seconds() << std::endl;


//            std::cout << "preprocessTime: " << retrArg->preprocessTime / 1E9 << std::endl;
//            std::cout << "ipTime: " << retrArg->ipTime / 1E9 << std::endl;
//            std::cout << "boundsTime: " << retrArg->boundsTime / 1E9 << std::endl;
//            std::cout << "scanTime: " << retrArg->scanTime / 1E9 << std::endl;
//            std::cout << "filterTime: " << retrArg->filterTime / 1E9 << std::endl;


            logging << "\t" << args.theta << "\t" << retrArg->comparisons << "\t" << getResultSetSize() << "\t";
            printTimes(t);

            if (args.resultsFile != "") {
                std::vector< std::vector<MatItem >* > resultsForWriting;
                resultsForWriting.push_back(thetaResults);
                writeResults(resultsForWriting, args.resultsFile);
            }

            logging.close();
        }

        inline void runTopkPerUser() {
            args.comparisons = 0;
            std::cout << "Multiplication starts! k: " << retrArg->k << std::endl;
            t.start();
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));
            retrArg->state->initializeForNewBucket(invLists);

            /////////////////////
            retrArg->allocTopkResults();

            for (row_type i = 0; i < queryMatrix.rowNum; i++) {

                double* query = queryMatrix.getMatrixRowPtr(i);
                retrArg->queryId = i;

                for (row_type j = 0; j < args.k; j++) {
                    double ip = queryMatrix.innerProduct(i, probeMatrix.getMatrixRowPtr(j));
                    retrArg->comparisons++;
                    retrArg->heap[j] = QueueElement(ip, j);
                }

                std::make_heap(retrArg->heap.begin(), retrArg->heap.end(), std::greater<QueueElement>()); //make the heap;
                probeBucket.ptrRetriever->runTopK(query, probeBucket, retrArg);
                retrArg->writeHeapToTopk(i);

            }

            topkResults = &(retrArg->topkResults);
            t.stop();

            std::cout << "Time for retrieval: " << t << std::endl;
            std::cout << "Comparisons: " << retrArg->comparisons << std::endl;
            std::cout << "Size of result: " << getResultSetSize() << std::endl;
            std::cout << "Preprocessing time: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Total time: " << (dataManipulationTime / 1E9) + t.elapsedTime().seconds() << std::endl;

//            std::cout << "preprocessTime: " << retrArg->preprocessTime / 1E9 << std::endl;
//            std::cout << "ipTime: " << retrArg->ipTime / 1E9 << std::endl;
//            std::cout << "boundsTime: " << retrArg->boundsTime / 1E9 << std::endl;
//            std::cout << "scanTime: " << retrArg->scanTime / 1E9 << std::endl;
//            std::cout << "filterTime: " << retrArg->filterTime / 1E9 << std::endl;


            logging << "\t" << args.k << "\t" << retrArg->comparisons << "\t" << getResultSetSize() << "\t";
            printTimes(t);
            logging.close();
        }

        inline void runTopkPerUser2() { // parallel version
            comp_type comparisons = 0;


            std::cout << "Multiplication starts! k: " << retrArg->k << std::endl;
            t.start();


            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));
            retrArg->state->initializeForNewBucket(invLists);
            omp_set_num_threads(args.threads);
            /////////////////////

            retrArg->topkResults.resize(queryMatrix.rowNum * args.k);

#pragma omp parallel reduction(+ : comparisons)
            {

                RetrievalArguments arg(probeMatrix.colNum, queryMatrix, probeMatrix, LEMP_TA, false);

                arg.k = args.k;
                arg.heap.resize(args.k);
                arg.init(0);

                arg.state->initializeForNewBucket(invLists);


#pragma omp for schedule(dynamic,10)
                for (row_type i = 0; i < queryMatrix.rowNum; i++) {

                    double* query = queryMatrix.getMatrixRowPtr(i);
                    arg.queryId = i;

                    for (row_type j = 0; j < args.k; j++) {
                        double ip = queryMatrix.innerProduct(i, probeMatrix.getMatrixRowPtr(j));
                        arg.comparisons++;
                        arg.heap[j] = QueueElement(ip, j);
                    }

                    std::make_heap(arg.heap.begin(), arg.heap.end(), std::greater<QueueElement>()); //make the heap;


                    probeBucket.ptrRetriever->runTopK(query, probeBucket, &arg);

                    row_type p = i * args.k;
                    std::copy(arg.heap.begin(), arg.heap.end(), retrArg->topkResults.begin() + p);

                }
                comparisons += arg.comparisons;

            }
            topkResults = &(retrArg->topkResults);
            t.stop();

            std::cout << "Time for retrieval: " << t << std::endl;
            std::cout << "Comparisons: " << comparisons << std::endl;
            std::cout << "Preprocessing time: " << dataManipulationTime / 1E9 << std::endl;
            std::cout << "Total time: " << (dataManipulationTime / 1E9) + t.elapsedTime().seconds() << std::endl;


            logging << "\t" << args.k << "\t" << retrArg->comparisons << "\t" << getResultSetSize() << "\t";
            printTimes(t);

            if (args.resultsFile != "") {
                std::vector<MatItem> results;
                localToGlobalIds(*topkResults, args.k, results, queryMatrix);
                std::vector< std::vector<MatItem >* > resultsForWriting;
                resultsForWriting.push_back(&results);
                writeResults(resultsForWriting, args.resultsFile);
            }

            logging.close();
        }

        inline void printTimes(rg::Timer& tAll) {

            logging << (dataManipulationTime / 1E9) << "\t" << tAll.elapsedTime().seconds() << "\t" << tAll.elapsedTime().seconds()+(dataManipulationTime / 1E9) << "\n";
        }

        std::vector<QueueElement> * getResultsTopk() {
            return topkResults;
        }

        std::vector<MatItem >* getResultsTheta() {
            return thetaResults;
        }

        inline int getResultSetSize() {
            if (args.k > 0) {
                return topkResults->size();
            } else {
                return thetaResults->size();
            }
        }





    };




}
#endif /* TA_ALL_H_ */
