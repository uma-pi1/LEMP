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

namespace ta {

    class TANRA_all {
    public:

        rg::Timer t;
        VectorMatrix queryMatrix, probeMatrix;
        ProbeBucket probeBucket, probeBucketK;
        std::unique_ptr<Retriever>  ptrMain;

        std::vector< MatItem >* thetaResults; // for a specific query holds the itemIDs + the score
        std::vector<QueueElement> * topkResults;


        LEMPArg args;
        RetrievalArguments* retrArg;
        std::ofstream logging;

        double dataManipulationTime;

        inline TANRA_all(LEMPArg& args) : args(args), dataManipulationTime(0),  ptrMain(new tanraRetriever()) {

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

                std::cout << "Problem with opening log-file. No log-file will be created" << std::endl;
            }

            
             if (args.isTARR) {
                std::cout << "ALGO: TANRA_all with RR" << std::endl;
            } else {
                std::cout << "ALGO: TANRA_all with MaxPiQi" << std::endl;
            }
            
            
            std::cout << "Threads: " << args.threads << std::endl;

            retrArg = new RetrievalArguments(probeMatrix.colNum, &queryMatrix, &probeMatrix, LEMP_TANRA, false, args.isTARR);
            retrArg->k = args.k;
            retrArg->theta = args.theta;

            retrArg->init(probeMatrix.rowNum);

            t.start();
            probeBucket.init(probeMatrix, 0, probeMatrix.rowNum, args); // initialize
            probeBucket.bucketScanThreshold = args.theta / probeBucket.normL2.second;

           
            
            if (probeBucket.ptrIndexes[SL] == 0)
                probeBucket.ptrIndexes[SL] = new QueueElementLists();

            static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL])->initializeLists(probeMatrix, 0, probeMatrix.rowNum);

            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.getIndex(SL));


            t.stop();
            dataManipulationTime += t.elapsedTime().nanos();

            logging << "TANRA_all";

            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";


        };

        inline ~TANRA_all() {
        }

        inline void multiply() {
            args.comparisons = 0;



            std::cout << "Multiplication starts! theta: " << retrArg->theta << std::endl;

            t.start();
            QueueElementLists* invLists = static_cast<QueueElementLists*> (probeBucket.ptrIndexes[SL]);
            retrArg->tanraState->initializeForNewBucket(invLists);


            for (row_type i = 0; i < queryMatrix.rowNum; i++) {
                const double* query = queryMatrix.getMatrixRowPtr(i);
                retrArg->queryId = i;
                ptrMain->run(query, probeBucket, retrArg); 
            }



            thetaResults = &(retrArg->results);
            int totalSize = getResultSetSize();
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

#endif	/* TANRA_ALL_H */

