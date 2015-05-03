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



namespace ta {

    class Naive {
        VectorMatrix queryMatrix;
        VectorMatrix probeMatrix;
        std::vector< MatItem >* thetaResults; // for a specific query holds the itemIDs + the score
        std::vector<QueueElement> * topkResults;

        LEMPArg& args;
        RetrievalArguments* retrArg;
        std::ofstream logging;

    public:

        inline Naive(LEMPArg& args) : args(args) {
            std::vector<std::vector<double> > leftMatrix, rightMatrix;


            if (args.querySideLeft) {
                leftMatrix = ta::readMatrix(args.usersFile, true);
                rightMatrix = ta::readMatrix(args.itemsFile, false);
            } else {
                leftMatrix = ta::readMatrix(args.itemsFile, false);
                rightMatrix = ta::readMatrix(args.usersFile, true);
            }

            queryMatrix.init(leftMatrix);
            probeMatrix.init(rightMatrix);


            //            queryMatrix.initForTopKPerUser(leftMatrix); // just to check
            //            probeMatrix.init(rightMatrix, true, true); // just to check
//                        probeMatrix.init(rightMatrix, true, false); // just to check


            retrArg = new RetrievalArguments(probeMatrix.colNum, queryMatrix, probeMatrix, args.method);
            retrArg->theta = args.theta;
            retrArg->k = args.k;

            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (logging.is_open()) {

                std::cout << "Writing output to " << args.logFile << std::endl;
            } else {

                std::cout << "Problem with opening log-file. No log-file will be created" << std::endl;
            }

            logging << "Naive";
            logging << "\t \"" << args.usersFile << "\"" << "\t" << args.threads << "\t";


        }

        inline ~Naive() {
            if (!retrArg)
                delete retrArg;
        }

               void multiply() {
                   rg::Timer t;
                   std::pair<bool, double> p;
                   retrArg->comparisons = 0;
                   LengthRetriever plainRetriever;
        
                   t.start();
                   for (row_type i = 0; i < queryMatrix.rowNum; i++) {
                       double* query = queryMatrix.getMatrixRowPtr(i);
                       retrArg->queryId = i;
                       plainRetriever.naive(query, 0, probeMatrix.rowNum, retrArg);
                   }
                   thetaResults = &(retrArg->results);
                   t.stop();
        
                   std::cout << "TIME for 2 sided: " << t << std::endl;
                   std::cout << "Result Size: " << getResultSetSize() <<" "<< (retrArg->results).size()<< std::endl;
                   std::cout << "Number of comparisons: " << retrArg->comparisons << std::endl;
		   
		   if(args.resultsFile != ""){
			std::vector< std::vector<MatItem >* > resultsForWriting;
			resultsForWriting.push_back(thetaResults);
			writeResults(resultsForWriting, args.resultsFile);
		   }
		   
		   
               }

//         void multiply() {// for parallel execution
//             rg::Timer t;
//             std::pair<bool, double> p;
//             retrArg->comparisons = 0;
//             LengthRetriever plainRetriever;
//             comp_type comparisons = 0;
//             t.start();
// 
//             omp_set_num_threads(args.threads);
// 
// #pragma omp parallel reduction(+ : comparisons)
//             {
// 
//                 RetrievalArguments arg(probeMatrix.colNum, queryMatrix, probeMatrix, args.method);
//                 arg.theta = args.theta;
// 
// #pragma omp for schedule(dynamic,10)  
//                 for (row_type i = 0; i < queryMatrix.rowNum; i++) {
//                     double* query = queryMatrix.getMatrixRowPtr(i);
//                     arg.queryId = i;
//                     plainRetriever.naive(query, 0, probeMatrix.rowNum, &arg);
//                 }
// 
//                 comparisons += arg.comparisons;
// 
//             }
// 
//             thetaResults = &(retrArg->results);
//             t.stop();
// 
//             std::cout << "TIME: " << t << std::endl;
// 	    std::cout << "Result Size: " << getResultSetSize() << std::endl;
//             std::cout << "Number of comparisons: " << comparisons << std::endl;
// 	    
// 	    if(args.resultsFile != ""){
// 		std::vector< std::vector<MatItem >* > resultsForWriting;
// 		resultsForWriting.push_back(thetaResults);
// 		writeResults(resultsForWriting, args.resultsFile);
// 	    }
// 	    
// 	    
//         }



//                 void topKperUser() {
//                     rg::Timer timer;
//                     double minScore = 0;
//                     retrArg->comparisons = 0;
//                     LengthRetriever plainRetriever;
//         
//                     timer.start();
//         
//                     retrArg->topkResults.resize(queryMatrix.rowNum * args.k);
//                     retrArg->heap.resize(args.k);
//         
//                     for (row_type i = 0; i < queryMatrix.rowNum; i++) {
//         
//                         double* query = queryMatrix.getMatrixRowPtr(i);
//                         retrArg->queryId = i;
//         
//                         for (row_type j = 0; j < args.k; j++) {
//                             double ip = queryMatrix.innerProduct(i, probeMatrix.getMatrixRowPtr(j));
//                             retrArg->comparisons++;
//                             retrArg->heap[j] = QueueElement(ip, j);                    
//                         }
//         
//                         std::make_heap(retrArg->heap.begin(), retrArg->heap.end(), std::greater<QueueElement>()); //make the heap;
//         
//                         plainRetriever.naiveTopk(query, retrArg->k, probeMatrix.rowNum, retrArg);
//         
//                         retrArg->writeHeapToTopk(i);
//         
//                     }
//         
//         
//                     topkResults = &(retrArg->topkResults);
//                     timer.stop();
//         
//                     std::cout << "TIME for 2 sided: " << timer << std::endl;
//                     std::cout << "Result Size: " << getResultSetSize() << std::endl;
//                     std::cout << "Number of comparisons: " << retrArg->comparisons << std::endl;
//         
// 		    
// 		    if(args.resultsFile != ""){
// 			  std::vector<MatItem> results;
// 			  localToGlobalIds(*topkResults, args.k, results, queryMatrix);
// 			  std::vector< std::vector<MatItem >* > resultsForWriting;
// 			  resultsForWriting.push_back(&results);
// 			  writeResults(resultsForWriting, args.resultsFile);
// 		    }
//         
//                 }


        // this is for the parallel experiment

       void topKperUser() {
           rg::Timer timer;
           double minScore = 0;
           retrArg->comparisons = 0;


           timer.start();

           retrArg->topkResults.resize(queryMatrix.rowNum * args.k);

           comp_type comparisons = 0;
           omp_set_num_threads(args.threads);

#pragma omp parallel reduction(+ : comparisons)
           {

               RetrievalArguments arg(probeMatrix.colNum, queryMatrix, probeMatrix, args.method);
               arg.theta = args.theta;
               arg.k = args.k;
               arg.heap.resize(args.k);

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


                   double minScore = arg.heap.front().data;

                   for (row_type j = args.k; j < probeMatrix.rowNum; j++) {
                       arg.comparisons++;
                       double ip = arg.probeMatrix->innerProduct(j, query);

                       if (ip > minScore) {
                           std::pop_heap(arg.heap.begin(), arg.heap.end(), std::greater<QueueElement>());
                           arg.heap.pop_back();
                           arg.heap.push_back(QueueElement(ip, arg.probeMatrix->getId(j)));
                           std::push_heap(arg.heap.begin(), arg.heap.end(), std::greater<QueueElement>());
                           minScore = arg.heap.front().data;
                       }
                   }

                   row_type p = i * args.k;
                   std::copy(arg.heap.begin(), arg.heap.end(), retrArg->topkResults.begin() + p);
               }

               comparisons += arg.comparisons;
           }




           topkResults = &(retrArg->topkResults);
           timer.stop();

           std::cout << "TIME for 2 sided: " << timer << std::endl;
           std::cout << "Result Size: " << getResultSetSize() << std::endl;
           std::cout << "Number of comparisons: " << comparisons << std::endl;

	   if(args.resultsFile != ""){
		std::vector<MatItem> results;
		localToGlobalIds(*topkResults, args.k, results, queryMatrix);
		std::vector< std::vector<MatItem >* > resultsForWriting;
		resultsForWriting.push_back(&results);
		writeResults(resultsForWriting, args.resultsFile);
	   }
	   
	   
           logging << "\t" << args.k << "\t" << comparisons << "\t" << getResultSetSize() << "\t" << timer << "\n";

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



#endif /* NAIVE_H_ */
