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
 * File:   PCATree.h
 * Author: chteflio
 *
 * Created on February 9, 2015, 10:30 AM
 */

#ifndef PCATREE_H
#define	PCATREE_H

namespace ta {

    comp_type GenerateID() {
        static comp_type nNextID = 0;
        return nNextID++;
    }

    struct PCA_tree {
        double median;
        col_type col;
        bool isLeaf;
        comp_type leaf_counter;

        PCA_tree* leftChild;
        PCA_tree* rightChild;
        VectorMatrix* leafData;


        // if depth=2 I will have levels 2(root) 1, 0(leaves)

        PCA_tree(const VectorMatrix& matrix,  int depth, col_type col, std::vector<row_type>& ids, 
        int k, std::vector<VectorMatrix*>& matricesInLeaves) :
        col(col), isLeaf(false), leafData(0) {



            std::vector<row_type> idsLeft, idsRight;

            if (ids.size() == 0 && col == 0) { // this is the root node
                ids.reserve(matrix.rowNum);
                for (int i = k; i < matrix.rowNum; i++)
                    ids.push_back(i);

		matricesInLeaves.resize(pow(2, depth));
            }

            if (depth == 0) {// this is a leaf            
                isLeaf = true;

                leaf_counter = GenerateID();               
                leafData = new VectorMatrix();
                leafData->addVectors(matrix, ids);
                matricesInLeaves[leaf_counter] = leafData;

            } else { // internal node

                if (ids.size() > 0) {

                    std::vector<double> values;
                    values.reserve(ids.size());

                    for (int i = 0; i < ids.size(); i++) {
                        double val = matrix.getMatrixRowPtr(ids[i])[col];
                        values.push_back(val);
                    }

                    std::sort(values.begin(), values.end(), std::greater<double>());


                    row_type pos = ids.size() / 2;
		    

                    median = values[pos];
		    if(ids.size()%2 == 0){
		       median += values[pos-1];
		       median /= 2;
		    }



                    for (int i = 0; i < ids.size(); i++) {

                        double val = matrix.getMatrixRowPtr(ids[i])[col];

                        if (val <= median) {
                            idsLeft.push_back(ids[i]);
                        } else {
                            idsRight.push_back(ids[i]);
                        }

                    }

                }


                leftChild = new PCA_tree(matrix, depth - 1, col + 1, idsLeft, k, matricesInLeaves);
                rightChild = new PCA_tree(matrix, depth - 1, col + 1, idsRight, k, matricesInLeaves);
            }


        }

        inline comp_type findBucketForQuery(const double* query) {

            if (isLeaf)
                return leaf_counter;

            if (query[col] <= median) {
                leftChild->findBucketForQuery(query);
            } else {
                rightChild->findBucketForQuery(query);
            }

        }



    };

    class SearchWithPCATree {
        PCA_tree* root;

        rg::Timer t;
        VectorMatrix probeFirstk, queryMatrix;

        std::vector<QueueElement> topkResults;
        std::vector<MatItem > results;


        std::vector<VectorMatrix*> matricesInLeaves;


        LEMPArg args;
        std::ofstream logging;
        std::vector<QueueElement> heap;
        comp_type comparisons;
        double minScore;

        double preprocessingTime, retrievalTime;

        inline void verify(const double* query, comp_type bucket) {
          
            for (row_type j = 0; j < matricesInLeaves[bucket]->rowNum; j++) {
                row_type probeId = matricesInLeaves[bucket]->lengthInfo[j].id;
                double euclideanDistance = matricesInLeaves[bucket]->L2Distance2(j, query);
                comparisons++;

                if (euclideanDistance < minScore) {
                    std::pop_heap(heap.begin(), heap.end(), std::less<QueueElement>());
                    heap.pop_back();
                    heap.push_back(QueueElement(euclideanDistance, probeId));
                    std::push_heap(heap.begin(), heap.end(), std::less<QueueElement>());
                    minScore = heap.front().data;
                }
                
            }
            
        }

    public:

        inline SearchWithPCATree(LEMPArg& args) : args(args), preprocessingTime(0), comparisons(0), minScore(0), retrievalTime(0), root(0) {
	    VectorMatrix probeMatrix;
            if (args.querySideLeft) {
                queryMatrix.readFromFile(args.usersFile, true);
                probeMatrix.readFromFile(args.itemsFile, false);
            } else {
                queryMatrix.readFromFile(args.itemsFile, false);
                probeMatrix.readFromFile(args.usersFile, true);
            }



            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (logging.is_open()) {

                std::cout << "Writing output to " << args.logFile << std::endl;
            } else {

                std::cout << "Problem with opening output file" << std::endl;
            }

            logging << "PCA_Tree" << "\t" << args.depth;
            std::cout << "PCA Tree with Depth: " << args.depth << std::endl;

            //create the tree
            std::vector<row_type> ids;
            t.start();
            root = new PCA_tree(probeMatrix,  args.depth, 0, ids, args.k, matricesInLeaves);
	    
	    ids.clear();
	    for(int i=0; i<args.k; i++){
	      ids.push_back(i);
	    }
	    probeFirstk.addVectors(probeMatrix, ids);
	    
            t.stop();
            preprocessingTime = t.elapsedTime().nanos();
            std::cout << "Time for tree construction: " << (preprocessingTime / 1E9) << std::endl;
        }

        ~SearchWithPCATree() {


        }

        inline void topKperUser() {

            std::cout << "Retrieval starts. k = " << args.k << std::endl;
            t.start();
            topkResults.resize(queryMatrix.rowNum * args.k);
            heap.resize(args.k);

            for (row_type i = 0; i < queryMatrix.rowNum; i++) {

                const double* query = queryMatrix.getMatrixRowPtr(i);

                for (row_type j = 0; j < args.k; j++) {
                    comparisons++;
                    double euclideanDistance = probeFirstk.L2Distance2(j, query);
                    heap[j] = QueueElement(euclideanDistance, j);
                }


                std::make_heap(heap.begin(), heap.end(), std::less<QueueElement>()); //make the heap;
                minScore = heap.front().data;

                // find bucket of query
                comp_type queryBucket = root->findBucketForQuery(query);
                verify(query, queryBucket);

                // now scan also buckets in hamming distance=1
                for (comp_type b = 0; b < matricesInLeaves.size(); b++) {

		    
		

                    if (queryBucket == b ||  matricesInLeaves[b]->rowNum == 0)
                        continue;

                    comp_type res = b ^ queryBucket;
                    int hammingDistance = __builtin_popcountll(res);

                    if (hammingDistance == 1) { // these guys are candidates
                        verify(query, b);
                    }
                }// for each bucket

                // write the results back
                row_type p = i * args.k;
                std::copy(heap.begin(), heap.end(), topkResults.begin() + p);

            }// for each query
            t.stop();
            retrievalTime = t.elapsedTime().nanos();


            localToGlobalIds(topkResults, args.k, results, queryMatrix);


            std::cout << "Time for Retrieval: " << (retrievalTime / 1E9) << std::endl;
            std::cout << "Size of result: " << results.size() << std::endl;
            std::cout << "Comparisons: " << comparisons << std::endl;

            logging << "\t" << args.k << "\t" << comparisons << "\t" << results.size() << "\t";
            double timeAll = (preprocessingTime / 1E9) + (retrievalTime / 1E9);
            std::cout << "Total Time: " << timeAll << std::endl;
            logging << preprocessingTime / 1E9 << "\t" << 0 << "\t" << retrievalTime / 1E9 << "\t" << timeAll << "\n";
            std::cout << preprocessingTime / 1E9 << "\t" << 0 << "\t" << retrievalTime / 1E9 << "\t" << timeAll << std::endl;



            std::vector< std::vector<MatItem >* > resultsForWriting;
            resultsForWriting.push_back(&results);


            if (args.resultsFile != "") {
                writeResults(resultsForWriting, args.resultsFile);
            }


            logging.close();
        }





    };



}

#endif	/* PCATREE_H */

