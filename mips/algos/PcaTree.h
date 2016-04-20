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

namespace mips {

    comp_type GenerateID() {
        static comp_type nNextID = 0;
        return nNextID++;
    }

    /* this files are for running the PCA-tree method AFTER the transformation of the input vectors. I.e. the datasets read need to be the transformed ones.
     */

    struct PCA_tree {
        double median;
        col_type col;
        bool isLeaf;
        comp_type leaf_counter;

        PCA_tree* leftChild;
        PCA_tree* rightChild;
        VectorMatrix* leafData;


        // if depth=2 I will have levels 2(root) 1, 0(leaves)

        PCA_tree(const VectorMatrix& matrix, int depth, col_type col, std::vector<row_type>& ids,
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
                    if (ids.size() % 2 == 0) {
                        median += values[pos - 1];
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

    class PcaTree : public Mip {
        PCA_tree* root;
        VectorMatrix probeFirstk;
        VectorMatrix probeMatrix;
        std::vector<RetrievalArguments> retrArg;
        std::vector<VectorMatrix*> matricesInLeaves;
        PcaTreeArguments args;
        bool isTransformed;

        arma::mat U;
        arma::vec mu;

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

            arma::mat A(rightMatrix.colNum + 1, rightMatrix.rowNum); //cols x rows

            mu.zeros(rightMatrix.colNum + 1);

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < rightMatrix.rowNum; i++) {
                double* dProbe = rightMatrix.getMatrixRowPtr(i);
                A(0, i) = 0;
                for (col_type j = 0; j < rightMatrix.colNum; ++j) {
                    A(j + 1, i) = dProbe[j];
                }
                double len = arma::norm(A.unsafe_col(i));
                A(0, i) = maxLen * maxLen - len * len;
                A(0, i) = (A(0, i) < 0) ? 0 : sqrt(A(0, i));
                mu += A.unsafe_col(i);
            }
            mu /= rightMatrix.rowNum;

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < rightMatrix.rowNum; i++) {
                A.unsafe_col(i) -= mu;
            }


            arma::vec s;
            arma::mat V;

            bool isOk = arma::svd_econ(U, s, V, A, "left");

            if (!isOk) {
                std::cout << "[ERROR] SVD failed " << std::endl;
                exit(1);
            }

            arma::inplace_trans(U);

            A = U*A;


            probeMatrix.rowNum = rightMatrix.rowNum;
            probeMatrix.colNum = rightMatrix.colNum + 1;

            probeMatrix.initializeBasics(probeMatrix.colNum, probeMatrix.rowNum, false);

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < probeMatrix.rowNum; i++) {
                double* dProbe = probeMatrix.getMatrixRowPtr(i);
                for (col_type j = 0; j < probeMatrix.colNum; ++j) {
                    dProbe[j] = A(j, i);
                }
                probeMatrix.setLengthInData(i, 1); // set to 1 by the transformation
            }
        }

        inline void transformQueryMatrix(const VectorMatrix& leftMatrix, VectorMatrix& queryMatrix) {
            // transform queryMatrix (transform ||q|| to 1 and q = [0;q])

            arma::mat A(leftMatrix.colNum + 1, leftMatrix.rowNum); //cols x rows

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < leftMatrix.rowNum; i++) {
                double* dProbe = leftMatrix.getMatrixRowPtr(i);
                A(0, i) = 0;
                for (col_type j = 0; j < leftMatrix.colNum; ++j) {
                    A(j + 1, i) = dProbe[j];
                }

                A.unsafe_col(i) -= mu;
            }

            A = U * A;

            queryMatrix.rowNum = leftMatrix.rowNum;
            queryMatrix.colNum = leftMatrix.colNum + 1;
            queryMatrix.initializeBasics(queryMatrix.colNum, queryMatrix.rowNum, false);

#pragma omp parallel for schedule(static,1000)
            for (row_type i = 0; i < queryMatrix.rowNum; i++) {
                double* dQuery = queryMatrix.getMatrixRowPtr(i);
                for (col_type j = 0; j < queryMatrix.colNum; ++j) {
                    dQuery[j] = A(j, i);
                }
                queryMatrix.setLengthInData(i, 1); // ||q|| = 1
            }
        }

        inline void verify(const double* query, comp_type bucket, row_type tid, double& minScore) {

            for (row_type j = 0; j < matricesInLeaves[bucket]->rowNum; j++) {
                row_type probeId = matricesInLeaves[bucket]->lengthInfo[j].id;
                double euclideanDistance = matricesInLeaves[bucket]->L2Distance2(j, query);
                //                std::cout<<"checking: "<<j<<" probe id: "<<probeId<<" result: "<<euclideanDistance<<std::endl;
                retrArg[tid].comparisons++;

                if (euclideanDistance < minScore) {
                    std::pop_heap(retrArg[tid].heap.begin(), retrArg[tid].heap.end(), std::less<QueueElement>());
                    retrArg[tid].heap.pop_back();
                    retrArg[tid].heap.push_back(QueueElement(euclideanDistance, probeId));
                    std::push_heap(retrArg[tid].heap.begin(), retrArg[tid].heap.end(), std::less<QueueElement>());
                    minScore = retrArg[tid].heap.front().data;
                }

            }

        }

        inline void printAlgoName(const VectorMatrix& leftMatrix) {
            logging << "PCA_TREE" << "\t" << args.threads << "\t d(" << args.depth << ")\t";
            std::cout << "[ALGORITHM] PCA_TREE with " << args.threads << " thread(s) and depth " << args.depth << std::endl;

            logging << "P(" << probeMatrix.rowNum << "x" << (0 + probeMatrix.colNum) << ")\t";
            logging << "Q^T(" << leftMatrix.rowNum << "x" << (0 + leftMatrix.colNum) << ")\t";
        }

        inline void initializeInternal(std::vector<VectorMatrix>& queryMatrices, const VectorMatrix& leftMatrix) {

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
                transformQueryMatrix(leftMatrix, queryMatrix);
                splitMatrices(queryMatrix, queryMatrices);
            } else {
                splitMatrices(leftMatrix, queryMatrices);
            }
            timer.stop();
            dataPreprocessingTimeLeft += timer.elapsedTime().nanos();


            retrArg.resize(myNumThreads);

            for (row_type i = 0; i < retrArg.size(); i++) {
                retrArg[i].initializeBasics(queryMatrices[i], probeMatrix, LEMP_L, args.theta, args.k, args.threads, 0, 0, 0, 0);
                retrArg[i].clear();
            }
        }

    public:

        inline PcaTree(InputArguments& input, int depth, bool isTransformed) : root(nullptr), isTransformed(isTransformed) {
            args.copyInputArguments(input);
            args.depth = depth;

            // now do the logging
            logging.open(args.logFile.c_str(), std::ios_base::app);

            if (!logging.is_open()) {
                std::cout << "[WARNING] No log will be created!" << std::endl;
            } else {
                std::cout << "[INFO] Logging in " << args.logFile << std::endl;
            }

            omp_set_num_threads(args.threads);
        }

        inline ~PcaTree() {
            logging.close();
        }

        inline void initialize(VectorMatrix& rightMatrix) {
            std::cout << "[INIT] ProbeMatrix contains " << rightMatrix.rowNum << " vectors with dimensionality " << (0 + rightMatrix.colNum) << std::endl;


            if (!isTransformed) {
                std::cout << "[INIT] ProbeMatrix will be transformed" << std::endl;
                timer.start();
                transformProbeMatrix(rightMatrix);
                timer.stop();
                dataPreprocessingTimeRight += timer.elapsedTime().nanos();
            } else {
                probeMatrix = rightMatrix;
            }

            //create the tree
            std::vector<row_type> ids;
            timer.start();
            root = new PCA_tree(probeMatrix, args.depth, 0, ids, args.k, matricesInLeaves);

            ids.clear();
            for (int i = 0; i < args.k; i++) {
                ids.push_back(i);
            }
            probeFirstk.addVectors(probeMatrix, ids);

            timer.stop();
            dataPreprocessingTimeRight = timer.elapsedTime().nanos();
        }

        inline void runTopK(VectorMatrix& leftMatrix, Results& results) {
            printAlgoName(leftMatrix);

            std::vector<VectorMatrix> queryMatrices;

            initializeInternal(queryMatrices, leftMatrix);
            results.resultsVector.resize(args.threads);

            std::cout << "[RETRIEVAL] Retrieval (k = " << args.k << ") starts ..." << std::endl;
            logging << "k(" << args.k << ")\t";

            timer.start();

            for (row_type i = 0; i < retrArg.size(); i++)
                retrArg[i].allocTopkResults();

            comp_type comparisons = 0;
#pragma omp parallel reduction(+ : comparisons)
            {
                row_type tid = omp_get_thread_num();
                double minScore = 0;

                for (row_type i = 0; i < queryMatrices[tid].rowNum; i++) {
                    const double* query = queryMatrices[tid].getMatrixRowPtr(i);

                    for (row_type j = 0; j < args.k; j++) {
                        retrArg[tid].comparisons++;
                        double euclideanDistance = probeFirstk.L2Distance2(j, query);
                        //                         std::cout<<"checking: "<<j<<" probe id: "<<0<<" result: "<<euclideanDistance<<std::endl;
                        retrArg[tid].heap[j] = QueueElement(euclideanDistance, j);
                    }

                    std::make_heap(retrArg[tid].heap.begin(), retrArg[tid].heap.end(), std::less<QueueElement>()); //make the heap;
                    minScore = retrArg[tid].heap.front().data;

                    // find bucket of query
                    comp_type queryBucket = root->findBucketForQuery(query);
                    verify(query, queryBucket, tid, minScore);

                    // now scan also buckets in hamming distance=1
                    for (comp_type b = 0; b < matricesInLeaves.size(); b++) {

                        if (queryBucket == b || matricesInLeaves[b]->rowNum == 0)
                            continue;

                        comp_type res = b ^ queryBucket;
                        int hammingDistance = __builtin_popcountll(res);

                        if (hammingDistance == 1) { // these guys are candidates
                            verify(query, b, tid, minScore);
                        }
                    }// for each bucket

                    // write the results back
                    retrArg[tid].writeHeapToTopk(i);

                    //                    std::cout << retrArg[tid].heap.front().data << " " << retrArg[tid].heap.front().id << std::endl;


                }// for each query

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

#endif	/* PCATREE_H */

