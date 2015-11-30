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
 * singleTree.h
 *
 *  Created on: Jul 21, 2014
 *      Author: chteflio
 */

#ifndef SINGLETREE_H_
#define SINGLETREE_H_



#include <taLib/my_mlpack/core.hpp>
#include <taLib/my_mlpack/core/metrics/ip_metric.hpp>
#include <taLib/my_mlpack/fastmks/fastmks_stat.hpp>
#include <taLib/my_mlpack/core/tree/cover_tree.hpp>
#include <taLib/my_mlpack/fastmks/fastmks.hpp>
#include <taLib/my_mlpack/fastmks/fastmks_rules.hpp>
#include <queue>
#include <taLib/my_mlpack/core/tree/cover_tree/single_tree_traverser.hpp>


using namespace ta::fastmks;
using namespace ta::kernel;
using namespace ta::metric;

namespace ta {

    class SingleTree : public Retriever {
        FastMKS<TreeType>* fastmks;

    public:

        inline SingleTree() {
            fastmks = new FastMKS<TreeType>(true, false);
        }

        inline ~SingleTree() {
            if (!fastmks)
                delete fastmks;
        }

        // arg should point to the correct tree

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            fastmks->SearchForTheta(arg->theta, arg->tree->tree, arg->probeMatrix, arg->queryMatrix, arg->results, arg->queryPos, arg->comparisons, arg->threads);
        }

        inline virtual void run(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) const {
            fastmks->Search(arg->k, arg->tree->tree, arg->probeMatrix, arg->queryMatrix, arg->heap, arg->queryPos, arg->comparisons, arg->threads);
        }

        // this is to be used only by the 1st bucket in Row-Top-k. It just QueryBatchializes the top-k elements

        inline virtual void runTopK(QueryBatch& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) const {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            row_type sampleSize = probeBucket.xValues->size();
            if (sampleSize > 0) {
                sampleTimes.resize(sampleSize);

                TreeIndex * index = static_cast<TreeIndex*> (probeBucket.ptrIndexes[TREE]);

                for (row_type i = 0; i < sampleSize; ++i) {

                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;

                    retrArg[t].tunerTimer.start();
                    fastmks->SearchForTheta(retrArg[t].theta, index->tree, retrArg[t].probeMatrix, retrArg[t].queryMatrix, retrArg[t].results, ind, retrArg[t].comparisons, 1);
                    retrArg[t].tunerTimer.stop();
                    sampleTimes[i] = retrArg[t].tunerTimer.elapsedTime().nanos();
                }

            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {
            row_type sampleSize = (probeBucket.xValues != nullptr ? probeBucket.xValues->size() : 0);
            sampleTimes.resize(sampleSize);
            TreeIndex * index = static_cast<TreeIndex*> (probeBucket.ptrIndexes[TREE]);

            for (row_type i = 0; i < sampleSize; ++i) {
                int t = probeBucket.xValues->at(i).i;
                int ind = probeBucket.xValues->at(i).j;

                const std::vector<QueueElement>& prevResults = prevBucket.sampleThetas->at(t).at(ind).results;

                retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                retrArg[0].tunerTimer.start();
                fastmks->Search(retrArg[0].k, index->tree, retrArg[0].probeMatrix, retrArg[t].queryMatrix, retrArg[0].heap, ind, retrArg[0].comparisons, 1);
                retrArg[0].tunerTimer.stop();
                sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
            }

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            TreeIndex * index = static_cast<TreeIndex*> (probeBucket.ptrIndexes[TREE]);


            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.isWorkDone())
                    continue;

                if (!index->isInitialized()) {
                    index->initializeTree(*(arg->probeMatrix), arg->threads, probeBucket.startPos, probeBucket.endPos);
                }

                row_type user = queryBatch.startPos;
                int start = queryBatch.startPos * arg->k;
                int end = queryBatch.endPos * arg->k;
                for (row_type i = start; i < end; i += arg->k) {

                    if (queryBatch.isQueryInactive(user)) {
                        user++;
                        continue;
                    }

                    double minScore = arg->topkResults[i].data;

                    if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                        queryBatch.inactivateQuery(user);
                        user++;
                        continue;
                    }

                    arg->moveTopkToHeap(i);
                    arg->queryId = arg->queryMatrix->getId(user);
                    fastmks->Search(arg->k, index->tree, arg->probeMatrix, arg->queryMatrix, arg->heap, user, arg->comparisons, arg->threads);

                    arg->writeHeapToTopk(user);
                    user++;
                }
            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            TreeIndex * index = static_cast<TreeIndex*> (probeBucket.ptrIndexes[TREE]);

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }

                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;
                    arg->queryId = arg->queryMatrix->getId(i);
                    fastmks->SearchForTheta(arg->theta, index->tree, arg->probeMatrix, arg->queryMatrix, arg->results, i, arg->comparisons, arg->threads);
                }
            }
        }

        inline virtual void cleanupAfterTuning() {

        }


    };




}

#endif /* SINGLETREE_H_ */
