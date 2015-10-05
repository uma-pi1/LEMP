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
 * File:   lshRetriever.h
 * Author: chteflio
 *
 * Created on June 26, 2015, 2:32 PM
 */

#ifndef LSHRETRIEVER_H
#define	LSHRETRIEVER_H

namespace ta {

    class LshRetriever : public Retriever {
        LengthRetriever plain;
        row_type t_b_indx;

        inline void processIndexesTopk(const double * query, row_type queryId,
                LshIndex* index, LshIndex* queryIndex, int activeBlocks,
                ProbeBucket& probeBucket, RetrievalArguments* arg) {

            row_type numCandidatesToVerify = 0;

#ifdef TIME_IT 
            arg->t.start();
#endif
                index->lshBins->getCandidates(queryIndex->cosSketches->sketches, queryId, arg->candidatesToVerify, numCandidatesToVerify,
                        arg->done, activeBlocks, probeBucket.startPos);

#ifdef TIME_IT
            arg->t.stop();
            arg->scanTime += arg->t.elapsedTime().nanos();
            arg->t.start();
#endif

            verifyCandidatesTopK_lengthTest(query, numCandidatesToVerify, arg);
#ifdef TIME_IT
            arg->t.stop();
            arg->ipTime += arg->t.elapsedTime().nanos();
#endif
        }

        inline void processIndexes(const double * query, row_type queryId,
                LshIndex* index, LshIndex* queryIndex, int activeBlocks,
                ProbeBucket& probeBucket, RetrievalArguments* arg) {


            row_type numCandidatesToVerify = 0;
#ifdef TIME_IT 
            arg->t.start();
#endif

                index->lshBins->getCandidates(queryIndex->cosSketches->sketches, queryId, arg->candidatesToVerify, numCandidatesToVerify,
                        arg->done, activeBlocks, probeBucket.startPos);
  
#ifdef TIME_IT
            arg->t.stop();
            arg->scanTime += arg->t.elapsedTime().nanos();
            arg->t.start();
#endif
            verifyCandidates_noLengthTest(query, numCandidatesToVerify, arg);

#ifdef TIME_IT
            arg->t.stop();
            arg->ipTime += arg->t.elapsedTime().nanos();
#endif

        }


    public:
        double t_b_LSH;

        LshRetriever() : Retriever(LEMP_LSH), t_b_LSH(1) {
        }

        ~LshRetriever() {
        }

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        // todo see that one from skratch

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {

            if (xValues->size() > 0) {
                // measure length
                plain.xValues = xValues;

                // access plain properly

                if (retrArg[0].competitorMethod != NULL) {
                    plain.sampleTimes = *(retrArg[0].competitorMethod);
                } else {
                    plain.tune(probeBucket, retrArg);
                }

                double avgLengthTime = plain.sampleTotalTime / xValues->size();


                // measure LSH
                sampleTimes.resize(xValues->size());
                LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

                // create lsh signatures for sample queries
                LshIndex queryIndex;
                queryIndex.initializeLists(xValues->size());


                double avgBlockTime = 0;
                double min_t_b = -1;
                int breakpoint = xValues->size() - 1;

                //                retrArg[0].tunerTimer.start();
                //                index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, LSH_SIGNATURES,
                //                        retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, retrArg[0].rig);
                //                
                //                
                //                retrArg[0].tunerTimer.stop();
                //                double blockOverhead = retrArg[0].tunerTimer.elapsedTime().nanos() / (LSH_SIGNATURES * probeBucket.activeQueries);
                //                double extraOverhead;

           

                for (int i = xValues->size() - 1; i >= 0; i--) {
                    int t = xValues->at(i).i;
                    int ind = xValues->at(i).j;

                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    double localTheta = probeBucket.bucketScanThreshold / query[-1];
                    int activeBlocks = retrArg[0].findActiveBlocks(localTheta);


                    if (activeBlocks > LSH_SIGNATURES || activeBlocks <= 0) {
                        sampleTimes[i] = plain.sampleTimes[i];

                    } else {
                        if (activeBlocks > index->initializedSketchesForIndex && i < xValues->size() - 1) {

                            double expectedTime = avgBlockTime * activeBlocks / (index->initializedSketchesForIndex * (xValues->size() - i));

                            if (expectedTime > avgLengthTime) {
                                min_t_b = localTheta;
                                breakpoint = i;
                                break;
                            }
                        }


                        index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, retrArg[0].rig);

                        retrArg[0].tunerTimer.start();
                        queryIndex.checkAndReallocateSingle(retrArg[t].queryMatrix, ind, i, activeBlocks, retrArg[t].sums, retrArg[t].rig);
                        processIndexes(query, i, index, &queryIndex, activeBlocks, probeBucket, &retrArg[0]);
                        retrArg[0].tunerTimer.stop();
                        sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                        avgBlockTime += sampleTimes[i] / activeBlocks;
                    }
                }

                if (breakpoint == xValues->size() - 1) {
                    // compare
                    findCutOffPoint(sampleTimes, plain.sampleTimes, sampleTotalTime, t_b_indx);

                    if (plain.sampleTotalTime < sampleTotalTime) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        double value = (t_b_indx == 0 ? -1 : xValues->at(t_b_indx).result);
                        probeBucket.setAfterTuning(1, value);
                    }

                } else {

                    if (breakpoint == xValues->size() - 2) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        probeBucket.setAfterTuning(1, min_t_b);
                    }
                }
            }

        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            row_type b = retrArg[0].bucketInd;


            if (xValues->size() > 0) {
                plain.xValues = xValues;

                for (row_type i = 0; i < xValues->size(); i++) {
                    int t = xValues->at(i).i;
                    int ind = xValues->at(i).j;

                    plain.sampleTimes.push_back(retrArg[0].globalData[b][t][ind].lengthTime);
                    plain.sampleTotalTime += retrArg[0].globalData[b][t][ind].lengthTime;
                }

                double avgLengthTime = plain.sampleTotalTime / xValues->size();

                // measure LSH
                sampleTimes.resize(xValues->size());
                LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

                if (!index->initialized) {
#if defined(TIME_IT)
                    retrArg[0].t.start();
#endif           
                    index->initializeLists(*(retrArg[0].probeMatrix), true, probeBucket.startPos, probeBucket.endPos);

#if defined(TIME_IT)
                    retrArg[0].t.stop();
                    retrArg[0].initializeListsTime += retrArg[0].t.elapsedTime().nanos();
#endif
                }
                // create lsh signatures for sample queries
                LshIndex queryIndex;
                queryIndex.initializeLists(xValues->size());



                double avgBlockTime = 0;
                double min_t_b = -1;
                int breakpoint = xValues->size() - 1;
              

                for (int i = xValues->size() - 1; i >= 0; i--) {

                    int t = xValues->at(i).i;
                    int ind = xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    std::vector<QueueElement>& prevResults = retrArg[0].globalData[b - 1][t][ind].results; //just reading


                    retrArg[0].heap.assign(prevResults.begin(), prevResults.end());
                    std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                    double minScore = retrArg[0].heap.front().data;
                    double localTheta = minScore * (minScore > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);
                    int activeBlocks = retrArg[0].findActiveBlocks(localTheta);

                    if (activeBlocks > LSH_SIGNATURES || activeBlocks <= 0) {
                        sampleTimes[i] = plain.sampleTimes[i];
                    } else {

                        if (activeBlocks > index->initializedSketchesForIndex && i < xValues->size() - 1) {

                            double expectedTime = avgBlockTime * activeBlocks / (index->initializedSketchesForIndex * (xValues->size() - i));

                            if (expectedTime > avgLengthTime) {
                                min_t_b = localTheta;
                                breakpoint = i;
                                break;
                            }
                        }
                        index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, retrArg[0].rig);


                        retrArg[0].tunerTimer.start();
                        queryIndex.checkAndReallocateSingle(retrArg[t].queryMatrix, ind, i, activeBlocks, retrArg[0].sums, retrArg[0].rig);

                        processIndexesTopk(query, i, index, &queryIndex, activeBlocks, probeBucket, &retrArg[0]);

                        retrArg[0].tunerTimer.stop();
                        sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                        avgBlockTime += sampleTimes[i] / activeBlocks;
                    }
                }


                if (breakpoint == xValues->size() - 1) { // no break at all
                    // compare
                    findCutOffPoint(sampleTimes, plain.sampleTimes, sampleTotalTime, t_b_indx);

                    if (plain.sampleTotalTime < sampleTotalTime) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        double value = (t_b_indx == 0 ? -1 : xValues->at(t_b_indx).result);
                        probeBucket.setAfterTuning(1, value);
                    }

                } else {

                    if (breakpoint == xValues->size() - 2) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        probeBucket.setAfterTuning(1, min_t_b);
                    }
                }
            }
     

        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {

            if (probeBucket.t_b == 1) {
                plain.runTopK(probeBucket, arg);
            } else { // do it per query

                LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

                if (!index->initialized) {
#ifdef TIME_IT
                    arg->t.start();
#endif           
                    index->initializeLists(*(arg->probeMatrix), true, probeBucket.startPos, probeBucket.endPos);

#ifdef TIME_IT
                    arg->t.stop();
                    arg->initializeListsTime += arg->t.elapsedTime().nanos();
#endif
                }

                for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                    if (arg->queryBatches[q].inactiveCounter == arg->queryBatches[q].rowNum)
                        continue;

                    QueryBucket_withTuning& queryBatch = arg->queryBatches[q];
#ifdef TIME_IT
                    arg->t.start();
#endif 

                    if (queryBatch.lshIndex == NULL) {
                        queryBatch.createLshIndex(*(arg->queryMatrix));
                    }

#ifdef TIME_IT
                    arg->t.stop();
                    arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif

                    LshIndex* queryIndex = queryBatch.lshIndex;

                    ///////////////////////
                    row_type user = queryBatch.startPos;
                    int start = queryBatch.startPos * arg->k;
                    int end = queryBatch.endPos * arg->k;
                    for (row_type i = start; i < end; i += arg->k) {
                        if (queryBatch.inactiveQueries[user - queryBatch.startPos]) {
                            user++;
                            continue;
                        }


                        const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                        double minScore = arg->topkResults[i].data;

                        if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets /////////////////
                            queryBatch.inactiveQueries[user - queryBatch.startPos] = true;
                            queryBatch.inactiveCounter++;
                            user++;
                            continue;
                        }

                        arg->moveTopkToHeap(i);


                        arg->queryId = arg->queryMatrix->getId(user);

                        if (probeBucket.t_b * probeBucket.normL2.second > minScore) {

                            plain.runTopK(query, probeBucket, arg);

                        } else {
#ifdef TIME_IT
                            arg->t.start();
#endif

                            double localTheta = minScore;

                            localTheta = localTheta * (localTheta > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);
                            int activeBlocks = arg->findActiveBlocks(localTheta);



#ifdef TIME_IT
                            arg->t.stop();
                            arg->boundsTime += arg->t.elapsedTime().nanos();
#endif
                            if (activeBlocks > LSH_SIGNATURES) {
                                plain.runTopK(query, probeBucket, arg);
                            } else {

#ifdef TIME_IT
                                arg->t.start();
#endif 
                                if (activeBlocks > index->initializedSketchesForIndex) {
                                    omp_set_lock(&(index->writelock));
                                    index->checkAndReallocateAll(arg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                            arg->sums, arg->countsOfBlockValues, arg->sketches, arg->rig); // need for lock here
                                    omp_unset_lock(&(index->writelock));
                                }                             

                                queryIndex->checkAndReallocateSingle(arg->queryMatrix, user, user - queryBatch.startPos, activeBlocks, arg->sums, arg->rig);
#ifdef TIME_IT
                                arg->t.stop();
                                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif 

                                processIndexesTopk(query, user - queryBatch.startPos, index, queryIndex, activeBlocks, probeBucket, arg);


                            }

                        }


                        arg->writeHeapToTopk(user);
                        user++;
                       
                    }
                }
            }
        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {

            LshIndex* index = static_cast<LshIndex*> (probeBucket.getIndex(LSH));

            for (row_type q = 0; q < arg->queryBatches.size(); q++) {

                if (arg->queryBatches[q].normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }

                QueryBucket_withTuning& queryBatch = arg->queryBatches[q];


                if (probeBucket.t_b == 1 || (probeBucket.t_b * queryBatch.normL2.first > probeBucket.bucketScanThreshold)) {
                    plain.run(queryBatch, probeBucket, arg);

                } else { // do it per query


#ifdef TIME_IT
                    arg->t.start();
#endif

                    if (queryBatch.lshIndex == NULL) {
                        queryBatch.createLshIndex(*(arg->queryMatrix));
                    }

#ifdef TIME_IT
                    arg->t.stop();
                    arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif


                    LshIndex* queryIndex = queryBatch.lshIndex;

                    for (row_type i = queryBatch.startPos; i < queryBatch.endPos; i++) {

                        const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                        if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                            break;

                        arg->queryId = arg->queryMatrix->getId(i);
                        if (probeBucket.t_b * query[-1] > probeBucket.bucketScanThreshold) {
                            plain.run(query, probeBucket, arg);
                        } else {

#ifdef TIME_IT
                            arg->t.start();
#endif
                            double localTheta = probeBucket.bucketScanThreshold / query[-1];
                            int activeBlocks = arg->findActiveBlocks(localTheta);

#ifdef TIME_IT
                            arg->t.stop();
                            arg->boundsTime += arg->t.elapsedTime().nanos();
#endif

                            if (activeBlocks > LSH_SIGNATURES || activeBlocks <= 0) {
                                plain.run(query, probeBucket, arg);
                            } else {

#ifdef TIME_IT
                                arg->t.start();
#endif  
                                if (activeBlocks > index->initializedSketchesForIndex) {
                                    omp_set_lock(&(index->writelock));
                                    index->checkAndReallocateAll(arg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                            arg->sums, arg->countsOfBlockValues, arg->sketches, arg->rig); // need for lock here
                                    omp_unset_lock(&(index->writelock));
                                }



                                queryIndex->checkAndReallocateSingle(arg->queryMatrix, i, i - queryBatch.startPos, activeBlocks, arg->sums, arg->rig);
#ifdef TIME_IT
                                arg->t.stop();
                                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif                    
                                processIndexes(query, i - queryBatch.startPos, index, queryIndex, activeBlocks, probeBucket, arg);
                            }
                        }
                    }
                }
            }

        }



    };




}

#endif	/* LSHRETRIEVER_H */

