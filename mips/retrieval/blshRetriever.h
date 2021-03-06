/* 
 * File:   blshRetriever.h
 * Author: chteflio
 *
 * Created on March 4, 2016, 3:39 PM
 */

#ifndef BLSHRETRIEVER_H
#define	BLSHRETRIEVER_H

namespace mips {

    class BlshRetriever : public LshRetriever {

        inline void processIndexesTopk(const double * query, row_type queryId,
                BlshIndex* index, LshIndex* queryIndex, int activeBlocks,
                ProbeBucket& probeBucket, RetrievalArguments* arg)const {

            row_type numCandidatesToVerify = 0;

#ifdef TIME_IT 
            arg->t.start();
#endif
            index->lshBins->getCandidatesWithBlshPruning(queryIndex->cosSketches->sketches, queryId, index->cosSketches->sketches,
                    index->minMatches, arg->candidatesToVerify, numCandidatesToVerify, arg->done, activeBlocks, probeBucket.startPos);

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
                BlshIndex* index, LshIndex* queryIndex, int activeBlocks,
                ProbeBucket& probeBucket, RetrievalArguments* arg)const {

            row_type numCandidatesToVerify = 0;

#ifdef TIME_IT 
            arg->t.start();
#endif

            index->lshBins->getCandidatesWithBlshPruning(queryIndex->cosSketches->sketches, queryId, index->cosSketches->sketches,
                    index->minMatches, arg->candidatesToVerify, numCandidatesToVerify, arg->done, activeBlocks, probeBucket.startPos);

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
        BlshRetriever() = default;
        ~BlshRetriever() = default;

        inline virtual void tune(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type sampleSize = probeBucket.xValues->size();

            if (sampleSize > 0) {

                // measure length
                // access plain properly
                if (retrArg[0].competitorMethod != NULL) {
                    plain.sampleTimes = *(retrArg[0].competitorMethod);
                } else {
                    plain.tune(probeBucket, prevBucket, retrArg);
                }

                double avgLengthTime = plain.sampleTotalTime / sampleSize;


                // measure LSH
                sampleTimes.resize(sampleSize);
                BlshIndex* index = static_cast<BlshIndex*> (probeBucket.getIndex(BLSH));

                // create lsh signatures for sample queries
                LshIndex queryIndex;
                queryIndex.initializeLists(probeBucket.xValues->size());

                double avgBlockTime = 0;
                double min_t_b = -1;
                int breakpoint = sampleSize - 1;

                for (int i = sampleSize - 1; i >= 0; i--) {
                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;

                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    double localTheta = probeBucket.bucketScanThreshold / query[-1];
                    int activeBlocks = repetitionsForTheta.findActiveBlocks(localTheta);


                    if (activeBlocks > LSH_SIGNATURES || activeBlocks <= 0) {
                        sampleTimes[i] = plain.sampleTimes[i];

                    } else {
                        if (activeBlocks > index->initializedSketchesForIndex && i < sampleSize - 1) {

                            double expectedTime = avgBlockTime * activeBlocks / (index->initializedSketchesForIndex * (probeBucket.xValues->size() - i));

                            if (expectedTime > avgLengthTime) {
                                min_t_b = localTheta;
                                breakpoint = i;
                                break;
                            }
                        }


                        index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, true);

                        retrArg[0].tunerTimer.start();
                        queryIndex.checkAndReallocateSingle(retrArg[t].queryMatrix, ind, i, activeBlocks, retrArg[t].sums);

                        LshRetriever::processIndexes(query, i, index, &queryIndex, activeBlocks, probeBucket, &retrArg[0]);
                        retrArg[0].tunerTimer.stop();
                        sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                        avgBlockTime += sampleTimes[i] / activeBlocks;
                    }
                }

                if (breakpoint == sampleSize - 1) {
                    // compare
                    findCutOffPoint(sampleTimes, plain.sampleTimes, sampleTotalTime, t_b_indx);

                    if (plain.sampleTotalTime < sampleTotalTime) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        double value = (t_b_indx == 0 ? -1 : probeBucket.xValues->at(t_b_indx).result);
                        probeBucket.setAfterTuning(1, value);
                    }

                } else {

                    if (breakpoint == sampleSize - 2) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        probeBucket.setAfterTuning(1, min_t_b);
                    }
                }
            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            BlshIndex* index = static_cast<BlshIndex*> (probeBucket.getIndex(BLSH));

            for (auto& queryBatch : arg->queryBatches) {

                if (queryBatch.maxLength() < probeBucket.bucketScanThreshold) {
                    break;
                }

                if (probeBucket.t_b == 1 || (probeBucket.t_b * queryBatch.minLength() > probeBucket.bucketScanThreshold)) {
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

                    for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {

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
                            int activeBlocks = repetitionsForTheta.findActiveBlocks(localTheta);

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
//                                 std::cout<<"worst: "<<index->worst<<std::endl;
//                                std::cout<<"About to process my first vector with BLSH"<<std::endl; 
//                                std::cout<<"activeBlocks: "<<activeBlocks<<" localTheta: "<<localTheta<<
//                                         " initializedSketchesForIndex: "<<index->initializedSketchesForIndex<<std::endl;
                                if (activeBlocks > index->initializedSketchesForIndex) {
                                    index->lockIndex();
                                    index->checkAndReallocateAll(arg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                            arg->sums, arg->countsOfBlockValues, arg->sketches, true); // need for lock here////////////////////////////
                                    index->unlockIndex();
                                }



                                queryIndex->checkAndReallocateSingle(arg->queryMatrix, i, i - queryBatch.startPos, activeBlocks, arg->sums);
#ifdef TIME_IT
                                arg->t.stop();
                                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif                    
                                
                                
                                processIndexes(query, i - queryBatch.startPos, index, queryIndex, activeBlocks, probeBucket, arg);
//                                std::cout<<"Done to process my first vector with LSH"<<std::endl;
//                                exit(1);
                            }
                        }
                    }
                }
            }

        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, const ProbeBucket& prevBucket, std::vector<RetrievalArguments>& retrArg) {

            row_type sampleSize = (probeBucket.xValues != nullptr ? probeBucket.xValues->size() : 0);

            if (sampleSize > 0) {
                plain.sampleTimes.reserve(sampleSize);

                for (row_type i = 0; i < sampleSize; ++i) {
                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;

                    plain.sampleTimes.emplace_back(probeBucket.sampleThetas[t][ind].lengthTime);
                    plain.sampleTotalTime += probeBucket.sampleThetas[t][ind].lengthTime;
                }

                double avgLengthTime = plain.sampleTotalTime / sampleSize;

                // measure LSH
                sampleTimes.resize(sampleSize);
                BlshIndex* index = static_cast<BlshIndex*> (probeBucket.getIndex(BLSH));

                if (!index->isInitialized()) {
#if defined(TIME_IT)
                    retrArg[0].t.start();
#endif           
                    index->initializeLists(*(retrArg[0].probeMatrix), std::numeric_limits<double>::max(), true, retrArg[0].R, probeBucket.startPos, probeBucket.endPos);

#if defined(TIME_IT)
                    retrArg[0].t.stop();
                    retrArg[0].initializeListsTime += retrArg[0].t.elapsedTime().nanos();
#endif
                }
                // create lsh signatures for sample queries
                LshIndex queryIndex;
                queryIndex.initializeLists(sampleSize);



                double avgBlockTime = 0;
                double min_t_b = -1;
                int breakpoint = sampleSize - 1;




                for (int i = sampleSize - 1; i >= 0; i--) {

                    int t = probeBucket.xValues->at(i).i;
                    int ind = probeBucket.xValues->at(i).j;
                    const double* query = retrArg[t].queryMatrix->getMatrixRowPtr(ind);

                    const std::vector<QueueElement>& prevResults = prevBucket.sampleThetas[t].at(ind).results; //just reading

                    std::copy(prevResults.begin(), prevResults.end(), retrArg[0].heap.begin());
                    std::make_heap(retrArg[0].heap.begin(), retrArg[0].heap.end(), std::greater<QueueElement>());

                    double minScore = retrArg[0].heap.front().data;
                    double localTheta = minScore * (minScore > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);
                    int activeBlocks = repetitionsForTheta.findActiveBlocks(localTheta);

                    if (activeBlocks > LSH_SIGNATURES || activeBlocks <= 0) {
                        sampleTimes[i] = plain.sampleTimes[i];
                    } else {

                        if (activeBlocks > index->initializedSketchesForIndex && i < sampleSize - 1) {

                            double expectedTime = avgBlockTime * activeBlocks / (index->initializedSketchesForIndex * (sampleSize - i));

                            if (expectedTime > avgLengthTime) {
                                min_t_b = localTheta;
                                breakpoint = i;
                                break;
                            }
                        }
                        index->checkAndReallocateAll(retrArg[0].probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                retrArg[0].sums, retrArg[0].countsOfBlockValues, retrArg[0].sketches, true);


                        retrArg[0].tunerTimer.start();
                        queryIndex.checkAndReallocateSingle(retrArg[t].queryMatrix, ind, i, activeBlocks, retrArg[0].sums);

                        LshRetriever::processIndexesTopk(query, i, index, &queryIndex, activeBlocks, probeBucket, &retrArg[0]);

                        retrArg[0].tunerTimer.stop();
                        sampleTimes[i] = retrArg[0].tunerTimer.elapsedTime().nanos();
                        avgBlockTime += sampleTimes[i] / activeBlocks;
                    }
                }


                if (breakpoint == sampleSize - 1) { // no break at all
                    // compare
                    findCutOffPoint(sampleTimes, plain.sampleTimes, sampleTotalTime, t_b_indx);

                    if (plain.sampleTotalTime < sampleTotalTime) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        double value = (t_b_indx == 0 ? -1 : probeBucket.xValues->at(t_b_indx).result);
                        probeBucket.setAfterTuning(1, value);
                    }

                } else {

                    if (breakpoint == sampleSize - 2) {
                        probeBucket.setAfterTuning(1, 1);
                    } else {
                        probeBucket.setAfterTuning(1, min_t_b);
                    }
                }
            } else {
                probeBucket.setAfterTuning(prevBucket.numLists, prevBucket.t_b);
            }


        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) const {

            if (probeBucket.t_b == 1) {
                plain.runTopK(probeBucket, arg);
            } else { // do it per query

                BlshIndex* index = static_cast<BlshIndex*> (probeBucket.getIndex(BLSH));

                if (!index->isInitialized()) {
#ifdef TIME_IT
                    arg->t.start();
#endif           
                    index->initializeLists(*(arg->probeMatrix), arg->worstMinScore, true, arg->R, probeBucket.startPos, probeBucket.endPos);

#ifdef TIME_IT
                    arg->t.stop();
                    arg->initializeListsTime += arg->t.elapsedTime().nanos();
#endif
                }
                if (index->minMatches == nullptr) {
                    double worstCaseTheta = arg->worstMinScore;

                    if (worstCaseTheta > 0) {
                        worstCaseTheta *= probeBucket.invNormL2.second;
                        worstCaseTheta = (worstCaseTheta > 1 ? 1 : worstCaseTheta); // it will break in the loop afterwards
                    } else {
                        worstCaseTheta *= probeBucket.invNormL2.first;
                        worstCaseTheta = (worstCaseTheta < -1 ? -1 : worstCaseTheta); // it will have to check everything
                    }
                    arg->worstMinScore = std::numeric_limits<double>::max();
                    index->lockIndex();                    
                    index->allocateBayesLSHMemory(worstCaseTheta);
                    index->unlockIndex();
                }



                for (auto& queryBatch : arg->queryBatches) {

                    if (queryBatch.isWorkDone())
                        continue;

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
                        if (queryBatch.isQueryInactive(user)) {
                            user++;
                            continue;
                        }


                        const double* query = arg->queryMatrix->getMatrixRowPtr(user);

                        double minScore = arg->topkResults[i].data;
                        double minScoreAppr = minScore;

#ifdef         HYBRID_APPROX             
                        //                        minScoreAppr += arg->queryMatrix->epsilonEquivalents[user];
                        //                        arg->currEpsilonAppr = arg->queryMatrix->epsilonEquivalents[user];

                        if (minScoreAppr >= 0) {
                            minScoreAppr *= (1 + arg->epsilon);
                            arg->currEpsilonAppr = (1 + arg->epsilon);
                        } else {
                            //                            minScoreAppr *= (1 - arg->epsilon);
                            arg->currEpsilonAppr = 1; //(1 - arg->epsilon);
                        }

#endif



                        if (probeBucket.normL2.second < minScoreAppr) {// skip this bucket and all other buckets
                            queryBatch.inactivateQuery(user);
                            user++;
#ifdef         HYBRID_APPROX                     

                            arg->moveTopkToHeap(i);
                            arg->totalErrorAfterResults += approximateHybridError(*arg, probeBucket.normL2.second, probeBucket.invNormL2.second);
#endif                   


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

#if defined(HYBRID_APPROX) && (defined(RELATIVE_APPROX)   || defined(ABS_APPROX))
                            double localTheta = minScoreAppr;
#else
                            double localTheta = minScore;
#endif

                            localTheta = localTheta * (localTheta > 0 ? probeBucket.invNormL2.second : probeBucket.invNormL2.first);
                            int activeBlocks = repetitionsForTheta.findActiveBlocks(localTheta);



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
                                    index->lockIndex();
                                    index->checkAndReallocateAll(arg->probeMatrix, true, probeBucket.startPos, probeBucket.endPos, activeBlocks,
                                            arg->sums, arg->countsOfBlockValues, arg->sketches, true); // need for lock here
                                    index->unlockIndex();
                                }

                                queryIndex->checkAndReallocateSingle(arg->queryMatrix, user, user - queryBatch.startPos, activeBlocks, arg->sums);
#ifdef TIME_IT
                                arg->t.stop();
                                arg->preprocessTime += arg->t.elapsedTime().nanos();
#endif 

                                processIndexesTopk(query, user - queryBatch.startPos, index, queryIndex, activeBlocks, probeBucket, arg);


                            }

                        }

                        minScore = arg->heap.front().data;

                        if (arg->worstMinScore > minScore) {
                            arg->worstMinScore = minScore;
                        }

                        arg->writeHeapToTopk(user);
                        user++;

                    }
                }
            }
        }

    };



}


#endif	/* BLSHRETRIEVER_H */

