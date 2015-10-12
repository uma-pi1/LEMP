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
 * File:   apRetriever.h
 * Author: chteflio
 * code based on L2AP
 * Created on September 24, 2014, 10:24 AM
 */

#ifndef APRETRIEVER_H
#define	APRETRIEVER_H

namespace ta {

    class apRetriever : public Retriever {

        void l2apProcessCandidates(const double* query, row_type nnzQuery, double maxQueryCoord,
                row_type numCandidates, double localTheta, L2apIndex* index,
                ProbeBucket& probeBucket, RetrievalArguments* arg) {


            double queryLength = query[-1];
            ssize_t i, j, p, r, sz, cid;
            val_t rwgt;
            double cval, simT, simTe, ip;
            ptr_t *rowptr;
            idx_t *rowind;
            val_t *rowval;


            row_type posInProbeMatrix;

            rowptr = index->params->docs->rowptr; // index in rowind/rowval where each row starts
            rowind = index->params->docs->rowind; // colid for each nnz
            rowval = index->params->docs->rowval; // val for each nnz


            simT = localTheta;
#if defined(L2CV)
            simTe = localTheta - LENBNDEPS;
#endif
#if defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
            rwgt = maxQueryCoord;
#endif

#if defined(DP5) || defined(DP6)
            sz = nnzQuery;
#endif

            for (i = 0; i < numCandidates; ++i) {

                cid = arg->candidatesToVerify[i];

                if (arg->accum[cid] == -2) {
                    arg->accum[cid] = -1;
                    continue;
                }


                cval = arg->accum[cid];
                arg->accum[cid] = -1;


#ifdef PSCV
                if (cval + ((double) index->ps[cid]) < simT) {
                    continue;
                }
#endif


#ifdef DP5
                if (cval + ((double) gk_min(index->endptr[cid] - rowptr[cid], sz) * index->rfiwgts[cid] * rwgt) < simT) {
                    continue;
                }
#endif

                r = rowptr[cid];

#if defined(DP2) || defined(DP3) || defined(DP4) \
		|| defined(DP6) || defined(DP7) || defined(DP8)


                col_type col = 0;
                // advance to next term in common
                // check again to ensure correct dotp hashVal contain query stuff


                for (j = index->endptr[cid] - 1, p = index->endptr[cid] - r; j >= r; j--, p--) {

                    if (arg->hashval[rowind[j]] != 0.0)
                        break;

                }

                if (j < r)
                    goto check;


#if defined(DP6)
                col = index->params->cpermProbeToOrig[rowind[j]];

                if (cval + ((double) gk_min(p, sz) * index->rfiwgts[cid] * arg->hashwgt[col]) < simT) {
                    continue;
                }
#endif

#endif


#ifdef L2CV

                for (; j >= r; j--)

                    if (arg->hashval[rowind[j]] != 0.0) {
                        cval += (double) arg->hashval[rowind[j]] * rowval[j];

                        if (j > r && cval + ((double) arg->hashlen[rowind[j]] * index->lengths[j]) < simTe) {
                            goto nexcid;
                        }
                    }
#endif

check: // check similarity value
                arg->comparisons++;



                if (cval >= simT) { //simT
                    posInProbeMatrix = index->params->rperm[cid] + probeBucket.startPos; // do I also need the startPos?                  

                    ip = cval * queryLength * arg->probeMatrix->getVectorLength(posInProbeMatrix);

                    if (arg->k == 0) {

                        if (ip >= arg->theta) { //simT                              
//                            arg->results.push_back(MatItem(ip, arg->queryId, arg->probeMatrix->getId(posInProbeMatrix)));
                            arg->results.emplace_back(ip, arg->queryId, arg->probeMatrix->getId(posInProbeMatrix));
                        }
                    }
                }
nexcid:
                continue;

            }//for all candidates


        }

        void l2apProcessCandidatesTopk(const double* query, row_type nnzQuery, double maxQueryCoord,
                row_type numCandidates, double localTheta, L2apIndex* index,
                ProbeBucket& probeBucket, RetrievalArguments* arg) {


            double queryLength = query[-1];
            double minScore = 0;
            minScore = arg->heap.front().data;

            ssize_t i, j, p, r, sz, cid;
            val_t rwgt;
            double cval, simT, simTe, ip;
            ptr_t *rowptr;
            idx_t *rowind;
            val_t *rowval;


            row_type posInProbeMatrix;

            rowptr = index->params->docs->rowptr; // index in rowind/rowval where each row starts
            rowind = index->params->docs->rowind; // colid for each nnz
            rowval = index->params->docs->rowval; // val for each nnz


            simT = localTheta;
#if defined(L2CV)
            simTe = localTheta - LENBNDEPS;
#endif
#if defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
            rwgt = maxQueryCoord;
#endif

#if defined(DP5) || defined(DP6)
            sz = nnzQuery;
#endif

            for (i = 0; i < numCandidates; ++i) {

                cid = arg->candidatesToVerify[i];

                if (arg->accum[cid] == -2) {
                    arg->accum[cid] = -1;
                    continue;
                }

                cval = arg->accum[cid];
                arg->accum[cid] = -1;


#ifdef PSCV
                if (cval + ((double) index->ps[cid]) < simT) {
                    continue;
                }
#endif



#ifdef DP5
                if (cval + ((double) gk_min(index->endptr[cid] - rowptr[cid], sz) * index->rfiwgts[cid] * rwgt) < simT) {
                    continue;
                }
#endif

                r = rowptr[cid];

#if defined(DP2) || defined(DP3) || defined(DP4) \
		|| defined(DP6) || defined(DP7) || defined(DP8)


                col_type col = 0;
                // advance to next term in common
                // check again to ensure correct dotp hashVal contain query stuff


                for (j = index->endptr[cid] - 1, p = index->endptr[cid] - r; j >= r; j--, p--) {

                    if (arg->hashval[rowind[j]] != 0.0)
                        break;

                }

                if (j < r)
                    goto check;


#if defined(DP6)
                col = index->params->cpermProbeToOrig[rowind[j]];

                if (cval + ((double) gk_min(p, sz) * index->rfiwgts[cid] * arg->hashwgt[col]) < simT) {
                    continue;
                }
#endif

#endif

#ifdef L2CV

                for (; j >= r; j--)

                    if (arg->hashval[rowind[j]] != 0.0) {
                        cval += (double) arg->hashval[rowind[j]] * rowval[j];

                        if (j > r && cval + ((double) arg->hashlen[rowind[j]] * index->lengths[j]) < simTe) {
                            goto nexcid;
                        }
                    }
#endif

check: // check similarity value
                arg->comparisons++;



                if (cval >= simT) { //simT
                    posInProbeMatrix = index->params->rperm[cid] + probeBucket.startPos; // do I also need the startPos?                  

                    ip = cval * queryLength * arg->probeMatrix->getVectorLength(posInProbeMatrix);

                    if (ip > minScore) {

                        std::pop_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                        arg->heap.pop_back();
//                        arg->heap.push_back(QueueElement(ip, arg->probeMatrix->getId(posInProbeMatrix)));
                        arg->heap.emplace_back(ip, arg->probeMatrix->getId(posInProbeMatrix));
                        std::push_heap(arg->heap.begin(), arg->heap.end(), std::greater<QueueElement>());
                        minScore = arg->heap.front().data;

                        if (minScore > 0) {
                            simT = minScore * probeBucket.invNormL2.second;
                        } else {
                            simT = minScore * probeBucket.invNormL2.first;
                        }
#if defined(L2CV)
                        simTe = simT - LENBNDEPS;
#endif
                    }
                }
nexcid:
                continue;

            }


        }

        void l2apFindMatches(const double * query, row_type nnzQuery, double maxQueryCoord,
                L2apIndex* index, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            double queryLength = query[-1];

            double localTheta;

            localTheta = probeBucket.bucketScanThreshold / queryLength;


            row_type numCandidatesToVerify = 0;

            ssize_t i, j;
            idx_t k, ncands, cid;
            ptr_t *rowptr, *starts, *ends;
            idx_t *rowind, *idxids;
            val_t v, myval;
            val_t *rowval, *idxvals, *idxlens;

            double simT, simTe, bnd, rs2, bsq, acc;

            rowptr = index->params->docs->rowptr; // index in rowind/rowval where each row starts
            rowind = index->params->docs->rowind; // col id for each nnz
            rowval = index->params->docs->rowval; // val for each nnz
            idxids = index->invIdx->ids;
            idxvals = index->invIdx->vals;
            starts = index->invIdx->starts; // starting points for lists in the inv index
            ends = index->invIdx->ends; // ending points for lists in the inv index

            ncands = 0; // number of candidates for this doc
            simT = localTheta;
#if defined(IDXL2)
            rs2 = 1.0;
            bsq = 1.0;
#endif
#if defined(L2CG) || defined(LENCG)
            idxlens = index->invIdx->lens;
#endif
#if defined(L2CG)
            simTe = localTheta - LENBNDEPS;
#endif

            // define vars for hashing DP data
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
            v = 0.0; /////////// to work with svd
#endif

            // compute rs1, the remaining score
#if defined(RS1) || defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4) \
		|| defined(DP6) || defined(DP7) || defined(DP8)

            //here access query in the correct order(see cperm)

            for (i = 0; i < probeBucket.colNum; i++) {

                double qi = query[index->params->cpermProbeToOrig[i]];
                if (qi == 0)
                    continue;


#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)

                if (fabs(qi) > v) {
                    v = fabs(qi);
                }
                arg->hashwgt[index->params->cpermProbeToOrig[i]] = v; // prefix scan of max row values, i.e. remaining max weight in row
#endif

            }
#endif


            //here access query in the correct order(see cperm)
            for (i = probeBucket.colNum - 1; i >= 0; i--) {
                k = index->params->cpermProbeToOrig[i];


                myval = query[k];

                if (myval == 0)
                    continue;

                arg->hashval[i] = myval;

#if defined(RS4)
                bnd = bsq;
#endif


#if defined(IDXL2)
                rs2 -= (double) myval * myval; //rst
                bsq = sqrt(rs2); //rs4
#endif


                arg->hashlen[i] = bsq;

                for (j = starts[i]; j < ends[i]; ++j) {


                    cid = idxids[j]; // potential candidate from the inv index
                    acc = arg->accum[cid];

                    if (acc > -1) {
                        // accumulate
                        acc = arg->accum[cid] += myval * idxvals[j];

                        // check l2-norm bound
#if defined(L2CG)
                        if (acc + bsq * idxlens[j] < simTe) {
                            arg->accum[cid] = -2;
                        }
#endif
                    } else if (bnd >= simT && acc > -2) {
                        // accumulate
                        if (acc == -1) {
                            acc = arg->accum[cid] = myval * idxvals[j];

                            arg->candidatesToVerify[numCandidatesToVerify] = cid;
                            numCandidatesToVerify++;

                        } else {
                            acc = arg->accum[cid] += myval * idxvals[j];
                        }


                        // check l2-norm bound
#if defined(L2CG)
                        if (acc + bsq * idxlens[j] < simTe) {
                            arg->accum[cid] = -2;
                        }
#endif
                    }

                }

            }


            l2apProcessCandidates(query, nnzQuery, maxQueryCoord, numCandidatesToVerify, localTheta, index, probeBucket, arg);

            // reset hashval
            std::fill(arg->hashval.begin(), arg->hashval.end(), 0);

        }

        void l2apFindMatchesTopk(const double * query, row_type nnzQuery, double maxQueryCoord,
                L2apIndex* index, ProbeBucket& probeBucket, RetrievalArguments* arg) {

            double localTheta;


            if (arg->heap.front().data > 0) {
                localTheta = arg->heap.front().data * probeBucket.invNormL2.second;
            } else {
                localTheta = arg->heap.front().data * probeBucket.invNormL2.first;
            }


            row_type numCandidatesToVerify = 0;

            ssize_t i, j;
            idx_t k, ncands, cid;
            ptr_t *rowptr, *starts, *ends;
            idx_t *rowind, *idxids;
            val_t v, myval;
            val_t *rowval, *idxvals, *idxlens;

            double simT, simTe, bnd, rs2, bsq, acc;

            rowptr = index->params->docs->rowptr; // index in rowind/rowval where each row starts
            rowind = index->params->docs->rowind; // col id for each nnz
            rowval = index->params->docs->rowval; // val for each nnz
            idxids = index->invIdx->ids;
            idxvals = index->invIdx->vals;
            starts = index->invIdx->starts; // starting points for lists in the inv index
            ends = index->invIdx->ends; // ending points for lists in the inv index

            ncands = 0; // number of candidates for this doc
            simT = localTheta;
#if defined(IDXL2)
            rs2 = 1.0;
            bsq = 1.0;
#endif
#if defined(L2CG) || defined(LENCG)
            idxlens = index->invIdx->lens;
#endif
#if defined(L2CG)
            simTe = localTheta - LENBNDEPS;
#endif

            // define vars for hashing DP data
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
            v = 0.0; /////////// to work with svd
#endif

            // compute rs1, the remaining score
#if defined(RS1) || defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4) \
		|| defined(DP6) || defined(DP7) || defined(DP8)

            //here access query in the correct order(see cperm)

            for (i = 0; i < probeBucket.colNum; i++) {

                double qi = query[index->params->cpermProbeToOrig[i]];
                if (qi == 0)
                    continue;


#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)

                if (fabs(qi) > v) {
                    v = fabs(qi);
                }
                arg->hashwgt[index->params->cpermProbeToOrig[i]] = v; // prefix scan of max row values, i.e. remaining max weight in row
#endif

            }
#endif


            //here access query in the correct order(see cperm)
            for (i = probeBucket.colNum - 1; i >= 0; i--) {
                k = index->params->cpermProbeToOrig[i];


                myval = query[k];

                if (myval == 0)
                    continue;

                arg->hashval[i] = myval;

#if defined(RS4)
                bnd = bsq;
#endif


#if defined(IDXL2)
                rs2 -= (double) myval * myval; //rst
                bsq = sqrt(rs2); //rs4
#endif


                arg->hashlen[i] = bsq;

                for (j = starts[i]; j < ends[i]; ++j) {


                    cid = idxids[j]; // potential candidate from the inv index
                    acc = arg->accum[cid];

                    if (acc > -1) {
                        // accumulate
                        acc = arg->accum[cid] += myval * idxvals[j];

                        // check l2-norm bound
#if defined(L2CG)
                        if (acc + bsq * idxlens[j] < simTe) {
                            arg->accum[cid] = -2;
                        }
#endif
                    } else if (bnd >= simT && acc > -2) {
                        // accumulate
                        if (acc == -1) {
                            acc = arg->accum[cid] = myval * idxvals[j];

                            arg->candidatesToVerify[numCandidatesToVerify] = cid;
                            numCandidatesToVerify++;

                        } else {
                            acc = arg->accum[cid] += myval * idxvals[j];
                        }


                        // check l2-norm bound
#if defined(L2CG)
                        if (acc + bsq * idxlens[j] < simTe) {
                            arg->accum[cid] = -2;
                        }
#endif
                    }

                }

            }


            l2apProcessCandidatesTopk(query, nnzQuery, maxQueryCoord, numCandidatesToVerify, localTheta, index, probeBucket, arg);

            // reset hashval
            std::fill(arg->hashval.begin(), arg->hashval.end(), 0);

        }




    public:

        apRetriever() = default;

        ~apRetriever() = default;

        inline virtual void run(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void run(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline void runTopK(const double* query, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(QueryBucket_withTuning& queryBatch, ProbeBucket& probeBucket, RetrievalArguments* arg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tune(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void tuneTopk(ProbeBucket& probeBucket, std::vector<RetrievalArguments>& retrArg) {
            std::cerr << "Error! You shouldn't have called that" << std::endl;
            exit(1);
        }

        inline virtual void runTopK(ProbeBucket& probeBucket, RetrievalArguments* arg) {
            L2apIndex* index = static_cast<L2apIndex*> (probeBucket.getIndex(AP));


            for(auto& queryBatch: arg->queryBatches){

                if (queryBatch.inactiveCounter == queryBatch.rowNum)
                    continue;


                if (!index->initialized) {
#if defined(TIME_IT)
                    arg->t.start();
#endif
                    index->initializeLists(*(arg->probeMatrix), arg->worstMinScore, arg->queryMatrix->cweights, probeBucket.startPos, probeBucket.endPos);
#if defined(TIME_IT)
                    arg->t.stop();
                    arg->initializeListsTime += arg->t.elapsedTime().nanos();
#endif
                    arg->worstMinScore = std::numeric_limits<double>::max();
                }


                //////////////////////////
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

                    if (probeBucket.normL2.second < minScore) {// skip this bucket and all other buckets
                        queryBatch.inactiveQueries[user - queryBatch.startPos] = true;
                        queryBatch.inactiveCounter++;
                        user++;
                        continue;
                    }

                    arg->moveTopkToHeap(i);

                    arg->queryId = arg->queryMatrix->getId(user);
                    row_type nnzQuery = arg->queryMatrix->vectorNNZ[user];
                    double maxQueryCoord = arg->queryMatrix->maxVectorCoord[user];

                    l2apFindMatchesTopk(query, nnzQuery, maxQueryCoord, index, probeBucket, arg);

                    minScore = arg->heap.front().data;

                    if (arg->worstMinScore > minScore) {
                        arg->worstMinScore = minScore;
                    }

                    arg->writeHeapToTopk(user);
                    user++;
                }
            }

        }

        inline virtual void run(ProbeBucket& probeBucket, RetrievalArguments* arg) {

            L2apIndex* index = static_cast<L2apIndex*> (probeBucket.getIndex(AP));

               for(auto& queryBatch: arg->queryBatches){

                if (queryBatch.normL2.second < probeBucket.bucketScanThreshold) {
                    break;
                }

                for (row_type i = queryBatch.startPos; i < queryBatch.endPos; ++i) {
                    const double* query = arg->queryMatrix->getMatrixRowPtr(i);

                    if (query[-1] < probeBucket.bucketScanThreshold)// skip all users from this point on for this bucket
                        break;

                    row_type nnzQuery = arg->queryMatrix->vectorNNZ[i];
                    double maxQueryCoord = arg->queryMatrix->maxVectorCoord[i];
                    arg->queryId = arg->queryMatrix->getId(i);
                    l2apFindMatches(query, nnzQuery, maxQueryCoord, index, probeBucket, arg);

                }
            }

        }


    };




}

#endif	/* APRETRIEVER_H */

