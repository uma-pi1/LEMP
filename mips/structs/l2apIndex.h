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
 * File:   l2apIndex.h
 * Author: chteflio
 * Code based on L2AP
 * Created on September 23, 2014, 2:10 PM
 */

#ifndef L2APINDEX_H
#define	L2APINDEX_H


namespace mips {

    class L2apIndex : public Index {

        /**
         * Index part of the vector after its search is complete
         */
        void l2apIndexRow(idx_t rid, da_csr_t *docs, std::vector<double>& cweights, val_t *lengths) { //the cweights should be from the query side

            std::vector<double> hashwgt(ncols);
            ssize_t j, k;
            val_t myval, maxrval;
            double b1, b2, b2l, sqb2, pscore, simT;
            ptr_t *rowptr = docs->rowptr; // index in rowind/rowval where each row starts
            idx_t *rowind = docs->rowind; // colid for each nnz
            val_t *rowval = docs->rowval; // val for each nnz
            simT = params->simT;
#ifdef PSCV
            pscore = 0.0;
#endif
#if defined(L2PS) || defined (IDXL2)
            b2 = sqb2 = b2l = 0.0;
#elif defined(LENPS) || defined (IDXLEN)
            b2 = 0.5;
#endif
#if ( defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7) ) && \
	!( defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8) )
            maxrval = 0.0;
#endif

            // this is my part
            double rs2 = 1;
            for (int i = rowptr[rid + 1] - 1; i >= rowptr[rid]; i--) {
                myval = rowval[i];
                rs2 -= (double) myval * myval; //rst               
                lengths[i] = sqrt(rs2); //rs4
            }


            double v = 0;
            // index terms in row i
            myval = b1 = 0.0;
            for (endptr[rid] = rowptr[rid + 1], j = rowptr[rid]; j < endptr[rid]; j++) {
                myval = rowval[j];

                /////////////////// now that is my block
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)

                if (fabs(myval) > v)
                    v = fabs(myval);
                hashwgt[rowind[j]] = v; // prefix scan of max row values, i.e. remaining max weight in row
#endif
                /////////////////////////////////////////////////
                col_type col = rowind[j];
                // now find the correct column wrt the permutation
                col = params->cpermProbeToOrig[col];
                b1 += (double) fabs(myval * cweights[col]); ////////////////// not sure if rwgts[rid] is relevant for me || into fabs for negatives
                ////////////////////////////////////////////////


                // adjust b2 estimate of prefix similarity with all other vectors
#if defined(L2PS) || defined (IDXL2)
                b2 += (double) myval * myval;
                sqb2 = sqrt(b2);
#elif defined(LENPS) || defined (IDXLEN)
                b2 += (double) 0.5 * myval * myval;
#endif

                // check whether to start indexing
#ifdef L2PS
                if (gk_min(b1, sqb2) >= simT)
                    break;
#elif defined(LENPS)
                if (gk_min(b1, b2) >= simT)
                    break;
#else
                if (b1 >= simT)
                    break;
#endif

#ifdef PSCV
                //update pscore
#ifdef L2PS
                pscore = gk_min(b1, sqb2);
#elif defined(LENPS)
                pscore = gk_min(b1, b2);
#else
                pscore = b1;
#endif
#endif

#if defined (L2PS)
                b2l = b2;
#endif

#if ( defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7) ) && \
	!( defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8) )
                // new f.i. max value
                if (myval > maxrval)
                    maxrval = myval;
#endif

#if defined(RS1) && defined(RS3)
                // adjust inverted index max col weights
                if (myval > criwgts[rowind[j]])
                    criwgts[rowind[j]] = myval;
#endif
            }
            // truncate the rest of the vector, since we're going to
            // index it.
            endptr[rid] = j;

            // update max forward index row val for row rid
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
            rfiwgts[rid] = j > rowptr[rid] ? hashwgt[rowind[j - 1]] : 0.0; //////////////// how is this supposed to work? hashwgt are for queries

#elif defined(DP1) || defined(DP3) || defined(DP5) || defined(DP7)
            rfiwgts[rid] = maxrval;
#endif

            //store the pscore for this row
#ifdef PSCV
            ps[rid] = pscore;
#endif

            // set b2 to the prefix norm up to current index but not including
#ifdef L2PS
            b2 = b2l;
            sqb2 = sqrt(b2);
#elif defined(LENPS)
            b2 -= 0.5;
#endif




            // index the remaining part of the vector
            for (; j < rowptr[rid + 1]; j++) {
                indexSize++;
                k = rowind[j];
#if defined (IDXL2)
                invIdx->lens[ invIdx->ends[k] ] = sqb2;
                b2 += (double) rowval[j] * rowval[j];
                sqb2 = sqrt(b2);
#elif defined(IDXLEN)
                invIdx->lens[ invIdx->ends[k] ] = b2;
                b2 += (double) 0.5 * rowval[j] * rowval[j];
#endif
                invIdx->ids[ invIdx->ends[k] ] = rid;
                invIdx->vals[ invIdx->ends[k] ] = rowval[j];
                invIdx->ends[k]++;
#if defined(RS1) && defined(RS3)
                if (rowval[j] > criwgts[k])
                    criwgts[k] = rowval[j];
#endif
            }

        }

        void freeMyParams(my_params_t** params) {
            da_csr_FreeAll(&(*params)->docs, LTERM); // &(*params)->neighbors,       
            gk_free((void**) &(*params)->rperm, &(*params)->cpermOrigToProbe, &(*params)->cpermProbeToOrig, LTERM); //&(*params)->iFile, &(*params)->oFile, &(*params)->dataset,&(*params)->filename, 
            gk_free((void**) params, LTERM);
        }

        void l2apReorderDocs(my_params_t *params, da_csr_t **docs, idx_t *rsizes, idx_t *csizes, val_t *rwgts,
                val_t *cwgts, idx_t *rperm, idx_t *cperm) {
            ssize_t i, j, k, nnz;
            idx_t nrows, ncols;
            da_iikv_t *orderColNnz = nullptr;
            da_ivkv_t *orderMaxVal = nullptr;
            char pv = 0;
            ptr_t *rowptr, *nrowptr;
            idx_t *rowind, *nrowind, *tperm;
            val_t v, *rowval, *nrowval;

            nrows = (*docs)->nrows;
            ncols = (*docs)->ncols;
            rowptr = (*docs)->rowptr;
            rowind = (*docs)->rowind;
            rowval = (*docs)->rowval;
            nnz = rowptr[nrows];

            //record col nnzs
            da_iset(ncols, 0, csizes);
            for (i = 0; i < nnz; i++)
                csizes[rowind[i]]++;

            //assign memory for ordering
            orderMaxVal = da_ivkvmalloc(nrows, nullptr); //"l2apReorderDocs: orderMaxVal"
            nrowptr = da_pmalloc(nrows + 1, nullptr); //"l2apReorderDocs: nrowptr"
            nrowind = da_imalloc(nnz, nullptr); // "l2apReorderDocs: nrowind")
            nrowval = da_vmalloc(nnz, nullptr); //"l2apReorderDocs: nrowval"

            orderColNnz = da_iikvmalloc(ncols, nullptr); //"l2apReorderDocs: orderColNnz"
            //get new column order


            for (j = 0; j < ncols; j++) {
                orderColNnz[j].key = j;
                orderColNnz[j].val = csizes[j];
            }

            da_iikvsortd(ncols, orderColNnz); // sort columns
            for (j = 0; j < ncols; j++) {
                cperm[orderColNnz[j].key] = j; // store column permutation
                csizes[j] = orderColNnz[j].val; // new column size after permuting

                params->cpermProbeToOrig[j] = orderColNnz[j].key;
            }

            // set new column ids in docs
            for (i = 0; i < nnz; i++)
                rowind[i] = cperm[rowind[i]];
            // get max weights for rows and columns + row sizes
            for (i = 0; i < nrows; i++) {
                rsizes[i] = (idx_t) (rowptr[i + 1] - rowptr[i]); //record row nnzs
                for (j = rowptr[i]; j < rowptr[i + 1]; j++) {

                    v = fabs(rowval[j]);

                    if (v > rwgts[i])
                        rwgts[i] = v; // store max val for row
                    if (v > cwgts[rowind[j]])
                        cwgts[rowind[j]] = v; //store max val for col
                }
                //set up kv array to get new row order
                orderMaxVal[i].key = i;
                orderMaxVal[i].val = rwgts[i];
            }
            da_ivkvsortd(nrows, orderMaxVal); // sort rows
            for (i = 0; i < nrows; i++) {
                rperm[i] = orderMaxVal[i].key; // store new order
                rwgts[i] = orderMaxVal[i].val; // new row weight after permuting
            }
            // permute rows in matrix
            for (nrowptr[0] = 0, nnz = 0, j = 0, i = 0; i < nrows; i++) {
                k = rperm[i];
                da_icopy(rowptr[k + 1] - rowptr[k], rowind + rowptr[k], nrowind + nnz);
                da_vcopy(rowptr[k + 1] - rowptr[k], rowval + rowptr[k], nrowval + nnz);
                rsizes[i] = rowptr[k + 1] - rowptr[k];
                nnz += rsizes[i];
                nrowptr[++j] = nnz;
            }
            gk_free((void**) &(*docs)->rowptr, &(*docs)->rowind, &(*docs)->rowval, LTERM);
            (*docs)->rowptr = nrowptr;
            (*docs)->rowind = nrowind;
            (*docs)->rowval = nrowval;

            da_csr_SortIndices((*docs), DA_ROW);

            //inverse the row permutation
            tperm = rperm;
            rperm = da_imalloc(nrows, nullptr); //"l2apReorderDocs: rperm"
            for (i = 0; i < nrows; i++)
                rperm[i] = tperm[i];

            gk_free((void**) &orderColNnz, &orderMaxVal, &tperm, LTERM); //&cperm, // keep the column permutations

            params->rperm = rperm;
            params->cpermOrigToProbe = cperm;

        }

    public:
        my_params_t *params;
        ssize_t nnz;
        idx_t nrows, ncols, mrlen, mclen;
        ptr_t *rowptr, *endptr;
        idx_t *cands, *rsizes, *csizes, *rperm, *cperm, *hashsz;
        val_t maxrval, myval;
        val_t *rwgts, *rfiwgts, *cwgts, *criwgts, *hashsum, *ps, *lengths, *sums; //*hashval, *hashlen, *hashwgt,
        da_invIdxJ_t *invIdx; //alternative inverted index for value-driven similarities//

        comp_type indexSize;

        inline L2apIndex() : params(nullptr), rowptr(nullptr), endptr(nullptr), cands(nullptr), rsizes(nullptr), 
        csizes(nullptr), rperm(nullptr), cperm(nullptr), hashsz(nullptr), rwgts(nullptr), 
        rfiwgts(nullptr), cwgts(nullptr), criwgts(nullptr), hashsum(nullptr), ps(nullptr), lengths(nullptr), sums(nullptr),
        invIdx(nullptr), indexSize(0) {//, hashval(NULL), hashlen(NULL),hashwgt(NULL),

        }

        inline ~L2apIndex() {

            if (initialized) {
                gk_free((void**) &rsizes, &csizes, &rwgts, &rfiwgts, &cwgts,
                        &cands, &endptr, //&hashval, &hashlen, &hashwgt, 
                        &lengths, &ps, //  &sums,
                        &invIdx->ids, &invIdx->vals, &invIdx->ends, &invIdx->starts,
                        &invIdx->lens, &invIdx, LTERM); //&invIdx->accum,

                if (params != nullptr)
                    freeMyParams(&params);
            }
        }

        inline void initializeLists(const VectorMatrix& matrix, double worstCaseTheta, std::vector<double>& cweights,
                ta_size_type start = 0, ta_size_type end = 0) {
            omp_set_lock(&writelock);


            if (!initialized) {
                params = (my_params_t *) gk_malloc(sizeof (my_params_t), nullptr);


                rg::Timer t;

                t.start();
                readInputData(params, matrix, start, end);
                t.stop();
                matrixToMatrixTime += t.elapsedTime().nanos();
                //            std::cout<<"read data"<<std::endl;

                ssize_t i, j, k;
                da_csr_t *docs = params->docs;
                nrows = docs->nrows; // num rows
                ncols = docs->ncols; // num cols


                if (worstCaseTheta > 0) {

                    worstCaseTheta /= matrix.getVectorLength(start);
                } else {
                    worstCaseTheta /= matrix.getVectorLength(end - 1);
                }


                params->simT = worstCaseTheta;
                params->simT2 = params->simT * params->simT;


                rsizes = da_imalloc(nrows, nullptr); // "l2apFindNeighbors: rsizes"
                csizes = da_ismalloc(ncols, 0.0, nullptr); //"l2apFindNeighbors: csizes"
                rperm = da_imalloc(nrows, nullptr); //"l2apFindNeighbors: rperm"
                cperm = da_imalloc(ncols, nullptr); //"l2apFindNeighbors: cperm"
                params->cpermProbeToOrig = da_imalloc(ncols, nullptr); //"l2apFindNeighbors: cperm"
                rwgts = da_vsmalloc(nrows, 0.0, nullptr); //"l2apFindNeighbors: rwgts"
                rfiwgts = da_vmalloc(nrows, nullptr); //"l2apFindNeighbors: rfiwgts"
                cwgts = da_vsmalloc(ncols, 0.0, nullptr); //"l2apFindNeighbors: cwgts"
#if defined(RS1) && defined(RS3)
                criwgts = da_vsmalloc(ncols, 0.0, "l2apFindNeighbors: criwgts");
#endif
                //            hashval = da_vsmalloc(ncols, 0, NULL); // hash values of x //"l2apFindNeighbors: hashval"
#if defined(DP7) || defined(DP8)
                hashsz = da_imalloc(ncols, "l2apFindNeighbors: hashsz");
#endif
#if defined(L2CV) || defined(LENCV)
                //            hashlen = da_vsmalloc(ncols, 0, NULL); // hash suffix lengths of x //"l2apFindNeighbors: hashlen"
#endif
#if defined(DP3) || defined(DP4)
                hashsum = da_vsmalloc(ncols, 0, "l2apFindNeighbors: hashsum"); // hash suffix sums of x
#endif
#if defined(DP2) || defined(DP4) || defined(DP6) || defined(DP8)
                //hashwgt = da_vmalloc(ncols, NULL); //"l2apFindNeighbors: hashwgt"
#endif
#ifdef PSCV
                ps = da_vmalloc(nrows, nullptr); // score prior to indexing threshold for each row //"l2ap2FindNeighbors: ps"
#endif

                // reorder matrix
                l2apReorderDocs(params, &docs, rsizes, csizes, rwgts, cwgts, rperm, cperm); /////////////////////////////////

                //            std::cout << "ok with REORDER" << std::endl;

                nnz = docs->rowptr[nrows]; // nnz in docs
                rowptr = docs->rowptr; // index in rowind/rowval where each row starts

                // allocate memory for candidates and the inverted indexes -- needs nnz & col counts
                cands = da_imalloc(nrows, nullptr); //"l2apFindNeighbors: candidates"
#if defined(DP1) || defined(DP2) || defined(DP3) || defined(DP4)
                sums = da_vsmalloc(nnz, 0.0, "l2apFindNeighbors: sums"); // keep track of suffix sums
#endif
                lengths = da_vsmalloc(nnz, 0.0, nullptr); // keep track of suffix lengths "l2apFindNeighbors: lengths"
                endptr = da_pmalloc(nrows, nullptr); //"l2apFindNeighbors: endptr"
                invIdx = (da_invIdxJ_t*) gk_malloc(sizeof (da_invIdxJ_t), nullptr); //"l2apFindNeighbors: invIdx"
                invIdx->ids = da_ismalloc(nnz, 0, nullptr); //"l2apFindNeighbors: invIdx->ids"
                invIdx->vals = da_vsmalloc(nnz, 0.0, nullptr); //"l2apFindNeighbors: invIdx->ids"
                invIdx->lens = da_vsmalloc(nnz, 0.0, nullptr); //"l2apFindNeighbors: invIdx->ids"
                invIdx->starts = da_pmalloc(ncols, nullptr); //"l2apFindNeighbors: invIdx->starts"
                invIdx->ends = da_pmalloc(ncols, nullptr); //"l2apFindNeighbors: invIdx->ends"
                for (j = 0, k = 0; j < ncols; j++) {
                    invIdx->starts[j] = invIdx->ends[j] = k;
                    k += csizes[j];
                }

                for (i = 0; i < nrows; i++) {
                    endptr[i] = rowptr[i + 1];
                    l2apIndexRow(i, docs, cweights, lengths);
                }

                initialized = true;

            }



            omp_unset_lock(&writelock);


        }

    };





}

#endif	/* L2APINDEX_H */

