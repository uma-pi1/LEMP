/* 
 * File:   connect.h
 * Author: chteflio
 *
 * Created on September 22, 2014, 4:33 PM
 */

#ifndef CONNECT_H
#define	CONNECT_H

#include <taLib/ap/includes.h>
#include <taLib/structs/VectorMatrix.h>

namespace ta {

    typedef struct {
        double simT; /* Similarity threshold */
        double simT2; /* Similarity threshold squared */
        int32_t seed; /* Seed for the random number generator */
        double epsilon; /* Recall parameter for BayesLSH methods */
        int nSketches; /* Number of NBITS_SKETCH bit sketches to check */


        da_csr_t *docs; /* Documents structure */


        /* bayeslsh related */
        double minMathesT; // threshold for matching, r = c2r(simT)
        int nHashBits; // nSketches * NBITS_SKETCH;
        uint32_t *sketches; // sketches
        int32_t *minMatches; // minimum number of hashes that should be observed to meet simT
        ssize_t *numPruned; // number pruned by each sketch limit

        /* row and col permutations */
        idx_t *rperm; // row permutation
        idx_t *cpermOrigToProbe; // col permutation
        idx_t *cpermProbeToOrig; // col permutation

    } my_params_t;

    void freeParams(my_params_t** params) {
        da_csr_Free(&(*params)->docs);
        gk_free((void**) &(*params)->rperm, LTERM);
        gk_free((void**) params, LTERM);
    }

    inline val_t dotProduct(idx_t *bInd, val_t *bVal, idx_t bLen, val_t *marker) {

        idx_t i;
        val_t result = 0.0;

        for (i = 0; i < bLen; i++)
            if (marker[bInd[i]] != 0.0) // marker is actually the hashval
                result += bVal[i] * marker[bInd[i]];

        return result;
    }

    inline static val_t dotProduct(idx_t *aInd, val_t *aVal, idx_t aLen,
            idx_t *bInd, val_t *bVal, idx_t bLen) {

        idx_t i, j;
        val_t result = 0;

        for (i = 0, j = 0; i < aLen && j < bLen;) {
            if (aInd[i] == bInd[j]) {
                result += aVal[i] * bVal[j];
                i++;
                j++;
            } else if (aInd[i] > bInd[j]) {
                j++;
            } else
                i++;
        }
        return result;
    }


    // this is my function
    inline da_csr_t *da_csr_Read(const VectorMatrix& matrix, row_type start, row_type end) {

        //        std::cout << "da_csr_Read" << std::endl;
        row_type k;
        //	ssize_t i, j, k, l;
        size_t nrows, ncols, nnz; //, nfields, fmt, ncon, lnlen, read, size;
        ptr_t *rowptr;
        idx_t *rowind;
        //int32_t *iinds, *jinds, *bptr = NULL, *bind = NULL, rid, len, nnz2, nr, nc, nz, rv;
        double *rowval = NULL; //, *vals, fval;
        //float *bval = NULL;
        //char readsizes, readwgts;
        //char *line = NULL, *head, *tail, fmtstr[256];
        //FILE *fpin;
        da_csr_t *mat = NULL;



        nrows = end - start;
        ncols = matrix.colNum;


        // also needs nnz

        //////////////////
        //		readsizes = 0;
        //		readwgts  = 0;
        //		readvals  = 1;
        //		numbering = 1;
        ////////////////////////	

        mat = da_csr_Create();
        mat->nrows = nrows;


        rowptr = mat->rowptr = da_pmalloc(nrows + 1, NULL);

        // find nnz
        nnz = 0;
        for (row_type i = start; i < end; i++) {
            const double* vector = matrix.getMatrixRowPtr(i);

            for (row_type j = 0; j < ncols; j++) {

                if (vector[j] != 0) {
                    //                if (fabs(vector[j]) >= 0.0000001) {/////////////////////////////   
                    nnz++;
                }
            }
        }

        rowind = mat->rowind = da_imalloc(nnz, NULL);

        rowval = mat->rowval = da_vsmalloc(nnz, 1.0, NULL);
        //rowval = mat->rowval = da_dsmalloc(nnz, 1.0, NULL);


        /*----------------------------------------------------------------------
         * Read the sparse matrix file
         *---------------------------------------------------------------------*/

        rowptr[0] = 0;
        k = 0;


        for (row_type i = start; i < end; i++) {

            const double* vector = matrix.getMatrixRowPtr(i);

            for (row_type j = 0; j < ncols; j++) {

                if (vector[j] != 0) {/////////////////////////////
                    //                if (fabs(vector[j]) >= 0.0000001) {/////////////////////////////   

                    rowind[k] = j; // column Ids of individual coordinates
                    ncols = gk_max(rowind[k], ncols); //////////////////////////////////////////////
                    rowval[k] = vector[j];
                    k++;
                }

            }
            rowptr[i - start + 1] = k;

        }
        mat->ncols = ncols; ////////////////// +1?     

        //        std::cout << "NROWS: " << mat->nrows << " NCOLS: " << mat->ncols << " NNZ: " << nnz << std::endl;

        return mat;
    }

    // my function

    inline void readInputData(my_params_t *params, const VectorMatrix& matrix, row_type start, row_type end) {
        da_csr_t *docs;

        //        std::cout << "In readInputData" << std::endl;

        docs = da_csr_Read(matrix, start, end);

        params->docs = docs;
    }



}

#endif	/* CONNECT_H */

