

/* 
 * File:   CosineSketches.h
 * Author: chteflio
 *
 * Created on June 25, 2015, 1:53 PM
 */

#ifndef COSINESKETCHES_H
#define	COSINESKETCHES_H

#include <boost/dynamic_bitset.hpp>

namespace ta {

    struct CosineSketches {
        row_type hashCodeLength; // must be a multiple of 8.
        row_type bytesPerCode;
        row_type numHashTables;
        row_type sketchSize;
        long long userOffset;
        row_type nVectors;
        uint8_t *sketches;

        CosineSketches(row_type nVectors, row_type hashCodeLength, row_type numHashTables) : numHashTables(numHashTables),
        hashCodeLength(hashCodeLength), sketchSize(numHashTables * hashCodeLength), nVectors(nVectors), bytesPerCode(hashCodeLength / 8), sketches(0) {
            assert(hashCodeLength % 8 == 0);
            userOffset = bytesPerCode * ((long)numHashTables);

        }

        inline void alloc() {
            long long  s = (long long)(bytesPerCode * numHashTables) * (long long)nVectors;	    
            sketches = new uint8_t[s](); /////////////
        }

        ~CosineSketches() {
            delete[] sketches;
        }

        void buildSingle(double* vec, row_type posInBucket, col_type colNum, std::vector<float>& sums, RandomIntGaussians* rig,
                uint8_t* sketch, row_type startBlock, row_type endBlock, row_type startHashBit, row_type endHashBit, row_type startOffset) {
            float val;
            long dimensionOffset;

            for (row_type j = 0; j < colNum; j++) { // dot product with nHashBits different random vectors
                val = (float) vec[j];
                if (val != 0) {
                    dimensionOffset = j * sketchSize + startHashBit;
                    for (long k = dimensionOffset, l = startHashBit; l < endHashBit; k++, l++) {
                        sums[l] += val * rig->intToFloatCache[rig->intGaussians[k]];
                    }
                }
            }

            long long t = posInBucket * userOffset + startOffset;
            row_type l = startHashBit; //0

            for (row_type k = startBlock; k < endBlock; k++) {

                for (row_type b = 0; b < bytesPerCode; b++) {
                    for (row_type n = 0; n < 8; n++) {
                        sketch[t] <<= 1;
                        if (sums[l] > 0)
                            sketch[t]++;
                        sums[l] = 0; // clear out sums
                        l++;
                    }
                    t++; // move byte by byte
                }
            }
        }

        void buildSingle2(double* vec, row_type posInBucket, col_type colNum, std::vector<float>& sums, RandomIntGaussians* rig,
                uint8_t* sketch, row_type startBlock, row_type endBlock, row_type startHashBit, row_type endHashBit, row_type startOffset) {
            float val;
            long dimensionOffset;

            for (row_type j = 0; j < colNum; j++) { // dot product with nHashBits different random vectors
                val = (float) vec[j];
                if (val != 0) {
                    dimensionOffset = j * sketchSize + startHashBit;
                    for (long k = dimensionOffset, l = startHashBit; l < endHashBit; k++, l++) {
                        sums[l] += val * rig->intToFloatCache[rig->intGaussians[k]];
                    }
                }
            }
            

            row_type t = startOffset;
            row_type l = startHashBit; //0

            for (row_type k = startBlock; k < endBlock; k++) {

                for (row_type b = 0; b < bytesPerCode; b++) {
                    for (row_type n = 0; n < 8; n++) {
                        sketch[t] <<= 1;
                        if (sums[l] > 0)
                            sketch[t]++;
                        sums[l] = 0; // clear out sums
                        l++;
                    }
                    t++; // move byte by byte
                }
            }
            
        }

        void buildBatch(const VectorMatrix* matrix, row_type start, row_type end, std::vector<float>& sums, uint8_t* sketch, RandomIntGaussians* rig, row_type startBlock, row_type endBlock) {

            row_type startOffset = bytesPerCode * startBlock;
            row_type startHashBit = startBlock * hashCodeLength;
            row_type endHashBit = endBlock * hashCodeLength;

            for (row_type i = start; i < end; i++) {
                // sums must be completely zero at this point.
                double* vec = matrix->getMatrixRowPtr(i);
                buildSingle(vec, (i - start), matrix->colNum, sums, rig, sketch, startBlock, endBlock, startHashBit, endHashBit, startOffset);
            }
            
//            std::cout<<"build batch finished"<<std::endl;
        }

        void printSketch(row_type posInBucket) {
            int k = 0;
            for (int i = 0; i < numHashTables; i++) {
                for (int j = 0; j < bytesPerCode; j++) {
                    std::cout << (int) sketches[posInBucket * numHashTables * bytesPerCode + k] << " ";
                    k++;
                }
                std::cout << " | ";
            }
            std::cout << std::endl;
        }
    };
}

#endif	/* COSINESKETCHES_H */

