

/* 
 * File:   CosineSketches.h
 * Author: chteflio
 *
 * Created on June 25, 2015, 1:53 PM
 */

#ifndef COSINESKETCHES_H
#define	COSINESKETCHES_H

#include <boost/dynamic_bitset.hpp>

#include "RandomIntGaussians.h"

namespace ta {
    


    class CosineSketches {


    public:
        uint8_t *sketches;

        CosineSketches() :  sketches(nullptr) {          
           
        }

        inline void alloc(row_type nVectors) {
            long long s = (long long) (coreLshInfo.bytesPerCode * coreLshInfo.numHashTables) * (long long) nVectors;
            sketches = new uint8_t[s](); /////////////
        }

        ~CosineSketches() {
            if (sketches != nullptr)
                delete[] sketches;
        }

        void buildSingle(const double* vec, row_type posInBucket, col_type colNum, std::vector<float>& sums,
                uint8_t* sketch, row_type startBlock, row_type endBlock, row_type startHashBit, row_type endHashBit, row_type startOffset) {
            float val;
            long dimensionOffset;

            RandomIntGaussians& rig = RandomIntGaussians::getInstance(colNum);

            for (row_type j = 0; j < colNum; ++j) { // dot product with nHashBits different random vectors
                val = (float) vec[j];
                if (val != 0) {
                    dimensionOffset = j * coreLshInfo.sketchSize + startHashBit;
                    for (long k = dimensionOffset, l = startHashBit; l < endHashBit; ++k, ++l) {
                        sums[l] += val * rig.intToFloatCache[rig.intGaussians[k]];
                    }
                }
            }

            long long t = posInBucket * coreLshInfo.userOffset + startOffset;
            row_type l = startHashBit; //0

            for (row_type k = startBlock; k < endBlock; ++k) {

                for (row_type b = 0; b < coreLshInfo.bytesPerCode; ++b) {
                    for (row_type n = 0; n < 8; ++n) {
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

        void buildBatch(const VectorMatrix* matrix, row_type start, row_type end, std::vector<float>& sums, uint8_t* sketch, row_type startBlock, row_type endBlock) {

            row_type startOffset = coreLshInfo.bytesPerCode * startBlock;
            row_type startHashBit = startBlock * coreLshInfo.hashCodeLength;
            row_type endHashBit = endBlock * coreLshInfo.hashCodeLength;

            for (row_type i = start; i < end; ++i) {
                // sums must be completely zero at this point.
                const double* vec = matrix->getMatrixRowPtr(i);
                buildSingle(vec, (i - start), matrix->colNum, sums, sketch, startBlock, endBlock, startHashBit, endHashBit, startOffset);
            }

        }

        void printSketch(row_type posInBucket) const {
            int k = 0;
            for (int i = 0; i < coreLshInfo.numHashTables; ++i) {
                for (int j = 0; j < coreLshInfo.bytesPerCode; ++j) {
                    std::cout << (int) sketches[posInBucket * coreLshInfo.numHashTables * coreLshInfo.bytesPerCode + k] << " ";
                    k++;
                }
                std::cout << " | ";
            }
            std::cout << std::endl;
        }
    };
}

#endif	/* COSINESKETCHES_H */

