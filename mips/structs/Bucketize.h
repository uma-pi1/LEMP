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
 * Bucketize.h
 *
 *  Created on: Mar 13, 2014
 *      Author: chteflio
 */

#ifndef BUCKETIZE_H_
#define BUCKETIZE_H_

#include <cmath>

namespace mips {

    /*
     * cacheSizeInKB : per processor
     */
    inline row_type computeBlockOffsetsByFactorCacheFittingForItems(const std::vector<QueueElement>& elements, row_type size,
            std::vector<row_type>& blockOffsets, double factor, row_type minItems, row_type cacheSizeInKB, row_type rank,
            const LempArguments& args) {

        int doubleSize = 8; //bytes
        int row_typeSize = sizeof (row_type);
        double t = 0.8;

        int singleVectorSpace = (rank + 1) * doubleSize + (doubleSize + row_typeSize); // the basic thing (coordinates+length) + lengthInfo

        cacheSizeInKB = t * cacheSizeInKB * 1024;

        switch (args.method) {
            case LEMP_I:
            case LEMP_LI:
                singleVectorSpace += row_typeSize; // for the candidatesToVerify
                singleVectorSpace += (doubleSize + doubleSize); // for the ext_cp_array
                singleVectorSpace += (doubleSize + row_typeSize) * rank; // index space

                break;

            case LEMP_C:
            case LEMP_LC:
                singleVectorSpace += row_typeSize; // for the candidatesToVerify
                singleVectorSpace += (doubleSize + row_typeSize); // for the cp_array
                singleVectorSpace += (doubleSize + row_typeSize) * rank; // index space

                break;

            case LEMP_TA:
                singleVectorSpace += (doubleSize + row_typeSize) * rank; // index space

                break;

            case LEMP_TREE:
                //singleVectorSpace *= 2; // in reality I have no idea how much space the tree will take but here I try to double the space
                //                singleVectorSpace += 2 * row_typeSize + 2 * row_typeSize + 3 * doubleSize + 3 * comp_typeSize;   

                singleVectorSpace += 5 * row_typeSize; // scale, start, end, point, numDescendants
                singleVectorSpace += 3 * row_typeSize; // & dataset, parent, child 
                singleVectorSpace += 3 * doubleSize; // base, parentDistance, furtherDistance
                singleVectorSpace += args.threads * (3 * doubleSize + 2 * row_typeSize); // stats


                break;

            case LEMP_AP: // change that later
                singleVectorSpace += row_typeSize; // for the candidatesToVerify
                singleVectorSpace += (2 * doubleSize + row_typeSize) * rank; // index space (max)
                singleVectorSpace += doubleSize; // for the accum

                break;

            case LEMP_BLSH:
	      singleVectorSpace +=  LSH_SIGNATURES; // for the sketches
	      
            case LEMP_LSH:

                singleVectorSpace += row_typeSize * LSH_SIGNATURES; // for the data
                cacheSizeInKB -= LSH_SIGNATURES * 257 * row_typeSize;
                singleVectorSpace += row_typeSize; // for the candidatesToVerify

		break;


        }

        blockOffsets.clear();
        blockOffsets.push_back(0);


        if (args.k > 0) {
            blockOffsets.push_back(args.k);
        }

        row_type maxItems = cacheSizeInKB / singleVectorSpace;


        // 	maxItems = 1000000;/////////////////////

        row_type ind;

        row_type maxBlockSize = 0;

        while (blockOffsets[blockOffsets.size() - 1] < size) {
            double value = elements[blockOffsets[blockOffsets.size() - 1]].data * factor;


            auto up = std::upper_bound(elements.begin(), elements.begin() + size, QueueElement(value, 0), std::greater<QueueElement>());
            ind = up - elements.begin();
            if (ind - blockOffsets[blockOffsets.size() - 1] < minItems) {
                ind = blockOffsets[blockOffsets.size() - 1] + minItems;
                ind = (ind > size ? size : ind);
            } else if (ind - blockOffsets[blockOffsets.size() - 1] > maxItems) {
                ind = blockOffsets[blockOffsets.size() - 1] + maxItems;
            }

            blockOffsets.push_back(ind);

            int blockSize = blockOffsets[blockOffsets.size() - 1] - blockOffsets[blockOffsets.size() - 2];
            if (blockSize > maxBlockSize)
                maxBlockSize = blockSize;
        }

        if (blockOffsets.size() > 1)
            blockOffsets.pop_back();


        return maxBlockSize * singleVectorSpace;
    };

    /*
     * cacheSizeInKB : per processor
     */
    inline void computeBlockOffsetsForUsersFixed(row_type size, std::vector<row_type>& blockOffsets,
            row_type cacheSizeInKB, row_type rank, const LempArguments& args, row_type maxBlockSize) {

        int doubleSize = 8; //bytes
        int row_typeSize = sizeof (row_type);
        int col_typeSize = sizeof (col_type);

        double t = 0.8;

        int singleVectorSpace = (rank + 1) * doubleSize + (doubleSize + row_typeSize); // the basic thing (coordinates+length)
        singleVectorSpace += args.k * (doubleSize + row_typeSize); // resultSet

        cacheSizeInKB = (cacheSizeInKB * 1024 - maxBlockSize) * t;


        switch (args.method) {
            case LEMP_I:
            case LEMP_LI:
            case LEMP_C:
            case LEMP_LC:
                singleVectorSpace += rank * col_typeSize; // queue
                break;
            case LEMP_LSH:
            case LEMP_BLSH:
                singleVectorSpace += row_typeSize * LSH_SIGNATURES; // for the sketches
                break;

        }

        row_type maxItems = cacheSizeInKB / singleVectorSpace;



        int blocks = ceil((double) size / maxItems);

        // 	blocks = 1;///////////////////////

        computeDefaultBlockOffsets(size, blocks, blockOffsets);
    }

    template<typename T>
    inline void bucketize(std::vector<T>& buckets, const VectorMatrix& matrix,
            const std::vector<row_type>& blockOffsets, const LempArguments& args) {

        row_type start, end;
        buckets.clear();
        buckets.resize(blockOffsets.size());

        for (row_type i = 0; i < blockOffsets.size(); ++i) {
            start = blockOffsets[i];
            end = (i == (blockOffsets.size() - 1) ? matrix.rowNum : blockOffsets[i + 1]);
            buckets[i].init(matrix, start, end, args);
        }
    }




}
#endif /* BUCKETIZE_H_ */
