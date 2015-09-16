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
 * VectorMatrix.h
 *
 *  Created on: Oct 10, 2013
 *      Author: chteflio
 */

#ifndef VECTORMATRIX_H_
#define VECTORMATRIX_H_

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <taLib/structs/read.cpp>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <util/exception.h>
#include <util/io.h>

#include <string>
#include <ostream>
#include <iomanip>
#include <boost/unordered_map.hpp>



#include <pmmintrin.h> //SSE3




using boost::unordered_map;

namespace ta {

    inline void computeDefaultBlockOffsets(row_type size, row_type blocks, std::vector<row_type>& blockOffsets, row_type start = 0) {

        blockOffsets.resize(blocks);
        row_type minSize = size / blocks;
        row_type remainder = size % blocks;
        for (row_type i = 0; i < blocks; i++) {
            if (i == 0) {
                blockOffsets[i] = start;
            } else {
                blockOffsets[i] = minSize + blockOffsets[i - 1];
                if (remainder > 0) {
                    ++blockOffsets[i];
                    --remainder;
                }
            }
        }
    };


    // typedef

    struct VectorMatrix {
        double* data;
        bool shuffled, normalized, extraMult;

        std::vector<double> cweights;
        std::vector<double> maxVectorCoord;
        std::vector<row_type> vectorNNZ;
        std::vector<QueueElement> colFrequencies;
        
        row_type colNum, offset, lengthOffset;
        row_type rowNum;
        std::vector<QueueElement> lengthInfo; // data: length id: vectorId


        std::vector<double> gammaEquivalents;

        // for simd instruction  
        int sizeDiv2;




        inline VectorMatrix() : data(0), shuffled(false), normalized(false), lengthOffset(1){////////////////////// 1 is for padding
        }

        inline ~VectorMatrix() {        
            free(data);
        }

        inline void readFromFile(std::string& fileName, bool left = true) {
            std::ifstream file(fileName.c_str(), std::ios_base::in);
            while (file.peek() == '%') {
                skipLineFromFile(file);
            }

            ta_size_type col; // columns
            ta_size_type row; // rows
            file >> row >> col;


            std::cout << "Reading file: " << fileName << " -- columns: " << col << "  rows: " << row << std::endl;

            if (left) {
                if (pow(2, sizeof (col_type) * 8) - 1 < col) {
                    std::cerr << "Your vectors have dimensionality " << col << " which is more than what lemp is compiled to store. Change the col_type in BasicStructs.h and recompile." << std::endl;
                    exit(1);
                }

                if (pow(2, sizeof (row_type) * 8) - 1 < row) {
                    std::cerr << "Your dataset has " << row << " vectors which is more than what lemp is compiled to store. Change the row_type in BasicStructs.h and recompile." << std::endl;
                    exit(1);
                }

            } else {
                if (pow(2, sizeof (col_type) * 8) - 1 < row) {
                    std::cerr << "Your vectors have dimensionality " << row << " which is more than what lemp is compiled to store. Change the col_type in BasicStructs.h and recompile." << std::endl;
                    exit(1);
                }

                if (pow(2, sizeof (row_type) * 8) - 1 < col) {
                    std::cerr << "Your dataset has " << col << " vectors which is more than what lemp is compiled to store. Change the row_type in BasicStructs.h and recompile." << std::endl;
                    exit(1);
                }
            }

            rowNum = (left ? row : col);
            colNum = (left ? col : row);
            offset = colNum + 2;
            sizeDiv2 = colNum & (-2);

            if (colNum < NUM_LISTS) {
                std::cout << "WARNING: Your vectors have dimensionality" << colNum << " and the tuner will try to search among " << NUM_LISTS <<
                        ". Perhaps you want to change the parameter NUM_LISTS in Definitions.h and recompile" << std::endl;
            }
            if (LOWER_LIMIT_PER_BUCKET >= rowNum) {
                std::cout << "WARNING: You have " << rowNum << " vectors and the tuner will try to take a sample of at least  " << LOWER_LIMIT_PER_BUCKET <<
                        " vectors per probe bucket. Perhaps you want to change the parameter LOWER_LIMIT_PER_BUCKET in Definitions.h and recompile" << std::endl;
            }


            extraMult = (sizeDiv2 < colNum);

            if (extraMult)
                offset++;


            //std::cout<<"rowNum: "<<rowNum<<" colNum: "<<(int)colNum<<std::endl;
            normalized = false;
            //            data = new double[offset * rowNum];            
            int res = posix_memalign((void **) &data, 16, sizeof (double)* offset * rowNum);
            //            std::cout << "res: " << res << " null? " << (data == NULL) << std::endl;

            if (res != 0) {
                std::cout << "Problem with allocating memory" << std::endl;
                exit(1);
            }

            for (int i = 0; i < rowNum; i++) {
                setLengthInData(i, 1);
            }

            if (extraMult) {
                for (int i = 0; i < rowNum; i++) {
                    data[(i + 1) * offset - 3] = 0; // zero-out the last padding
                }
            }


            if (left) {
                if (file) {
                    for (ta_size_type i = 0; i < col; i++) {// read one column
                        for (ta_size_type j = 0; j < row; j++) {
                            double f;
                            file >> f;

                            double* d = getMatrixRowPtr(j);
                            d[i] = f;
                        }
                    }
                }
                file.close();

            } else {
                if (file) {
                    for (ta_size_type i = 0; i < col; i++) {// read one column
                        for (ta_size_type j = 0; j < row; j++) {
                            double f;
                            file >> f;

                            double* d = getMatrixRowPtr(i);
                            d[j] = f;
                        }
                    }
                }
                file.close();
            }
        }

        inline void init(VectorMatrix& matrix, bool sort, bool ignoreLength) {
            rowNum = matrix.rowNum;
            colNum = matrix.colNum;
            offset = colNum + 2;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            if (extraMult)
                offset++;

            normalized = true;

            //            data = new double[offset * rowNum];
            int res = posix_memalign((void **) &data, 16, sizeof (double)* offset * rowNum);
            //            std::cout << "res: " << res << " null? " << (data == NULL) << std::endl;

            if (res != 0) {
                std::cout << "Problem with allocating memory" << std::endl;
                exit(1);
            }

            if (extraMult) {
                for (int i = 0; i < rowNum; i++) {
                    data[(i + 1) * offset - 3] = 0; // zero-out the last padding
                }
            }

            lengthInfo.resize(rowNum);

            if (ignoreLength) {
#pragma omp parallel for  schedule(static, 1000)
                // get lengths
                for (int i = 0; i < rowNum; i++) {
                    double len = 0;
                    const double* vec = matrix.getMatrixRowPtr(i);
                    for (int j = 0; j < colNum; j++) {
                        len += vec[j] * vec[j];
                    }
                    len = sqrt(len);

                    lengthInfo[i] = QueueElement(1, i);

                    setLengthInData(i, 1); ////////////
                    double x = 1 / len;

                    double * d1 = getMatrixRowPtr(i);
                    double * d2 = matrix.getMatrixRowPtr(i);
                    for (int j = 0; j < colNum; j++) {
                        d1[j] = d2[j] * x;
                    }
                }


            } else {
#pragma omp parallel for schedule(static,1000)
                for (int i = 0; i < rowNum; i++) {
                    double len = 0;
                    const double* vec = matrix.getMatrixRowPtr(i);
                    for (int j = 0; j < colNum; j++) {
                        len += vec[j] * vec[j];
                    }
                    len = sqrt(len);

                    lengthInfo[i] = QueueElement(len, i);
                }

                if (sort) {
                    shuffled = true;
                    std::sort(lengthInfo.begin(), lengthInfo.end(), std::greater<QueueElement>());
                }



#pragma omp parallel for schedule(static,1000)
                for (int i = 0; i < rowNum; i++) {
                    setLengthInData(i, lengthInfo[i].data);
                    double x = 1 / lengthInfo[i].data;


                    double * d1 = getMatrixRowPtr(i);
                    double * d2 = matrix.getMatrixRowPtr(lengthInfo[i].id);

                    for (int j = 0; j < colNum; j++) {
                        d1[j] = d2[j] * x;
                    }
                }

            }
        }

        inline void addVectors(const VectorMatrix& matrix, std::vector<row_type>& dataIds) {
            rowNum = dataIds.size();
            colNum = matrix.colNum;
            offset = colNum + 2;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            if (extraMult)
                offset++;

            normalized = false;

            //            data = new double[offset * rowNum];
            int res = posix_memalign((void **) &data, 16, sizeof (double)* offset * rowNum);
            //            std::cout << "res: " << res << " null? " << (data == NULL) << std::endl;

            if (res != 0) {
                std::cout << "Problem with allocating memory" << std::endl;
                exit(1);
            }

            if (extraMult) {
                for (int i = 0; i < rowNum; i++) {
                    data[(i + 1) * offset - 3] = 0; // zero-out the last padding
                }
            }

            lengthInfo.resize(rowNum);


            for (int i = 0; i < rowNum; i++) {

                const double* vec = matrix.getMatrixRowPtr(dataIds[i]);

                lengthInfo[i] = QueueElement(1, dataIds[i]);

                double * d1 = getMatrixRowPtr(i);
                for (int j = 0; j < colNum; j++) {
                    d1[j] = vec[j];
                }
            }
        }


        // please only assign to const double *

        inline double* getMatrixRowPtr(row_type row) const {// the row starts from pos 1. Do ptr[-1] to get the length
            return &data[row * offset + 1 + lengthOffset];
        }

        inline double getVectorLength(row_type row) const {
            return data[row * offset + lengthOffset];
        }

        inline double setLengthInData(row_type row, double len) {
            return data[row * offset + lengthOffset] = len;
        }


        inline row_type getId(row_type row) const {
            return (normalized ? lengthInfo[row].id : row);
        }

        inline double cosine(row_type row, const double* query) {

            const double* d_ptr = getMatrixRowPtr(row);
            double cosine = 0;

#ifdef  WITH_SIMD

            //            __m128d sum = _mm_set1_pd(0.0);
            //
            //            for (int i = 0; i < sizeDiv2; i += 2) {
            //                __m128d xmm1 = _mm_loadu_pd(d_ptr + i);
            //                __m128d xmm2 = _mm_loadu_pd(query + i);
            //                sum = _mm_add_pd(sum, _mm_mul_pd(xmm1, xmm2));
            //            }
            //
            //            cosine = _mm_cvtsd_f64(_mm_hadd_pd(sum, sum));
            //
            //            cosine += extraMult * d_ptr[sizeDiv2] * query[sizeDiv2];
            //
            //            //            if(sizeDiv2 < colNum)
            //            //                cosine += d_ptr[sizeDiv2] * query[sizeDiv2];

            __m128d sum = _mm_set1_pd(0.0);

            int size = colNum + extraMult;

            for (int i = 0; i < size; i += 2) {
                //                __m128d xmm1 = _mm_load_pd(d_ptr + i);
                //                __m128d xmm2 = _mm_load_pd(query + i);
                sum = _mm_add_pd(sum, _mm_mul_pd(_mm_load_pd(d_ptr + i), _mm_load_pd(query + i)));
            }

            cosine = _mm_cvtsd_f64(_mm_hadd_pd(sum, sum));

            //            cosine += extraMult * d_ptr[sizeDiv2] * query[sizeDiv2];

            //            if(sizeDiv2 < colNum)
            //                cosine += d_ptr[sizeDiv2] * query[sizeDiv2];

            return cosine;
#else


            for (int i = 0; i < colNum; i++) {
                cosine += query[i] * d_ptr[i];
            }
            return cosine;

#endif
        }

        inline double L2Distance(row_type row, const double* query)const {

            const double* d_ptr = getMatrixRowPtr(row);

            double dist = 0;

            if (normalized) {

                for (int i = 0; i < colNum; i++) {
                    double value = query[i] * query[-1] - d_ptr[i] * d_ptr[-1]; // unnormalize
                    dist += value * value;
                }

                dist = sqrt(dist);
            } else {
                for (int i = 0; i < colNum; i++) {
                    dist += (query[i] - d_ptr[i]) * (query[i] - d_ptr[i]);
                }
                dist = sqrt(dist);

            }

            return dist;
        }
        // I assume non normalized case as needed in PCA trees

        inline double L2Distance2(row_type row, const double* query)const {
            const double* d_ptr = getMatrixRowPtr(row);

            double dist = 0;
            for (int i = 0; i < colNum; i++) {
                dist += (query[i] - d_ptr[i]) * (query[i] - d_ptr[i]);
            }
            return dist;
        }

        inline double innerProduct(row_type row, const double* query) {
            const double ip = query[-1] * getVectorLength(row) * cosine(row, query);
            return ip;
        }

        inline std::pair<bool, double> passesThreshold(row_type row, const double* query, double theta) {

            std::pair<bool, double> p;
            double ip = 1;

            if (normalized) {
                ip = query[-1] * getVectorLength(row);

                if (ip < theta) {
                    p.first = false;
                    return p;
                }
            }

            ip *= cosine(row, query);
            p.second = ip;

            if (ip < theta) {
                p.first = false;
                return p;
            } else {
                p.first = true;
                return p;
            }
        }


    };

    /*  map: id: original matrix id, first: thread second: posInMatrix
     */
    inline void initializeMatrices(VectorMatrix& originalMatrix, std::vector<VectorMatrix>& matrices, bool sort, bool ignoreLengths, double gamma = 0) {

        row_type threads = matrices.size();

        if (threads == 1) {

            matrices[0].rowNum = originalMatrix.rowNum;
            matrices[0].colNum = originalMatrix.colNum;
            matrices[0].offset = matrices[0].colNum + 2;
            matrices[0].sizeDiv2 = matrices[0].colNum & (-2);
            matrices[0].extraMult = (matrices[0].sizeDiv2 < matrices[0].colNum);
            if (matrices[0].extraMult)
                matrices[0].offset++;

            matrices[0].normalized = true;
            //            matrices[0].data = new double[(matrices[0].offset) * matrices[0].rowNum];
            int res = posix_memalign((void **) &(matrices[0].data), 16, sizeof (double)* matrices[0].offset * matrices[0].rowNum);
            //            std::cout << "res: " << res << " null? " << (matrices[0].data == NULL) << std::endl;

            if (res != 0) {
                std::cout << "Problem with allocating memory" << std::endl;
                exit(1);
            }

            matrices[0].lengthInfo.resize(matrices[0].rowNum);


            if (matrices[0].extraMult) {
                for (int i = 0; i < matrices[0].rowNum; i++) {
                    matrices[0].data[(i + 1) * matrices[0].offset - 3] = 0; // zero-out the last padding
                }
            }

            if (ignoreLengths) {

#ifdef ABS_APPROX
                matrices[0].gammaEquivalents.resize(matrices[0].rowNum, gamma);
#endif


                for (int i = 0; i < matrices[0].rowNum; i++) {
                    double len = 0;
                    const double* vec = originalMatrix.getMatrixRowPtr(i);
                    for (int j = 0; j < matrices[0].colNum; j++) {
                        len += vec[j] * vec[j];
                    }
                    len = sqrt(len);

                    matrices[0].lengthInfo[i] = QueueElement(1, i);
                    matrices[0].setLengthInData(i, 1);
                    double x = 1 / len;
#ifdef ABS_APPROX
                    matrices[0].gammaEquivalents[i] *= x;
                    //                    std::cout<<"len: "<<len<<" 1/len: "<<x<<" "<<matrices[0].gammaEquivalents[i]<<std::endl;
#endif

                    double * d1 = matrices[0].getMatrixRowPtr(i);
                    double * d2 = originalMatrix.getMatrixRowPtr(i);

                    for (int j = 0; j < originalMatrix.colNum; j++) {
                        d1[j] = d2[j] * x;
                    }

                }

            } else {
                for (int i = 0; i < matrices[0].rowNum; i++) {
                    double len = 0;
                    const double* vec = originalMatrix.getMatrixRowPtr(i);
                    for (int j = 0; j < matrices[0].colNum; j++) {
                        len += vec[j] * vec[j];
                    }
                    len = sqrt(len);

                    matrices[0].lengthInfo[i] = QueueElement(len, i);
                }

                if (sort) {
                    matrices[0].shuffled = true;
                    std::sort(matrices[0].lengthInfo.begin(), matrices[0].lengthInfo.end(), std::greater<QueueElement>());
                }

                for (int i = 0; i < matrices[0].rowNum; i++) {
                    matrices[0].setLengthInData(i, matrices[0].lengthInfo[i].data);

                    double x = 1 / matrices[0].lengthInfo[i].data;


                    double * d1 = matrices[0].getMatrixRowPtr(i);
                    double * d2 = originalMatrix.getMatrixRowPtr(matrices[0].lengthInfo[i].id);

                    for (int j = 0; j < originalMatrix.colNum; j++) {
                        d1[j] = d2[j] * x;
                    }

                }

            }

        } else { // multiple threads


            std::vector<row_type> permuteVector(originalMatrix.rowNum);
            for (int i = 0; i < permuteVector.size(); i++) {
                permuteVector[i] = i;
            }
            rg::Random32 random(123);
            rg::shuffle(permuteVector.begin(), permuteVector.end(), random);
            std::vector<row_type> blockOffsets;
            computeDefaultBlockOffsets(permuteVector.size(), threads, blockOffsets);


#pragma omp parallel
            {
                row_type tid = omp_get_thread_num();

                matrices[tid].colNum = originalMatrix.colNum;
                matrices[tid].offset = matrices[tid].colNum + 2;
                matrices[tid].sizeDiv2 = matrices[tid].colNum & (-2);
                matrices[tid].extraMult = (matrices[tid].sizeDiv2 < matrices[tid].colNum);
                if (matrices[tid].extraMult)
                    matrices[tid].offset++;

                row_type start = blockOffsets[tid];
                row_type end = (tid == blockOffsets.size() - 1 ? originalMatrix.rowNum : blockOffsets[tid + 1]);

                matrices[tid].rowNum = end - start;

                matrices[tid].normalized = true;
                matrices[tid].lengthInfo.resize(matrices[tid].rowNum);

                //                matrices[tid].data = new double[(matrices[tid].offset) * matrices[tid].rowNum];
                int res = posix_memalign((void **) &(matrices[tid].data), 16, sizeof (double)* matrices[tid].offset * matrices[tid].rowNum);
                //                std::cout << "res: " << res << " null? " << (matrices[tid].data == NULL) << std::endl;
                if (res != 0) {
                    std::cout << "Problem with allocating memory" << std::endl;
                    exit(1);
                }

                if (matrices[tid].extraMult) {
                    for (int i = 0; i < matrices[tid].rowNum; i++) {
                        matrices[tid].data[(i + 1) * matrices[tid].offset - 3] = 0; // zero-out the last padding
                    }
                }



                if (ignoreLengths) {

#ifdef ABS_APPROX
                    matrices[tid].gammaEquivalents.resize(matrices[tid].rowNum, gamma);
#endif

                    for (int i = start; i < end; i++) {
                        double len = 0;

                        row_type ind = permuteVector[i];

                        const double* vec = originalMatrix.getMatrixRowPtr(ind);
                        for (int j = 0; j < matrices[tid].colNum; j++) {
                            len += vec[j] * vec[j];
                        }
                        len = sqrt(len);


                        matrices[tid].lengthInfo[i - start] = QueueElement(1, i - start);


                        matrices[tid].setLengthInData(i - start, 1);

                        double x = 1 / len;
#ifdef ABS_APPROX
                        matrices[tid].gammaEquivalents[i] *= x;
#endif

                        double * d1 = matrices[tid].getMatrixRowPtr(i - start);
                        double * d2 = originalMatrix.getMatrixRowPtr(ind);

                        for (int j = 0; j < originalMatrix.colNum; j++) {
                            d1[j] = d2[j] * x;
                        }

                        matrices[tid].lengthInfo[i - start].id = ind; // the original id

                    }
                } else {
                    for (int i = start; i < end; i++) {
                        double len = 0;
                        row_type ind = permuteVector[i];
                        const double* vec = originalMatrix.getMatrixRowPtr(ind);
                        for (int j = 0; j < matrices[tid].colNum; j++) {
                            len += vec[j] * vec[j];
                        }
                        len = sqrt(len);

                        matrices[tid].lengthInfo[i - start] = QueueElement(len, i - start);
                    }

                    if (sort) {
                        matrices[tid].shuffled = true;
                        std::sort(matrices[tid].lengthInfo.begin(), matrices[tid].lengthInfo.end(), std::greater<QueueElement>());
                    }


                    for (int i = 0; i < matrices[tid].rowNum; i++) {
                        matrices[tid].setLengthInData(i, matrices[tid].lengthInfo[i].data);

                        double x = 1 / matrices[tid].lengthInfo[i].data;

                        row_type ind = permuteVector[matrices[tid].lengthInfo[i].id + start];

                        double * d1 = matrices[tid].getMatrixRowPtr(i);
                        double * d2 = originalMatrix.getMatrixRowPtr(ind);

                        for (int j = 0; j < originalMatrix.colNum; j++) {
                            d1[j] = d2[j] * x;
                        }


                        matrices[tid].lengthInfo[i].id = ind; // the original id
                    }
                }


            }

        }
    }

    inline void localToGlobalIds(const std::vector<std::vector<QueueElement> >& localResults, std::vector<MatItem>& globalResults,
            const VectorMatrix& userMatrix) {
        globalResults.clear();
        globalResults.reserve(localResults.size() * localResults[0].size()); // queries*k
        for (long i = 0; i < localResults.size(); i++) {
            row_type queryId = userMatrix.lengthInfo[i].id;
            for (long j = 0; j < localResults[i].size(); j++) {
                globalResults.push_back(MatItem(localResults[i][j].data, queryId, localResults[i][j].id));
            }
        }
    }

   
    
    inline void localToGlobalIds(const std::vector<QueueElement>& localResults, int k, std::vector<MatItem>& globalResults,
            const VectorMatrix& userMatrix) {
        globalResults.clear();
        globalResults.reserve(localResults.size()); // queries*k
        row_type user = 0;
        if (userMatrix.normalized) {
            for (long i = 0; i < localResults.size(); i += k) {
                row_type queryId = userMatrix.lengthInfo[user].id;
                for (int j = 0; j < k; j++) {
                    globalResults.push_back(MatItem(localResults[i + j].data, queryId, localResults[i + j].id));
                }
                user++;
            }
        } else {
            for (long i = 0; i < localResults.size(); i += k) {
                row_type queryId = user;
                for (int j = 0; j < k; j++) {
                    globalResults.push_back(MatItem(localResults[i + j].data, queryId, localResults[i + j].id));
                }
                user++;
            }
        }
    }

    void calculateAPneededForQuery(std::vector<VectorMatrix>& matrices, double thres, int k, std::vector<double>& global_cweights) {

        global_cweights.resize(matrices[0].colNum, 0);

#pragma omp parallel
        {

            row_type tid = omp_get_thread_num();
            col_type colNum = matrices[tid].colNum;
            row_type rowNum = matrices[tid].rowNum;
            row_type endUser = rowNum;

            if (k == 0) {
                std::vector<QueueElement>::const_iterator up = std::lower_bound(matrices[tid].lengthInfo.begin(), matrices[tid].lengthInfo.end(), QueueElement(thres, 0), std::greater<QueueElement>());
                endUser = up - matrices[tid].lengthInfo.begin();
            }

            matrices[tid].cweights.resize(colNum);
            matrices[tid].maxVectorCoord.resize(rowNum);
            matrices[tid].vectorNNZ.resize(rowNum, 0);

            for (int i = 0; i < endUser; i++) {

                double * d = matrices[tid].getMatrixRowPtr(i);

                for (int j = 0; j < colNum; j++) {

                    if (matrices[tid].cweights[j] < fabs(d[j]))
                        matrices[tid].cweights[j] = fabs(d[j]);

                    if (d[j] != 0)
                        matrices[tid].vectorNNZ[i]++;
                    if (fabs(d[j]) > matrices[tid].maxVectorCoord[i])
                        matrices[tid].maxVectorCoord[i] = fabs(d[j]);
                }


            }

#pragma omp critical
            {
                for (int i = 0; i < colNum; i++) {
                    if (global_cweights[i] < matrices[tid].cweights[i])
                        global_cweights[i] = matrices[tid].cweights[i];
                }
            }
        }
    }
}
#endif /* VECTORMATRIX_H_ */
