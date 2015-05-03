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
        col_type* queues;
        col_type maxLists;
        row_type sizeOfQueue;

        bool shuffled, normalized, extraMult;

        std::vector<double> cweights;
        std::vector<double> maxVectorCoord;
        std::vector<row_type> vectorNNZ;
        std::vector<QueueElement> colFrequencies;
        row_type colNum, offset;
        row_type rowNum;
        std::vector<QueueElement> lengthInfo; // data: length id: vectorId

        std::vector<double> clusterTheta;
        std::vector<row_type> id2Pos;

        // for simd instruction  
        int sizeDiv2;

        inline void initQueues(col_type _maxLists, double thres = 0) {
            maxLists = _maxLists;

            if (thres == 0) {
                queues = new col_type[rowNum * maxLists];
                sizeOfQueue = rowNum;
            } else {
                std::vector<QueueElement>::const_iterator up = std::lower_bound(lengthInfo.begin(), lengthInfo.end(), QueueElement(thres, 0), std::greater<QueueElement>());
                row_type activeQueries = up - lengthInfo.begin();

                sizeOfQueue = activeQueries;
                queues = new col_type[activeQueries * maxLists];
            }
        }

        inline void preprocessQueues(row_type start, row_type end) {
            std::vector<QueueElement> tmp(maxLists);

            end = (end > sizeOfQueue ? sizeOfQueue : end);

            for (row_type j = start; j < end; j++) {
                const double* query = getMatrixRowPtr(j);

                for (col_type i = 0; i < maxLists; i++) {
                    double value = fabs(query[i]);
                    tmp[i] = QueueElement(value, i);
                }
                std::make_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());

                for (col_type i = maxLists; i < colNum; i++) {
                    double value = fabs(query[i]);
                    if (value > tmp.front().data) {
                        std::pop_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                        tmp.pop_back();
                        tmp.push_back(QueueElement(value, i));
                        std::push_heap(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                    }
                }

                std::sort(tmp.begin(), tmp.end(), std::greater<QueueElement>());
                for (col_type i = 0; i < maxLists; i++) {
                    queues[j * maxLists + i] = tmp[i].id;
                }
            }
        }

        inline col_type* getQueue(row_type user) const {
            return &queues[user * maxLists];
        }








        //	inline friend void mapToItemIds(std::vector<MatItem>& results, const VectorMatrix& first, const VectorMatrix& second);
        //	inline friend void mapToItemIds(std::vector<MatItem>& results, const VectorMatrix& itemMatrix);
        inline friend void writeResults(std::vector<std::vector<QueueElement> >& results, const VectorMatrix& userMatrix, const VectorMatrix& itemMatrix, std::string& fileName);

        inline VectorMatrix() : data(0), queues(0), shuffled(false), normalized(false) {
        }

        inline ~VectorMatrix() {
            delete[] data;
            if (!queues)
                delete[] queues;
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
	    
	    if(left){
		if(pow(2, sizeof(col_type) * 8)-1 < col){
		    std::cerr<<"Your vectors have dimensionality "<<col<<" which is more than what lemp is compiled to store. Change the col_type in BasicStructs.h and recompile."<<std::endl;
		    exit(1);
		}
		
		if(pow(2, sizeof(row_type) * 8)-1 < row){
		    std::cerr<<"Your dataset has "<<row<<" vectors which is more than what lemp is compiled to store. Change the row_type in BasicStructs.h and recompile."<<std::endl;
		    exit(1);
		}
	      
	    }else{
		if(pow(2, sizeof(col_type) * 8)-1 < row){
		    std::cerr<<"Your vectors have dimensionality "<<row<<" which is more than what lemp is compiled to store. Change the col_type in BasicStructs.h and recompile."<<std::endl;
		    exit(1);
		}
		
		if(pow(2, sizeof(row_type) * 8)-1 < col){
		    std::cerr<<"Your dataset has "<<col<<" vectors which is more than what lemp is compiled to store. Change the row_type in BasicStructs.h and recompile."<<std::endl;
		    exit(1);
		}	      
	    }
	    
	    rowNum = (left ? row : col);
            colNum = (left ? col : row);
            offset = colNum + 1;
            sizeDiv2 = colNum & (-2);
	    
	    if(colNum < NUM_LISTS){
		std::cout<<"WARNING: Your vectors have dimensionality"<<colNum<<" and the tuner will try to search among "<<NUM_LISTS<<
		". Perhaps you want to change the parameter NUM_LISTS in Definitions.h and recompile"<<std::endl;	      
	    }
	    if(LOWER_LIMIT_PER_BUCKET >= rowNum){
	        std::cout<<"WARNING: You have "<<rowNum<<" vectors and the tuner will try to take a sample of at least  "<<LOWER_LIMIT_PER_BUCKET<<
		" vectors per probe bucket. Perhaps you want to change the parameter LOWER_LIMIT_PER_BUCKET in Definitions.h and recompile"<<std::endl;	      
	    }
            
            
            extraMult = (sizeDiv2 < colNum);
                
                

            //std::cout<<"rowNum: "<<rowNum<<" colNum: "<<(int)colNum<<std::endl;
            normalized = false;
            data = new double[offset * rowNum];

            for (int i = 0; i < rowNum; i++) {
                data[i * offset] = 1;
            }

            if (left) {
                if (file) {
                    for (ta_size_type i = 0; i < col; i++) {// read one column
                        for (ta_size_type j = 0; j < row; j++) {
                            double f;
                            file >> f;
                            data[j * offset + i + 1] = f;
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
                            data[i * offset + j + 1] = f;
                        }
                    }
                }
                file.close();
            }
        }

        inline void init(std::vector<std::vector<double> >& matrix, bool norm = false, bool sort = false) {
            rowNum = matrix.size();
            colNum = matrix[0].size();
            offset = colNum + 1;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            normalized = norm;

            data = new double[offset * rowNum];

            if (normalized) {
                lengthInfo.resize(rowNum);
#pragma omp parallel
                {

                    // get lengths
#pragma omp for schedule(static,1000)
                    for (int i = 0; i < rowNum; i++) {
                        double len = 0;
                        for (int j = 0; j < colNum; j++) {
                            len += matrix[i][j] * matrix[i][j];
                        }
                        len = sqrt(len);

                        lengthInfo[i] = QueueElement(len, i);
                    }

#pragma omp single
                    {
                        if (sort) {
                            shuffled = true;
                            std::sort(lengthInfo.begin(), lengthInfo.end(), std::greater<QueueElement>());
                        }
                    }


#pragma omp for nowait schedule(static,1000)
                    for (int i = 0; i < rowNum; i++) {
                        data[i * offset] = lengthInfo[i].data;
                        double x = 1 / lengthInfo[i].data;
                        for (int j = 1; j < offset; j++) {
                            data[i * offset + j] = matrix[lengthInfo[i].id][j - 1] * x;
                        }
                    }
                }


            } else {
#pragma omp parallel for schedule(static,1000)
                for (int i = 0; i < rowNum; i++) {
                    data[i * offset] = 1;
                    for (int j = 1; j < offset; j++) {
                        data[i * offset + j] = matrix[i][j - 1];
                    }
                }
            }
        }

        inline void init(VectorMatrix& matrix, bool sort = false) {
            rowNum = matrix.rowNum;
            colNum = matrix.colNum;
            offset = colNum + 1;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            normalized = true;

            data = new double[offset * rowNum];


            lengthInfo.resize(rowNum);

            //get lengths     

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
                data[i * offset] = lengthInfo[i].data;
                double x = 1 / lengthInfo[i].data;

                for (int j = 1; j < offset; j++) {
                    data[i * offset+ j] = matrix.data[lengthInfo[i].id * offset + j] * x;
                }
            }


        }

        void calculateColFrequencies(row_type end) {

            colFrequencies.resize(colNum);

            for (int i = 0; i < end; i++) {
                for (int j = 1; j < offset; j++) {

                    if (data[i * offset + j] != 0) {
                        colFrequencies[j - 1].id = j - 1;
                        colFrequencies[j - 1].data++;
                    }
                }
            }
            std::sort(colFrequencies.begin(), colFrequencies.end(), std::greater<QueueElement>());
        }

        void calculateAPneededForQuery(row_type end) {
            cweights.resize(colNum);
            maxVectorCoord.resize(rowNum);
            vectorNNZ.resize(rowNum, 0);

            for (int i = 0; i < end; i++) {

                for (int j = 1; j < offset; j++) {

                    if (cweights[j - 1] < fabs(data[i * offset + j]))
                        cweights[j - 1] = fabs(data[i * offset + j]);

                    if (data[i * offset + j] != 0)
                        vectorNNZ[i]++;
                    if (fabs(data[i * offset + j]) > maxVectorCoord[i])
                        maxVectorCoord[i] = fabs(data[i * offset + j]);


                }
            }
        }

        inline void initForTopKPerUser(std::vector<std::vector<double> >& matrix) {
            rowNum = matrix.size();
            colNum = matrix[0].size();
            offset = colNum + 1;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            normalized = true;

            data = new double[offset * rowNum];

            lengthInfo.resize(rowNum);

            // get lengths
            for (int i = 0; i < rowNum; i++) {
                double len = 0;
                for (int j = 0; j < colNum; j++) {
                    len += matrix[i][j] * matrix[i][j];
                }
                len = sqrt(len);

                lengthInfo[i] = QueueElement(len, i);
            }



            for (int i = 0; i < rowNum; i++) {
                data[i * offset] = 1;
                double x = 1 / lengthInfo[i].data;
                for (int j = 1; j < offset; j++) {
                    data[i * offset + j] = matrix[lengthInfo[i].id][j - 1] * x;
                }
                lengthInfo[i].data = 1; // no need to make all 1
            }
        }

        inline void initForTopKPerUser(VectorMatrix& matrix) {
            rowNum = matrix.rowNum;
            colNum = matrix.colNum;
            offset = colNum + 1;
            sizeDiv2 = colNum & (-2);
            extraMult = (sizeDiv2 < colNum);
            normalized = true;

            data = new double[offset * rowNum];

            lengthInfo.resize(rowNum);


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

                data[i * offset] = 1;
                double x = 1 / len;
                for (int j = 1; j < offset; j++) {
                    data[i * offset + j] = matrix.data[i * offset + j] * x;

                }
            }

        }


        // please only assign to const double *

        inline double* getMatrixRowPtr(row_type row) const {// the row starts from pos 1. Do ptr[-1] to get the length
            return &data[row * offset + 1];
        }

        inline double getVectorLength(row_type row) const {
            return data[row * offset];
        }

        inline double calcNormL2(row_type row) {
            const double* d_ptr = getMatrixRowPtr(row);
            double len = 0;
            for (int i = 0; i < colNum; i++) {
                len += d_ptr[i] * d_ptr[i];
            }
            return sqrt(len);
        }

        inline row_type getId(row_type row) const {
            return (normalized ? lengthInfo[row].id : row);
        }

        //         inline double cosine(row_type row, const double* query){
        // 
        //             const double* d_ptr = getMatrixRowPtr(row);
        //             double cosine = 0;
        //             for (int i = 0; i < colNum; i++) {
        //                 cosine += query[i] * d_ptr[i];
        //             }
        //             return cosine;
        //         }

        inline double cosine(row_type row, const double* query) {

            const double* d_ptr = getMatrixRowPtr(row);
            double cosine = 0;

#ifdef  WITH_SIMD
            
            __m128d sum = _mm_set1_pd(0.0); 
           
            for (int i = 0; i < sizeDiv2; i += 2) {
                __m128d xmm1 = _mm_loadu_pd(d_ptr + i);
                __m128d xmm2 = _mm_loadu_pd(query + i);
                sum = _mm_add_pd(sum, _mm_mul_pd(xmm1, xmm2));
            }

            cosine = _mm_cvtsd_f64(_mm_hadd_pd(sum, sum));            
            
            cosine += extraMult * d_ptr[sizeDiv2] * query[sizeDiv2];
            
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

        inline void printVector(row_type row)const {
            std::cout << "L2norm: " << data[row * (colNum + 1)] << " |";
            for (int j = 1; j < colNum + 1; j++) {

                if (data[row * (colNum + 1) + j] < 0.001 && data[row * (colNum + 1) + j]>-0.001)
                    std::cout << 0 << " ";
                else
                    std::cout << data[row * (colNum + 1) + j] << " ";
            }
            std::cout << std::endl;
        }

        inline double getSumPi2(row_type row, std::vector<col_type>& coord) {
            const double* vec = getMatrixRowPtr(row);
            double sum = 0;

            for (col_type i = 0; i < coord.size(); i++) {
                sum += vec[coord[i]] * vec[coord[i]];
            }
            return sum;
        }



    };


    //inline void mapToItemIds(std::vector<MatItem>& results, const VectorMatrix& first, const VectorMatrix& second){
    //	if (first.normalized || second.normalized){
    //		for(long i=0; i<results.size(); i++){
    //			if (first.normalized)
    //				results[i].i = first.lengthInfo[results[i].i].id;
    //			if (second.normalized)
    //				results[i].j = second.lengthInfo[results[i].j].id;
    //		}
    //	}//if not normalized nothing to be done here
    //}
    //
    ///*
    // * Assumes that MatItem.i is the queryId already
    // */
    //inline void mapToItemIds(std::vector<MatItem>& results, const VectorMatrix& itemMatrix){
    //	if (itemMatrix.normalized){
    //		for(long i=0; i<results.size(); i++){
    //			results[i].j = itemMatrix.lengthInfo[results[i].j].id;
    //		}
    //	}//if not normalized nothing to be done here
    //}
    //
    //inline void mapToItemIds(std::vector<std::vector<MatItem> >& results, const VectorMatrix& itemMatrix){
    //	if (itemMatrix.normalized){
    //		for(long i=0; i<results.size(); i++){
    //			for(long j=0; j<results[i].size(); j++){
    //				results[i][j].j = itemMatrix.lengthInfo[results[i][j].j].id;
    //			}
    //
    //		}
    //	}//if not normalized nothing to be done here
    //}

    inline void writeResults(std::vector<std::vector<QueueElement> >& results, const VectorMatrix& userMatrix, const VectorMatrix& itemMatrix, std::string& fileName) {


        std::cout << "Writing results to: " << fileName << std::endl;

        std::ofstream out(fileName.c_str());
        if (!out.is_open())
            RG_THROW(rg::IOException, std::string("Cannot open file ") + fileName);


        if (!userMatrix.shuffled && !itemMatrix.shuffled) {// ids correspond to original ids
            for (long i = 0; i < results.size(); i++) {
                out << i << ": ";
                for (long j = 0; j < results[i].size(); j++) {
                    out << "(" << results[i][j].data << ", " << results[i][j].id << ")\t";
                }
                out << std::endl;
            }
        } else if (userMatrix.shuffled && !itemMatrix.shuffled) {// topk case
            for (long i = 0; i < results.size(); i++) {
                out << i << ": ";
                for (long j = 0; j < results[i].size(); j++) {
                    out << "(" << results[i][j].data << ", " << itemMatrix.lengthInfo[results[i][j].id].id << ")\t";
                }
                out << std::endl;
            }
        } else if (userMatrix.shuffled && itemMatrix.shuffled) { //above-theta case
            for (long i = 0; i < results.size(); i++) {
                out << userMatrix.lengthInfo[i].id << ": ";
                for (long j = 0; j < results[i].size(); j++) {
                    out << "(" << results[i][j].data << ", " << itemMatrix.lengthInfo[results[i][j].id].id << ")\t";
                }
                out << std::endl;
            }
        }
        // done
        out.close();

        std::cout << "Done with writing results" << std::endl;
    }

    // for Above-theta

    inline void mapToGlobalIds(const std::vector<std::vector<QueueElement> >& localResults, int totalSize, std::vector<MatItem>& globalResults,
            const VectorMatrix& userMatrix, const VectorMatrix& itemMatrix) {
        globalResults.clear();
        globalResults.reserve(totalSize);
        for (long i = 0; i < localResults.size(); i++) {
            for (long j = 0; j < localResults[i].size(); j++) {
                globalResults.push_back(MatItem(localResults[i][j].data, userMatrix.lengthInfo[i].id, itemMatrix.lengthInfo[localResults[i][j].id].id));
            }
        }
    }

    // for RowTopk

    inline void mapToGlobalIds(std::vector<std::vector<QueueElement> >& localResults, const VectorMatrix& itemMatrix) {

        // first update itemIDs and sort in place
#pragma omp  parallel for schedule(static,100)
        for (long i = 0; i < localResults.size(); i++) {

            for (long j = 0; j < localResults[i].size(); j++) {
                localResults[i][j].id = itemMatrix.lengthInfo[localResults[i][j].id].id;
            }

            std::sort(localResults[i].begin(), localResults[i].end(), std::greater<QueueElement>());
        }
    }

    /*  map: id: original matrix id, first: thread second: posInMatrix
     */
    inline void initializeMatrices(VectorMatrix& originalMatrix, std::vector<VectorMatrix>& matrices, bool sort = false) {

        row_type threads = matrices.size();

        if (threads == 1) {

            matrices[0].rowNum = originalMatrix.rowNum;
            matrices[0].colNum = originalMatrix.colNum;
            matrices[0].offset = matrices[0].colNum + 1;
            matrices[0].sizeDiv2 = matrices[0].colNum & (-2);
            matrices[0].extraMult = (matrices[0].sizeDiv2 < matrices[0].colNum);
            matrices[0].normalized = true;
            matrices[0].data = new double[(matrices[0].offset) * matrices[0].rowNum];
            matrices[0].lengthInfo.resize(matrices[0].rowNum);

            //std::cout << matrices[0].lengthInfo.size() << std::endl;


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
                matrices[0].data[i * (matrices[0].offset)] = matrices[0].lengthInfo[i].data;
                double x = 1 / matrices[0].lengthInfo[i].data;

                for (int j = 1; j < matrices[0].offset; j++) {
                    matrices[0].data[i * (matrices[0].offset) + j] = originalMatrix.data[matrices[0].lengthInfo[i].id * (matrices[0].offset) + j] * x;
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


            //std::cout << "ok with blockoffsets: " << blockOffsets.size() << std::endl;
#pragma omp parallel
            {


                row_type tid = omp_get_thread_num();

                matrices[tid].colNum = originalMatrix.colNum;
                matrices[tid].offset = matrices[tid].colNum + 1;
                matrices[tid].sizeDiv2 = matrices[tid].colNum & (-2);
                matrices[tid].extraMult = (matrices[tid].sizeDiv2 < matrices[tid].colNum);
                row_type start = blockOffsets[tid];
                row_type end = (tid == blockOffsets.size() - 1 ? originalMatrix.rowNum : blockOffsets[tid + 1]);

                matrices[tid].rowNum = end - start;

                matrices[tid].normalized = true;
                matrices[tid].lengthInfo.resize(matrices[tid].rowNum);

                matrices[tid].data = new double[(matrices[tid].offset) * matrices[tid].rowNum];



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
                    matrices[tid].data[i * (matrices[tid].offset)] = matrices[tid].lengthInfo[i].data;
                    double x = 1 / matrices[tid].lengthInfo[i].data;

                    row_type ind = permuteVector[matrices[tid].lengthInfo[i].id + start];


                    for (int j = 1; j < matrices[tid].offset; j++) {
                        matrices[tid].data[i * (matrices[tid].offset) + j] = originalMatrix.data[ind * (matrices[tid].offset) + j] * x;

                    }
                    matrices[tid].lengthInfo[i].id = ind; // the original id
                }

            }

        }


    }

    inline void initializeMatricesForTopKPerUser(VectorMatrix& originalMatrix, std::vector<VectorMatrix>& matrices) {

        row_type threads = matrices.size();

        if (threads == 1) {
            matrices[0].rowNum = originalMatrix.rowNum;
            matrices[0].colNum = originalMatrix.colNum;
            matrices[0].offset = matrices[0].colNum + 1;
            matrices[0].sizeDiv2 = matrices[0].colNum & (-2);
            matrices[0].extraMult = (matrices[0].sizeDiv2 < matrices[0].colNum);
            matrices[0].normalized = true;
            matrices[0].data = new double[(matrices[0].offset) * matrices[0].rowNum];
            matrices[0].lengthInfo.resize(matrices[0].rowNum);

            // get lengths
            for (int i = 0; i < matrices[0].rowNum; i++) {
                double len = 0;
                const double* vec = originalMatrix.getMatrixRowPtr(i);
                for (int j = 0; j < matrices[0].colNum; j++) {
                    len += vec[j] * vec[j];
                }
                len = sqrt(len);

                matrices[0].lengthInfo[i] = QueueElement(1, i);

                matrices[0].data[i * (matrices[0].offset)] = 1;
                double x = 1 / len;
                for (int j = 1; j < matrices[0].offset; j++) {
                    matrices[0].data[i * (matrices[0].offset) + j] = originalMatrix.data[i * (originalMatrix.offset) + j] * x;
                }
            }



        } else {

            std::vector<row_type> permuteVector(originalMatrix.rowNum);
            for (int i = 0; i < permuteVector.size(); i++) {
                permuteVector[i] = i;
            }

            rg::Random32 random(123);

            rg::shuffle(permuteVector.begin(), permuteVector.end(), random);

            std::vector<row_type> blockOffsets;

            computeDefaultBlockOffsets(permuteVector.size(), threads, blockOffsets);


            //std::cout << "ok with blockoffsets: " << blockOffsets.size() << std::endl;
#pragma omp parallel
            {


                row_type tid = omp_get_thread_num();

                matrices[tid].colNum = originalMatrix.colNum;
                matrices[tid].offset = matrices[tid].colNum + 1;
                matrices[tid].sizeDiv2 = matrices[tid].colNum & (-2);
                matrices[tid].extraMult = (matrices[tid].sizeDiv2 < matrices[tid].colNum);
                row_type start = blockOffsets[tid];
                row_type end = (tid == blockOffsets.size() - 1 ? originalMatrix.rowNum : blockOffsets[tid + 1]);

                matrices[tid].rowNum = end - start;

                matrices[tid].normalized = true;
                matrices[tid].lengthInfo.resize(matrices[tid].rowNum);

                matrices[tid].data = new double[(matrices[tid].offset) * matrices[tid].rowNum];

                for (int i = start; i < end; i++) {
                    double len = 0;

                    row_type ind = permuteVector[i];

                    const double* vec = originalMatrix.getMatrixRowPtr(ind);
                    for (int j = 0; j < matrices[tid].colNum; j++) {
                        len += vec[j] * vec[j];
                    }
                    len = sqrt(len);


                    matrices[tid].lengthInfo[i - start] = QueueElement(1, i - start);

                    matrices[tid].data[(i - start) * (matrices[tid].offset)] = 1;

                    double x = 1 / len;
                    for (int j = 1; j < matrices[tid].offset; j++) {
                        matrices[tid].data[(i - start) * (matrices[tid].offset) + j] = originalMatrix.data[ind * (originalMatrix.offset) + j] * x;

                    }

                    matrices[tid].lengthInfo[i - start].id = ind; // the original id

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

        //        for (long i=0; i<userMatrix.rowNum; i++){
        //            row_type queryId = userMatrix.lengthInfo[i].id;
        //            
        //            int offset = i*k;
        //            for(int j = 0; j<k; j++){
        //                globalResults.push_back(MatItem(localResults[offset+j].data, queryId, localResults[offset+j].id));
        //            
        //            }
        //        }

        row_type user = 0;

//        double global = 0;

        if (userMatrix.normalized) {
            for (long i = 0; i < localResults.size(); i += k) {
                row_type queryId = userMatrix.lengthInfo[user].id;
                
//                global += localResults[i].data;

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
        
//        std::cout<<"global: "<<global/user<<std::endl;

    }

    void calculateAPneededForQuery(std::vector<VectorMatrix>& matrices, double thres, int k, std::vector<double>& global_cweights) {

        global_cweights.resize(matrices[0].colNum, 0);

#pragma omp parallel
        {

            row_type tid = omp_get_thread_num();
            col_type colNum = matrices[tid].colNum;
            col_type offset = colNum + 1;
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

                for (int j = 1; j < offset; j++) {

                    if (matrices[tid].cweights[j - 1] < fabs(matrices[tid].data[i * offset + j]))
                        matrices[tid].cweights[j - 1] = fabs(matrices[tid].data[i * offset + j]);

                    if (matrices[tid].data[i * offset + j] != 0)
                        matrices[tid].vectorNNZ[i]++;
                    if (fabs(matrices[tid].data[i * offset + j]) > matrices[tid].maxVectorCoord[i])
                        matrices[tid].maxVectorCoord[i] = fabs(matrices[tid].data[i * offset + j]);


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
