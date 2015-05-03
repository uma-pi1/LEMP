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
 * read.h
 *
 *  Created on: Sep 13, 2012
 *      Author: chteflio
 */


#include <taLib/structs/read.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>




namespace ta {

    inline void skipLineFromFile(std::ifstream & file) {
        char c = 0;
        while (c != '\n' && !file.eof() && file.good()) {
            file >> std::noskipws >> c;
        }
        file >> std::skipws;
    }


    // assumes that 0 column of the matrix is spared for storing the length

    void normalizeMatrix(std::vector<std::vector<double> > & matrix) {
        assert(matrix.size() > 0);
        ta_size_type colnum = matrix[0].size();
        ta_size_type rownum = matrix.size();

        std::cout << "normalizing... " << std::endl;
        for (ta_size_type j = 0; j < rownum; j++) {
            double length = 0.0;
            for (ta_size_type i = 1; i < colnum; i++) {
                length += matrix[j][i] * matrix[j][i];
            }
            length = sqrt(length);
            matrix[j][0] = length; //store the length on the 0-column
            for (ta_size_type i = 1; i < colnum; i++) {
                matrix[j][i] /= length;
            }
        }
    }

    std::vector<std::vector<double> > readMatrix(std::string& fileName, bool left) {

        std::vector<std::vector<double> > matrix;
        matrix.clear();
        std::ifstream file(fileName.c_str(), std::ios_base::in);
        while (file.peek() == '%') {
            skipLineFromFile(file);
        }

        ta_size_type colnum; // columns
        ta_size_type rownum; // rows
        file >> rownum >> colnum;
#ifdef DEBUG
        std::cout << "Reading file: " << fileName << " -- columns: " << colnum << "  rows: " << rownum << std::endl;
#endif
        // left matrix:  vector of lines
        // right matrix: vector of columns
        matrix.resize(left ? rownum : colnum);
        if (file) {
            for (ta_size_type i = 0; i < colnum; i++) {// read one column
                for (ta_size_type j = 0; j < rownum; j++) {
                    double f;
                    file >> f;
                    //				if(fabs(f)<=1E-6){///////////////////
                    //					f = 0;
                    //				}
                    matrix[left ? j : i].push_back(f);
                }
            }
        }
        file.close();
        //	std::cout << "Done reading matrix " << fileName << std::endl;



        return matrix;
    }

//    std::vector< std::vector<row_type> > readClusters(std::string& fileName, int numOfClusters, row_type numOfItems) {
//        std::ifstream file(fileName.c_str(), std::ios_base::in);
//
//        std::vector<std::vector<row_type> > clusters;
//        clusters.resize(numOfClusters);
//
//
//        int c;
//
//        if (file) {
//
//            for (int i = 0; i < numOfItems; i++) {
//                file >> c;
//                clusters[c].push_back(i);
//            }
//
//
//
//        }
//        file.close();
//
//        return clusters;
//    }
    
        void readClusters(std::string& fileName, std::vector<int> & clusterIds, row_type numOfItems) {
        std::ifstream file(fileName.c_str(), std::ios_base::in);
        
        clusterIds.reserve(numOfItems);

        int c;
        
         std::cout << "Reading clustering file: " << fileName << std::endl;

        if (file) {

            for (int i = 0; i < numOfItems; i++) {
                file >> c;
                clusterIds.push_back(c);
            }
        }
        file.close();

    }
    
    
    //void writeResults(std::string& fileName, std::vector<MatItem>& results){
    //	std::sort (results.begin(), results.end(), std::greater<MatItem>());
    //
    //	// open file
    //	std::ofstream out(fileName.c_str());
    //	if (!out.is_open())
    //		RG_THROW(rg::IOException, std::string("Cannot open file ") + fileName);
    //
    //	for (long i=0; i<results.size(); i++){
    //		out<<results[i].result<<"\t"<<results[i].i<<"\t"<<results[i].j<<std::endl;
    //	}
    //	// done
    //	out.close();
    //}

    inline void writeResults(std::vector< std::vector<MatItem >* >& globalResults, std::string& file) {

        std::ofstream out(file.c_str());
        
         if (!out.is_open())
             std::cout<<"file "<<file<<" not open"<<std::endl;
         else
            std::cout<<"file "<<file<<" open for writing"<<std::endl;

        for (int t = 0; t < globalResults.size(); t++) {//each thread

            std::sort(globalResults[t]->begin(), globalResults[t]->end(), std::less<MatItem>());

            for (int j = 0; j < globalResults[t]->size(); j++) {// each result            
                out << globalResults[t]->at(j).i << " " << globalResults[t]->at(j).j << " " << globalResults[t]->at(j).result << std::endl;
            }

        }
        out.close();
    }
    
    

    void readAndCompareResults(std::string& fileName1, std::string& fileName2) {
        std::ifstream file1(fileName1.c_str(), std::ios_base::in);
        std::ifstream file2(fileName2.c_str(), std::ios_base::in);

        std::string line1, line2;
        std::vector<std::string> res1, res2;

        if (!file1.is_open())
            std::cout<<"file "<<fileName1<<" not open"<<std::endl;
        if (!file2.is_open())
            std::cout<<"file "<<fileName2<<" not open"<<std::endl;
        

        while (std::getline(file1, line1) && std::getline(file2, line2)) {

            res1.push_back(line1);
            res2.push_back(line2);
        }

        file1.close();
        file2.close();

        std::cout << "res1 size: " << res1.size() << std::endl;
        std::cout << "res2 size: " << res2.size() << std::endl;

        row_type nonSame = 0;
        if (res1.size() != res2.size()) {
            std::cout << "of non equal size" << std::endl;
            exit(1);
        }


        for (int i = 0; i < res1.size(); i++) {
            bool thisOk = false;

            for (int j = 0; j < res2.size(); j++) {
                if (res1[i] == res2[j]) {
                    thisOk = true;
                    break;
                }
            }
            if (!thisOk) {
                nonSame++;
            }
        }

        std::cout << "nonSame: " << nonSame << " out of " << res1.size() << std::endl;

    }

    bool compareResults(std::string& fileName1, std::string& fileName2) {
        std::ifstream file1(fileName1.c_str(), std::ios_base::in);
        std::ifstream file2(fileName2.c_str(), std::ios_base::in);

        std::string line1, line2;
        ta_size_type i = 0;
        while (std::getline(file1, line1) && std::getline(file2, line2)) {
            i++;
            if (line1 != line2) {
                file1.close();
                file2.close();
                std::cout << "The results are NOT EQUAL" << std::endl;
                std::cout << "First inconsistency in line " << i << std::endl;
                return false;
            }

        }

        file1.close();
        file2.close();
        std::cout << "The results are EQUAL" << std::endl;
        return true;
    }








}





