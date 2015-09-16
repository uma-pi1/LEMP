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
        return matrix;
    }

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
}





