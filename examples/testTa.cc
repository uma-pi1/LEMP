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

#include <iostream>
#include <mips/mips.h>

#include <sys/time.h>
#include <fstream>
#include <climits>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace mips;
using namespace rg;

// EXAMPLE OF TA
int main() {

    // first we have to define the input arguments for the algorithm
    InputArguments args;
    args.theta = 0.035; // theta for the Above-theta problem (will be used when running .runAboveTheta())
//    args.k = 10; // k for the Row-Top-k problem (will be used when running .runTopK())
    args.threads = 1; // number of threads the algorithm will use for execution
    args.logFile = "../../results/log.txt"; // log file which will contain the timings

    // TA will also need some special arguments
    bool isTARR = true; // true will use ROUND ROBIN, false will use a MAX HEAP
    
    // we will need some input matrices for the algorithm (create them like (1), (2) or (3))
    VectorMatrix leftMatrix, rightMatrix;

    // (1) create a VectorMatrix by using .mma files
    string queryFile = "../../dataset/w-gnmf-sample.mma";
    string probeFile = "../../dataset/h-gnmf-sample.mma";
    leftMatrix.readFromFile(queryFile, 0,0, true);
    rightMatrix.readFromFile(probeFile, 0,0, false);

    // (2) create a VectorMatrix by using .csv files
//    int r = 50; // dimensionality of vectors
//    int m = 1000; // number of query-vectors
//    int n = 1000; // number of probe-vectors
//    string queryFile = "../../dataset/w-gnmf-sample.csv";
//    string probeFile = "../../dataset/h-gnmf-sample.csv";
//    leftMatrix.readFromFile(queryFile, r, m, true);
//    rightMatrix.readFromFile(probeFile, r, n, false);
    
    // (3) create a VectorMatrix by using several vectors of double values 
    std::vector<std::vector<double>> m(1); // a matrix containing only 1 vector
    m[0].resize(50); // this vector will have dimensionality 50
    std::iota(m[0].begin(), m[0].end(), 0); // we will just fill in some values as an example
    VectorMatrix leftMatrixStreaming(m);

    // create an instance of the algorithm you wish to use
    mips::Ta algo(args, isTARR);

    // initialize the algorithm with a probe matrix (algorithm will store a copy of it inside)
    algo.initialize(rightMatrix);

    // create a structure for storing the results
    Results results;
    
    // decide between the Row-Top-k problem or the Above-theta problem
    if (args.k > 0) { // Row-Top-k
        
        // solves the Row-Top-k problem and stores the results in the specified results structure 
        algo.runTopK(leftMatrix, results);
        
    } else { // Above-theta
         
        // solves the Above-theta problem and stores the results in the specified results structure
        algo.runAboveTheta(leftMatrix, results);

//        algo.setTheta(1.68); // we might also change the theta if we want so
        
        // we might also run another matrix which will append the results to the specified results structure (otherwise you will have to create another results instance)
        algo.runAboveTheta(leftMatrixStreaming, results);
        
    }
    
    // write the results to a file
    string resultsFile = "../../results/output.data";
    results.writeToFile(resultsFile);

    return 0;

}
