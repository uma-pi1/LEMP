/*
 * runAlgoWithTuner.cc
 *
 *  Created on: Feb 10, 2014
 *      Author: chteflio
 */
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


//#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <mips/mips.h>


using namespace std;
using namespace mips;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    double theta, R, epsilon;
    string usersFile;
    string itemsFile;
    string logFile, resultsFile;

    bool querySideLeft = true;
    bool isTARR = true;
    int k, cacheSizeinKB, threads, r, m, n;
    std::string methodStr;
    LEMP_Method method;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&usersFile), "file containing the query matrix (left side)")
            ("P", value<string>(&itemsFile), "file containing the probe matrix (right side)")
            ("theta", value<double>(&theta), "theta value")
            ("R", value<double>(&R)->default_value(0.97), "recall parameter for LSH")
	    ("epsilon", value<double>(&epsilon)->default_value(0.0), "epsilon value for LEMP-LI with Absolute or Relative Approximation")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("isTARR", value<bool>(&isTARR)->default_value(true), "for LEMP-TA. If 1 Round Robin schedule is used (default). Otherwise Max PiQi")
            ("method", value<string>(&methodStr), "LEMP_X where X: L, LI, LC, I, C, TA, TREE, AP, LSH")
            ("k", value<int>(&k)->default_value(0), "top k (default 0). If 0 Above-theta will run")
            ("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
	    ("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
            ("cacheSizeinKB", value<int>(&cacheSizeinKB)->default_value(8192), "cache size in KB")
            ("t", value<int>(&threads)->default_value(1), "num of threads (default 1)")
            ("r", value<int>(&r)->default_value(0), "num of coordinates in each vector (needed when reading from csv files)")
            ("m", value<int>(&m)->default_value(0), "num of vectors in Q^T (needed when reading from csv files)")
            ("n", value<int>(&n)->default_value(0), "num of vectors in P (needed when reading from csv files)")
            ;

    positional_options_description pdesc;
    pdesc.add("Q^T", 1);
    pdesc.add("P", 2);

    variables_map vm;
    store(command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
    notify(vm);

    if (vm.count("help") || vm.count("Q^T") == 0 || vm.count("P") == 0) {
        cout << "runLemp [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }
        
    InputArguments args;
    args.logFile = logFile;
    args.theta = theta;
    args.k = k;
    args.threads = threads;

    if (methodStr.compare("LEMP_LI") == 0) {
        method = LEMP_LI;
    } else if (methodStr.compare("LEMP_LC") == 0) {
        method = LEMP_LC;
    } else if (methodStr.compare("LEMP_L") == 0) {
        method = LEMP_L;
    } else if (methodStr.compare("LEMP_I") == 0) {
        method = LEMP_I;
    } else if (methodStr.compare("LEMP_C") == 0) {
        method = LEMP_C;
    } else if (methodStr.compare("LEMP_TA") == 0) {
        method = LEMP_TA;
    } else if (methodStr.compare("LEMP_TREE") == 0) {
        method = LEMP_TREE;
    } else if (methodStr.compare("LEMP_AP") == 0) {
        method = LEMP_AP;
    } else if (methodStr.compare("LEMP_LSH") == 0) {
        method = LEMP_LSH;
    } else if (methodStr.compare("LEMP_BLSH") == 0) {
        method = LEMP_BLSH;
    } 
    else {
        cout << "[ERROR] This method is not possible. Please try {LEMP_L, LEMP_LI, LEMP_LC, LEMP_I, LEMP_C, LEMP_TA, LEMP_TREE, LEMP_AP, LEMP_LSH, LEMP_BLSH}" << endl << endl;
        cout << desc << endl;
        return 1;
    }

    VectorMatrix leftMatrix, rightMatrix;

    if (querySideLeft) {
        leftMatrix.readFromFile(usersFile, r, m, true);
        rightMatrix.readFromFile(itemsFile, r, n, false);
    } else {
        leftMatrix.readFromFile(itemsFile, r, n, false);
        rightMatrix.readFromFile(usersFile, r, m, true);
    }

    mips::Lemp algo(args, cacheSizeinKB, method, isTARR, R, epsilon);
    
    algo.initialize(rightMatrix);

    Results results;
    if (args.k > 0) {
        algo.runTopK(leftMatrix, results);
    } else {
        algo.runAboveTheta(leftMatrix, results);
    }
    
    if (resultsFile != "") {
        results.writeToFile(resultsFile);
    }

    return 0;
}


