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
 * runTA_all.cc
 *
 *  Created on: Mar 11, 2014
 *      Author: chteflio
 */


#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <mips/mips.h>


using namespace std;
using namespace mips;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    double theta;
    string usersFile;
    string itemsFile;
    string logFile, resultsFile;

    bool querySideLeft = true;
    bool isTARR = true;
    int k, r, m, n, threads;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&usersFile), "file containing Q^T")
            ("P", value<string>(&itemsFile), "file containing the P")
            ("theta", value<double>(&theta), "theta value")
            ("k", value<int>(&k)->default_value(0), "top k (default 0). If 0 Above-theta will run")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("isTARR", value<bool>(&isTARR)->default_value(true), "If 1 Round Robin schedule is used (default). Otherwise Max PiQi")
            ("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
            ("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
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
        cout << "runTa [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }
    
    InputArguments args;
    args.theta = theta;
    args.threads = threads;
//    args.usersFile = leftMatrix;
//    args.itemsFile = rightMatrix;
    args.logFile = logFile;
//    args.resultsFile = resultsFile;
    args.k = k;
//    args.querySideLeft = querySideLeft;
//    args.method = LEMP_TA;
//    args.isTARR = isTARR;
//    args.r = r;
//    args.m = m;
//    args.n = n;

    VectorMatrix leftMatrix, rightMatrix;

    if (querySideLeft) {
        leftMatrix.readFromFile(usersFile, r, m, true);
        rightMatrix.readFromFile(itemsFile, r, n, false);
    } else {
        leftMatrix.readFromFile(itemsFile, r, n, false);
        rightMatrix.readFromFile(usersFile, r, m, true);
    }
    
    mips::TaNra algo(args, isTARR);

    algo.initialize(rightMatrix);
    
    Results results;
    if (args.k > 0) {
//         algo.runTopK(leftMatrix, results);
	    std::cout<<"top-k for TA-NRA not yet implemented"<<std::endl;
    } else {
        algo.runAboveTheta(leftMatrix, results);
    }
    
    if (resultsFile != "") {
        results.writeToFile(resultsFile);
    }

    return 0;
}

