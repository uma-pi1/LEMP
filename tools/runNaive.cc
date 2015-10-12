
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
 * runNaive.cc
 *
 *  Created on: Sep 21, 2012
 *      Author: chteflio
 */

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <taLib/taLib.h>


using namespace std;
using namespace ta;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    double theta;
    int k, r, m, n;
    rg::Timer t;
    string leftMatrix;
    string rightMatrix;
    string logFile, resultsFile;
    string traceVar;
    bool querySideLeft = true;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&leftMatrix), "file containing Q^T")
            ("P", value<string>(&rightMatrix), "file containing the P")
            ("theta", value<double>(&theta), "theta")
            ("k", value<int>(&k)->default_value(0), "top k (default 0). If 0 Above-theta will run")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
            ("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
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
        cout << "runNaive [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }
    std::cout << "Running NAIVE " << std::endl;

    LEMPArg args;
    args.theta = theta;
    args.usersFile = leftMatrix;
    args.itemsFile = rightMatrix;
    args.logFile = logFile;
    args.resultsFile = resultsFile;
    args.k = k;
    args.querySideLeft = querySideLeft;
    args.r = r;
    args.m = m;
    args.n = n;

    Naive naive(args);

    if (k > 0) {
        std::cout << "k: " << k << std::endl;
        naive.topKperUser();
    } else {
        std::cout << "theta: " << theta << std::endl;
        naive.multiply();
    }


    return 0;
}


