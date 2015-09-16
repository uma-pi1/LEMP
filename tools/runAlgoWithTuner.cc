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
#include <taLib/taLib.h>


using namespace std;
using namespace ta;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    double theta, epsilon, gamma;
    string leftMatrix;
    string rightMatrix;
    string logFile, resultsFile;


    bool querySideLeft = true;
    int k, cacheSizeinKB, threads;
    std::string methodStr;
    LEMP_Method method;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&leftMatrix), "file containing the query matrix (left side)")
            ("P", value<string>(&rightMatrix), "file containing the probe matrix (right side)")
            ("theta", value<double>(&theta), "theta value")
            ("eps", value<double>(&epsilon)->default_value(0.03), "epsilon value (experimental)")
	    ("gamma", value<double>(&gamma)->default_value(0.0), "gamma value (experimental)")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("method", value<string>(&methodStr), "LEMP_X where X: L, LI, LC, I, C, TA, TREE, AP, BLSH")
            ("k", value<int>(&k)->default_value(0), "top k (default 0). If 0 Above-theta will run")
            ("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
	    ("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
            ("cacheSizeinKB", value<int>(&cacheSizeinKB)->default_value(8192), "cache size in KB")
            ("t", value<int>(&threads)->default_value(1), "num of threads (default 1)")
            ;




    positional_options_description pdesc;
    pdesc.add("Q^T", 1);
    pdesc.add("P", 2);

    variables_map vm;
    store(command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
    notify(vm);

    if (vm.count("help") || vm.count("Q^T") == 0 || vm.count("P") == 0) {
        cout << "run Algo with Tuner [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }
    
    
    LEMPArg args;
    args.theta = theta;
    args.usersFile = leftMatrix;
    args.itemsFile = rightMatrix;
    args.logFile = logFile;
    args.resultsFile = resultsFile;
    args.k = k;
    args.querySideLeft = querySideLeft;
    args.cacheSizeinKB = cacheSizeinKB;
    args.threads = threads;
    args.epsilon = epsilon;
    args.gamma = gamma;
    
   

    if (methodStr.compare("LEMP_LI") == 0) {
        args.method = LEMP_LI;
    } else if (methodStr.compare("LEMP_LC") == 0) {
        args.method = LEMP_LC;
    } else if (methodStr.compare("LEMP_L") == 0) {
        args.method = LEMP_L;
    } else if (methodStr.compare("LEMP_I") == 0) {
        args.method = LEMP_I;
    } else if (methodStr.compare("LEMP_C") == 0) {
        args.method = LEMP_C;
    } else if (methodStr.compare("LEMP_TA") == 0) {
        args.method = LEMP_TA;
    } else if (methodStr.compare("LEMP_TREE") == 0) {
        args.method = LEMP_TREE;
    } else if (methodStr.compare("LEMP_AP") == 0) {
        args.method = LEMP_AP;
    } else if (methodStr.compare("LEMP_LSH") == 0) {
        args.method = LEMP_LSH;
    } 
    else {
        cout << "I do not know this method. Try {LEMP_LI, LEMP_LC, LEMP_L, LEMP_I, LEMP_C, LEMP_TA, LEMP_TREE, LEMP_AP, LEMP_BLSH, LEMP_LSH}" << endl << endl;
        cout << desc << endl;
        return 1;
    }



    ta::Algo_withTuning algo(args);

    if (args.k > 0) {
        algo.runTopkPerUser();
    } else {

        algo.multiply();
    }

    return 0;
}


