/* 
 * File:   runLSH_all.cc
 * Author: chteflio
 *
 * Created on March 30, 2015, 6:12 PM
 */

#include <boost/program_options.hpp>

#include <iostream>
#include <taLib/taLib.h>


using namespace std;
using namespace ta;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    double epsilon;
    string leftMatrix;
    string rightMatrix;
    string logFile, resultsFile;

    bool querySideLeft = true;
    int k, r, m, n;


    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&leftMatrix), "file containing the query matrix (left side)")
            ("P", value<string>(&rightMatrix), "file containing the probe matrix (right side)")
//            ("eps", value<double>(&epsilon)->default_value(0.03), "false negative rate (unused)")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("k", value<int>(&k), "top k ")
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
        cout << "run SimpleLSH [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }


    LEMPArg args;
    args.usersFile = leftMatrix;
    args.itemsFile = rightMatrix;
    args.logFile = logFile;
    args.k = k;
    args.querySideLeft = querySideLeft;
//    args.epsilon = epsilon;
    args.resultsFile = resultsFile;
    args.r = r;
    args.m = m;
    args.n = n;

    ta::SimpleLSH algo(args);
    algo.runTopkPerUser();


    return 0;
}

