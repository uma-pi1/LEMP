/* 
 * File:   runLSH_all.cc
 * Author: chteflio
 *
 * Created on March 30, 2015, 6:12 PM
 */

#include <boost/program_options.hpp>

#include <iostream>
#include <mips/mips.h>

using namespace std;
using namespace mips;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    string usersFile;
    string itemsFile;
    string logFile, resultsFile;

    bool querySideLeft = true;
    bool isTransformed = false;
    int k, r, m, n, threads;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&usersFile), "file containing the query matrix (left side)")
            ("P", value<string>(&itemsFile), "file containing the probe matrix (right side)")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("isTransformed", value<bool>(&isTransformed)->default_value(false), "1 if the input files contain the tranformed vectors. If 0 the program will transform them")
            ("k", value<int>(&k), "top k ")
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
        cout << "runSimpleLsh [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }

    InputArguments args;
    args.logFile = logFile;
    args.k = k;
    args.threads = threads;

    VectorMatrix leftMatrix, rightMatrix;

    if (querySideLeft) {
        leftMatrix.readFromFile(usersFile, r, m, true);
        rightMatrix.readFromFile(itemsFile, r, n, false);
    } else {
        leftMatrix.readFromFile(itemsFile, r, n, false);
        rightMatrix.readFromFile(usersFile, r, m, true);
    }
    
    mips::SimpleLsh algo(args, isTransformed);
    
    algo.initialize(rightMatrix);

    Results results;
    algo.runTopK(leftMatrix, results);
    
    if (resultsFile != "") {
        results.writeToFile(resultsFile);
    }

    return 0;
}

