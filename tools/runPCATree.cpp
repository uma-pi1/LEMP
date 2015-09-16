#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <taLib/taLib.h>


using namespace std;
using namespace ta;
using namespace boost::program_options;

int main(int argc, char *argv[]) {
    int k, depth;
    string leftMatrix;
    string rightMatrix;
    string logFile, resultsFile;
    bool querySideLeft = true;

    // read command line
    options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("Q^T", value<string>(&leftMatrix), "file containing the query matrix (left side)")
            ("P", value<string>(&rightMatrix), "file containing the probe matrix (right side)")
            ("k", value<int>(&k)->default_value(5), "top k (default=5)")
            ("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k") 
            ("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
	    ("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
            ("d", value<int>(&depth)->default_value(0), "depth of the tree (default=0)")
            ;

    positional_options_description pdesc;
    pdesc.add("Q^T", 1);
    pdesc.add("P", 2);

    variables_map vm;
    store(command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
    notify(vm);

    if (vm.count("help") || vm.count("Q^T") == 0 || vm.count("P") == 0) {
        cout << "run PCA Tree [options] <Q^T> <P>" << endl << endl;
        cout << desc << endl;
        return 1;
    }
    std::cout << "Running PCA tree " << std::endl;

    LEMPArg args;
    args.usersFile = leftMatrix;
    args.itemsFile = rightMatrix;
    args.resultsFile = resultsFile;
    args.logFile = logFile;
    args.k = k;
    args.depth = depth;
    args.querySideLeft = querySideLeft;


    ta::SearchWithPCATree algo(args);


    std::cout << "k: " << k << std::endl;
    algo.topKperUser();



    return 0;
}
