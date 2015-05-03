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
#include <taLib/taLib.h>


using namespace std;
using namespace ta;
using namespace boost::program_options;





int main(int argc, char *argv[]) {
	double theta;
	string leftMatrix;
	string rightMatrix;
	string logFile, resultsFile;

	bool querySideLeft = true;
	int k;

	// read command line
	options_description desc("Options");
	desc.add_options()
				("help", "produce help message")
				("Q^T", value<string>(&leftMatrix), "file containing Q^T")
				("P", value<string>(&rightMatrix), "file containing the P")
				("theta", value<double>(&theta), "theta value")
				("k", value<int>(&k)->default_value(0), "top k (default 0). If 0 Above-theta will run")
				("querySideLeft", value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")	   
				("logFile", value<string>(&logFile)->default_value(""), "output File (contains runtime information)")
				("resultsFile", value<string>(&resultsFile)->default_value(""), "output File (contains the results)")
                               
				
	;

	positional_options_description pdesc;
	pdesc.add("Q^T", 1);
	pdesc.add("P", 2);

	variables_map vm;
	store(command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
	notify(vm);

	if (vm.count("help") || vm.count("Q^T")==0 || vm.count("P")==0) {
		cout << "run TA all [options] <Q^T> <P>" << endl << endl;
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
        args.method = LEMP_TA;



	ta::TA_all algo(args);

	if(args.k > 0){
		algo.runTopkPerUser();
	}else{
		algo.multiply();
	}

	return 0;
}

