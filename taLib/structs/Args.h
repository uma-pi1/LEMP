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
 * Args.h
 *
 *  Created on: Mar 13, 2014
 *      Author: chteflio
 */

#ifndef ARGS_H_
#define ARGS_H_

namespace ta{

enum LEMP_Method{
        LEMP_L = 0,
	LEMP_I = 1,
        LEMP_TREE = 2,
        LEMP_LSH = 3,
        LEMP_TA = 4,        
	LEMP_LI = 5,
	LEMP_LC = 6,	
	LEMP_C = 7,	
        LEMP_AP = 8,
        LEMP_TANRA = 9
    
};

struct LEMPArg{
        double theta;
	int k;
	std::string usersFile;
	std::string itemsFile;
	bool querySideLeft;
        bool isTARR;
        std::string logFile, resultsFile;
        comp_type comparisons;
    
	LEMP_Method method;
	int cacheSizeinKB;
	int threads;
        int depth;// for PCA Trees
        double R, epsilon; // for LSH R:recall
	col_type listsForExplorationMode;
	double dataManipulationTime, tuningTime;
        
        int r, m, n;// rank, numUsers, numItems: necessary only when reading a csv file

	double boundsTime, ipTime, scanTime, preprocessTime, filterTime, initializeListsTime, speedyTime;


	LEMPArg():k(0), querySideLeft(true), cacheSizeinKB (8192), method(LEMP_LI), boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0), speedyTime(0), initializeListsTime(0),
			dataManipulationTime(0), tuningTime(0), threads(1), listsForExplorationMode(1), depth(5), R(1.0), epsilon(0), isTARR (false), r(0), m(0), n(0){}

	void printTimes(){
		std::cout<< "Preprocessing Time: "<<dataManipulationTime/1E9 <<std::endl;
		std::cout<< "Tuning Time: "<<tuningTime/1E9 <<std::endl;
	}


};






}
#endif /* ARGS_H_ */
