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
        LEMP_BLSH = 9,        
        LEMP_NB = 10 // uses many and does naive bayes with costs
};



struct Arg{
	double theta;
	int k;
	std::string usersFile;
	std::string itemsFile;
	bool querySideLeft;
	
        std::string clusterFile;

	// logging
	comp_type comparisons;
	std::string logFile, resultsFile;

	Arg(): k(0), querySideLeft(true){}

};


struct LEMPArg: public Arg{
	LEMP_Method method;
	int cacheSizeinKB;
	int threads;
        int depth;// for PCA Trees
        int numClusters;
        double epsilon; // for LSH
	col_type listsForExplorationMode;
	double dataManipulationTime, tuningTime;

	double boundsTime, ipTime, scanTime, preprocessTime, filterTime, initializeListsTime, speedyTime;
	double stepTime, queueTime, thresTime, updateStateTime, preTime; // todo change names and add to constructor

	LEMPArg():Arg(), cacheSizeinKB (8192), method(LEMP_LI), boundsTime(0), ipTime(0), scanTime(0), preprocessTime(0), filterTime(0), speedyTime(0), initializeListsTime(0),
			dataManipulationTime(0), tuningTime(0), stepTime(0), queueTime(0), thresTime(0), updateStateTime(0), preTime(0), threads(1),
			listsForExplorationMode(1), depth(5), epsilon(0.03), numClusters(0){}

	void printTimes(){
		

		std::cout<< "Preprocessing Time: "<<dataManipulationTime/1E9 <<std::endl;
		std::cout<< "Tuning Time: "<<tuningTime/1E9 <<std::endl;
		
#ifdef TIME_IT
		std::cout<< "-------------"<<std::endl;
		if (method != LEMP_TA){

			//LOG4CXX_INFO(logger, "initializeListsTime: "<<initializeListsTime/1E9 );
		}else{
			std::cout<<"preTime: "<<preTime/1E9 <<std::endl;
			std::cout<< "stepTime: "<<stepTime/1E9<<std::endl;
			std::cout<< "queueTime: "<<queueTime/1E9<<std::endl;
			std::cout<< "thresTime: "<<thresTime/1E9 <<std::endl;
			std::cout<< "updateStateTime: "<<updateStateTime/1E9 <<std::endl;
			std::cout<< "ipTime: "<<ipTime/1E9 <<std::endl;
		}
		std::cout<< "-------------" <<std::endl;
#endif
	}


};

struct TreeArg: public Arg{
	int N0;
	TreeArg(): Arg(), N0(2){}

};





}
#endif /* ARGS_H_ */
