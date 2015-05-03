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
 * read.h
 *
 *  Created on: Sep 14, 2012
 *      Author: chteflio
 */

#ifndef TALIB_STRUCTS_READ_H_
#define TALIB_STRUCTS_READ_H_

#include <string>
#include <vector>
#include <ostream>
#include <iomanip>
#include <taLib/structs/BasicStructs.h>
#include <boost/numeric/interval.hpp>


#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <boost/tokenizer.hpp>

#include <boost/algorithm/string.hpp>



namespace ta{

std::vector< std::vector<double> > readMatrix(std::string& fileName, bool left=true);

//std::vector< std::vector<row_type> > readClusters(std::string& fileName, int numOfClusters, row_type numOfItems);
void readClusters(std::string& fileName, std::vector<int> & clusterIds, row_type numOfItems) ;

//void writeResults(std::string& fileName, std::vector<MatItem>& results);



bool compareResults(std::string& fileName1, std::string& fileName2);

template <typename T>
std::ostream & operator<<(std::ostream & os, const std::vector<std::vector<T> > & m)
{
	os << "matrix: " << m.size() << " x " << (m.size() == 0 ? 0 : m[0].size()) << std::endl;
	os << std::fixed;
	for (ta_size_type i = 0; i < m.size(); i++) {
		for (ta_size_type j = 0; j < m[0].size(); j++) {
			os << std::setprecision(5) << m[i][j] << "  ";
		}
		os << std::endl;
	}
	if (m.size() > 0) {
		for (ta_size_type j = 0; j < m[0].size(); j++) {
			os << "-------- ";
		}
		os << std::endl;
	}
	return os;
}

template <typename T>
std::ostream & operator<<(std::ostream & os, const std::vector<T> & m)
{
	os << "Vector: " << m.size()  << std::endl;
	os << std::fixed;
	for (ta_size_type i = 0; i < m.size(); i++) {
		os << std::setprecision(5) << m[i] << "  ";
	}
	os << std::endl;
	return os;
}


//int readBucketBestResults(std::string& fileName, std::vector<row_type>& ids, std::vector<double>& localThetas, std::vector<row_type>& labels){
//
//	std::ifstream file(fileName.c_str(), std::ios_base::in);
//
//	int queryId, bucketId, method, numLists;
//	double localTheta, time;
//
//	while (file) {
//		file >> queryId >> bucketId >> localTheta >> method >> numLists >> time;
//		ids.push_back(queryId);
//		localThetas.push_back(localTheta);
//
//		if(method == 5){
//			labels.push_back(method*10+numLists);
//		}else if(method == 4){
//			labels.push_back(method*100+numLists);
//		}else{
//			labels.push_back(method);
//		}
//
//		//labels.push_back(method);
//
//
//	}
//	file.close();
//
//	// for some reason the last line is read twice
//	ids.pop_back();
//	localThetas.pop_back();
//	labels.pop_back();
//
//	std::cout<<"Observations: "<<ids.size()<<" "<<localThetas.size()<<" "<<labels.size()<<std::endl;
//
//	return bucketId;
//
//}
int readBucketBestResults(std::string& fileName, std::vector<row_type>& ids, std::vector<double>& localThetas, std::vector<row_type>& labels,
        std::vector<double>& timings){

	std::ifstream file(fileName.c_str(), std::ios_base::in);

	int queryId, bucketId, method, numLists;
	double localTheta, time;

	while (file) {
		file >> queryId >> bucketId >> localTheta >> method >> numLists >> time;
		ids.push_back(queryId);
		localThetas.push_back(localTheta);

		if(method == 5){
			labels.push_back(method*10+numLists);
		}else if(method == 4){
			labels.push_back(method*1000+numLists);
		}else{
			labels.push_back(method);
		}

		//labels.push_back(method);
                timings.push_back(time);

	}
	file.close();

	// for some reason the last line is read twice
	ids.pop_back();
	localThetas.pop_back();
	labels.pop_back();

	std::cout<<"Observations: "<<ids.size()<<" "<<localThetas.size()<<" "<<labels.size()<<std::endl;

	return bucketId;

}




}


#endif /* READ_H_ */
