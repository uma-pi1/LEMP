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
 * AuxStructs.h
 *
 *  Created on: Oct 17, 2013
 *      Author: chteflio
 */

#ifndef AUXSTRUCTS_H_
#define AUXSTRUCTS_H_

#include<taLib/structs/BasicStructs.h>
#include<utility>

namespace ta{

struct NormVector{
	std::vector<double> coordinates; // the coordinates of the vector
	double length;
	row_type id;
	col_type nonZero;




	NormVector(){};
	NormVector(std::vector<double>& vec, row_type vectorId=0){
		id = vectorId;
		col_type len = vec.size();
		coordinates.reserve(len);

		length=0;
		for (col_type i = 0; i < len; i++)
		{
			length += vec[i] * vec[i];
		}
		length = sqrt(length);

		nonZero = 0;

		for (col_type i = 0; i < len; i++) {
			coordinates.push_back(vec[i]/length);
			if (coordinates[i]>0.001)
				nonZero++;
		}
	}

	// GETTERS
	inline double getLengthsProduct(NormVector& q){return length * q.getLength();}
	inline double getLength()const{return length;}
	inline std::vector<double>& getCoordinates(){return coordinates;}
	inline double getCoordInPos(col_type i){return coordinates[i];}
	inline row_type getId(){return id;}
	inline col_type getNonZeros(){return nonZero;}



	// SETTERS
	inline void setCoordInPos(col_type i, double value){coordinates[i] = value;}
	inline void setLength(double value){length = value;}

	inline void clear(col_type size=0){
		coordinates.clear();
		length = 0;
		if (size != 0){
			coordinates.resize(size);
		}
	};

	inline std::string print(){
		std::ostringstream  s1;
		s1<<length;
		s1<<coordinates;
		return s1.str();
	}

	/*
	 * The following functions order TAVectors according to their lengths
	 */
	inline bool operator<(const NormVector & other) const {return length < other.getLength();}
	inline bool operator<=(const NormVector & other) const {return length <= other.getLength();}
	inline bool operator>(const NormVector & other) const {return length > other.getLength();}
	inline bool operator>=(const NormVector & other) const {return length >= other.getLength();}

	inline double cosine(NormVector& q){
		double sum = 0;

		for (col_type i=0, len=coordinates.size(); i<len; i++){
			sum += coordinates[i] * q.coordinates[i];
		}
		return sum;
	}

	inline double innerProduct(NormVector& q){return cosine(q) * getLengthsProduct(q);}

	// bool: passes the theta or not
	// double the actual inner product
	inline std::pair<bool, double> passesThreshold(NormVector& q, double theta){
		std::pair<bool, double> p;
		p.second = 0;
		double sum = 0;

		double ip = getLengthsProduct(q);
		if (ip < theta){ //prune quickly
			p.first = false;
			return p;
		}

		p.second = ip * cosine(q);
		if (p.second >= theta)
			p.first = true;
		else
			p.first = false;

		return p;
	}

	inline double L1norm(){
		double sum = 0;
		for (col_type i=0, len=coordinates.size(); i<len; i++){
			sum += fabs(coordinates[i]);
		}
		sum *= length;
		return sum;
	}
};

struct NormMatrix{
	std::vector<NormVector> items;
public:

	inline NormMatrix(){};
	row_type rowNum;
	col_type colNum;
	col_type maxNonZeros;
	double avgNonZeros;

	inline void init(std::vector<std::vector<double> >& matrix){

		rowNum = matrix.size();
		colNum = matrix[0].size();
		items.clear();
		items.reserve(rowNum);

		double avgSimNonZeros = 0;

		maxNonZeros = 0;
		avgNonZeros = 0;
		//		rg::Timer t;
		//		t.start();
		for (row_type i=0; i<rowNum; i++){
			items.push_back(NormVector(matrix[i], i));

			avgNonZeros += items[i].getNonZeros();

			double x =  items[i].L1norm() / items[i].getLength();
			x *= x;
			avgSimNonZeros += x;

			if (items[i].getNonZeros()>maxNonZeros || maxNonZeros == 0)
				maxNonZeros = items[i].getNonZeros();
		}

		avgNonZeros /= rowNum;
		avgSimNonZeros /= rowNum;

		std::cout<<"Avg NonZeros: "<<avgNonZeros<<" avgSimNonZeros: "<<avgSimNonZeros<<std::endl;
	}

	inline NormVector& getVecInPos(row_type row){return items[row];}

	inline std::pair<bool, double> passesThreshold(NormVector& q, row_type rowId, double theta){return items[rowId].passesThreshold(q, theta);}

	inline double innerProduct(NormVector& q, row_type rowId){return items[rowId].innerProduct(q);}

	inline void print(){
		std::cout<< "length |"<<" coordinates"<<std::endl;
		for (row_type i=0; i<rowNum; i++){
			items[i].print();
		}
	}
	inline double getValue(row_type rowId, col_type coord){return items[rowId].coordinates[coord];}
	inline double getLength(row_type rowId){return items[rowId].length;}
	inline void sort(){
		std::sort (items.begin(), items.end(), std::greater<NormVector>());
	}
};


}


#endif /* AUXSTRUCTS_H_ */
