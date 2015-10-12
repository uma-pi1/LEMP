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
 * Definitions.h
 *
 *  Created on: Mar 18, 2014
 *      Author: chteflio
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <util/evaluation.h>


// #define TUNE
// #define DEBUG
// #define TIME_IT
// #define WITH_SIMD
// #define WITH_SIMD_INCR


// for bucketizing
#define FACTOR 0.9
#define ITEMS_PER_BLOCK  30

// for tuning
#define UPPER_LIMIT_PER_BUCKET 250
#define LOWER_LIMIT_PER_BUCKET  20
#define NUM_LISTS    10

// #define RELATIVE_APPROX
// #define ABS_APPROX


// for probeBuckets
#define NUM_INDEXES 6 // 0: no index 1: sorted list 2: int sorted list 3: tree 4: AP 5:LSH 
#define LSH_SIGNATURES 200 //so that if I am about to get more than 80% of the probe vectors I will run length || only multiples of 4 please
#define LSH_CODE_LENGTH 8


#define INVPI  1 / PI

// for parallelizing the Trees
#define PADDING 8 // padding for doubles so that there will be no false sharing




double matrixToMatrixTime = 0;



#endif /* DEFINITIONS_H_ */
