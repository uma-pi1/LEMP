/*
 * testAlgoTuning.cc
 *
 *  Created on: Feb 5, 2014
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



#include <iostream>
#include <taLib/taLib.h>

#include <sys/time.h>
#include <fstream>
#include <climits>
#include <cmath>
#include <algorithm>


using namespace std;
using namespace ta;
using namespace rg;

int main() {

    LEMPArg args;
    args.theta = 0.035;
    //       args.k = 10;
    //        args.querySideLeft = false;
    args.method = LEMP_LI;
    args.cacheSizeinKB = 2048;
    args.threads = 1;
    args.R = 0.5;
    args.isTARR = true;
    args.r = 50;
    args.m = 771611;
    args.n = 132209;
           


//    args.usersFile = "/home/chteflio/data/ta-data/w-50-gnmf-zero.mma";
//    args.itemsFile = "/home/chteflio/data/ta-data/h-50-gnmf-zero.mma";

    args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w-50-gnmf-zero.csv";
    args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h-50-gnmf-zero.csv";

    //     args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w50-net-noav.mma";
    //     args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h50-net-noav.mma";


    args.logFile = "/home/chteflio/workspace/Log.txt";
    //    args.resultsFile = "/home/chteflio/data/test.data";

    ta::Algo_withTuning algo(args);


    if (args.k == 0) {
        algo.multiply();

    } else {

        algo.runTopkPerUser();

        std::vector<QueueElement> rs1 = algo.getResultsTopk();
    }


    return 0;
}
