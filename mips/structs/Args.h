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

namespace mips {

    enum LEMP_Method {
        LEMP_L = 0,
        LEMP_I = 1,
        LEMP_TREE = 2,
        LEMP_LSH = 3,
        LEMP_TA = 4,
        LEMP_LI = 5,
        LEMP_LC = 6,
        LEMP_C = 7,
        LEMP_AP = 8,
        LEMP_TANRA = 9,
        LEMP_BLSH = 10

    };

    struct InputArguments {
        double theta;
        int k;
        int threads;
        std::string logFile;
        
        InputArguments():k(0),threads(1){}

        void copyInputArguments(InputArguments& input) {
            theta = input.theta;
            k = input.k;
            threads = input.threads;
            logFile = input.logFile;
        }
    };

    struct LempArguments : public InputArguments {
        bool isTARR;
        comp_type comparisons;
        LEMP_Method method;
        int cacheSizeinKB;
        double R, epsilon; // for LSH R:recall
        int numTrees;
        int search_k;

        LempArguments() : cacheSizeinKB(sysconf(_SC_LEVEL2_CACHE_SIZE) / pow(2, 10)),
        method(LEMP_LI),  R(1.0), epsilon(0), isTARR(false), numTrees(1), search_k(1000) {
        }
    };

    struct PcaTreeArguments : public InputArguments {
        int depth = 5;
        PcaTreeArguments():InputArguments(){}
    };





}
#endif /* ARGS_H_ */
