// Copyright 2015 Christina Teflioudi
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
/*
* runNaive.cc
*
* Created on: Sep 21, 2012
* Author: chteflio
*/
/* 
 * File:   testTANRA_all.cpp
 * Author: chteflio
 *
 * Created on August 24, 2015, 10:04 AM
 */

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
    args.theta = 1.74;

    args.usersFile = "/home/chteflio/data/ta-data/so_x_p_proj5_1-w-50-svd.mma";
    args.itemsFile = "/home/chteflio/data/ta-data/so_x_p_proj5_1-h-50-svd.mma";
    
    
    

    args.logFile = "/home/chteflio/workspace/Log.txt";

    ta::TANRA_all algo(args);
    algo.multiply();

    return 0;
}


