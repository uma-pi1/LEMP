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
    args.logFile = "/home/chteflio/workspace/Log.txt";
    args.k = 10;   
    args.depth = 6;

    args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w50-net-reduced-noav.mma";
    args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h50-net-reduced-noav.mma";
    args.resultsFile= "/home/chteflio/newpca.data";
    
    
    ta::SearchWithPCATree algo(args);
    
    algo.topKperUser();
    



    return 0;
}
