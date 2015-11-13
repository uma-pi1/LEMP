/* 
 * File:   testLSH_all.cc
 * Author: chteflio
 *
 * Created on March 26, 2015, 12:31 PM
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
    args.k = 10;
    args.R = 0.9;
    args.isTransformed = false;

//    args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w50-net-srebro-noav.mma";
//    args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h50-net-srebro-noav.mma";

    args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w50-net-noav.mma";
    args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h50-net-noav.mma";


    args.logFile = "/home/chteflio/workspace/Log.txt";
//    args.resultsFile="/home/chteflio/data/test1.data";


    ta::SimpleLSH algo(args);


    algo.runTopkPerUser();

    return 0;
}








