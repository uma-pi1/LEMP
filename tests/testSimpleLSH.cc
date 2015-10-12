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



    args.usersFile = "/home/chteflio/chteflioInMPI/data/ta-data/w-net-srebro.mma";
    args.itemsFile = "/home/chteflio/chteflioInMPI/data/ta-data/h-net-srebro.mma";


    args.logFile = "/home/chteflio/workspace/Log.txt";



    ta::SimpleLSH algo(args);


    algo.runTopkPerUser();

    return 0;
}








