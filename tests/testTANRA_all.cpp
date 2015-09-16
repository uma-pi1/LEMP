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


