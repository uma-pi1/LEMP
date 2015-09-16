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
