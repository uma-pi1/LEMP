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

#ifndef TA_TALIB_H
#define TA_TALIB_H


#include <omp.h>

#include <util/random.h>

#include <mips/structs/Definitions.h>
#include <mips/structs/BasicStructs.h>
#include <mips/structs/Args.h>


#include <mips/structs/VectorMatrix.h>
#include <mips/structs/VectorMatrix2.h>
#include <mips/structs/Results.h>
#include <mips/structs/Output.h>
#include <mips/structs/Lists.h>

////////////////////////////// mlpack stuff//////////////////////
#include <mips/my_mlpack/core.hpp>
#include <mips/my_mlpack/fastmks/fastmks_stat.hpp>
#include <mips/my_mlpack/fastmks/fastmks.hpp>
#include <mips/my_mlpack/fastmks/fastmks_impl.hpp>
#include <mips/my_mlpack/fastmks/fastmks_rules.hpp>
#include <mips/my_mlpack/fastmks/fastmks_rules_impl.hpp>
///////////////////////////////////////////////////////////////////
/////////// l2ap stuff ///////////////
#include <mips/ap/includes.h>
#include <mips/ap/connect.h>
/////////////////////////////////////

#include <mips/lsh/RandomIntGaussians.h>
#include <mips/lsh/CosineSketches.h>
#include <mips/lsh/LshBins.h>

#include <mips/structs/l2apIndex.h>
#include <mips/structs/LshIndex.h>
#include <mips/structs/BlshIndex.h>


#include <mips/structs/TreeIndex.h>
#include <mips/structs/Candidates.h>
#include <mips/structs/TAState.h>
#include <mips/structs/TANRAState.h>
#include <mips/structs/RetrievalArguments.h>//////////////////////
#include <mips/structs/QueryBatch.h>
#include <mips/structs/ProbeBucket.h>
#include <mips/structs/Bucketize.h>
#include <mips/structs/CandidateVerification.h>

#include <mips/retrieval/Retriever.h>
#include <mips/retrieval/ListsTuneData.h>
#include <mips/retrieval/TuneTopk.h>
#include <mips/retrieval/coord.h>
#include <mips/retrieval/icoord.h>
#include <mips/retrieval/ta.h>
#include <mips/retrieval/tanra.h>
#include <mips/retrieval/mixed.h>
#include <mips/retrieval/singleTree.h>
#include <mips/retrieval/apRetriever.h>
#include <mips/retrieval/lshRetriever.h> 
#include <mips/retrieval/blshRetriever.h> 

#include <mips/algos/Mip.h>
#include <mips/algos/Naive.h>
#include <mips/algos/PcaTree.h>
#include <mips/algos/Ta.h>
#include <mips/algos/TaNra.h>
#include <mips/algos/SimpleLsh.h>
#include <mips/algos/Lemp.h>


#endif
