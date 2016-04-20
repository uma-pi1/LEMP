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
 * TreeIndex.h
 *
 *  Created on: Jul 21, 2014
 *      Author: chteflio
 */

#ifndef TREEINDEX_H_
#define TREEINDEX_H_


#include <mips/my_mlpack/core.hpp>
#include <mips/my_mlpack/core/metrics/ip_metric.hpp>
#include <mips/my_mlpack/fastmks/fastmks.hpp>

using namespace std;
using namespace mips::fastmks;
using namespace mips::kernel;
using namespace mips::metric;

namespace mips {

    typedef tree::MyCoverTree<IPMetric<LinearKernel>, tree::FirstPointIsRoot, FastMKSStat> TreeType; // Convenience typedef.

    class TreeIndex : public Index {
        double base;
        const arma::mat dataset; //dummy

    public:

        TreeType* tree;

        inline TreeIndex() : base(1.3), tree(0) {
        }

        inline ~TreeIndex() {
            if (!tree) {
                delete tree;
            }
        }

        inline void initializeTree(VectorMatrix& matrix, int threads, ta_size_type start = 0, ta_size_type end = 0) {
            omp_set_lock(&writelock);

            if (!initialized) {
                if (start == end) {
                    start = 0;
                    end = matrix.rowNum;
                }

                // now build the tree
                tree = new TreeType(dataset, threads, base, &matrix, start, end);

                initialized = true;
            }
            omp_unset_lock(&writelock);
        }

    };





}


#endif /* TREEINDEX_H_ */
