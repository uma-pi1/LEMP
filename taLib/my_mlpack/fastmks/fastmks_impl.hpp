/**
 * @file fastmks_impl.hpp
 * @author Ryan Curtin
 *
 * Implementation of the FastMKS class (fast max-kernel search).
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_METHODS_FASTMKS_FASTMKS_IMPL_HPP
#define __MY_MLPACK_METHODS_FASTMKS_FASTMKS_IMPL_HPP

// In case it hasn't yet been included.
#include "fastmks.hpp"

#include "fastmks_rules.hpp"
#include <queue>
#include <taLib/structs/BasicStructs.h>


namespace ta {
    namespace fastmks {


        // my stuff

        template<typename TreeType>
        FastMKS<TreeType>::FastMKS(
                const bool single,
                const bool naive) :
        treeOwner(false),
        single(single),
        naive(naive),
        metric(IPMetric<LinearKernel>()) //direct hack
        {
        }

        template<typename TreeType>
        FastMKS<TreeType>::~FastMKS() {
        }

        template<typename TreeType>
        void FastMKS<TreeType>::Search(const size_t k,
                TreeType* referenceTree, VectorMatrix* probeMatrix, VectorMatrix* queryMatrix,
                std::vector<QueueElement> & finalResults,
                size_t queryInd, comp_type& comparisons, int threads) {

            typedef FastMKSRules<TreeType> RuleType;
            std::vector<RuleType> rules(threads, RuleType(metric.Kernel(), probeMatrix, queryMatrix));

            typename TreeType::template SingleTreeTraverser<RuleType> traverser(rules);
            traverser.Traverse(queryInd, *referenceTree, finalResults);

            for (int i=0; i<threads; i++)
                comparisons += rules[i].BaseCases();

            return;
        }




        // my code just for single

        template<typename TreeType>
        void FastMKS<TreeType>::SearchForTheta(const double theta,
                TreeType* referenceTree, VectorMatrix* probeMatrix, VectorMatrix* queryMatrix,
                std::vector<MatItem> & finalResults,
                size_t queryInd, comp_type& comparisons, int threads) {

            typedef FastMKSRules<TreeType> RuleType;
            //RuleType rules(metric.Kernel(), probeMatrix, queryMatrix);
            std::vector<RuleType> rules(threads, RuleType(metric.Kernel(), probeMatrix, queryMatrix));


            typename TreeType::template SingleTreeTraverser<RuleType> traverser(rules);
            traverser.TraverseForTheta(queryInd, *referenceTree, finalResults, theta);

             for (int i=0; i<threads; i++)
                comparisons += rules[i].BaseCases();

            return;

        }




    }; // namespace fastmks
}; // namespace mlpack

#endif
