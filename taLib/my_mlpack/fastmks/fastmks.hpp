/**
 * @file fastmks.hpp
 * @author Ryan Curtin
 *
 * Definition of the FastMKS class, which implements fast exact max-kernel
 * search.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_METHODS_FASTMKS_FASTMKS_HPP
#define __MY_MLPACK_METHODS_FASTMKS_FASTMKS_HPP

#include <taLib/my_mlpack/core.hpp>
#include <taLib/my_mlpack/core/metrics/ip_metric.hpp>
#include <taLib/my_mlpack/core/kernels/linear_kernel.hpp>
#include "fastmks_stat.hpp"
#include <taLib/my_mlpack/core/tree/cover_tree.hpp>

#include <taLib/structs/BasicStructs.h>

using namespace ta::kernel;
using namespace ta::tree;
using namespace ta::metric;

namespace ta {
namespace fastmks /** Fast max-kernel search. */ {

typedef LinearKernel KernelType;

/**
 * An implementation of fast exact max-kernel search.  Given a query dataset and
 * a reference dataset (or optionally just a reference dataset which is also
 * used as the query dataset), fast exact max-kernel search finds, for each
 * point in the query dataset, the k points in the reference set with maximum
 * kernel value K(p_q, p_r), where k is a specified parameter and K() is a
 * Mercer kernel.
 *
 * For more information, see the following paper.
 *
 * @code
 * @inproceedings{curtin2013fast,
 *   title={Fast Exact Max-Kernel Search},
 *   author={Curtin, Ryan R. and Ram, Parikshit and Gray, Alexander G.},
 *   booktitle={Proceedings of the 2013 SIAM International Conference on Data
 *       Mining (SDM 13)},
 *   year={2013}
 * }
 * @endcode
 *
 * This class allows specification of the type of kernel and also of the type of
 * tree.  FastMKS can be run on kernels that work on arbitrary objects --
 * however, this only works with cover trees and other trees that are built only
 * on points in the dataset (and not centroids of regions or anything like
 * that).
 *
 * @tparam KernelType Type of kernel to run FastMKS with.
 * @tparam TreeType Type of tree to run FastMKS with; it must have metric
 *     IPMetric<KernelType>.
 */




template<
    typename TreeType = tree::MyCoverTree<metric::IPMetric<KernelType>,
        tree::FirstPointIsRoot, FastMKSStat>
>
class FastMKS
{
 public:



  FastMKS(const bool single = false,
          const bool naive = false);



  //! Destructor for the FastMKS object.
  ~FastMKS();



  void Search(const size_t k,
  		TreeType* referenceTree, VectorMatrix* probeMatrix, VectorMatrix* queryMatrix,
  		std::vector<QueueElement> & finalResults,
  		size_t queryInd, comp_type& comparisons, int threads);
  
  void SearchForTheta(const double theta, 
  		TreeType* referenceTree, VectorMatrix* probeMatrix, VectorMatrix* queryMatrix,
  		std::vector<MatItem> & finalResults,
  		size_t queryInd, comp_type& comparisons, int threads);


  //! Get the inner-product metric induced by the given kernel.
  const metric::IPMetric<KernelType>& Metric() const { return metric; }
  //! Modify the inner-product metric induced by the given kernel.
  metric::IPMetric<KernelType>& Metric() { return metric; }

  //TreeType* referenceTree;
  //! The instantiated inner-product metric induced by the given kernel.
  metric::IPMetric<KernelType> metric;

 private:

  //! If true, this object created the trees and is responsible for them.
  bool treeOwner;

  //! If true, single-tree search is used.
  bool single;
  //! If true, naive (brute-force) search is used.
  bool naive;

};

}; // namespace fastmks
}; // namespace mlpack

// Include implementation.
#include "fastmks_impl.hpp"

#endif
