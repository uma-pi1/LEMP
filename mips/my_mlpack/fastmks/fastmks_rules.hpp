/**
 * @file fastmks_rules.hpp
 * @author Ryan Curtin
 *
 * Rules for the single or dual tree traversal for fast max-kernel search.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_METHODS_FASTMKS_FASTMKS_RULES_HPP
#define __MY_MLPACK_METHODS_FASTMKS_FASTMKS_RULES_HPP

#include <mips/my_mlpack/core.hpp>
#include <mips/my_mlpack/core/tree/cover_tree/cover_tree.hpp>
#include <mips/structs/BasicStructs.h>


namespace mips {
namespace fastmks {

/**
 * The base case and pruning rules for FastMKS (fast max-kernel search).
 */
template<typename TreeType>
class FastMKSRules
{
 public:
  FastMKSRules(KernelType& kernel,
               VectorMatrix* probeMatrix=NULL, VectorMatrix* queryMatrix=NULL);



  //! Compute the base case (kernel value) between two points.
  double BaseCase(const size_t queryIndex, const size_t referenceIndex, std::vector<QueueElement>& results);

  double BaseCaseForTheta(const size_t queryIndex, const size_t referenceIndex, std::vector<MatItem>& results, 
          const double theta);

  /**
   * Get the score for recursion order.  A low score indicates priority for
   * recursion, while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).
   *
   * @param queryIndex Index of query point.
   * @param referenceNode Candidate to be recursed into.
   */
  double Score(const size_t queryIndex, TreeType& referenceNode, std::vector<QueueElement>& results);
  double ScoreForTheta(const size_t queryIndex, TreeType& referenceNode, std::vector<MatItem>& results, 
          const double theta);

  /**
   * Get the score for recursion order.  A low score indicates priority for
   * recursion, while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).
   *
   * @param queryNode Candidate query node to be recursed into.
   * @param referenceNode Candidate reference node to be recursed into.
   */
  double Score(TreeType& queryNode, TreeType& referenceNode, std::vector<QueueElement>& results);


  double ScoreForTheta(TreeType& queryNode, TreeType& referenceNode, std::vector<MatItem>& results, double theta);

  /**
   * Re-evaluate the score for recursion order.  A low score indicates priority
   * for recursion, while DBL_MAX indicates that a node should not be recursed
   * into at all (it should be pruned).  This is used when the score has already
   * been calculated, but another recursion may have modified the bounds for
   * pruning.  So the old score is checked against the new pruning bound.
   *
   * @param queryIndex Index of query point.
   * @param referenceNode Candidate node to be recursed into.
   * @param oldScore Old score produced by Score() (or Rescore()).
   */
  double Rescore(const size_t queryIndex,
                 TreeType& referenceNode,
                 const double oldScore, std::vector<QueueElement>& results) const;
  double RescoreForTheta(TreeType& /*referenceNode*/, const double oldScore, const double theta) const;


  //! Get the number of times BaseCase() was called.
  size_t BaseCases() const { return baseCases; }
  //! Modify the number of times BaseCase() was called.
  size_t& BaseCases() { return baseCases; }

  //! Get the number of times Score() was called.
  size_t Scores() const { return scores; }
  //! Modify the number of times Score() was called.
  size_t& Scores() { return scores; }

  VectorMatrix* probeMatrix;
  VectorMatrix* queryMatrix;
  //std::vector<std::vector<QueueElement> >* finalResults; // my code

 private:


  //! The instantiated kernel.
  KernelType& kernel;

  //! The last query index BaseCase() was called with.
  size_t lastQueryIndex;
  //! The last reference index BaseCase() was called with.
  size_t lastReferenceIndex;
  //! The last kernel evaluation resulting from BaseCase().
  double lastKernel;

  //! Calculate the bound for a given query node.
  double CalculateBound(TreeType& queryNode) const;

  //! For benchmarking.
  size_t baseCases;
  //! For benchmarking.
  size_t scores;
};

}; // namespace fastmks
}; // namespace mlpack

// Include implementation.
#include "fastmks_rules_impl.hpp"

#endif
