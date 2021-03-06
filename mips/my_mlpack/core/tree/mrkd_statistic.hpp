/**
 * @file mrkd_statistic.hpp
 * @author James Cline
 *
 * Definition of the statistic for multi-resolution kd-trees.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_CORE_TREE_MRKD_STATISTIC_HPP
#define __MY_MLPACK_CORE_TREE_MRKD_STATISTIC_HPP

#include <mips/my_mlpack/core.hpp>

namespace mips {
namespace tree {

/**
 * Statistic for multi-resolution kd-trees.
 */
class MRKDStatistic
{
 public:
  //! Initialize an empty statistic.
  MRKDStatistic();

  /**
   * This constructor is called when a node is finished initializing.
   *
   * @param node The node that has been finished.
   */
  template<typename TreeType>
  MRKDStatistic(const TreeType& /* node */);

  /**
   * Returns a string representation of this object.
   */
  std::string ToString() const;

  //! Get the index of the initial item in the dataset.
  size_t Begin() const { return begin; }
  //! Modify the index of the initial item in the dataset.
  size_t& Begin() { return begin; }

  //! Get the number of items in the dataset.
  size_t Count() const { return count; }
  //! Modify the number of items in the dataset.
  size_t& Count() { return count; }

  //! Get the center of mass.
  const arma::colvec& CenterOfMass() const { return centerOfMass; }
  //! Modify the center of mass.
  arma::colvec& CenterOfMass() { return centerOfMass; }

  //! Get the index of the dominating centroid.
  size_t DominatingCentroid() const { return dominatingCentroid; }
  //! Modify the index of the dominating centroid.
  size_t& DominatingCentroid() { return dominatingCentroid; }

  //! Access the whitelist.
  const std::vector<size_t>& Whitelist() const { return whitelist; }
  //! Modify the whitelist.
  std::vector<size_t>& Whitelist() { return whitelist; }

 private:
  //! The data points this object contains.
  const arma::mat* dataset;
  //! The initial item in the dataset, so we don't have to make a copy.
  size_t begin;
  //! The number of items in the dataset.
  size_t count;
  //! The left child.
  const MRKDStatistic* leftStat;
  //! The right child.
  const MRKDStatistic* rightStat;
  //! A link to the parent node; NULL if this is the root.
  const MRKDStatistic* parentStat;

  // Computed statistics.
  //! The center of mass for this dataset.
  arma::colvec centerOfMass;
  //! The sum of the squared Euclidean norms for this dataset.
  double sumOfSquaredNorms;

  // There may be a better place to store this -- HRectBound?
  //! The index of the dominating centroid of the associated hyperrectangle.
  size_t dominatingCentroid;

  //! The list of centroids that cannot own this hyperrectangle.
  std::vector<size_t> whitelist;
  //! Whether or not the whitelist is valid.
  bool isWhitelistValid;
};

}; // namespace tree
}; // namespace mlpack

// Include implementation.
#include "mrkd_statistic_impl.hpp"

#endif // __MLPACK_CORE_TREE_MRKD_STATISTIC_HPP
