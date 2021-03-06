/**
 * @file cosine_distance_impl.hpp
 * @author Ryan Curtin
 *
 * This implements the cosine distance.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_CORE_KERNELS_COSINE_DISTANCE_IMPL_HPP
#define __MY_MLPACK_CORE_KERNELS_COSINE_DISTANCE_IMPL_HPP

#include "cosine_distance.hpp"

namespace mips {
namespace kernel {

template<typename VecType>
double CosineDistance::Evaluate(const VecType& a, const VecType& b)
{
  // Since we are using the L2 inner product, this is easy.  But we have to make
  // sure we aren't dividing by zero (if we are, then the cosine similarity is
  // 0: we reason this value because the cosine distance is just a normalized
  // dot product; take away the normalization, and if ||a|| or ||b|| is equal to
  // 0, then a^T b is zero too).
  const double denominator = norm(a, 2) * norm(b, 2);
  if (denominator == 0.0)
    return 0;
  else
    return dot(a, b) / denominator;
}

}; // namespace kernel
}; // namespace mlpack

#endif
