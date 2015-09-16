/**
 * @file fastmks_rules_impl.hpp
 * @author Ryan Curtin
 *
 * Implementation of FastMKSRules for cover tree search.
 *
 * This file is part of mlpack 1.0.12.
 *
 * mlpack is free software; you may redstribute it and/or modify it under the
 * terms of the 3-clause BSD license.  You should have received a copy of the
 * 3-clause BSD license along with mlpack.  If not, see
 * http://www.opensource.org/licenses/BSD-3-Clause for more information.
 */
#ifndef __MY_MLPACK_METHODS_FASTMKS_FASTMKS_RULES_IMPL_HPP
#define __MY_MLPACK_METHODS_FASTMKS_FASTMKS_RULES_IMPL_HPP

// In case it hasn't already been included.
#include "fastmks_rules.hpp"
#include <taLib/structs/BasicStructs.h>


namespace ta {
namespace fastmks {

template<typename TreeType>
FastMKSRules<TreeType>::FastMKSRules(
		KernelType& kernel,
		VectorMatrix* probeMatrix, VectorMatrix* queryMatrix) :
		kernel(kernel),
		lastQueryIndex(-1),
		lastReferenceIndex(-1),
		lastKernel(0.0),
		baseCases(0),
		scores(0),
		probeMatrix(probeMatrix), queryMatrix(queryMatrix)
		{}



template<typename TreeType>
inline force_inline
double FastMKSRules<TreeType>::BaseCase(
		const size_t queryIndex,
		const size_t referenceIndex,
                std::vector<QueueElement>& results)
		{

	// Score() always happens before BaseCase() for a given node combination.  For
	// cover trees, the kernel evaluation between the two centroid points already
	// happened.  So we don't need to do it.  Note that this optimizes out if the
	// first conditional is false (its result is known at compile time).

	if ((queryIndex == lastQueryIndex) &&
			(referenceIndex == lastReferenceIndex))
		return lastKernel;

	// Store new values.
	lastQueryIndex = queryIndex;
	lastReferenceIndex = referenceIndex;


	++baseCases;

	const double* query = queryMatrix->getMatrixRowPtr(queryIndex);

	double kernelEval = probeMatrix->innerProduct(referenceIndex, query);

	lastKernel = kernelEval;


	// If this is a better candidate, insert it into the list.
	if (kernelEval < results.front().data)
		return kernelEval;

	std::pop_heap (results.begin(), results.end(), std::greater<QueueElement>());
	results.pop_back();
	results.push_back(QueueElement(kernelEval,  probeMatrix->getId(referenceIndex)));
	std::push_heap (results.begin(), results.end(), std::greater<QueueElement>());


	return kernelEval;


		}

template<typename TreeType>
inline force_inline
double FastMKSRules<TreeType>::BaseCaseForTheta(
		const size_t queryIndex,
		const size_t referenceIndex,
        std::vector<MatItem>& results,
		const double theta)
		{


	// Score() always happens before BaseCase() for a given node combination.  For
	// cover trees, the kernel evaluation between the two centroid points already
	// happened.  So we don't need to do it.  Note that this optimizes out if the
	// first conditional is false (its result is known at compile time).
	if ((queryIndex == lastQueryIndex) &&
			(referenceIndex == lastReferenceIndex))
		return lastKernel;

	// Store new values.
	lastQueryIndex = queryIndex;
	lastReferenceIndex = referenceIndex;


	++baseCases;
	double kernelEval = probeMatrix->innerProduct(referenceIndex, queryMatrix->getMatrixRowPtr(queryIndex));

	// Update the last kernel value, if we need to.

	lastKernel = kernelEval;

	// If this is a better candidate, insert it into the list.
	if (kernelEval < theta)
		return kernelEval;

	results.push_back(MatItem(kernelEval, queryMatrix->getId(queryIndex), probeMatrix->getId(referenceIndex)));

	return kernelEval;


		}

template<typename TreeType>
double FastMKSRules<TreeType>::Score(const size_t queryIndex,
		TreeType& referenceNode, std::vector<QueueElement>& results)
		{

	const double bestKernel = results.front().data;

	// See if we can perform a parent-child prune.
	const double furthestDist = referenceNode.FurthestDescendantDistance();
	if (referenceNode.Parent() != NULL)
	{
            
		double maxKernelBound;
		const double parentDist = referenceNode.ParentDistance();
		const double combinedDistBound = parentDist + furthestDist;
		const double lastKernel = referenceNode.Parent()->Stat().LastKernel(omp_get_thread_num());

		maxKernelBound = lastKernel +
				combinedDistBound * queryMatrix->getVectorLength(queryIndex);

		if (maxKernelBound < bestKernel)
			return DBL_MAX;
	}

	// Calculate the maximum possible kernel value, either by calculating the
	// centroid or, if the centroid is a point, use that.
	++scores;
	double kernelEval;

	// Could it be that this kernel evaluation has already been calculated?
	if (tree::TreeTraits<TreeType>::HasSelfChildren &&
			referenceNode.Parent() != NULL &&
			referenceNode.Point(0) == referenceNode.Parent()->Point(0))
	{
		kernelEval = referenceNode.Parent()->Stat().LastKernel(omp_get_thread_num());
	}
	else
	{
		kernelEval = BaseCase(queryIndex, referenceNode.Point(0), results);
	}


	referenceNode.Stat().LastKernel(omp_get_thread_num()) = kernelEval;

	double maxKernel;

	maxKernel = kernelEval + furthestDist * queryMatrix->getVectorLength(queryIndex);


	// We return the inverse of the maximum kernel so that larger kernels are
	// recursed into first.
	return (maxKernel > bestKernel) ? (1.0 / maxKernel) : DBL_MAX;
        
		}

//my code here
template<typename TreeType>
double FastMKSRules<TreeType>::ScoreForTheta(const size_t queryIndex,
		TreeType& referenceNode, std::vector<MatItem>& results, const double theta)
		{
 
	// See if we can perform a parent-child prune.
	const double furthestDist = referenceNode.FurthestDescendantDistance();
        
	if (referenceNode.Parent() != NULL)
	{
		const double parentDist = referenceNode.ParentDistance();
		const double combinedDistBound = parentDist + furthestDist;
		const double lastKernel = referenceNode.Parent()->Stat().LastKernel(omp_get_thread_num());

		double maxKernelBound = lastKernel +
				combinedDistBound * queryMatrix->getVectorLength(queryIndex);


		if (maxKernelBound < theta)
			return DBL_MAX;
	}
	// Calculate the maximum possible kernel value, either by calculating the
	// centroid or, if the centroid is a point, use that.
	++scores;
	double kernelEval;


	// Could it be that this kernel evaluation has already been calculated?
	if (tree::TreeTraits<TreeType>::HasSelfChildren &&
			referenceNode.Parent() != NULL &&
			referenceNode.Point(0) == referenceNode.Parent()->Point(0))
	{
		kernelEval = referenceNode.Parent()->Stat().LastKernel(omp_get_thread_num());
	}
	else
	{
		kernelEval = BaseCaseForTheta(queryIndex, referenceNode.Point(0), results, theta);
	}
      
	referenceNode.Stat().LastKernel(omp_get_thread_num()) = kernelEval;

	double maxKernel = kernelEval + furthestDist * queryMatrix->getVectorLength(queryIndex);
    

	// We return the inverse of the maximum kernel so that larger kernels are
	// recursed into first.
	return (maxKernel >= theta) ? (1.0 / maxKernel) : DBL_MAX;
		}


template<typename TreeType>
double FastMKSRules<TreeType>::Rescore(const size_t queryIndex,
		TreeType& /*referenceNode*/,
		const double oldScore, std::vector<QueueElement>& results) const
		{
	const double bestKernel = results.front().data;

	return ((1.0 / oldScore) >= bestKernel) ? oldScore : DBL_MAX;
		}

//my code here
template<typename TreeType>
double FastMKSRules<TreeType>::RescoreForTheta(TreeType& /*referenceNode*/,
		const double oldScore, const double theta) const
		{
	return ((1.0 / oldScore) >= theta) ? oldScore : DBL_MAX;
		}





}; // namespace fastmks
}; // namespace ta

#endif
