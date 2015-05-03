//    Copyright 2015 Rainer Gemulla
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

#ifndef RG_SAMPLING_AGGREGATOR_H
#define RG_SAMPLING_AGGREGATOR_H

#include <util/config.h>

namespace rg {

/** Interface for an aggregator. An aggregator accumulates its input values 
 * in order to compute aggregate information. Some aggregagtors also support 
 * removal of values. */
template <typename T>
struct Aggregator {
	/** Accumulate the given value into the aggregate. 
	 * 
	 * @param t a value to be accumulated 
	 */
	void add(const T& t);

	/** Remove the given value from the aggregate (optional). This method 
	 * assumes that the value to be removed has been inserted before.
	 * 
	 * @param t a value to be removed
	 * @throw UnsupportedOperationException when not implemented by the specific
	 * 				aggregator
	 */
	void remove(const T& t);
	
	/** Clear the aggregator. This method resets the aggregator to its initial
	 * state. */
	void clear();
	
	/** Type of the values that are aggregated */ 
	typedef T value_type;
};
	
} // namespace rg
	
#endif 
