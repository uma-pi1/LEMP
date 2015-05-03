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


#ifndef RG_SAMPLING_SHUFFLE_H
#define RG_SAMPLING_SHUFFLE_H

#include <iterator>
#include <util/io/utils.h>

namespace rg {

/** Shuffles data in place using Knuth shuffle */
template<typename RandomAccessIterator>
void shuffle(RandomAccessIterator begin, RandomAccessIterator end, rg::Random32& random) {
	typedef typename std::iterator_traits<RandomAccessIterator>::value_type T;
	typedef typename std::iterator_traits<RandomAccessIterator>::difference_type D;
	const D n = end - begin;
	if (n==0) return;
	T temp = *begin;
	for (D i=0; i<n; i++) {
		D j = random.nextInt(n-i);
		temp = *begin;
		*begin = *(begin+j);
		*(begin+j) = temp;
		begin++;
	}
};

}

#endif
