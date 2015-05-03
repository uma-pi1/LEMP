
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

#ifndef RG_TIMEVAL_H
#define RG_TIMEVAL_H

#include <util/config.h>

#include <util/portability/inttypes.h>
#include <util/portability/time.h>
#include <util/evaluation/duration.h>

#include <iostream>
#include <iomanip>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- timeval conversions
// ------------------------------------------------------------------------------------------------	
namespace detail {
	inline uint64_t timevalToNanos(const timeval& v) {
		return (uint64_t)v.tv_sec*1000000000U + (uint64_t)v.tv_usec*1000U;
	}
}

// ------------------------------------------------------------------------------------------------
// -- timeval operators 
// ------------------------------------------------------------------------------------------------	

/** Add two timeval's (without checking input). */
inline timeval operator+(const timeval& v1, const timeval& v2) {
	timeval result;
	if (v1.tv_usec + v2.tv_usec < 1000000) {
		result.tv_sec = v1.tv_sec + v2.tv_sec;
		result.tv_usec = v1.tv_usec + v2.tv_usec;
	} else {
		result.tv_sec = v1.tv_sec + v2.tv_sec + 1;
		result.tv_usec = v1.tv_usec + v2.tv_usec - 1000000;
	}
	return result;
}

/** Add a timeval (without checking input). */
inline timeval& operator+=(timeval& v1, const timeval& v2) {
	if (v1.tv_usec + v2.tv_usec < 1000000) {
		v1.tv_sec += v2.tv_sec;
		v1.tv_usec += v2.tv_usec;
	} else {
		v1.tv_sec += v2.tv_sec + 1;
		v1.tv_usec += v2.tv_usec - 1000000;
	}
	return v1;
}

/** Substract two timeval's (without checking input). */
inline timeval operator-(const timeval& v1, const timeval& v2) {
	timeval result;
	if (v1.tv_usec >= v2.tv_usec) {
		result.tv_sec = v1.tv_sec - v2.tv_sec;
		result.tv_usec = v1.tv_usec - v2.tv_usec;
	} else {
		result.tv_sec = v1.tv_sec - v2.tv_sec - 1;
		result.tv_usec = v1.tv_usec + (1000000 - v2.tv_usec);
	}
	return result;
}

// ------------------------------------------------------------------------------------------------
// -- timeval/duration operators 
// ------------------------------------------------------------------------------------------------	

/** Add a timeval to a duration. */
template<class NanoType>
inline NanoDuration<NanoType>& operator+=(NanoDuration<NanoType>& d, const timeval& v) {
	return d += NanoDuration<NanoType>(detail::timevalToNanos(v));
}

} // namespace rg

#endif
