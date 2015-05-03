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


#ifndef RG_MEASURE_H
#define RG_MEASURE_H

#include <util/config.h>

#include <util/evaluation/formatting.h>
#include <util/evaluation/Timer.h>

namespace rg {

// ------------------------------------------------------------------------------------------------
// -- BenchmarkResult class
// ------------------------------------------------------------------------------------------------	

/** Describes the result of a benchmark */
class BenchmarkResult {
	/** number of repeated calls of the benchmarked function */
	int m_repetitions;
	
	/** total time required for all these calls */
	NanoDuration<double> m_totalTime; 
	
public:
	/** Constructor. */
	explicit BenchmarkResult() : m_repetitions(0), m_totalTime(0) { };
					
	/** Constructor.
	 * 
	 * @param repetitions number of repeated calls of the benchmarked function
	 * @param totalTime total time required for all these calls
	 */
	explicit BenchmarkResult(int repetitions, const Duration& totalTime) 
				: m_repetitions(repetitions), m_totalTime(totalTime.nanos()) { };
	
	/** Returns mean execution time */
	NanoDuration<double> mean() const {
		return NanoDuration<double>(m_totalTime.nanos()/m_repetitions);
	}
	
	/** Returns total execution time */
	NanoDuration<double> total() const {
		return m_totalTime;
	}
	
	/** Returns the number of repetitions */
	int repetitions() const {
		return m_repetitions;
	}		
};


// ------------------------------------------------------------------------------------------------
// -- BenchmarkResult output
// ------------------------------------------------------------------------------------------------	

/** Formatted output of a benchmark result */
template<typename CharT, typename Traits> 
std::basic_ostream<CharT, Traits>& operator<<(
		std::basic_ostream<CharT, Traits>& out, const BenchmarkResult& r) 
{
	out << "mean: ";
	outputDuration<NanoDuration<double>, CharT, Traits >(out, r.mean());
	out << " ";
	outputRate(out, r.mean());
	return out;
};

// ------------------------------------------------------------------------------------------------
// -- benchmarking macros
// ------------------------------------------------------------------------------------------------	

/** Benchmark WHAT by benchmarking the total time of a loop containing WHAT.
 * 
 * @param BENCHMARK_RESULT variable of type BenchmarkResult (holds the result)
 * @param REPETITIONS number of times WHAT is executed
 * @param WHAT code to be benchmarked
 */
#define BENCHMARK_LOOP(BENCHMARK_RESULT, REPETITIONS, WHAT) { 				\
	Timer __t__; 															\
	__t__.start(); 															\
	for (int __i__=0; __i__<REPETITIONS; __i__++) { 						\
		WHAT; 																\
	} 																		\
	__t__.stop();															\
																			\
	BENCHMARK_RESULT = BenchmarkResult(REPETITIONS, __t__.elapsedTime()); 	\
}

} // namespace rg


#endif
