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

#ifndef RG_RANDOM_BOOSTRANDOMWRAPPER_H
#define RG_RANDOM_BOOSTRANDOMWRAPPER_H

#include <util/config.h>
#include <util/exception.h>

#include <limits>

#include <boost/random.hpp>

namespace rg {

/** 
 * This class wraps the pseudo-random number generators (PRNG) in order to add useful utility 
 * functions. The wrapper works only with (and is carefully optimized for) PRNGs that produce 
 * random numbers in the full 32 bit integer range. 
 *  
 * Many algorithms and classes in librg are parametrized with a random number generator.
 * This class defines the interface required for these generators. It also complies to 
 * RandomNumberGenerator concept of the C++ standard library.
 *
 * @tparam PRNG A 32-bit pseudo-random number generator. The underlying PRNG must have a default 
 * 			    constructor and a constructor that takes an initial seed as argument. The 32bit 
 *              random number are obtained using the PRNG::operator(). All PRNGs in the Boost library 
 *              comply to these requirements.
 * @tparam Seed type of seed (numeric)
 */
template<class PRNG, class Seed>
class RandomWrapper32 {
public:
	typedef PRNG Prng;

private:

	/** The underlying generator */
	PRNG m_prng;
	
	/** the default instance */
	static RandomWrapper32<PRNG, Seed> m_defaultInstance;
	
public:
	/** Creates a new instance of the pseudo-random number generator using the default seed. */
	RandomWrapper32() : m_prng() { };

	/** Creates a new instance of the pseudo-random number generator using the specified seed. 
	 * 
	 * @param seed the seed of the pseudo-random-number generator 
	 */
	explicit RandomWrapper32(Seed seed) : m_prng(seed) { }

	/** Returns a pseudorandom, uniformly distributed integer value between 0 (inclusive) and the 
	 * specified value (exclusive), drawn from this random number generator's sequence. 
	 * 
	 * @param n the bound on the random number to be returned. Must be non-negative.
	 * @return a pseudorandom, uniformly distributed integer value 
	 * 		   between 0 (inclusive) and n (exclusive).
	 * 
	 * @throws InvalidArgumentException when n is negative
	 */
	int32_t nextInt(int32_t n) {
		if (n < 0) RG_THROW(InvalidArgumentException, "n must be positive");
		return operator()(n); // conversions are ok
	}

	/** Returns a pseudorandom, uniformly distributed integer value between 0 (inclusive) and the
	 * specified value (exclusive), drawn from this random number generator's sequence.
	 *
	 * @param n the bound on the random number to be returned.
	 * @return a pseudorandom, uniformly distributed integer value between 0 (inclusive)
	 * 		   and n (exclusive).
	 *
	 * TODO: needs improvement
	 */
	int64_t nextInt(int64_t n) {
		return (int64_t)(nextDouble()*n);
	}

	/** Returns a pseudorandom, uniformly distributed integer value between 0 (inclusive) and the 
	 * specified value (exclusive), drawn from this random number generator's sequence. 
	 * 
	 * @param n the bound on the random number to be returned.
	 * @return a pseudorandom, uniformly distributed integer value between 0 (inclusive) 
	 * 		   and n (exclusive).		 
	 */
	uint32_t nextInt(uint32_t n) {
		return operator()(n);
	}

	/** Returns a pseudorandom, uniformly distributed integer value between 0 (inclusive) and the
	 * specified value (exclusive), drawn from this random number generator's sequence.
	 *
	 * @param n the bound on the random number to be returned.
	 * @return a pseudorandom, uniformly distributed integer value between 0 (inclusive)
	 * 		   and n (exclusive).
	 *
	 * TODO: needs improvement
	 */
	uint64_t nextInt(uint64_t n) {
		return (uint64_t)(nextDouble()*n);
	}

	/** Returns the next pseudorandom, uniformly distributed float value between 0.0 
	 * and 1.0 from this random number generator's sequence. */
	float nextFloat() {
		// compile time constants
		const uint32_t andf = (1 << std::numeric_limits<float>::digits) - 1;
		const float divf = 1 << std::numeric_limits<float>::digits;

		return (operator()() & andf) / divf;
	}
	
	/** Returns the next pseudorandom, uniformly distributed double value between 0.0 
	 * and 1.0 from this random number generator's sequence.
	 *
	 * @return the next pseudorandom, uniformly distributed double value 
	 *         between 0.0 and 1.0 from this random number generator's sequence 
	 */
	double nextDouble() {
		// compile time constants
		const uint64_t andf = ((uint64_t)1 << std::numeric_limits<double>::digits) - 1;
		const double divf = (uint64_t)1 << std::numeric_limits<double>::digits;
		
		uint64_t u = ((uint64_t)operator()() << 32) + operator()();
		return (u & andf) / divf;
	}
	
	/** Returns a pseudorandom, uniformly distributed 32bit unsigned integer value, 
	 * drawn from this random number generator's sequence. 
	 * 
	 * @return a pseudorandom, uniformly distributed 32bit unsigned integer value 
	 */
	uint32_t operator()() {
		return m_prng(); 
	};

	/** Returns a pseudorandom, uniformly distributed 32bit unsigned integer value 
	 * between 0 (inclusive) and the specified value (exclusive), drawn from this 
	 * random number generator's sequence. 
	 * 
	 * @param n the bound on the random number to be returned.
	 * @return a pseudorandom, uniformly distributed 32bit unsigned integer value 
	 * 		   between 0 (inclusive) and n (exclusive).		 
	 */
	uint32_t operator()(uint32_t n) {
		// adapted from Sun's java implementation (Random#nextInt(int))		
		if ((n & -n) == n)  // i.e., n is a power of 2
			return ((uint64_t)n * operator()()) >> 32;

		uint32_t bits, val;
		do {
			bits = operator()();
			val = bits % n;
		} while ((uint32_t)(bits - val + (n-1U)) < bits - val); // acceptance-rejection step  
		return val;
	};
	
	/** Returns the default instance of this pseudo-random number generator. This method always returns
	 * the same generator when called repeatedly. */
	static RandomWrapper32<PRNG, Seed>& defaultInstance() {
		return m_defaultInstance;
	}

	/** Returns the underlying random number generator */
	PRNG& prng() {
		return m_prng;
	}
};

/** Initialize default instance */
template<class PRNG, class Seed>
RandomWrapper32<PRNG, Seed> RandomWrapper32<PRNG, Seed>::m_defaultInstance
	= RandomWrapper32<PRNG, Seed>((Seed)time(NULL));

// for windows: replace boost::uint32_t of mt19937 by uint32_t to avoid (invalid) ambiguity warnings
typedef boost::random::mersenne_twister<uint32_t,32,624,397,31,0x9908b0df,11,
  7,0x9d2c5680,15,0xefc60000,18, 3346425566U> my_mt19937;

/** The default random type used within this library */   
typedef RandomWrapper32<my_mt19937, uint32_t> Random32; 

} // namespace rg

#endif
