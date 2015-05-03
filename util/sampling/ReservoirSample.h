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


#ifndef RG_SAMPLING_RESERVOIRSAMPLE_H
#define RG_SAMPLING_RESERVOIRSAMPLE_H

#include <util/config.h>
#include <util/sampling/Sample.h>

#include <util/exception.h>
#include <util/math.h>
#include <util/random.h>

#include <cmath>
#include <vector>

namespace rg {

using std::vector;

// ------------------------------------------------------------------------------------------------
// -- ReservoirSample class
// ------------------------------------------------------------------------------------------------

/** 
 * This class implements regular reservoir sampling. Note that reservoir 
 * sampling does not support deletions. 
 */ 
template <class T, bool useSkipCounter = true, class Random = Random32>
class ReservoirSample {
private:
	/** desired sample size / upper bound on sample size */
	int m_targetSize;
	
	/** pseudo-random number generator */
	Random &m_random;
	
	/** number of elements in the underlying data */
	unsigned int m_N;

	/** number of records to skip before next sample insertion */
	int m_skip; 

	/** the sample */
	vector<T> m_sample;	
	
	/** constant that indicates that the skip counter m_skip must be recomputed */
	static const int SKIP_INVALID = -1;
	
public:
	/** Default constructor.
	 * 
	 * @param targetSize desired sample size / upper bound on sample size
	 * @param random a pseudo-random number generator  
	 */
	ReservoirSample(int targetSize, Random& random = Random::defaultInstance()) 
			: m_targetSize(targetSize), m_random(random), m_N(0), m_skip(SKIP_INVALID) { 
		m_sample.reserve(m_targetSize);
	};

			
	// -- Sample interface --------------------------------------------------------------
			
	void add(const T& t) {
		++m_N;
		
		if (size() < m_targetSize) { 
			// initial fill
			m_sample.push_back(t);
		} else if (!useSkipCounter) {
			// simple algorithm without skip counters
			int index = m_random.nextInt(m_N);
			if (index < m_targetSize) {
				m_sample[index] = t;
			}
		} else { 
			// use skip counters
			if (m_skip == SKIP_INVALID) { // generate
				m_skip = rsSkipZ(m_random, m_targetSize, m_N);	
			}
			if (m_skip > 0) {
				--m_skip;
			} else { // accept
				int index = m_random.nextInt((uint32_t)m_targetSize);
				m_sample[index] = t;
				m_skip = SKIP_INVALID;
			}
		}
	}
	
	/** Not implemented.
	 * 
	 * @throw UnsupportedOperationException
	 */
	void remove(const T& t) {
		RG_THROW(UnsupportedOperationException, "plain reservoir sampling does not support deletions");
	}
	
	void clear() {
		m_sample.clear();
		m_N = 0;
		m_skip = SKIP_INVALID;		
	}
	
	int size() const {
		return m_sample.size();
	}

	// -- iterators -------------------------------------------------------------------------------
	
	/** type of constant iterator over the sample items */
	typedef typename vector<T>::const_iterator const_iterator;
	
	/** Returns an iterator referring to the first element in the sample.
	 * 
	 * @return iterator referring to the first element 
	 */
	const_iterator begin() const {
		return m_sample.begin();
	};
	
	/** Returns an iterator referring to the past-the-end element in the 
	 * sample.
	 * 
	 * @return iterator referring to the past-the-end element 
	 */
	const_iterator end() const {
		return m_sample.end();
	}

	
	
}; // class ReservoirSample


// ------------------------------------------------------------------------------------------------
// -- Skip counter computation for reservoir sampling 
// ------------------------------------------------------------------------------------------------

/** Returns the number of tuples which would by skipped by reservoir sampling before the next
 * update of the sample using algorithm R (=the naive algorithm). 
 * 
 * See algorithm R in "Random Sampling With A Reservoir" by J.S. Vitter (1985) for more details.
 * 
 * @param random a pseudo-random number generator
 * @param n sample size
 * @param N current population size + 1
 * @return the number of tuples skipped before the next update of the sample 
 */ 
template<class Random>
inline unsigned int rsSkipR(Random& random, unsigned int n, unsigned int N) {
	int S = 0;
	while (n <= random.nextInt(N+S)) {
		S++;
	}
	return S;		
}

/** Returns the number of tuples which would by skipped by reservoir sampling before the next
 * update of the sample using algorithm X (less calls to PRNG than rsSkipR). Note that for fast 
 * random number generators (such as the mersenne twister), rsSkipR is faster than rsSkipX. 
 * 
 * See algorithm X in "Random Sampling With A Reservoir" by J.S. Vitter (1985) for more details.
 * 
 * @param random a pseudo-random number generator
 * @param n sample size
 * @param N current population size + 1
 * @return the number of tuples skipped before the next update of the sample
 */
template<class Random>
inline unsigned int rsSkipX(Random& random, unsigned int n, unsigned int N) {
	// that's not in the original paper but speeds up things
	//if (n<<2 <= N) return skipR(random, n, N);
	
	// original skipX algorithm
	double V = random.nextDouble();
	int S = 0;
	int num = N-n;
	double quot = (double)num / N;
	while (quot > V) {
		++S;
		++N;
		++num;
		quot = (quot*num)/N;
	}
	return S;
}

/** Returns the number of tuples which would by skipped by reservoir sampling before the next
 * update of the sample using algorithm Z (advanced algorithm). This is the fastest algorithm 
 * and should be preferred over rsSkipR and rsSkipX.
 * 
 * See algorithm Z in "Random Sampling With A Reservoir" by J.S. Vitter (1985) for more details.
 * 
 * @param random a pseudo-random number generator
 * @param n sample size
 * @param N current population size + 1
 * @return the number of tuples skipped before the next update of the sample
 */
template<class Random>
inline unsigned int rsSkipZ(Random& random, unsigned int n, unsigned int N) {
	const double THRESHOLD = 22; // threshold is important for speed and to avoid overflows (!)
	const int t = N-1;
	if (t <= THRESHOLD*n) {
		//return skipX(random, n, N); 
		return rsSkipR(random, n, N); // superior with fast PRNGs
	} else {
		int S;
		int term = t-n+1;
		double W = exp(-log(random.nextDouble())/n);
		while (true) {
			double U = random.nextDouble();
			double H = t*(W-1.0);
			S = (int)floor(H);

			// test1
			double lhs = exp(log(((U*(square((t+1)/term)))*(term+S))/(t+H))/n);
			double rhs = (((t+H)/(term+S))*term)/t;
			if (lhs <= rhs) break;
			
			// test2
			double y = (((U*(t+1))/term)*(t+S+1))/(t+H);
			int denom;
			int numer_lim;  
			if (n<S) {
				denom = t;
				numer_lim = term+S;
			} else {
				denom = t-n+S;
				numer_lim = t+1;
			}
			for (int numer=t+S; numer>=numer_lim; numer--) {
				y = (y*numer)/denom;
				--denom;
			}
			if (exp(log(y)/n)<= (t+H)/t) break;
			W = exp(-log(random.nextDouble())/n);
		}
		return S;
	}
};


} // namespace

#endif
