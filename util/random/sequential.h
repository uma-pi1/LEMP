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


#ifndef RG_SAMPLING_SEQUENTIAL_H
#define RG_SAMPLING_SEQUENTIAL_H

#include <util/io/utils.h>

namespace rg {

/** Returns the number of tuples which would by skipped before the next update
 * of the sample using algorithm A.
 * See algorithm A in "Faster Methods for Random Sampling" by J.S. Vitter, 1984
 * for more details about record skipping.
 *
 * @param n Number of tuples remaining from the sample.
 * @param N Number of tuples remaining in the population
 * @return the number of tuples skipped before the next update of the sample
 */
template<typename UInt>
inline UInt skipSequentialA(rg::Random32& random, UInt n, UInt N) {
	UInt S;

	if(n >= 2) { // Algorithm A
		double V = random.nextDouble();
		S = 0;
		UInt top = N - n;
		double quot = (double)top / N;
		while (quot > V) {
			S++;
			top--;
			N--;
			quot = (quot * top) / N;
		}
	} else if (n == 1) { // select one at random (optimization)
		S = random.nextInt(N);
	} else { // the sample is complete
		S = -1;
	}

	return S;
}

/** Returns the number of tuples which would by skipped before the next update
 * of the sample using algorithm D.
 * See algorithm D in "Faster Methods for Random Sampling" by J.S. Vitter, 1984
 * for more details about record skipping.
 *
 * @param n Number of tuples remaining from the sample.
 * @param N Number of tuples remaining in the population
 * @return the number of tuples skipped before the next update of the sample
 */
template<typename UInt>
inline UInt skipSequentialD(rg::Random32& random, UInt n, UInt N) {
	if (n <= 1 || 13*n >= N) {
		return skipSequentialA(random, n, N);
	}

	double nreal = n;
	double ninv = 1.0/nreal;
	double Nreal = N;
	double Vprime = exp(log(random.nextDouble())*ninv);
	UInt qu1 = N -n + 1;
	double qu1real = qu1;
	double nmin1inv = 1.0/(-1.0+nreal);
	UInt S;
	while (true) {
		// D2
		double X;
		while (true) {
			X = Nreal*(-Vprime + 1.0);
			S = (UInt)floor(X);
			if (S < qu1) break;
			Vprime = exp(log(random.nextDouble())*ninv);
		}
		double U = random.nextDouble();
		double negSreal = -S;

		// D3
		double y1 = exp(log(U*Nreal/qu1real)*nmin1inv);
		Vprime = y1 * (-X/Nreal+1)*(qu1real/(negSreal+qu1real));
		if (Vprime <= 1.0) break;

		// D4
		double y2 = 1;
		double top = -1.0 + Nreal;
		double bottom;
		UInt limit;
		if (n-1 > S) {
			bottom = -nreal + Nreal;
			limit = -S+N;
		} else {
			bottom = -1.0 + negSreal + Nreal;
			limit = qu1;
		}
		for (UInt t=N-1; t>=limit; t--) {
			y2 = (y2*top)/bottom;
			top--;
			bottom--;
		}
		if (Nreal/(Nreal-X) >= y1*exp(log(y2)*nmin1inv)) {
			break;
		}
		Vprime = exp(log(random.nextDouble())*ninv);
	}

	// D5
	return S;
}

template<typename UInt>
inline UInt skipSequential(rg::Random32& random, UInt n, UInt N) {
	return skipSequentialD(random, n, N);
}

/** Samples n values from the interval [0,N) without replacement. The resulting vector is sorted.
 *
 * @param n sample size
 * @param N population size
 */
template<typename UInt>
std::vector<UInt> sample(Random32& random, UInt n, UInt N) {
	if (n>N) RG_THROW(rg::InvalidArgumentException,
			paste("n=", n, " must not be larger than N=", N));
	std::vector<UInt> result;

	// generate
	UInt current = 0;
	while (n > 0) {
		UInt skip = skipSequential(random, n, N-current);
		current += skip;
		result.push_back(current);
		current++;
		n--;
	}
	return result;
}

}

#endif
