//    Copyright 2015 Christina Teflioudi
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

/*
 * Candidates.h
 *
 *  Created on: Nov 7, 2013
 *      Author: chteflio
 */

#ifndef CANDIDATES_H_
#define CANDIDATES_H_


namespace ta {

    struct Candidate_incr {
#ifdef  WITH_SIMD_INCR
        double x[2] __attribute__((aligned(16))); //x[0]: sum pi*qi x[1]: sum pi^2
#else
        double x[2]; //x[0]: sum pi*qi x[1]: sum pi^2        
#endif

        inline Candidate_incr() {
        }

        inline void addFirst(double qi, double pi) {

#ifdef  WITH_SIMD_INCR
            __m128d xmm1 = _mm_set1_pd(pi);
            __m128d xmm2 = _mm_setr_pd(qi, pi);

            _mm_store_pd(x, _mm_mul_pd(xmm1, xmm2));

#else
            x[0] = pi*qi;
            x[1] = pi*pi;
#endif
        }

        inline void add(double qi, double pi) {

#ifdef  WITH_SIMD_INCR
            __m128d xmm1 = _mm_set1_pd(pi);
            __m128d xmm2 = _mm_setr_pd(qi, pi);
            __m128d sum = _mm_load_pd(x);

            sum = _mm_add_pd(sum, _mm_mul_pd(xmm1, xmm2));
            _mm_store_pd(x, sum);

#else
            x[0] += pi*qi;
            x[1] += pi*pi;
#endif

        }

        inline bool prune(double len, double theta, double seenQi2) {


            double len2 = len * len;
            double rightPart = theta - x[0] * len;

            if (rightPart < 0 || (1 - x[1]) * seenQi2 * len2 >= rightPart * rightPart)
                return false;

            return true;



//            x[0] = theta - x[0] * len;
//
//            if (x[0] < 0 || (1 - x[1]) * seenQi2 * len * len >= x[0] * x[0])
//                return false;
//
//            return true;



        }

    };

}


#endif /* CANDIDATES_H_ */
