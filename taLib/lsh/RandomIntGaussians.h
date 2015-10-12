
/* 
 * File:   lsh_structs.h
 * Author: chteflio
 *
 * Created on June 24, 2015, 7:42 PM
 */

#ifndef RANDOMINTGAUSSIANS_H
#define	RANDOMINTGAUSSIANS_H
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


namespace ta {

    struct RandomIntGaussians {
        //         gsl_rng *r;

        boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)

        uint16_t* intGaussians; // dimensions x nHashes 
        float *intToFloatCache;

        float leftLimit, rightLimit;
        float floatToIntFactor;
        float intToFloatFactor;

        int range; // if we only want gaussian samples between -8 and  +8, set range to 16.
        int nDimensions;
        int nHashes;

        inline RandomIntGaussians(int dimensions, int numHashes) : nDimensions(dimensions), nHashes(numHashes), range(16), //[-8,8]
        intGaussians(nullptr), intToFloatCache(nullptr) {
            assert(numHashes % 8 == 0);

            leftLimit = (0 - range) / 2.0;
            rightLimit = leftLimit + range;
            floatToIntFactor = (65535) * 1.0 / range;
            intToFloatFactor = range * 1.0 / (65535); // 65535 largest int in 16 bits

            long s = nHashes * nDimensions; //maxHashes            
            intGaussians = new uint16_t[s];

            //            r = gsl_rng_alloc(gsl_rng_taus2);
            //            gsl_rng_set(r, 123); // very much not-random ///////////////////////////////


            rng.seed(123);



            int cacheSize = 65536;
            intToFloatCache = new float[cacheSize];

            intToFloatCache[0] = leftLimit;    
            for (int i = 1; i < cacheSize; ++i) {
                intToFloatCache[i] = intToFloatCache[i - 1] + intToFloatFactor;
            }

            fill();
        }

        inline void fill() {
            float rf;

            boost::normal_distribution<> nd(0.0, 1.0);
            boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_normal(rng, nd);

            for (int b = 0; b < nHashes; ++b) {

                for (int d = 0; d < nDimensions; ++d) {
                    long long k = d * nHashes + b;

                    //                    rf = (float) gsl_ran_gaussian_ziggurat(r, 1.0);

                    rf = (float) var_normal();


                    if (rf < leftLimit) {
                        intGaussians[k] = (uint16_t) 0;
                        continue;
                    }
                    if (rf > rightLimit) {
                        intGaussians[k] = (uint16_t) 65535;
                        continue;
                    }
                    intGaussians[k] = floatToInt(rf);

                }
            }

        }

        inline float intToFloat(uint16_t a) {
            return (float) (a * intToFloatFactor + leftLimit);
        }

        inline uint16_t floatToInt(float a) {
            return (uint16_t) round((a - leftLimit) * floatToIntFactor);
        }

        ~RandomIntGaussians() {
            if (intGaussians != nullptr)
                delete[] intGaussians;
            //             gsl_rng_free(r);
            if (intToFloatCache != nullptr)
                delete[] intToFloatCache;
        }

    };




}

#endif	/* RANDOMINTGAUSSIANS_H */

