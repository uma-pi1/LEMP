
/* 
 * File:   lsh_structs.h
 * Author: chteflio
 *
 * Created on June 24, 2015, 7:42 PM
 */

#ifndef RANDOMINTGAUSSIANS_H
#define	RANDOMINTGAUSSIANS_H

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


namespace ta {


    // This is a singleton class

    class RandomIntGaussians {
        boost::mt19937 rng;

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
            rng.seed(123); // not that random
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

        RandomIntGaussians(const RandomIntGaussians& rhs) = delete;
        RandomIntGaussians& operator=(const RandomIntGaussians& rhs) = delete;


    public:

        float *intToFloatCache;
        uint16_t* intGaussians; // dimensions x nHashes

        static RandomIntGaussians& getInstance(int dimensions, int numHashes) {// call as RandomIntGaussians::getInstance() gives you functionality similar to a global variable
            static RandomIntGaussians rig(dimensions, numHashes); // static: protection from data races
            return rig;
        }

        ~RandomIntGaussians() {
            if (intGaussians != nullptr)
                delete[] intGaussians;
            if (intToFloatCache != nullptr)
                delete[] intToFloatCache;
        }

    };

    class ActiveLshRepetitionsForTheta {// somehow calling function on singleton takes too much time 
        std::vector<QueueElement> thetaRepetitions; // id: repetitions needed, double:theta


        ActiveLshRepetitionsForTheta& operator=(const ActiveLshRepetitionsForTheta& rhs) = delete;
        ActiveLshRepetitionsForTheta(const ActiveLshRepetitionsForTheta& rhs) = delete;

    public:

        ActiveLshRepetitionsForTheta() {
        }

        inline void init(double recallLevel) {
            thetaRepetitions.reserve(LSH_SIGNATURES);
            double logNom = log(1 - recallLevel);
            double exponent = 1 / ((double) LSH_CODE_LENGTH); //0.125; // 1/8 LSH_CODE_LENGTH

            for (int b = 0; b < LSH_SIGNATURES; ++b) {
                double logDenom = logNom / (b + 1);
                double denom = exp(logDenom);
                double thres8 = 1 - denom;
                double thres = pow(thres8, exponent);
                double acosThres = (1 - thres) * PI;
                double theta = cos(acosThres);

                thetaRepetitions.emplace_back(theta, b + 1);
            }
            std::sort(thetaRepetitions.begin(), thetaRepetitions.end(), std::less<QueueElement>());
        }

        inline int findActiveBlocks(double theta) {
            auto it = std::upper_bound(thetaRepetitions.begin(), thetaRepetitions.end(), QueueElement(theta, 0));

            int pos = it - thetaRepetitions.begin();
            if (pos >= thetaRepetitions.size()) {
                pos = thetaRepetitions.size() - 1;
                return thetaRepetitions[pos].id;
            }

            int b = thetaRepetitions[pos].id + 1;
            return b;
        }


    };


    ActiveLshRepetitionsForTheta repetitionsForTheta;


}

#endif	/* RANDOMINTGAUSSIANS_H */

