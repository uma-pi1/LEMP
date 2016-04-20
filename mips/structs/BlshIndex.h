/* 
 * File:   BlshIndex.h
 * Author: chteflio
 *
 * Created on March 4, 2016, 3:19 PM
 */

#ifndef BLSHINDEX_H
#define	BLSHINDEX_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_errno.h>

namespace mips {
    

    class BlshIndex : public LshIndex {

        int computeMinMatches(int nTrials) {
            double p;
            double epsilon = 1 - R;
            int guess = nTrials / 2;
            int start = 0;
            int end = nTrials;

            //             std::cout<<"threshold: "<<threshold<<" epsilon: "<<epsilon<<std::endl;

            int retValue = -1;
            do {
                p = 1.0 - posteriorCdf(guess, nTrials, minMathesT);
                if (p > epsilon) {
                    if (start == end) {
                        retValue = start - 1;
                        break;
                    }
                    end = guess - 1;
                    if (end < start) {
                        retValue = guess - 1;
                        break;
                    }
                } else if (p < epsilon) {
                    if (start == end) {
                        retValue = start;
                        break;
                    }
                    start = guess + 1;
                    if (start > end) {
                        retValue = guess;
                        break;
                    }
                } else {
                    retValue = guess;
                    break;
                }
                guess = (start + end) / 2;
                if (end < 0 || start > nTrials) {
                    std::cout << "Problem in computeMinMatches start " << start << " end " << end << std::endl;
                    exit(1);
                }

            } while (1);

            if (retValue < 0)
                retValue = 0;
            else if (retValue > nTrials)
                retValue = nTrials;

            return retValue;
        }


        inline double posteriorCdf(double s, double n, double x) {
            if (x >= 1.0)
                return 1.0;
            if (x <= 0.5)
                return 0;
            double b1 = 1.0;
            double bHalf = boost::math::ibeta(s + 1, n - s + 1, 0.5);
            double bx = boost::math::ibeta(s + 1, n - s + 1, x);
            double den = b1 - bHalf;
            if (den < 1.0e-15)
                return exp(log(bx - bHalf) - log(den));
            else
                return (bx - bHalf) / den;


        }

    public:
        double minMathesT; // threshold for matching, r = c2r(simT)
        long long hashGroups = 0;
        int32_t *minMatches; // minimum number of hashes that should be observed to meet simT   
        double R;
        double worst;

        BlshIndex() : LshIndex(), minMathesT(0), minMatches(nullptr) {
        }

        inline ~BlshIndex() {
            if (minMatches != nullptr)
                delete[] minMatches;
        }

        inline void allocateBayesLSHMemory(double worstCaseTheta) {
            if (minMatches == nullptr) {
                // set up cache space and pre-compute all minimum matches               
                minMatches = da_i32malloc(hashGroups, NULL); //"allocateBayesLSHMemory: minMatches"
                

                //                                std::cout << "worstCaseTheta: " << worstCaseTheta << std::endl;                 
                minMathesT = (1.0 - acos(worstCaseTheta) * INVPI); // min matches threshold
                for (int i = 0; i < hashGroups; i++) {
                    minMatches[i] = computeMinMatches((i + 1) * 32);
                    //                                        std::cout<<minMatches[i]<<" ";
                }
                //                                std::cout <<  std::endl;

            }

        }

        inline void initializeLists(const VectorMatrix& matrix, double worstCaseTheta, bool forProbeVectors, double recall,
                ta_size_type start = 0, ta_size_type end = 0) {

            omp_set_lock(&writelock);

            if (!initialized) {
                R = recall;
                end = (end == 0 ? matrix.rowNum : end);
                row_type nVectors = end - start;
                cosSketches = new CosineSketches(nVectors, LSH_CODE_LENGTH, LSH_SIGNATURES);
                hashGroups = LSH_CODE_LENGTH * LSH_SIGNATURES / 32;
                initializedSketches.resize(nVectors, 0);
                cosSketches->alloc();

                if (forProbeVectors) {
                    switch (LSH_CODE_LENGTH) {
                        case 8:
                            lshBins = new LshBinsDense();
                            break;
                        case 16:
                            lshBins = new LshBinsSparse<uint16_t>();
                            break;
                        case 24:
                        case 32:
                            lshBins = new LshBinsSparse<uint32_t>();
                            break;
                        default:
                            lshBins = new LshBinsSparse<uint64_t>();
                            break;
                    }

                    lshBins->init(cosSketches->bytesPerCode, cosSketches->numHashTables, cosSketches->nVectors);
                }

                if (worstCaseTheta < std::numeric_limits<double>::max()) {
                    if (worstCaseTheta > 0) {
                        worstCaseTheta /= matrix.getVectorLength(start);
                        worstCaseTheta = (worstCaseTheta > 1 ? 1 : worstCaseTheta); // it will break in the loop afterwards

                    } else {
                        worstCaseTheta /= matrix.getVectorLength(end - 1);
                        worstCaseTheta = (worstCaseTheta < -1 ? -1 : worstCaseTheta); // it will have to check everything
                    }

                    worst = worstCaseTheta;
                    allocateBayesLSHMemory(worstCaseTheta);
                }


                initialized = true;
            }
            omp_unset_lock(&writelock);
        }

    };

}


#endif	/* BLSHINDEX_H */

