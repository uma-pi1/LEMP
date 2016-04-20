/* 
 * File:   Mip.h
 * Author: chteflio
 *
 * Created on February 4, 2016, 7:25 PM
 */

#ifndef MIP_H
#define	MIP_H


namespace mips {

    class Mip {
    protected:
        VectorMatrix probeMatrix;
        std::shared_ptr<VectorMatrix> probeMatrix_ptr;
        double dataPreprocessingTimeLeft = 0;
        double dataPreprocessingTimeRight = 0;
        double tuningTime = 0;
        double retrievalTime = 0;
        rg::Timer timer;
        comp_type totalComparisons = 0;

        /**
         * Format allows to log multiple run-methods and multiple algorithms to a single logfile:
         * 
         * algorithm threads
         * P rows cols
         * Q^T/thetaOrK rows cols thetaOrK resultsSize comparisons preprocessingTime tuningTime retrievalTime totalTime
         */
        std::ofstream logging;

        inline void outputStats() {
            std::cout << "[STATS] comparisons = " << totalComparisons << std::endl;
            logging << totalComparisons << "\t";

            double dataPreprocessingTime = dataPreprocessingTimeLeft + dataPreprocessingTimeRight;

            std::cout << "[STATS] preprocessing time = " << (dataPreprocessingTime / 1E9) << "s" << std::endl;
            logging << dataPreprocessingTime / 1E9 << "\t";

            std::cout << "[STATS] tuning time = " << (tuningTime / 1E9) << "s" << std::endl;
            logging << tuningTime / 1E9 << "\t";

            std::cout << "[STATS] retrieval time = " << (retrievalTime / 1E9) << "s" << std::endl;
            logging << retrievalTime / 1E9 << "\t";

            std::cout << "[STATS] total time = " << ((dataPreprocessingTime + tuningTime + retrievalTime) / 1E9) << "s" << std::endl;
            logging << (dataPreprocessingTime + tuningTime + retrievalTime) / 1E9 << "\n";
        }
    public:

        Mip() {
        }

        virtual ~Mip() {

        }
        

        inline virtual void runAboveTheta(VectorMatrix& leftMatrix, Results& results) {
            std::cerr << "You should not have called that." << std::endl;
            exit(1);
        }

        inline virtual void runTopK(VectorMatrix& leftMatrix, Results& results) {
            std::cerr << "You should not have called that." << std::endl;
            exit(1);
        }

        inline void clearTimings() {
            dataPreprocessingTimeLeft = 0;
            tuningTime = 0;
            retrievalTime = 0;
            totalComparisons = 0;
        }

    };

}

#endif	/* MIP_H */
