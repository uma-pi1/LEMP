/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Results.h
 * Author: bertsch
 *
 * Created on February 12, 2016, 1:46 PM
 */

#ifndef RESULTS_H
#define RESULTS_H

namespace mips {

    class Results {
    public:

        std::vector< std::vector<MatItem > > resultsVector;

        Results() {

        }

        //resize according to number or threads

        //assign results from retrievalArguments

        inline ~Results() {

        }

        inline void writeToFile(std::string& file) {
            std::ofstream out(file.c_str());

            if (!out.is_open())
                std::cout << "[WARNING] File " << file << " is not open!" << std::endl;
            else
                std::cout << "[INFO] Results will be written to " << file << std::endl;

            for (auto& threadResult : resultsVector) {
                //            std::sort(threadResult.begin(), threadResult.end(), std::less<MatItem>());

                for (auto& result : threadResult) {
                    out << result.i << " " << result.j << " " << result.result << std::endl;
                }
            }

            out.close();
        }

        void moveAppend(std::vector<MatItem>& src, int tid) {
            if (resultsVector[tid].empty()) {
                resultsVector[tid] = std::move(src);
            } else {
                resultsVector[tid].reserve(resultsVector[tid].size() + src.size());
                std::move(std::begin(src), std::end(src), std::back_inserter(resultsVector[tid]));
                src.clear();
            }
        }

        void clear() {
            for (auto& threadResult : resultsVector) {
                threadResult.clear();
            }
        }

        void clearAll() {
            resultsVector.clear();
        }

        comp_type getResultSize() {
            comp_type counter = 0;
            for (int t = 0; t < resultsVector.size(); ++t) {
                counter += resultsVector[t].size();
            }
            return counter;
        }



    };

}

#endif /* RESULTS_H */

