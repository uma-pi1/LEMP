/* 
 * File:   Output.h
 * Author: chteflio
 *
 * Created on February 4, 2016, 12:32 PM
 */

#ifndef OUTPUT_H
#define	OUTPUT_H

#include <unordered_map>

namespace mips{
    
    inline void getCorrectResultsFromFile(std::vector<std::unordered_map<row_type, bool> >& correctResults,
            const std::string& filename, int numQueries, int k) {


        correctResults.resize(numQueries);

        std::ifstream in(filename.c_str(), std::ios_base::in);

        if (!in.is_open()) {
            std::cout << "[WARNING] File " << filename << " is not open!" << std::endl;
        } else {
            std::cout << "[INFO] Correct results will be read from " << filename << std::endl;
        }

        if (in.good()) {
            for (int i = 0; i < numQueries * k; i++) {
                row_type query;
                in >> query;

                row_type id;
                in >> id;

                // not used but need to read the value
                double value;
                in >> value;

                correctResults[query][id] = true;
            }
        }

    }

    inline double getRecall(const std::vector<std::unordered_map<row_type, bool> >& correctResults,
            const std::vector<std::vector<MatItem> >& results) {

        comp_type hits = 0;
        comp_type total = 0;

        for (int t = 0; t < results.size(); ++t) {
            for (MatItem result : results[t]) {
                total++;

                if (correctResults[result.i].count(result.j) > 0) {
                    hits++;
                }

            }
        }

        return ((double) hits) / total;

    }


}

#endif	/* OUTPUT_H */

