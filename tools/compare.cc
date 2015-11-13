#include <boost/program_options.hpp>

#include <taLib/taLib.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {
    // OPTIONS
    
    std::string filenameExactResults;
    std::string filenameApproximateResults;
    std::string filenameQueryMatrix;
    std::string filenameProbeMatrix;
    bool querySideLeft;
    int numElements, k;
    int r, m, n;
    
    boost::program_options::options_description desc("Options");
    desc.add_options()
            ("help", "produce help message")
            ("exact", boost::program_options::value<string>(&filenameExactResults), "file containing the exact results")
            ("approximate", boost::program_options::value<string>(&filenameApproximateResults), "file containing the approximate results")
            ("Q^T", boost::program_options::value<string>(&filenameQueryMatrix), "file containing the query matrix (left side)")
            ("P", boost::program_options::value<string>(&filenameProbeMatrix), "file containing the probe matrix (right side)")
            ("querySideLeft", boost::program_options::value<bool>(&querySideLeft)->default_value(true), "1 if Q^T contains the queries (default). Interesting for Row-Top-k")
            ("k", boost::program_options::value<int>(&k), "top k")
            ("r", boost::program_options::value<int>(&r)->default_value(0), "num of coordinates in each vector (needed when reading from csv files)")
            ("m", boost::program_options::value<int>(&m)->default_value(0), "num of vectors in Q^T (needed when reading from csv files)")
            ("n", boost::program_options::value<int>(&n)->default_value(0), "num of vectors in P (needed when reading from csv files)")
            ;
    
    boost::program_options::positional_options_description pdesc;
    pdesc.add("exact", 1);
    pdesc.add("approximate", 2);
    pdesc.add("Q^T", 3);
    pdesc.add("P", 4);
    pdesc.add("k", 5);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(pdesc).run(), vm);
    boost::program_options::notify(vm);

    if (vm.count("help") || vm.count("exact") == 0 || vm.count("approximate") == 0 || vm.count("Q^T") == 0 || vm.count("P") == 0 || vm.count("k") == 0) {
        std::cout << "run compare [options] <exact> <approximate> <Q^T> <P> <k>" << std::endl << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }
    
    // COMPARING
    
    // reading matrices
    std::cout << "reading query matrix and probe matrix" << std::endl;
    ta::VectorMatrix queryMatrix, probeMatrix;
    if (querySideLeft) {
        queryMatrix.readFromFile(filenameQueryMatrix, r, m, true);
        probeMatrix.readFromFile(filenameProbeMatrix, r, n, false);
    } else {
        queryMatrix.readFromFile(filenameProbeMatrix, r, n, false);
        probeMatrix.readFromFile(filenameQueryMatrix, r, m, true);
    }
    numElements = queryMatrix.rowNum;
    
    std::ifstream fileExact(filenameExactResults.c_str(), std::ios_base::in);
    std::vector< std::vector <std::pair<int, double> > > resultsExact(numElements);
    std::ifstream fileApproximate(filenameApproximateResults.c_str(), std::ios_base::in);
    std::vector< std::vector <std::pair<int, double> > > resultsApproximate(numElements);

    if (fileExact && fileApproximate) {
        // reading exact
	std::cout << "reading " << filenameExactResults << " (exact results)" << std::endl;
        if(fileExact.good()) {
            for(int i = 0; i < numElements * k; i++) {
                int query;
                fileExact >> query;

                int id;
                fileExact >> id;

                // not used but need to read the value
                double value;
                fileExact >> value;
                
                double ip =  probeMatrix.innerProduct(id, queryMatrix.getMatrixRowPtr(query));
                
                resultsExact[query].push_back(std::make_pair(id, ip));
            }
        }
        
        // sorting exact
        std::cout << "sorting exact results (descending)" << std::endl;
        for(int i = 0; i < numElements; i++) {
            // descending
            std::sort(resultsExact[i].begin(), resultsExact[i].end(), [](const std::pair<int, double> &left, const std::pair<int, double> &right) {return left.second > right.second;});
        }
        
        // reading approximate
        std::cout << "reading " << filenameApproximateResults << " (approximate results)" << std::endl;
        if(fileApproximate.good()) {
            for(int i = 0; i < numElements * k; i++) {
                int query;
                fileApproximate >> query;

                int id;
                fileApproximate >> id;

                // not used but need to read the value
                double value;
                fileApproximate >> value;
                
                double ip =  probeMatrix.innerProduct(id, queryMatrix.getMatrixRowPtr(query));
                
                resultsApproximate[query].push_back(std::make_pair(id, ip));
            }
        }
        
        // sorting approximate
        std::cout << "sorting approximate results (descending)" << std::endl;
        for(int i = 0; i < numElements; i++) {
            // descending
            std::sort(resultsApproximate[i].begin(), resultsApproximate[i].end(), [](const std::pair<int, double> &left, const std::pair<int, double> &right) {return left.second > right.second;});
        }
        
        // calculate recall
        std::cout << "calculating recall" << std::endl;
	int same = 0;
        for(int i = 0; i < numElements; i++) {
            for(int j = 0; j < k; j++) {
                for(int n = 0; n < k; n++) {
                    if(resultsExact[i][j].first == resultsApproximate[i][n].first) {
                        same++;
			break;
                    }
		}
            }
	}
        double recall = static_cast<double>(same)/(numElements * k);

        // calculate rmse
        std::cout << "calculating rmse" << std::endl;
        double rmse = 0;
        double absoluteMax = 0;
        for(int i = 0; i < numElements; i++) {
            double sum = 0;
            for(int j = 0; j < k; j++) {
                double absoluteError = abs(resultsExact[i][j].second - resultsApproximate[i][j].second);
                
                if(absoluteError > absoluteMax) {
                    absoluteMax = absoluteError;
                }
                
                sum += absoluteError * absoluteError;
            }
            rmse += sqrt(sum/k);
        }
        double rmseAverage = rmse/numElements;

        // calculate are
        std::cout << "calculating are" << std::endl;
        double are = 0;
        double relativeMax = 0;
        for(int i = 0; i < numElements; i++) {
            double sum = 0;
            for(int j = 0; j < k; j++) {
                double relativeError = abs((resultsExact[i][j].second - resultsApproximate[i][j].second)/resultsExact[i][j].second);
                
                if(relativeError > relativeMax) {
                    relativeMax = relativeError;
                }
                
                sum += relativeError;
            }
            are += sum/k;
        }
        double areAverage = are/numElements;

        // print results
	std::cout << std::endl << same << " out of " << (numElements * k) << " results are the same" << std::endl;
        std::cout << "recall = " << recall << std::endl;
        std::cout << "rmse (average) = " << rmseAverage << std::endl;
        std::cout << "absolute error (maximum) = " << absoluteMax << std::endl;
        std::cout << "are (average) = " << areAverage << std::endl;
        std::cout << "relative error (maximum) = " << relativeMax << std::endl;
    }

}
