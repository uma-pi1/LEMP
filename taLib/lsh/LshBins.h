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
 * File:   LshBins.h
 * Author: chteflio
 *
 * Created on June 27, 2015, 11:22 AM
 */

#ifndef LSHBINS_H
#define	LSHBINS_H
#include <boost/dynamic_bitset.hpp>

namespace ta {

    class LshBins {
    public:
        long long userOffset;
        row_type bytesPerCode;
        row_type numHashTables;
        row_type nVectors;
        long long numBinsPerHashTable;

        inline LshBins() {
        }

        inline void init(int bytesPerCode_, int numHashTables_, row_type nVectors_) {
            bytesPerCode = bytesPerCode_;
            numHashTables = numHashTables_;
            nVectors = nVectors_;
            userOffset = bytesPerCode * numHashTables;
            numBinsPerHashTable = pow(2, bytesPerCode * 8);
        }

        inline virtual void printBins() {
        }

        inline virtual void resizeBins(row_type size) {
        }

        inline virtual void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {
        }

        inline virtual void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
        }

        inline virtual void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {
        }

        inline virtual void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {
        }


    };

    class LshBins8 : public LshBins {
    public:
        //        std::vector<unordered_map<uint8_t, std::vector<row_type> > > bins;
        std::vector<row_type > binsOffsets;
        std::vector<row_type> data;

        inline LshBins8() {
        }

        inline virtual void printBins() {
            //            for (int block = 0; block < bins.size(); block++) {
            //                unordered_map<uint8_t, std::vector<row_type> >::iterator it;
            //                std::cout << "block: " << block << std::endl;
            //                for (it = bins[block].begin(); it != bins[block].end(); it++) {
            //                    std::cout << "hv: " << (int) it->first << " size: " << it->second.size() << std::endl;
            //                }
            //            }
        }

        inline virtual void resizeBins(row_type size) {
            //            bins.resize(size);
        }

        //        inline void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {
        //            bins.resize(endBlock);
        //            row_type startOffset = bytesPerCode * startBlock;
        //
        //            for (row_type i = 0; i < nVectors; i++) {
        //                int t = i * userOffset + startOffset;
        //
        //                for (int block = startBlock; block < endBlock; block++) {
        //                    uint8_t key = 0;
        //                    memcpy(&key, sketch + t, bytesPerCode);
        //                    memset(sketch + t, 0, bytesPerCode);
        //
        //                    t += bytesPerCode;
        //                    bins[block][key].push_back(i);
        //
        //                }
        //            }
        //        }
        //        inline void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
        //                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {
        //            done.reset();
        //            long long t = queryPos * userOffset;
        //            unordered_map<uint8_t, std::vector<row_type> >::iterator it;
        //
        //            for (row_type block = 0; block < activeBlocks; block++) {
        //
        //                uint8_t key = 0;
        //                memcpy(&key, querySketches + t, bytesPerCode);
        //
        //
        //                t += bytesPerCode;
        //                it = bins[block].find(key);
        //
        //                if (it != bins[block].end()) {
        //
        //                    for (row_type i = 0; i < it->second.size(); i++) {
        //                        row_type uId = it->second.at(i);
        //
        //                        if (done[uId])
        //                            continue;
        //
        //                        done[uId] = true;
        //
        //                        candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
        //                        numCandidatesToVerify++;
        //                    }
        //                }
        //
        //            }
        //        }
        //        inline void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
        //            if (bins.size() < endBlock)
        //                bins.resize(endBlock);
        //
        //            row_type startOffset = bytesPerCode * startBlock;
        //            int t = startOffset;
        //
        //            for (int block = startBlock; block < endBlock; block++) {
        //                uint8_t key = 0;
        //                memcpy(&key, sketch + t, bytesPerCode);
        //                memset(sketch + t, 0, bytesPerCode);
        //                t += bytesPerCode;
        //                bins[block][key].push_back(probeId);
        //            }
        //        }
        //        inline void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
        //                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {
        //
        //            unordered_map<uint8_t, std::vector<row_type> >::iterator it;
        //
        //            uint8_t key = 0;
        //            memcpy(&key, querySketches + t, bytesPerCode);
        //
        //
        //            it = bins[binNum].find(key);
        //
        //            if (it != bins[binNum].end()) {
        //
        //                for (row_type i = 0; i < it->second.size(); i++) {
        //                    row_type uId = it->second.at(i);
        //
        //                    if (done[uId])
        //                        continue;
        //
        //                    done[uId] = true;
        //
        //                    candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
        //                    numCandidatesToVerify++;
        //                }
        //            }
        //        }

        inline void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {

            data.resize(endBlock * nVectors);
            row_type numValuesPerBlock = 256;
            binsOffsets.resize(endBlock * (numValuesPerBlock + 1));


            // update hashBuckets
            for (int block = startBlock; block < endBlock; block++) {
                row_type startOffset = bytesPerCode * block;
                int blockOffsetOnData = block * nVectors;

                for (row_type i = 0; i < nVectors; i++) {
                    long long t = i * userOffset + startOffset;
                    uint8_t key = sketch[t];
                    //                    std::cout<<" hv: "<<(int)key<<" t: "<<t<<std::endl;
                    countsOfBlockValues[key]++;
                }

                binsOffsets[block * (numValuesPerBlock + 1) ] = 0;
                row_type startOffset2 = block * (numValuesPerBlock + 1);
                for (int value = 1; value <= numValuesPerBlock; value++) {
                    binsOffsets[startOffset2 + value] = binsOffsets[startOffset2 + value - 1] + countsOfBlockValues[value - 1];
                }

                for (row_type i = 0; i < nVectors; i++) {
                    long long t = i * userOffset + startOffset;
                    uint8_t key = sketch[t];
                    sketch[t] = 0;
                    int j = binsOffsets[startOffset2 + key + 1] - countsOfBlockValues[key];
                    countsOfBlockValues[key]--;
                    data[blockOffsetOnData + j] = i;
                }
            }

        }

        inline void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
//             if (bins.size() < endBlock)
//                 bins.resize(endBlock);
// 
//             row_type startOffset = bytesPerCode * startBlock;
//             int t = startOffset;
// 
//             for (int block = startBlock; block < endBlock; block++) {
//                 uint8_t key = 0;
//                 memcpy(&key, sketch + t, bytesPerCode);
//                 memset(sketch + t, 0, bytesPerCode);
//                 t += bytesPerCode;
//                 bins[block][key].push_back(probeId);
//             }
        }

        inline void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {

            done.reset();
            long long t = queryPos * userOffset;
            row_type offset2 = 256 + 1;
            row_type offset = 0;
	    

            for (row_type block = 0; block < activeBlocks; block++) {

                uint8_t key = querySketches[t];
                t++;

                row_type start = binsOffsets[block * offset2 + key];
                row_type end = binsOffsets[block * offset2 + key + 1];
                //                row_type offset = block * nVectors;
		

                for (row_type j = start; j < end; j++) { // find all with same hash value in probe

                    row_type uId = data[offset + j];

                    if (done[uId])
                        continue;

                    done[uId] = true;

                    candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                    numCandidatesToVerify++;
                }
                offset += nVectors;

            }
        }

        inline void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {

            //             long long t = queryPos * userOffset + binNum;
            row_type offset2 = 256 + 1;
            row_type offset = nVectors * binNum;

            uint8_t key = querySketches[t];

            row_type start = binsOffsets[binNum * offset2 + key];
            row_type end = binsOffsets[binNum * offset2 + key + 1];

            for (row_type j = start; j < end; j++) { // find all with same hash value in probe

                row_type uId = data[offset + j];

                if (done[uId])
                    continue;

                done[uId] = true;

                candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                numCandidatesToVerify++;
            }
        }

    };

    class LshBins16 : public LshBins {
    public:
        std::vector<unordered_map<uint16_t, std::vector<row_type> > > bins;

        inline LshBins16() {
        }

        inline void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {
            bins.resize(endBlock);
            row_type startOffset = bytesPerCode * startBlock;


            for (row_type i = 0; i < nVectors; i++) {
                long long t = i * userOffset + startOffset;

                for (int block = startBlock; block < endBlock; block++) {
                    uint16_t key = 0;
                    memcpy(&key, sketch + t, bytesPerCode);
                    memset(sketch + t, 0, bytesPerCode);

                    t += bytesPerCode;
                    bins[block][key].push_back(i);
                }
            }
        }

        inline void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
            //            std::cout << "in populateBinsSingle" << std::endl;

            row_type startOffset = bytesPerCode * startBlock;
            int t = startOffset;

            for (int block = startBlock; block < endBlock; block++) {
                uint16_t key = 0;
                memcpy(&key, sketch + t, bytesPerCode);
                memset(sketch + t, 0, bytesPerCode);

                t += bytesPerCode;
                bins[block][key].push_back(probeId);
            }
        }

        inline virtual void printBins() {
            for (int block = 0; block < bins.size(); block++) {
                unordered_map<uint16_t, std::vector<row_type> >::iterator it;
                std::cout << "block: " << block << std::endl;
                for (it = bins[block].begin(); it != bins[block].end(); it++) {
                    std::cout << "hv: " << (int) it->first << " size: " << it->second.size() << std::endl;
                }
            }
        }

        inline virtual void resizeBins(row_type size) {
            bins.resize(size);
        }

        inline void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {

            done.reset();
            long long t = queryPos * userOffset;
            unordered_map<uint16_t, std::vector<row_type> >::iterator it;
            for (row_type block = 0; block < activeBlocks; block++) {

                uint16_t key = 0;
                memcpy(&key, querySketches + t, bytesPerCode);
                t += bytesPerCode;

                it = bins[block].find(key);

                if (it != bins[block].end()) {

                    for (row_type i = 0; i < it->second.size(); i++) {
                        row_type uId = it->second.at(i);

                        if (done[uId])
                            continue;

                        done[uId] = true;

                        candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                        numCandidatesToVerify++;
                    }
                }
            }
        }

        inline void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {


            //             long long t = queryPos * userOffset + bytesPerCode*binNum;
            unordered_map<uint16_t, std::vector<row_type> >::iterator it;

            uint16_t key = 0;
            memcpy(&key, querySketches + t, bytesPerCode);


            it = bins[binNum].find(key);

            if (it != bins[binNum].end()) {

                for (row_type i = 0; i < it->second.size(); i++) {
                    row_type uId = it->second.at(i);

                    if (done[uId])
                        continue;

                    done[uId] = true;

                    candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                    numCandidatesToVerify++;
                }
            }
        }


    };

    class LshBins32 : public LshBins {
    public:
        std::vector<unordered_map<uint32_t, std::vector<row_type> > > bins;

        inline LshBins32() {
        }

        inline virtual void printBins() {
            for (int block = 0; block < bins.size(); block++) {
                unordered_map<uint32_t, std::vector<row_type> >::iterator it;
                std::cout << "block: " << block << std::endl;
                for (it = bins[block].begin(); it != bins[block].end(); it++) {
                    std::cout << "hv: " << (int) it->first << " size: " << it->second.size() << std::endl;

                }
            }
        }

        inline void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {
            bins.resize(endBlock);
            row_type startOffset = bytesPerCode * startBlock;

            std::cout << "bytesPerCode: " << bytesPerCode << std::endl;


            for (row_type i = 0; i < nVectors; i++) {
                long long t = i * userOffset + startOffset;

                for (int block = startBlock; block < endBlock; block++) {
                    uint32_t key = 0;
                    memcpy(&key, sketch + t, bytesPerCode);
                    memset(sketch + t, 0, bytesPerCode);
                    t += bytesPerCode;
                    bins[block][key].push_back(i);

                }
            }
        }

        inline virtual void resizeBins(row_type size) {
            bins.resize(size);
        }

        inline void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
            if (bins.size() < endBlock)
                bins.resize(endBlock);

            row_type startOffset = bytesPerCode * startBlock;
            int t = startOffset;

            for (int block = startBlock; block < endBlock; block++) {
                uint32_t key = 0;
                memcpy(&key, sketch + t, bytesPerCode);
                memset(sketch + t, 0, bytesPerCode);
                t += bytesPerCode;
                bins[block][key].push_back(probeId);
            }
        }

        inline void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {

            done.reset();
            long long t = queryPos * userOffset;
            unordered_map<uint32_t, std::vector<row_type> >::iterator it;

            for (row_type block = 0; block < activeBlocks; block++) {

                uint32_t key = 0;
                memcpy(&key, querySketches + t, bytesPerCode);


                t += bytesPerCode;
                it = bins[block].find(key);

                if (it != bins[block].end()) {

                    for (row_type i = 0; i < it->second.size(); i++) {
                        row_type uId = it->second.at(i);

                        if (done[uId])
                            continue;

                        done[uId] = true;

                        candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                        numCandidatesToVerify++;
                    }
                }

            }
        }

        inline void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {

            //             long long t = queryPos * userOffset + bytesPerCode*binNum;
            unordered_map<uint32_t, std::vector<row_type> >::iterator it;

            uint32_t key = 0;
            memcpy(&key, querySketches + t, bytesPerCode);


            it = bins[binNum].find(key);

            if (it != bins[binNum].end()) {

                for (row_type i = 0; i < it->second.size(); i++) {
                    row_type uId = it->second.at(i);

                    if (done[uId])
                        continue;

                    done[uId] = true;

                    candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                    numCandidatesToVerify++;
                }
            }
        }
    };

    class LshBins64 : public LshBins {
    public:
        std::vector<unordered_map<uint64_t, std::vector<row_type> > > bins;

        inline LshBins64() {
        }

        inline virtual void printBins() {
            for (int block = 0; block < bins.size(); block++) {
                unordered_map<uint64_t, std::vector<row_type> >::iterator it;
                std::cout << "block: " << block << std::endl;
                for (it = bins[block].begin(); it != bins[block].end(); it++) {
                    std::cout << "hv: " << (int) it->first << " size: " << it->second.size() << std::endl;

                }
            }
        }

        inline void populateBins(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues) {
            bins.resize(endBlock);
            row_type startOffset = bytesPerCode * startBlock;


            for (row_type i = 0; i < nVectors; i++) {
                long long t = i * userOffset + startOffset;

                for (int block = startBlock; block < endBlock; block++) {
                    uint64_t key = 0;
                    memcpy(&key, sketch + t, bytesPerCode);
                    memset(sketch + t, 0, bytesPerCode);
                    t += bytesPerCode;
                    bins[block][key].push_back(i);

                }
            }
        }

        inline virtual void resizeBins(row_type size) {
            bins.resize(size);
        }

        inline void populateBinsSingle(uint8_t* sketch, row_type startBlock, row_type endBlock, std::vector<row_type>& countsOfBlockValues, row_type probeId) {
            if (bins.size() < endBlock)
                bins.resize(endBlock);

            row_type startOffset = bytesPerCode * startBlock;
            int t = startOffset;

            for (int block = startBlock; block < endBlock; block++) {
                uint64_t key = 0;
                memcpy(&key, sketch + t, bytesPerCode);
                memset(sketch + t, 0, bytesPerCode);
                t += bytesPerCode;
                bins[block][key].push_back(probeId);
            }
        }

        inline void getCandidates(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type activeBlocks, row_type probeBucketStartPos) {

            done.reset();
            long long t = queryPos * userOffset;
            unordered_map<uint64_t, std::vector<row_type> >::iterator it;

            for (row_type block = 0; block < activeBlocks; block++) {

                uint64_t key = 0;
                memcpy(&key, querySketches + t, bytesPerCode);


                t += bytesPerCode;
                it = bins[block].find(key);

                if (it != bins[block].end()) {

                    for (row_type i = 0; i < it->second.size(); i++) {
                        row_type uId = it->second.at(i);

                        if (done[uId])
                            continue;

                        done[uId] = true;

                        candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                        numCandidatesToVerify++;
                    }
                }

            }
        }

        inline void getCandidatesOfSingleBin(uint8_t* querySketches, row_type queryPos, std::vector<row_type>& candidatesToVerify, row_type& numCandidatesToVerify,
                boost::dynamic_bitset<>& done, row_type binNum, row_type probeBucketStartPos, long long t) {

            //             long long t = queryPos * userOffset + bytesPerCode*binNum;
            unordered_map<uint64_t, std::vector<row_type> >::iterator it;

            uint64_t key = 0;
            memcpy(&key, querySketches + t, bytesPerCode);


            it = bins[binNum].find(key);

            if (it != bins[binNum].end()) {

                for (row_type i = 0; i < it->second.size(); i++) {
                    row_type uId = it->second.at(i);

                    if (done[uId])
                        continue;

                    done[uId] = true;

                    candidatesToVerify[numCandidatesToVerify] = uId + probeBucketStartPos;
                    numCandidatesToVerify++;
                }
            }
        }
    };

}

#endif	/* LSHBINS_H */

