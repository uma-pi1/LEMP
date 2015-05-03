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
 * MatItem.h
 *
 *  Created on: Apr 1, 2012
 *      Author: olka
 */

#ifndef TA_TALIB_STRUCTS_BASICSTRUCTS_H_
#define TA_TALIB_STRUCTS_BASICSTRUCTS_H_

#include <iostream>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <cstring>



namespace ta {

    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::row_major> DenseMatrix; // for vector matrix
    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> DenseMatrixCM; // for lists
    typedef boost::numeric::ublas::matrix_row<DenseMatrix> MatrixRow;
    typedef DenseMatrix::size_type ta_size_type;

    typedef uint8_t col_type;
    typedef uint32_t row_type;
    typedef uint64_t comp_type;
    typedef uint8_t code_type; // this might need to be changed

    typedef col_type SimpleCandidate;

    comp_type scanned = 0;

    struct BestMethodTriplet {
        double time;
        row_type method;
        col_type numLists;

        BestMethodTriplet() : time(std::numeric_limits<double>::max()), method(0), numLists(0) {
        };

        BestMethodTriplet(double time, row_type method, col_type numLists) : time(time), method(method), numLists(numLists) {
        };

        bool operator==(const BestMethodTriplet & other) const {
            return time == other.time;
        }

        bool equals(const BestMethodTriplet & other)const {
            return time == other.time && method == other.method && numLists == other.numLists;
        }

        inline bool operator!=(const BestMethodTriplet & other) const {
            return !operator==(other);
        }

        bool operator<(const BestMethodTriplet & other) const {
            if (time < other.time) {
                return true;
            } else if (time == other.time) {
                if (method < other.method) {
                    return true;
                } else if (method == other.method) {
                    if (numLists < other.numLists)
                        return true;
                }
                return false;
            }
            return false;
        }

        inline bool operator<=(const BestMethodTriplet & other) const {
            return operator<(other) || operator==(other);
        }

        bool operator>(const BestMethodTriplet & other) const {
            return !operator<=(other);
        }

        inline bool operator>=(const BestMethodTriplet & other) const {
            return !operator<(other);
        }
    };

    struct MatItem {
        double result;
        ta_size_type i;
        ta_size_type j;

        MatItem(){}
        MatItem(double result, ta_size_type i, ta_size_type j) : result(result), i(i), j(j) {
        };

        bool operator==(const MatItem & other) const {
            return i == other.i && j == other.j;
        }

        bool equals(const MatItem & other)const {
            return result == other.result && i == other.i && j == other.j;
        }

        inline bool operator!=(const MatItem & other) const {
            return !operator==(other);
        }

        bool operator<(const MatItem & other) const {
            if (result < other.result) {
                return true;
            } else if (result == other.result) {
                if (j < other.j) {
                    return true;
                } else if (j == other.j) {
                    if (i < other.i)
                        return true;
                }
                return false;
            }
            return false;
        }

        inline bool operator<=(const MatItem & other) const {
            return operator<(other) || operator==(other);
        }

        bool operator>(const MatItem & other) const {
            return !operator<=(other);
        }

        inline bool operator>=(const MatItem & other) const {
            return !operator<(other);
        }
    };

    struct ResultSetItem {
        ta_size_type userId;
        ta_size_type itemId;

        inline ResultSetItem(ta_size_type userId, ta_size_type itemId) : itemId(itemId), userId(userId) {
        }

        inline ~ResultSetItem() {
        }
    };

    inline std::ostream & operator<<(std::ostream & os, const ResultSetItem & v) {
        os << "[" << v.userId << ", " << v.itemId << "] ";
        return os;
    }

    struct ResultSetMatrix {
        // 1D queryId 2D: list of itemIds for this query
        std::vector<std::vector<ta_size_type> > result;

        inline ResultSetMatrix(ta_size_type users) {
            result.resize(users);
        }

        inline ~ResultSetMatrix() {
        }

        inline std::vector<ta_size_type>& getVectorForQuery(ta_size_type queryId) {
            return result[queryId];
        }

        inline ta_size_type size() {
            ta_size_type size = 0;
            for (ta_size_type i = 0; i < result.size(); i++) {
                size += result[i].size();
            }
            return size;
        }

    };


    //std::ostream & operator<<(std::ostream & os, const MatItem & v);

    inline std::ostream & operator<<(std::ostream & os, const MatItem & v) {
        os << "(" << v.result << "  [" << v.i << ", " << v.j << "])";
        return os;
    }

    struct QueueElement {
        double data;
        ta_size_type id;

        inline QueueElement(double data, ta_size_type id) : data(data), id(id) {
        };

        inline QueueElement() : data(0.0), id(0) {
        };

        inline ~QueueElement() {
        };

        bool operator==(const QueueElement & other) const {
            if (this == &other)
                return true;
            return data == other.data && id == other.id;
        };

        inline bool operator!=(const QueueElement & other) const {
            if (this == &other)
                return false;
            return data != other.data || id != other.id;
        };

        inline bool operator<(const QueueElement & other) const {
            return data < other.data;
        }

        inline bool operator<=(const QueueElement & other) const {
            return data <= other.data;
        }

        inline bool operator>(const QueueElement & other) const {
            return data > other.data;
        }

        inline bool operator>=(const QueueElement & other) const {
            return data >= other.data;
        }



    };


    //std::ostream & operator<<(std::ostream & os, const QueueElement & v);

    inline std::ostream & operator<<(std::ostream & os, const QueueElement & v) {
        os << "(" << v.data << ", " << v.id << ")";
        return os;
    }


    //max heap for Items

    class maxHeap {
    private:
        std::vector<QueueElement> heap;
        ta_size_type size;
        //long heapifies;

        /** Returns the Vector index of the left child */
        inline ta_size_type left(ta_size_type parent) {
            return ((parent + 1) << 1) - 1;
        }

        /** Returns the Vector index of the right child */
        inline ta_size_type right(ta_size_type parent) {
            return ((parent + 1) << 1);
        }

        /** Returns the Vector index of the parent */
        inline ta_size_type parent(ta_size_type child) {
            return ((child + 1) >> 1) - 1;
        }

        inline void exchange(ta_size_type i, ta_size_type j) {
            double tempData = heap[j].data;
            ta_size_type tempId = heap[j].id;

            heap[j].data = heap[i].data;
            heap[j].id = heap[i].id;

            heap[i].data = tempData;
            heap[i].id = tempId;
        }

        /** Also known as downheap, restores the heap condition
         * starting at node i and working its way down. */
        inline void heapify(ta_size_type i) {
            ta_size_type l = left(i);
            ta_size_type r = right(i);
            ta_size_type largest;

            //heapifies++;

            if (l < size && heap[l] > heap[i]) {
                largest = l;
                //__builtin_prefetch(&heap[((largest + 1) << 1) - 1], 1, 0);
            } else {
                largest = i;
            }


            if (r < size && heap[r] > heap[largest]) {
                largest = r;
                //__builtin_prefetch(&heap[((largest + 1) << 1)], 1, 0);
            }


            if (largest != i) {
                exchange(i, largest);
                heapify(largest);
            } else {
                //std::cout<<"heap[l]: "<<heap[l]<<" heap[r]: "<<heap[r]<<std::endl;
            }
        }

    public:

        maxHeap() {
            //heapifies = 0;
        };

        ~maxHeap() {
        };

        inline QueueElement* getRoot() {
            return &heap[0];
        }

        inline void updateRoot(double newData) {
            heap[0].data = newData;
            heapify(0);
        }

        inline double getMaxValue() {
            return heap[0].data;
        }

        inline ta_size_type getMaxId() {
            return heap[0].id;
        }

        inline void clear() {
            heap.clear();
        }

        /** Inserts key into the heap, and then upheaps that key to a
         * position where the heap property is satisfied. */
        inline bool add(QueueElement key) {
            size = heap.size();
            ta_size_type i = size;
            heap.resize(size + 1);

            // upheap if necessary
            while (i > 0 && heap[parent(i)] < key) {
                heap[i] = heap[parent(i)];
                i = parent(i);
            }
            heap[i] = key;
            return true;
        }


    };

    //max heap for Items

    class minHeap {
    private:
        std::vector<QueueElement> heap;
        ta_size_type size;
        //long heapifies;

        /** Returns the Vector index of the left child */
        inline ta_size_type left(ta_size_type parent) {
            return ((parent + 1) << 1) - 1;
        }

        /** Returns the Vector index of the right child */
        inline ta_size_type right(ta_size_type parent) {
            return ((parent + 1) << 1);
        }

        /** Returns the Vector index of the parent */
        inline ta_size_type parent(ta_size_type child) {
            return ((child + 1) >> 1) - 1;
        }

        inline void exchange(ta_size_type i, ta_size_type j) {
            double tempData = heap[j].data;
            ta_size_type tempId = heap[j].id;

            heap[j].data = heap[i].data;
            heap[j].id = heap[i].id;

            heap[i].data = tempData;
            heap[i].id = tempId;
        }

        /** Also known as downheap, restores the heap condition
         * starting at node i and working its way down. */
        inline void heapify(ta_size_type i) {
            ta_size_type l = left(i);
            ta_size_type r = right(i);
            ta_size_type largest;

            //heapifies++;

            if (l < size && heap[l] < heap[i]) {
                largest = l;
                //__builtin_prefetch(&heap[((largest + 1) << 1) - 1], 1, 0);
            } else {
                largest = i;
            }


            if (r < size && heap[r] < heap[largest]) {
                largest = r;
                //__builtin_prefetch(&heap[((largest + 1) << 1)], 1, 0);
            }


            if (largest != i) {
                exchange(i, largest);
                heapify(largest);
            } else {
                //std::cout<<"heap[l]: "<<heap[l]<<" heap[r]: "<<heap[r]<<std::endl;
            }
        }

    public:

        minHeap() {
            //heapifies = 0;
        };

        ~minHeap() {
        };

        inline QueueElement* getRoot() {
            return &heap[0];
        }

        inline void updateRoot(double newData) {
            heap[0].data = newData;
            heapify(0);
        }

        inline double getMaxValue() {
            return heap[0].data;
        }

        inline ta_size_type getMinId() {
            return heap[0].id;
        }

        inline void clear() {
            heap.clear();
        }

        /** Inserts key into the heap, and then upheaps that key to a
         * position where the heap property is satisfied. */
        inline bool add(QueueElement key) {
            size = heap.size();
            ta_size_type i = size;
            heap.resize(size + 1);

            // upheap if necessary
            while (i > 0 && heap[parent(i)] > key) {
                heap[i] = heap[parent(i)];
                i = parent(i);
            }
            heap[i] = key;
            return true;
        }


    };

    struct IntervalElement {
        row_type start;
        row_type end;
        col_type col;

        IntervalElement(col_type col, row_type start, row_type end) : col(col), start(start), end(end) {
        }

        IntervalElement() : col(0), start(0), end(0) {
        }

        inline bool operator<(const IntervalElement & other) const {
            return end - start < other.end - other.start;
        }

        inline bool operator<=(const IntervalElement & other) const {
            return end - start <= other.end - other.start;
        }

        inline bool operator>(const IntervalElement & other) const {
            return end - start > other.end - other.start;
        }

        inline bool operator>=(const IntervalElement & other) const {
            return end - start >= other.end - other.start;
        }



    };

//    struct RowTopkResults {
//        MatItem* results;
//        int k;
//        row_type rowNum;
//
//        RowTopkResults() : results(0) {
//        }
//
//        ~RowTopkResults() {
//            if (!results) {
//                delete[] results;
//            }
//        }
//
//        inline void init(int _k, row_type _rowNum) {
//            k = _k;
//            rowNum = _rowNum;
//            results = new MatItem[k * rowNum];
//        }
//
//        inline void insertElement(double ip, row_type queryIndex, row_type pid) {
//
//            row_type insertPosition = 0;
//            for (; insertPosition < k; ++insertPosition)
//                if (ip >= results[queryIndex*k + insertPosition].result)//products(insertPosition, queryIndex))
//                    break;
//
//            // We only memmove() if there is actually a need to shift something.
//            if (insertPosition < (k - 1)) {
//                int len = (k - 1) - insertPosition;
//                memmove(results + queryIndex*k + insertPosition + 1,  results + queryIndex*k + insertPosition, sizeof (MatItem) * len);
//          
//            }
//
//            // Now put the new information in the right index.
//            results[queryIndex*k + insertPosition].result = ip;
//            results[queryIndex*k + insertPosition].j = pid;
//
//        }
//        
//        inline void assignElement(double ip, row_type queryIndex, row_type pos, row_type pid) {
//            results[queryIndex*k + pos] = MatItem(ip, queryIndex, pid);
//        }
//        
//        inline void sort(row_type queryIndex){
//            std::sort(results + queryIndex * k, results + (queryIndex+1) * k, std::greater<MatItem>());       
//        }
//        inline double getMinScore(row_type queryIndex){
//            return results[queryIndex * k + k-1].result;
//        }
//
//
//    };



}
#endif /* MATITEM_H_ */
