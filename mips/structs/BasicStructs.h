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



namespace mips {

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

    struct MatItem {
        double result;
        ta_size_type i;
        ta_size_type j;

        inline ~MatItem() = default;
        inline MatItem() = default;

        inline MatItem(double result, ta_size_type i, ta_size_type j) : result(result), i(i), j(j) {
        };

        inline bool operator==(const MatItem & other) const {
            return i == other.i && j == other.j;
        }

        inline bool equals(const MatItem & other)const {
            return result == other.result && i == other.i && j == other.j;
        }

        inline bool operator!=(const MatItem & other) const {
            return !operator==(other);
        }

        inline bool operator<(const MatItem & other) const {
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

        inline bool operator>(const MatItem & other) const {
            return !operator<=(other);
        }

        inline bool operator>=(const MatItem & other) const {
            return !operator<(other);
        }
    };

    inline std::ostream & operator<<(std::ostream & os, const MatItem & v) {
        os << "(" << v.result << "  [" << v.i << ", " << v.j << "])";
        return os;
    }

    struct QueueElement {
        double data = 0.0;
        ta_size_type id = 0;

        inline QueueElement(double data, ta_size_type id) : data(data), id(id) {
        };

        inline QueueElement() = default;

        inline ~QueueElement() = default;

        inline bool operator==(const QueueElement & other) const {
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

            if (l < size && heap[l] > heap[i]) {
                largest = l;
            } else {
                largest = i;
            }
            if (r < size && heap[r] > heap[largest]) {
                largest = r;
            }
            if (largest != i) {
                exchange(i, largest);
                heapify(largest);
            }
        }

    public:

        maxHeap() = default;

        ~maxHeap() = default;

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

        inline ta_size_type getMaxId() const {
            return heap[0].id;
        }

        inline void clear() {
            heap.clear();
        }

        /** Inserts key into the heap, and then upheaps that key to a
         * position where the heap property is satisfied. */
        inline bool add(QueueElement key) { //QueueElement&& ?
            size = heap.size();
            ta_size_type i = size;
            heap.resize(size + 1);

            // upheap if necessary
            while (i > 0 && heap[parent(i)] < key) {
                heap[i] = heap[parent(i)];
                i = parent(i);
            }
            heap[i] = key;
            size = heap.size();
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

    
    
    /* These functions are needed for the tuning part
    */
    inline void calculateTimeInCutoff(const std::vector<double>& method1, const std::vector<double>& method2, row_type cutoffSample, double& time) {
        if (cutoffSample == 0) {
            time = 0;
            for (auto& element : method1)
                time += element;

        } else {
            time -= method1[cutoffSample - 1];
            time += method2[cutoffSample - 1]; // other can be LENGTH for example
        }
    }

    inline void findCutOffPoint(const std::vector<double>& method1, const std::vector<double>& method2, double& bestTime, row_type& best_t_b_ind) {
        double time;
        for (row_type i = 0; i < method1.size(); ++i) {
            calculateTimeInCutoff(method1, method2, i, time);

            if ((i == 0) || bestTime > time) {
                bestTime = time;
                best_t_b_ind = i;
            }

        }
    }

    inline void findCutOffPointInv(const std::vector<double>& method1, const std::vector<double>& method2, double& bestTime, row_type& best_t_b_ind) {
        double time = bestTime;
        for (int i = method1.size() - 1; i >= 0; i--) {
            time -= method1[i];
            time += method2[i];
            if (bestTime > time) {
                bestTime = time;
                best_t_b_ind = i;
            }
        }
    }


}
#endif /* MATITEM_H_ */
