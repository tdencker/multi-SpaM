/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */
//#ifndef MERGESORT_HPP_
//#define MERGESORT_HPP_
//
//#ifdef _OPENMP
//#include <algorithm>
//#include <omp.h>
//
//template<typename Iter>
//void mergeSortBody(Iter begin, Iter end, int n)
//{
//    auto len = std::distance(begin, end);
//    if (len <= 1024 || n < 2)
//    {
//        std::sort(begin,end);
//        return;
//    }
//    Iter mid = std::next(begin, len/2);
//    #pragma omp task
//    mergeSortBody(begin, mid, n-2);
//    #pragma omp task
//    mergeSortBody(mid, end, n-2);
//    #pragma omp taskwait
//    std::inplace_merge(begin, mid, end);
//}
//template<typename Iter>
//void mergeSort(Iter begin, Iter end)
//{
//    #pragma omp parallel
//    {
//        #pragma omp single
//        mergeSortBody(begin, end, omp_get_num_threads());
//    }
//}
//#else
//template<typename Iter>
//void mergeSort(Iter begin, Iter end)
//{
//    std::sort(begin, end);
//}
//#endif
//
//#endif
