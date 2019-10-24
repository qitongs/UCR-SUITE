/* 
 * File:   utils.h
 * Author: michele, Qitong Wang
 *
 * Created on October 25, 2014, 8:41 AM
 * Modified on October 24, 2019, 15:03 PM
 */

#ifndef UTILS_H
#define	UTILS_H

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))
#define dist(x,y) ((x-y)*(x-y))

/// Pseudo infinite number for this code
#define INF 1e20

using namespace std;

template <typename T>
void sort_indexes(const T * values, int * indexes, int length) {
    iota(indexes, indexes + length, 0);
    sort(indexes, indexes + length, [values](int i1, int i2) {return values[i1] - values[i2];});
}

/// Print function for debugging
template <typename T>
void print_array(T * x, int len) {
    for(int i = 0; i < len; ++ i)
        cout << x[i] << " ";
    cout << endl;
}

/// If expected error happens, teminated the program.
void error(int id) {
    switch (id) {
        case 1: cout << "ERROR: Memory can't be allocated" << endl; break;
        case 2: cout << "ERROR: File not Found" << endl; break;
        case 3: cout << "ERROR: Can't create Output File" << endl; break;
        case 4: cout << "ERROR: Invalid Arguments" << endl; break;
        case 5: cout << "ERROR: Input file broken" << endl; break;
        default: cout << "ERROR: Un-managed" << endl;
    }
    exit(1);
}

#endif	/* UTILS_H */
