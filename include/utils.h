/* 
 * File:   utils.h
 * Author: Qitong Wang
 *
 * Modified on October 24, 2019, 15:03 PM
 */

#ifndef UTILS_H
#define UTILS_H

#include <iostream>

#define min(x, y) ((x)<(y)?(x):(y))
#define max(x, y) ((x)>(y)?(x):(y))
#define dist(x, y) ((x-y)*(x-y))

#define INF 1e20

using namespace std;

void error(int id) {
    switch (id) {
        case 1:
            cout << "ERROR: Memory can't be allocated" << endl;
            break;
        case 2:
            cout << "ERROR: File not Found" << endl;
            break;
        case 3:
            cout << "ERROR: Can't create Output File" << endl;
            break;
        case 4:
            cout << "ERROR: Invalid Arguments" << endl;
            break;
        case 5:
            cout << "ERROR: Input file broken" << endl;
            break;
        case 6:
            cout << "ERROR: Couldn't get CPU info" << endl;
            break;
        default:
            cout << "ERROR: Un-managed" << endl;
    }
    exit(1);
}

#endif    /* UTILS_H */
