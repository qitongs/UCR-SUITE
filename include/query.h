//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_QUERY_H
#define UCR_SUITE_QUERY_H

#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>

#include "utils.h"
#include "bounds.h"

class Query {
public:
    int length, overlap_length, warping_window;
    int *sorted_indexes;
    double *values, *sorted_values, *sorted_upper_envelop, *sorted_lower_envelop;

    Query(ifstream &query_ifs, int length, int warping_window, int overlap_length);

    ~Query();

private:
    static void sort_indices(const double *values, int *indexes, int length);
};

Query::Query(ifstream &query_ifs, int length, int warping_window, int overlap_length) : length(length),
                                                                                        warping_window(warping_window),
                                                                                        overlap_length(overlap_length) {
    auto upper_envelop = (double *) malloc(sizeof(double) * length);
    auto lower_envelop = (double *) malloc(sizeof(double) * length);

    this->sorted_indexes = (int *) malloc(sizeof(int) * length);
    this->values = (double *) malloc(sizeof(double) * length);
    this->sorted_values = (double *) malloc(sizeof(double) * length);
    this->sorted_upper_envelop = (double *) malloc(sizeof(double) * length);
    this->sorted_lower_envelop = (double *) malloc(sizeof(double) * length);

    if (this->sorted_indexes == nullptr || this->values == nullptr || upper_envelop == nullptr ||
        lower_envelop == nullptr || this->sorted_values == nullptr || this->sorted_upper_envelop == nullptr ||
        this->sorted_lower_envelop == nullptr) {
        error(1);
    }

    for (int i = 0; i < length; ++i) {
        query_ifs >> values[i];
    }

    get_envelops_lemire(this->values, length, warping_window, lower_envelop, upper_envelop);
    sort_indices(this->values, this->sorted_indexes, length);

    for (int i = 0; i < length; ++i) {
        this->sorted_values[i] = this->values[this->sorted_indexes[i]];
        this->sorted_upper_envelop[i] = upper_envelop[this->sorted_indexes[i]];
        this->sorted_lower_envelop[i] = lower_envelop[this->sorted_indexes[i]];
    }

    free(upper_envelop);
    free(lower_envelop);
}

Query::~Query() {
    free(this->sorted_indexes);
    free(this->values);
    free(this->sorted_values);
    free(this->sorted_upper_envelop);
    free(this->sorted_lower_envelop);
}

void Query::sort_indices(const double *values, int *indexes, int length) {
    iota(indexes, indexes + length, 0);
    sort(indexes, indexes + length, [values](int i1, int i2) { return abs(values[i1]) > abs(values[i2]); });
}

#endif //UCR_SUITE_QUERY_H
