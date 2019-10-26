//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_QUERY_H
#define UCR_SUITE_QUERY_H

#include <cmath>
#include <cstdio>
#include <fstream>
#include "utils.h"
#include "bounds.h"

class Query {
public:
    int length;
    int *sorted_indexes;
    double *normalized_points, *sorted_normalized_points, *sorted_upper_envelop, *sorted_lower_envelop;

    Query(ifstream &query_ifs, int length, int warping_window);

    ~Query();
};

Query::Query(ifstream &query_ifs, int length, int warping_window) : length(length) {
    auto points = (double *) malloc(sizeof(double) * length);
    auto upper_envelop = (double *) malloc(sizeof(double) * length);
    auto lower_envelop = (double *) malloc(sizeof(double) * length);

    this->sorted_indexes = (int *) malloc(sizeof(int) * length);
    this->normalized_points = (double *) malloc(sizeof(double) * length);
    this->sorted_normalized_points = (double *) malloc(sizeof(double) * length);
    this->sorted_upper_envelop = (double *) malloc(sizeof(double) * length);
    this->sorted_lower_envelop = (double *) malloc(sizeof(double) * length);

    if (this->sorted_indexes == nullptr || points == nullptr || this->normalized_points == nullptr ||
        upper_envelop == nullptr || lower_envelop == nullptr || this->sorted_normalized_points == nullptr ||
        this->sorted_upper_envelop == nullptr || this->sorted_lower_envelop == nullptr) {
        error(1);
    }

    double sum = 0, squared_sum = 0, mean, std;

    for (int i = 0; i < length; ++i) {
        query_ifs >> points[i];
        sum += points[i];
        squared_sum += (points[i] * points[i]);
    }

    mean = sum / length;
    std = sqrt(squared_sum / length - mean * mean);
    for (int i = 0; i < length; ++i) {
        this->normalized_points[i] = (points[i] - mean) / std;
    }

    lower_upper_lemire(this->normalized_points, length, warping_window, lower_envelop, upper_envelop);
    sort_indexes(&this->normalized_points, this->sorted_indexes, length);

    for (int i = 0; i < length; ++i) {
        this->sorted_normalized_points[i] = this->normalized_points[this->sorted_indexes[i]];
        this->sorted_upper_envelop[i] = upper_envelop[this->sorted_indexes[i]];
        this->sorted_lower_envelop[i] = lower_envelop[this->sorted_indexes[i]];
    }

    free(points);
    free(upper_envelop);
    free(lower_envelop);
}

Query::~Query() {
    free(this->sorted_indexes);
    free(this->normalized_points);
    free(this->sorted_upper_envelop);
    free(this->sorted_lower_envelop);
}

#endif //UCR_SUITE_QUERY_H
