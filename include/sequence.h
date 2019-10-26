//
// Created by Qitong Wang on 2019/10/23.
//

#ifndef UCR_SUITE_SEQUENCE_H
#define UCR_SUITE_SEQUENCE_H

#include <fstream>

class Sequence {
public:
    int length;
    double *points;

    Sequence(ifstream &sequence_ifs, int length);

    ~Sequence();
};

Sequence::Sequence(ifstream &sequence_ifs, int length) {
    this->length = length;
    this->points = (double *) malloc(sizeof(double) * length);

    if (this->points == nullptr) {
        error(1);
    }

    for (int i = 0; i < length; ++i) {
        sequence_ifs >> points[i];
    }
}

Sequence::~Sequence() {
    free(this->points);
}

#endif //UCR_SUITE_SEQUENCE_H
