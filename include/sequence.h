//
// Created by Qitong Wang on 2019/10/23.
//

#ifndef UCR_SUITE_SEQUENCE_H
#define UCR_SUITE_SEQUENCE_H

#include <fstream>

class Sequence {
public:
    int length;
    double *values;

    Sequence(ifstream &sequence_ifs, int length);

    ~Sequence();
};

Sequence::Sequence(ifstream &sequence_ifs, int length) {
    this->length = length;
    this->values = (double *) malloc(sizeof(double) * length);

    if (this->values == nullptr) {
        error(1);
    }

    for (int i = 0; i < length; ++i) {
        sequence_ifs >> values[i];
    }
}

Sequence::~Sequence() {
    free(this->values);
}

#endif //UCR_SUITE_SEQUENCE_H
