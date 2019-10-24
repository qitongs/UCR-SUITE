//
// Created by Qitong Wang on 2019/10/23.
//

#ifndef UCR_SUITE_SEQUENCE_H
#define UCR_SUITE_SEQUENCE_H

class sequence {
public:
    int length;
    double *points;

    sequence(FILE *database_file, int length) {
        this->length = length;
        this->points = (double *) malloc(sizeof(double) * length);

        if (this->points == nullptr) {
            error(1);
        }

        double point;
        for (int i = 0; i < length; ++i) {
            if (fscanf(database_file, "%lf", &point) == EOF) {
                error(5);
            }
            points[i] = point;
        }

    }

    ~sequence() {
        free(this->points);
    }
};

#endif //UCR_SUITE_SEQUENCE_H
