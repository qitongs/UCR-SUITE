//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_DISTANCES_H
#define UCR_SUITE_DISTANCES_H

double dtw(const double *sequence, const double *query, const double *bounds_remaining, int length, int warping_window,
           double bsf) {

    double *distances, *distances_previous, *tmp;
    int i, sequence_iter, query_iter, distances_iter;
    double underneath, left, lower_left, partial_distance;

    distances = (double *) malloc(sizeof(double) * (2 * warping_window + 1));
    for (i = 0; i < 2 * warping_window + 1; i++) {
        distances[i] = INF;
    }

    distances_previous = (double *) malloc(sizeof(double) * (2 * warping_window + 1));
    for (i = 0; i < 2 * warping_window + 1; i++) {
        distances_previous[i] = INF;
    }

    for (sequence_iter = 0; sequence_iter < length; sequence_iter++) {
        distances_iter = max(0, warping_window - sequence_iter);
        partial_distance = INF;

        for (query_iter = max(0, sequence_iter - warping_window);
             query_iter <= min(length - 1, sequence_iter + warping_window);
             query_iter++, distances_iter++) {

            if ((sequence_iter == 0) && (query_iter == 0)) {
                distances[distances_iter] = dist(sequence[0], query[0]);
                partial_distance = distances[distances_iter];
                continue;
            }

            // query is horizontal
            if ((query_iter == 0) || (distances_iter == 0)) {
                left = INF;
            } else {
                left = distances[distances_iter - 1];
            }

            // sequence is vertical
            if ((sequence_iter == 0) || (distances_iter + 1 > 2 * warping_window)) {
                underneath = INF;
            } else {
                underneath = distances_previous[distances_iter + 1];
            }

            if ((sequence_iter == 0) || (query_iter == 0)) {
                lower_left = INF;
            } else {
                lower_left = distances_previous[distances_iter];
            }

            distances[distances_iter] = min(min(underneath, left), lower_left)
                                        + dist(sequence[sequence_iter], query[query_iter]);

            // Find minimum distances in row for early abandoning (possibly to use column instead of row)
            // TODO any hint for a tighter bound?
            if (distances[distances_iter] < partial_distance) {
                partial_distance = distances[distances_iter];
            }
        }

        if (partial_distance + bounds_remaining[sequence_iter + warping_window + 1] >= bsf &&
            sequence_iter + warping_window < length - 1) {
            free(distances);
            free(distances_previous);
            return partial_distance + bounds_remaining[sequence_iter + warping_window + 1];
        }

        tmp = distances;
        distances = distances_previous;
        distances_previous = tmp;
    }

    distances_iter -= 1;
    double distance = distances_previous[distances_iter];

    free(distances);
    free(distances_previous);

    return distance;
}

#endif //UCR_SUITE_DISTANCES_H
