//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_BOUNDS_H
#define UCR_SUITE_BOUNDS_H

#include "circular_array.h"
#include "utils.h"

// TODO check the potential length = 0, warping_window = 0 problem
/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is introduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void get_envelops_lemire(const double *sequence, const int length, const int warping_window, double *lower_envelop,
                        double *upper_envelop) {
    CircularArray upper, lower;

    init(&upper, 2 * warping_window + 2);
    init(&lower, 2 * warping_window + 2);

    push_back(&upper, 0);
    push_back(&lower, 0);

    for (int i = 1; i < length; i++) {
        if (i > warping_window) {
            upper_envelop[i - warping_window - 1] = sequence[front(&upper)];
            lower_envelop[i - warping_window - 1] = sequence[front(&lower)];
        }

        if (sequence[i] > sequence[i - 1]) {
            pop_back(&upper);

            while (!empty(&upper) && sequence[i] > sequence[back(&upper)]) {
                pop_back(&upper);
            }
        } else {
            pop_back(&lower);
            
            while (!empty(&lower) && sequence[i] < sequence[back(&lower)]) {
                pop_back(&lower);
            }
        }

        push_back(&upper, i);
        push_back(&lower, i);

        if (i == 2 * warping_window + 1 + front(&upper)) {
            pop_front(&upper);
        } else if (i == 2 * warping_window + 1 + front(&lower)) {
            pop_front(&lower);
        }
    }

    for (int i = length; i < length + warping_window + 1; i++) {
        upper_envelop[i - warping_window - 1] = sequence[front(&upper)];
        lower_envelop[i - warping_window - 1] = sequence[front(&lower)];

        if (i - front(&upper) >= 2 * warping_window + 1) {
            pop_front(&upper);
        }
        if (i - front(&lower) >= 2 * warping_window + 1) {
            pop_front(&lower);
        }
    }

    destroy(&upper);
    destroy(&lower);
}

/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give significant benefits.
/// And using the first and last points can be computed in constant time.
/// The pruning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double get_kim(const double *subsequence, const double *query_normalized, const int start, const int query_length,
               const double mean, const double std, const double bsf) {
    double distance, partial_bound;
    double head_0 = (subsequence[start] - mean) / std;
    double tail_0 = (subsequence[start + query_length - 1] - mean) / std;

    partial_bound = dist(head_0, query_normalized[0]) + dist(tail_0, query_normalized[query_length - 1]);
    if (partial_bound >= bsf) {
        return partial_bound;
    }

    double head_1 = (subsequence[(start + 1)] - mean) / std;
    distance = min(dist(head_1, query_normalized[0]), dist(head_0, query_normalized[1]));
    distance = min(distance, dist(head_1, query_normalized[1]));
    partial_bound += distance;
    if (partial_bound >= bsf) {
        return partial_bound;
    }

    double tail_1 = (subsequence[(start + query_length - 2)] - mean) / std;
    distance = min(dist(tail_1, query_normalized[query_length - 1]), dist(tail_0, query_normalized[query_length - 2]));
    distance = min(distance, dist(tail_1, query_normalized[query_length - 2]));
    partial_bound += distance;
    if (partial_bound >= bsf) {
        return partial_bound;
    }

    double head_2 = (subsequence[(start + 2)] - mean) / std;
    distance = min(dist(head_0, query_normalized[2]), dist(head_1, query_normalized[2]));
    distance = min(distance, dist(head_2, query_normalized[2]));
    distance = min(distance, dist(head_2, query_normalized[1]));
    distance = min(distance, dist(head_2, query_normalized[0]));
    partial_bound += distance;
    if (partial_bound >= bsf) {
        return partial_bound;
    }

    double tail_2 = (subsequence[(start + query_length - 3)] - mean) / std;
    distance = min(dist(tail_0, query_normalized[query_length - 3]), dist(tail_1, query_normalized[query_length - 3]));
    distance = min(distance, dist(tail_2, query_normalized[query_length - 3]));
    distance = min(distance, dist(tail_2, query_normalized[query_length - 2]));
    distance = min(distance, dist(tail_2, query_normalized[query_length - 1]));
    partial_bound += distance;

    return partial_bound;
}

double get_keogh(const int *query_sorted_indices, const double *subsequence, double *subsequence_normalized,
                 const double *query_upper_envelop, const double *query_lower_envelop, double *local_bounds,
                 const int start, const int query_length, const double mean, const double std, const double bsf) {
    double partial_bound = 0, local_distance;

    for (int i = 0; i < query_length && partial_bound < bsf; i++) {
        subsequence_normalized[query_sorted_indices[i]] = (subsequence[start + query_sorted_indices[i]] - mean) / std;

        local_distance = 0;
        if (subsequence_normalized[query_sorted_indices[i]] > query_upper_envelop[i]) {
            local_distance = dist(subsequence_normalized[query_sorted_indices[i]], query_upper_envelop[i]);
        } else if (subsequence_normalized[query_sorted_indices[i]] < query_lower_envelop[i]) {
            local_distance = dist(subsequence_normalized[query_sorted_indices[i]], query_lower_envelop[i]);
        }

        partial_bound += local_distance;
        local_bounds[query_sorted_indices[i]] = local_distance;
    }

    return partial_bound;
}

double get_keogh_converse(const int *query_sorted_indices, const double *query_sorted_normalized_values,
                          double *local_bounds, const double *subsequence_lower_envelop,
                          const double *subsequence_upper_envelop, const int query_length, const double mean,
                          const double std, const double bsf) {
    double partial_bound = 0, normalized_upper_envelop, normalized_lower_envelop, local_distance;

    for (int i = 0; i < query_length && partial_bound < bsf; i++) {
        normalized_upper_envelop = (subsequence_upper_envelop[query_sorted_indices[i]] - mean) / std;

        local_distance = 0;
        if (query_sorted_normalized_values[i] > normalized_upper_envelop) {
            local_distance = dist(query_sorted_normalized_values[i], normalized_upper_envelop);
        } else {
            normalized_lower_envelop = (subsequence_lower_envelop[query_sorted_indices[i]] - mean) / std;
            if (query_sorted_normalized_values[i] < normalized_lower_envelop) {
                local_distance = dist(query_sorted_normalized_values[i], normalized_lower_envelop);
            }
        }

        partial_bound += local_distance;
        local_bounds[query_sorted_indices[i]] = local_distance;
    }

    return partial_bound;
}

#endif //UCR_SUITE_BOUNDS_H
