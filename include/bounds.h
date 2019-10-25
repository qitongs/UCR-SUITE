//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_BOUNDS_H
#define UCR_SUITE_BOUNDS_H

#include "circular_array.h"
#include "utils.h"

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is introduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern Recognition 42(9), 2009.
void lower_upper_lemire(double *t, int len, int r, double *l, double *u) {
    struct CircularArray du, dl;

    init(&du, 2 * r + 2);
    init(&dl, 2 * r + 2);

    push_back(&du, 0);
    push_back(&dl, 0);

    for (int i = 1; i < len; i++) {
        if (i > r) {
            u[i - r - 1] = t[front(&du)];
            l[i - r - 1] = t[front(&dl)];
        }
        if (t[i] > t[i - 1]) {
            pop_back(&du);
            while (!empty(&du) && t[i] > t[back(&du)])
                pop_back(&du);
        } else {
            pop_back(&dl);
            while (!empty(&dl) && t[i] < t[back(&dl)])
                pop_back(&dl);
        }
        push_back(&du, i);
        push_back(&dl, i);
        if (i == 2 * r + 1 + front(&du))
            pop_front(&du);
        else if (i == 2 * r + 1 + front(&dl))
            pop_front(&dl);
    }
    for (int i = len; i < len + r + 1; i++) {
        u[i - r - 1] = t[front(&du)];
        l[i - r - 1] = t[front(&dl)];
        if (i - front(&du) >= 2 * r + 1)
            pop_front(&du);
        if (i - front(&dl) >= 2 * r + 1)
            pop_front(&dl);
    }
    destroy(&du);
    destroy(&dl);
}

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give significant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is not long, say in length 128.
double lb_kim_hierarchy(double *t, double *q, int j, int len, double mean, double std, double bsf = INF) {
    /// 1 point at front and back
    double d, lb;
    double x0 = (t[j] - mean) / std;
    double y0 = (t[(len - 1 + j)] - mean) / std;
    lb = dist(x0, q[0]) + dist(y0, q[len - 1]);
    if (lb >= bsf) return lb;

    /// 2 points at front
    double x1 = (t[(j + 1)] - mean) / std;
    d = min(dist(x1, q[0]), dist(x0, q[1]));
    d = min(d, dist(x1, q[1]));
    lb += d;
    if (lb >= bsf) return lb;

    /// 2 points at back
    double y1 = (t[(len - 2 + j)] - mean) / std;
    d = min(dist(y1, q[len - 1]), dist(y0, q[len - 2]));
    d = min(d, dist(y1, q[len - 2]));
    lb += d;
    if (lb >= bsf) return lb;

    /// 3 points at front
    double x2 = (t[(j + 2)] - mean) / std;
    d = min(dist(x0, q[2]), dist(x1, q[2]));
    d = min(d, dist(x2, q[2]));
    d = min(d, dist(x2, q[1]));
    d = min(d, dist(x2, q[0]));
    lb += d;
    if (lb >= bsf) return lb;

    /// 3 points at back
    double y2 = (t[(len - 3 + j)] - mean) / std;
    d = min(dist(y0, q[len - 3]), dist(y1, q[len - 3]));
    d = min(d, dist(y2, q[len - 3]));
    d = min(d, dist(y2, q[len - 2]));
    d = min(d, dist(y2, q[len - 1]));
    lb += d;

    return lb;
}

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for early abandoning in DTW.
double lb_keogh_cumulative(int *order, double *t, double *uo, double *lo, double *cb, int j, int len, double mean,
                           double std, double best_so_far = INF) {
    double lb = 0;
    double x, d;

    for (int i = 0; i < len && lb < best_so_far; i++) {
        x = (t[(order[i] + j)] - mean) / std;
        d = 0;
        if (x > uo[i])
            d = dist(x, uo[i]);
        else if (x < lo[i])
            d = dist(x, lo[i]);
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int *order, double *tz, double *qo, double *cb, double *l, double *u, int len,
                                double mean, double std, double best_so_far = INF) {
    double lb = 0;
    double uu, ll, d;

    for (int i = 0; i < len && lb < best_so_far; i++) {
        uu = (u[order[i]] - mean) / std;
        ll = (l[order[i]] - mean) / std;
        d = 0;
        if (qo[i] > uu)
            d = dist(qo[i], uu);
        else {
            if (qo[i] < ll)
                d = dist(qo[i], ll);
        }
        lb += d;
        cb[order[i]] = d;
    }
    return lb;
}

#endif //UCR_SUITE_BOUNDS_H
