//
// Created by Qitong Wang on 2019/10/22.
// Modified from UCR_SUITE
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iostream>
#include <queue>
#include "../include/bounds.h"
#include "../include/distances.h"
#include "../include/parameters.h"
#include "../include/query.h"
#include "../include/sequence.h"

using namespace std;

struct hit {
    long long location;
    double distance;
};

struct compare {
    bool operator()(const hit &a, const hit &b) {
        return a.distance > b.distance;
    }
};

int main(int argc, char *argv[]) {
    double *t;
    double *tz, *cb, *cb1, *cb2, *u_d, *l_d;

    double d;
    long long i, j;
    double ex, ex2, mean, std;
    double t1, t2;
    int kim = 0, keogh = 0, keogh2 = 0;
    double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
    double *buffer, *u_buff, *l_buff;

    /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum square), will be restarted for reducing the floating point error.
    int EPOCH = 100000;

    parameters parameters(argc, argv);
    FILE *database_fp = fopen(parameters.database_filename.c_str(), "r");
    FILE *queries_fp = fopen(parameters.queries_filename.c_str(), "r");
    if (database_fp == nullptr || queries_fp == nullptr) {
        error(2);
    }

    /// start the clock
    t1 = clock();

    cb = (double *) malloc(sizeof(double) * parameters.query_length);
    if (cb == nullptr)
        error(1);
    cb1 = (double *) malloc(sizeof(double) * parameters.query_length);
    if (cb1 == nullptr)
        error(1);
    cb2 = (double *) malloc(sizeof(double) * parameters.query_length);
    if (cb2 == nullptr)
        error(1);
    u_d = (double *) malloc(sizeof(double) * parameters.query_length);
    if (u_d == nullptr)
        error(1);
    l_d = (double *) malloc(sizeof(double) * parameters.query_length);
    if (l_d == nullptr)
        error(1);
    t = (double *) malloc(sizeof(double) * parameters.query_length * 2);
    if (t == nullptr)
        error(1);
    tz = (double *) malloc(sizeof(double) * parameters.query_length);
    if (tz == nullptr)
        error(1);
    buffer = (double *) malloc(sizeof(double) * EPOCH);
    if (buffer == nullptr)
        error(1);
    u_buff = (double *) malloc(sizeof(double) * EPOCH);
    if (u_buff == nullptr)
        error(1);
    l_buff = (double *) malloc(sizeof(double) * EPOCH);
    if (l_buff == nullptr)
        error(1);

    int sequence_length = 1000000;
    if (sequence_length < parameters.query_length) {
        cout << "queried sequence not long enough" << endl;
        exit(1);
    }
    query query(queries_fp, parameters.query_length, parameters.warping_window);
    sequence sequence(database_fp, sequence_length);
    // TODO reserved for multiple queries
    fclose(queries_fp);
    fclose(database_fp);

    /// Initial the cumulative lower bound
    for (i = 0; i < parameters.query_length; i++) {
        cb[i] = 0;
        cb1[i] = 0;
        cb2[i] = 0;
    }

    i = 0;          /// current position of the data in current chunk of size EPOCH
    j = 0;          /// the starting position of the data in the circular array, t
    ex = ex2 = 0;
    int it = 0, ep = 0, k = 0;
    long long I;    /// the starting position of the data in current chunk of size EPOCH
    long long last_location = -1;
    long long loc;

    // TODO overlap_ratio is set to prune contiguous matches. Replace with better methods
    float overlap_ratio = 0.5;

    double bsf = INF;
    priority_queue<hit, vector<hit>, compare> bsf_pq;

    int next = 0;
    while (true) {
        if (next + query.length > sequence.length) {
            break;
        }

        if (next + EPOCH > sequence.length) {
            copy(sequence.points + next, sequence.points + sequence.length, buffer);
            next = sequence.length;
            ep = sequence.length - next;
        } else {
            copy(sequence.points + next, sequence.points + next + EPOCH, buffer);
            next += (EPOCH - query.length + 1);
            ep = EPOCH;
        }

        lower_upper_lemire(buffer, ep, parameters.warping_window, l_buff, u_buff);

        /// Do main task here..
        ex = 0;
        ex2 = 0;
        for (i = 0; i < ep; i++) {
            /// A bunch of data has been read and pick one of them at a time to use
            d = buffer[i];

            /// Calcualte sum and sum square
            ex += d;
            ex2 += d * d;

            /// t is a circular array for keeping current data
            t[i % parameters.query_length] = d;

            /// Double the size for avoiding using modulo "%" operator
            t[(i % parameters.query_length) + parameters.query_length] = d;

            /// Start the task when there are more than parameters.query_length-1 points in the current chunk
            if (i >= parameters.query_length - 1) {
                mean = ex / parameters.query_length;
                std = ex2 / parameters.query_length;
                std = sqrt(std - mean * mean);

                /// compute the start location of the data in the current circular array, t
                j = (i + 1) % parameters.query_length;
                /// the start location of the data in the current chunk
                I = i - (parameters.query_length - 1);

                /// Use a constant lower bound to prune the obvious subsequence
                lb_kim = lb_kim_hierarchy(t, query.normalized_points, j, parameters.query_length, mean, std, bsf);

                if (lb_kim < bsf) {
                    /// Use a linear time lower bound to prune; z_normalization of t will be computed on the fly.
                    /// uo, lo are envelop of the query.
                    lb_k = lb_keogh_cumulative(query.sorted_indexes, t, query.sorted_upper_envelop,
                                               query.sorted_lower_envelop, cb1, j, parameters.query_length, mean,
                                               std, bsf);
                    if (lb_k < bsf) {
                        /// Take another linear time to compute z_normalization of t.
                        /// Note that for better optimization, this can merge to the previous function.
                        for (k = 0; k < parameters.query_length; k++) {
                            tz[k] = (t[(k + j)] - mean) / std;
                        }

                        /// Use another lb_keogh to prune
                        /// qo is the sorted query. tz is unsorted z_normalized data.
                        /// l_buff, u_buff are big envelop for all data in this chunk
                        lb_k2 = lb_keogh_data_cumulative(query.sorted_indexes, tz, query.sorted_normalized_points,
                                                         cb2, l_buff + I, u_buff + I, parameters.query_length, mean,
                                                         std, bsf);
                        if (lb_k2 < bsf) {
                            /// Choose better lower bound between lb_keogh and lb_keogh2 to be used in early abandoning DTW
                            /// Note that cb and cb2 will be cumulative summed here.
                            if (lb_k > lb_k2) {
                                cb[parameters.query_length - 1] = cb1[parameters.query_length - 1];
                                for (k = parameters.query_length - 2; k >= 0; k--)
                                    cb[k] = cb[k + 1] + cb1[k];
                            } else {
                                cb[parameters.query_length - 1] = cb2[parameters.query_length - 1];
                                for (k = parameters.query_length - 2; k >= 0; k--)
                                    cb[k] = cb[k + 1] + cb2[k];
                            }

                            /// Compute DTW and early abandoning if possible
                            dist = dtw(tz, query.normalized_points, cb, parameters.query_length,
                                       parameters.warping_window, bsf);

                            if (dist < bsf) {   /// Update bsf
                                /// loc is the real starting location of the nearest neighbor in the file
                                loc = (it) * (EPOCH - parameters.query_length + 1) + i - parameters.query_length +
                                      1;
                                if (loc > last_location + parameters.query_length * overlap_ratio) {
                                    hit current = {loc, dist};
                                    bsf_pq.push(current);
                                    if (bsf_pq.size() > parameters.num_neighbors) {
                                        bsf_pq.pop();
                                    }
                                    bsf = bsf_pq.top().distance;
                                    last_location = loc;
                                }
                            }
                        } else
                            keogh2++;
                    } else
                        keogh++;
                } else
                    kim++;

                /// Reduce obsolute points from sum and sum square
                ex -= t[j];
                ex2 -= t[j] * t[j];
            }
        }
    }

    free(cb);
    free(cb1);
    free(cb2);
    free(tz);
    free(t);
    free(l_d);
    free(u_d);
    free(l_buff);
    free(u_buff);

    t2 = clock();
    printf("\n");

    /// Note that loc and i are long long.
    for (int m = 0; m < parameters.num_neighbors; ++m) {
        hit current = bsf_pq.top();
        cout << "Location : " << current.location << endl;
        cout << "Distance : " << sqrt(current.distance) << endl;
        bsf_pq.pop();
    }

    cout << "Data Scanned : " << sequence.length << endl;
    cout << "Total Execution Time : " << (t2 - t1) / CLOCKS_PER_SEC << " sec" << endl;

    printf("\n");
    printf("Pruned by LB_Kim    : %6.2f%%\n", ((double) kim / sequence.length) * 100);
    printf("Pruned by LB_Keogh  : %6.2f%%\n", ((double) keogh / sequence.length) * 100);
    printf("Pruned by LB_Keogh2 : %6.2f%%\n", ((double) keogh2 / sequence.length) * 100);
    printf("DTW Calculation     : %6.2f%%\n", 100 - (((double) kim + keogh + keogh2) / sequence.length * 100));
    return 0;
}
