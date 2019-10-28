//
// Created by Qitong Wang on 2019/10/22.
// Modified from UCR_SUITE
//

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <queue>
#include <vector>

#include <boost/thread.hpp>

#include "bounds.h"
#include "distances.h"
#include "parameters.h"
#include "query.h"
#include "sequence.h"
#include "thread_pool.h"

using namespace std;

struct Hit {
    Hit(int location, double distance) : location(location), distance(distance) {}

    int location;
    double distance;
};

struct compare {
    bool operator()(const Hit *prior, const Hit *later) {
        return prior->distance < later->distance;
    }
};

void conduct_query(const Sequence *sequence, const Query *query, const Parameters *parameters, int sequence_id,
                   int query_id, boost::mutex *results_mutex, ofstream *results_ofs, clock_t start_time,
                   atomic<int> *processed) {
    int i, num_to_skip = 0, next = 0, epoch_size = parameters->epoch, num_epoches = 0, to_calculate, start, position_in_epoch, location;
    double sum, squared_sum, mean, std, bsf = INF, distance = 0, lb_kim = 0, lb_keogh_1 = 0, lb_keogh_2 = 0;
    priority_queue<Hit *, vector<Hit *>, compare> bsf_pq;
    Hit *hit;
    double *tighter_bounds_keogh;

    auto buffer = (double *) malloc(sizeof(double) * parameters->epoch);
    auto upper_envelop = (double *) malloc(sizeof(double) * parameters->epoch);
    auto lower_envelop = (double *) malloc(sizeof(double) * parameters->epoch);
    auto current = (double *) malloc(sizeof(double) * query->length * 2);
    auto current_normalized = (double *) malloc(sizeof(double) * query->length);
    auto bounds_keogh_left = (double *) malloc(sizeof(double) * query->length);
    auto bounds_keogh_1 = (double *) malloc(sizeof(double) * query->length);
    auto bounds_keogh_2 = (double *) malloc(sizeof(double) * query->length);

    if (buffer == nullptr || upper_envelop == nullptr || lower_envelop == nullptr || current == nullptr ||
        current_normalized == nullptr || bounds_keogh_left == nullptr || bounds_keogh_1 == nullptr ||
        bounds_keogh_2 == nullptr) {
        error(1);
    }

    while (next + query->length - 1 < sequence->length) {
        if (next + parameters->epoch > sequence->length) {
            copy(sequence->points + next, sequence->points + sequence->length, buffer);
            epoch_size = sequence->length - next;
            next = sequence->length;
        } else {
            copy(sequence->points + next, sequence->points + next + parameters->epoch, buffer);
            next += (parameters->epoch - query->length + 1);
        }

        sum = 0, squared_sum = 0;
        // TODO here buffer is not z-normalized, but in query it's z-normalized
        lower_upper_lemire(buffer, epoch_size, query->warping_window, lower_envelop, upper_envelop);

        for (to_calculate = 0; to_calculate < epoch_size; ++to_calculate) {
            current[to_calculate % query->length] = buffer[to_calculate];
            current[(to_calculate % query->length) + query->length] = buffer[to_calculate];

            sum += buffer[to_calculate];
            squared_sum += (buffer[to_calculate] * buffer[to_calculate]);

            if (to_calculate < query->length - 1) {
                continue;
            }

            start = (to_calculate + 1) % query->length;

            // TODO skip the next following overlap sub-sequences to check is an easy but not optimal implementation
            if (num_to_skip-- > 0) {
                goto UPDATE_STATISTICS;
            }

            mean = sum / query->length;
            std = sqrt(squared_sum / query->length - mean * mean);
            lb_kim = lb_kim_hierarchy(current, query->normalized_points, start, query->length, mean, std, bsf);

            if (lb_kim >= bsf) {
                goto UPDATE_STATISTICS;
            }

            position_in_epoch = to_calculate + 1 - query->length;
            lb_keogh_1 = lb_keogh_cumulative(query->sorted_indexes, current, query->sorted_upper_envelop,
                                             query->sorted_lower_envelop, bounds_keogh_1, start, query->length, mean,
                                             std, bsf);
            if (lb_keogh_1 >= bsf) {
                goto UPDATE_STATISTICS;
            }

            // TODO Note that for better optimization, this can be merged into the previous lb_keogh_cumulative
            for (i = 0; i < query->length; ++i) {
                current_normalized[i] = (current[(start + i)] - mean) / std;
            }

            lb_keogh_2 = lb_keogh_data_cumulative(query->sorted_indexes, current_normalized,
                                                  query->sorted_normalized_points, bounds_keogh_2,
                                                  lower_envelop + position_in_epoch, upper_envelop + position_in_epoch,
                                                  query->length, mean, std, bsf);

            if (lb_keogh_2 >= bsf) {
                goto UPDATE_STATISTICS;
            }

            tighter_bounds_keogh = lb_keogh_1 > lb_keogh_2 ? bounds_keogh_1 : bounds_keogh_2;
            copy(tighter_bounds_keogh, tighter_bounds_keogh + query->length, bounds_keogh_left);
            for (i = query->length - 2; i >= 0; --i) {
                bounds_keogh_left[i] += bounds_keogh_left[i + 1];
            }

            distance = dtw(current_normalized, query->normalized_points, bounds_keogh_left, query->length,
                           query->warping_window, bsf);

            if (distance < bsf) {
                location = num_epoches * (parameters->epoch - query->length + 1) + to_calculate + 1 - query->length;
                if (bsf_pq.size() == parameters->num_neighbors) {
                    hit = bsf_pq.top();
                    bsf_pq.pop();
                    hit->distance = distance;
                    hit->location = location;
                    bsf_pq.push(hit);
                } else {
                    bsf_pq.push(new Hit(location, distance));
                }
                bsf = bsf_pq.top()->distance;
                num_to_skip = query->overlap_length;
            }

            UPDATE_STATISTICS:
            sum -= current[start];
            squared_sum -= (current[start] * current[start]);
        }
        num_epoches += 1;
    }

    free(buffer);
    free(upper_envelop);
    free(lower_envelop);
    free(current);
    free(current_normalized);
    free(bounds_keogh_left);
    free(bounds_keogh_1);
    free(bounds_keogh_2);

    {
        // TODO asynchronously write to file (but not critical)
        boost::mutex::scoped_lock lock(*results_mutex);
        for (i = 0; i < parameters->num_neighbors; ++i) {
            *results_ofs << sequence_id << " " << query_id << " " << sqrt(bsf_pq.top()->distance) << " "
                         << bsf_pq.top()->location << " " << query->length << endl;
            bsf_pq.pop();
        }
    }

    while (!bsf_pq.empty()) {
        free(bsf_pq.top());
        bsf_pq.pop();
    }

    cout << ++*processed << " in " << (clock() - start_time) / CLOCKS_PER_SEC << "s" << endl;
}

int main(int argc, char *argv[]) {
    boost::mutex results_mutex;
    Parameters parameters(argc, argv);

    remove(parameters.results_path.c_str());
    ifstream database_ifs(parameters.database_path);
    ifstream queries_ifs(parameters.queries_path);
    ofstream results_ofs(parameters.results_path, std::ofstream::out | std::ofstream::app);

    if (!database_ifs || !queries_ifs || !results_ofs) {
        error(2);
    }

    int size, length, id;
    clock_t start_time = clock();

    database_ifs >> size;
    vector<const Sequence *> sequences;
    vector<int> sequence_ids;
    for (int i = 0; i < size; ++i) {
        database_ifs >> id >> length;
        if (length > 0) {
            sequence_ids.push_back(id);
            sequences.push_back(new Sequence(database_ifs, length));
        } else {
            cout << "Empty sequence " << id << endl;
        }
    }
    database_ifs.close();

    queries_ifs >> size;
    vector<const Query *> queries;
    vector<int> query_ids;
    for (int i = 0; i < size; ++i) {
        queries_ifs >> id >> length;
        if (length > 0) {
            query_ids.push_back(id);
            queries.push_back(new Query(queries_ifs, length, (int) (parameters.warping_window * (double) length),
                                        (int) (parameters.overlap_ratio * (double) length)));
        } else {
            cout << "Empty query " << id << endl;
        }
    }
    queries_ifs.close();

    cout << "load in " << (clock() - start_time) / CLOCKS_PER_SEC << "s" << endl;

    {
        ThreadPool threadPool;
        atomic<int> processed{0};
        start_time = clock();

        for (int i = 0; i < (int) sequences.size(); ++i) {
            for (int j = 0; j < (int) queries.size(); ++j) {
                if ((long) queries[j]->length * parameters.min_query_length_multiple < (long) sequences[i]->length) {
                    cout << "Too exhaustive for " << query_ids[j] << " on " << sequence_ids[i] << ": "
                         << queries[j]->length << " / " << sequences[i]->length << endl;
                    continue;
                }
                threadPool.enqueue(boost::bind(conduct_query, sequences[i], queries[j], &parameters, sequence_ids[i],
                                               query_ids[j], &results_mutex, &results_ofs, start_time, &processed));
            }
        }
    }

    results_ofs.close();
    return 0;
}
