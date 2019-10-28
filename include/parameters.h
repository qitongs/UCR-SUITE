//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_PARAMETERS_H
#define UCR_SUITE_PARAMETERS_H

#include <string>

#include <boost/program_options.hpp>

#include "utils.h"

// min_length_for_kim should prune more than (1 / query_length) to be worthwhile in complexity (but constants matter)
// roughly all applied as shown in the following examples, thus not deployed
// [length, skipped, by_kim, by_keogh, by_keogh_converse]
// 188 0.000875 0.223911 0.471361 0.281039
// 174 0.00150521 0.325654 0.514832 0.102458
// 459 0.00453125 0.153158 0.758717 0.0257066
// 151 0.00132813 0.312524 0.41979 0.248977
// 142 0.000850694 0.468589 0.463522 0.0386623
// 420 0.00601563 0.13398 0.329016 0.409307
// 390 0.00501042 0.253196 0.475569 0.162555
// 321 0.00258333 0.145105 0.438671 0.384721
// 594 0.00450694 0.103576 0.554319 0.255978
// 420 0.00441146 0.183122 0.368269 0.263959
// 289 0.00332986 0.281138 0.425273 0.198531
// 230 0.00261545 0.334196 0.417409 0.152005
// 222 0.00179514 0.304003 0.493819 0.15058
// 494 0.00510417 0.128636 0.676745 0.134022
// 470 0.00563021 0.149609 0.408197 0.250408
// 433 0.00515104 0.161164 0.354253 0.296301
// 446 0.00294097 0.113105 0.730171 0.0847083
// 253 0.00134549 0.194232 0.760005 0.0282569
// 373 0.0049783 0.161444 0.307638 0.397655
// 200 0.00237847 0.303751 0.375157 0.273594
// 328 0.00308333 0.169898 0.409855 0.338328
// 246 0.00191667 0.266346 0.498723 0.198414
// 286 0.00274653 0.258409 0.374628 0.317203
// 241 0.00216667 0.180495 0.500257 0.27237
// 459 0.00367187 0.144428 0.673668 0.0700313
// 549 0.00825 0.134356 0.50417 0.168525
// 490 0.00727344 0.160368 0.386308 0.182649
// 355 0.00388889 0.202942 0.415668 0.301661
// 240 0.00245833 0.23039 0.390745 0.320414
// 180 0.001375 0.256753 0.40511 0.304148

class Parameters {
public:
    string database_path, queries_path, results_path;
    // Cumulative values (sum, squared_sum, etc) will be restarted to reduce floating point errors for every EPOCH points
    int num_neighbors, epoch = 100000, min_length_for_kim = 500;
    long min_query_length_multiple = 100000;
    double overlap_ratio, warping_window = 1;

    Parameters(int argc, char *argv[]);
};

Parameters::Parameters(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description desc("Allowed options");
    desc.add_options()
            ("help", "Help message")
            ("db", po::value<string>(&this->database_path)->required(), "Sequences database path")
            ("q", po::value<string>(&this->queries_path)->required(), "Queries path")
            ("r", po::value<string>(&this->results_path)->required(), "Results path")
            ("w", po::value<double>(), "DTW warping ratio (default: 1)")
            ("k", po::value<int>(), "k in kNN (default: 1)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        error(4);
    }

    if (vm.count("k")) {
        this->num_neighbors = vm["k"].as<int>();
    } else {
        this->num_neighbors = 1;
    }

    if (vm.count("w")) {
        this->warping_window = vm["w"].as<double>();
    }

    // TODO currently overlap_ratio is set to warping_window to prevent points matching to themselves
    // But the behaviour is un-predicted as a result of z-normalization
    this->overlap_ratio = this->warping_window;
}

#endif //UCR_SUITE_PARAMETERS_H
