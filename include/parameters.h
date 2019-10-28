//
// Created by Qitong Wang on 2019/10/24.
//

#ifndef UCR_SUITE_PARAMETERS_H
#define UCR_SUITE_PARAMETERS_H

#include <string>

#include <boost/program_options.hpp>

#include "utils.h"

class Parameters {
public:
    string database_path, queries_path, results_path;
    // Cumulative values (sum, squared_sum, etc) will be restarted to reduce floating point errors for every EPOCH points
    int num_neighbors, epoch = 100000, min_length_for_kim = 128;
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
