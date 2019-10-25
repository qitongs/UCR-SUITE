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
    string database_path;
    string queries_path;
    string results_path;
    int query_length;
    int warping_window;
    int num_neighbors;
    /// For every EPOCH points, all cumulative values, such sum, squared_sum, will be restarted for reducing the floating point error.
    int epoch = 100000;
    // TODO overlap_ratio is set to prune contiguous matches. Replace with better methods
    float overlap_ratio = 0.5;

    Parameters(int argc, char *argv[]) {
        namespace bpo = boost::program_options;

        bpo::options_description desc("Allowed options");
        desc.add_options()
                ("help", "Help message")
                ("db", bpo::value<string>(&this->database_path)->required(), "Sequences database path")
                ("q", bpo::value<string>(&this->queries_path)->required(), "Queries path")
                ("r", bpo::value<string>(&this->results_path)->required(), "Results path")
                ("w", bpo::value<double>(), "DTW warping ratio")
                ("k", bpo::value<int>(), "k in kNN");

        bpo::variables_map vm;
        bpo::store(bpo::parse_command_line(argc, argv, desc), vm);
        bpo::notify(vm);

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
            this->warping_window = floor(vm["w"].as<double>() * this->query_length);
        } else {
            this->warping_window = this->query_length;
        }
    }
};

#endif //UCR_SUITE_PARAMETERS_H
