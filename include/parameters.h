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
// 174 0.00323177 0.218978 0.555281 0.141543
// 188 0.00373437 0.163182 0.606129 0.156316
// 420 0.00809375 0.00988108 0.0956615 0.560175
// 151 0.00371094 0.366174 0.477936 0.101372
// 142 0.00230903 0.293179 0.520549 0.113337
// 390 0.00832813 0.0105712 0.193819 0.425608
// 321 0.00730556 0.0257656 0.374333 0.393884
// 420 0.00853125 0.0279184 0.0571988 0.513456
// 459 0.00789063 0.0127917 0.251446 0.360835
// 470 0.00709896 0.0146007 0.0748707 0.510562
// 594 0.00993576 0.00533854 0.55368 0.14694
// 289 0.00619792 0.0924939 0.451851 0.260115
// 222 0.00439236 0.103855 0.546001 0.233701
// 230 0.00499132 0.10946 0.153451 0.514066
// 490 0.00880469 0.0141554 0.0937448 0.484572
// 446 0.00668403 0.043276 0.81939 0.0409488
// 373 0.00886458 0.088855 0.146733 0.461378
// 253 0.00423177 0.0323533 0.801076 0.096079
// 200 0.00465278 0.0544766 0.284578 0.47851
// 433 0.00589757 0.0153212 0.162497 0.425831
// 328 0.00619444 0.0812274 0.583511 0.17612
// 246 0.00470833 0.0643672 0.608461 0.224613
// 494 0.0114844 0.00332639 0.135052 0.46099
// 241 0.00475 0.0430408 0.533663 0.254803
// 286 0.0074375 0.0739731 0.153036 0.529484
// 459 0.0104688 0.0125052 0.425639 0.275077
// 549 0.00970312 0.00296528 0.0548168 0.519038
// 240 0.00552083 0.094099 0.482586 0.297999
// 355 0.00817274 0.0151293 0.20596 0.492408
// 180 0.0035 0.258672 0.579885 0.108595

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
