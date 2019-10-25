//
// Created by Qitong Wang on 2019/10/25.
//

#ifndef UCR_SUITE_THREAD_POOL_H
#define UCR_SUITE_THREAD_POOL_H

#include <iostream>
#include <thread>
#include <future>
#include <utility>
#include <chrono>

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/asio.hpp>

#include "utils.h"

class ThreadPool {
public:
    explicit ThreadPool() : work_(io_service_) {
        int cores = (int) std::thread::hardware_concurrency();
        if (cores == 0) {
            error(6);
        }
        cout << "Num of cores: " << cores << endl;
        this->create_threads(cores);
    }

    explicit ThreadPool(size_t size) : work_(io_service_) {
        this->create_threads(size);
    }

    ~ThreadPool() {
        io_service_.stop();
        workers_.join_all();
    }

    template<class F>
    void enqueue(F f) {
        io_service_.post(f);
    }

private:
    boost::thread_group workers_;
    boost::asio::io_service io_service_;
    boost::asio::io_service::work work_;

    void create_threads(int num_threads) {
        for (size_t i = 0; i < num_threads; ++i) {
            this->workers_.create_thread(boost::bind(&boost::asio::io_service::run, &this->io_service_));
        }
    }
};

#endif //UCR_SUITE_THREAD_POOL_H
