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

typedef std::unique_ptr<boost::asio::io_service::work> work_ptr;
typedef boost::asio::io_service io_service;

class ThreadPool {
public:
    explicit ThreadPool();

    explicit ThreadPool(size_t size);

    ~ThreadPool() {
        work_.reset();
        workers_.join_all();
        // io_service_.stop() interrupt all threads
        // might be commented out as threads already stops after workers_.join_all()
        // TODO figure out why this must be commented out
//        io_service_.stop();
    }

    template<class F>
    void enqueue(F f) {
        io_service_.post(f);
    }

private:
    boost::thread_group workers_;
    io_service io_service_;
    work_ptr work_;

    void create_threads(int num_threads);
};

ThreadPool::ThreadPool() : io_service_(), work_(new work_ptr::element_type(io_service_)) {
    int cores = (int) std::thread::hardware_concurrency();
    if (cores == 0) {
        error(6);
    }

    cout << "Num of cores: " << cores << endl;
    this->create_threads(cores);
}

ThreadPool::ThreadPool(size_t size) : io_service_(), work_(new work_ptr::element_type(io_service_)) {
    this->create_threads(size);
}

void ThreadPool::create_threads(int num_threads) {
    for (size_t i = 0; i < num_threads; ++i) {
        this->workers_.create_thread(boost::bind(&io_service::run, &this->io_service_));
    }
}

#endif //UCR_SUITE_THREAD_POOL_H
