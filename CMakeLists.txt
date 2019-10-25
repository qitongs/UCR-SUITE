cmake_minimum_required (VERSION 2.6)

project (UCR-SUITE CXX)
set (UCR-SUITE_VERSION_MAJOR 1)
set (UCR-SUITE_VERSION_MINOR 1)

set(CMAKE_CXX_STANDARD 11)

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_DEBUG_LIBS OFF)
set(Boost_USE_RELEASE_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.71.0 COMPONENTS program_options thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories(include)
set(SOURCES src/UCRDTW_kNN.cpp)
set(CXX_FLAGS -g -Wall -pthread)

add_executable(kNNDTW ${SOURCES})
target_link_libraries(kNNDTW ${Boost_LIBRARIES})