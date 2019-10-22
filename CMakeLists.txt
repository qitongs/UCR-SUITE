cmake_minimum_required (VERSION 2.6)

project (UCR-SUITE)
set (UCR-SUITE_VERSION_MAJOR 1)
set (UCR-SUITE_VERSION_MINOR 0)

#Bring the headers, such as Student.h into the project
include_directories(include)

#Can manually add the sources using the set command as follows:
set(SOURCES src/UCR_kNNDTW.cpp)

#However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")

add_executable(kNNDTW ${SOURCES})