#ifndef _BASE_H
#define _BASE_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <utility>
#include <ctime>

#define THREADS_NUM 48

namespace py = pybind11;
using std::cout;

inline double sqr(double x) {
    return x * x;
}

template <typename T>
inline T min(const T&a, const T& b) {
    return a < b ? a : b;
}

template <typename T>
inline T max(const T&a, const T& b) {
    return a > b ? a : b;
}

double getDist(double x1, double y1, double x2, double y2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

#endif