#ifndef _COMMON_H
#define _COMMON_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

const double GAMMA = 1.4;
typedef unsigned int uint;


#include "adolc/adolc.h"
#include "adolc/sparse/sparsedrivers.h"



#define ENABLE_EIGEN
//#define ENABLE_ARMA

#if defined(ENABLE_EIGEN) && defined(ENABLE_ARMA)
#error "Cannot use both the linear algebra library at the same time."
#endif

#if defined(ENABLE_ARMA)
#include <armadillo>
#elif defined(ENABLE_EIGEN)
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SuperLUSupport>
#endif


#include <chrono>
class Timer {
public:
    Timer() {
        reset();
    }
    void reset() {
        m_timestamp = std::chrono::high_resolution_clock::now();
    }
    float diff() {
        std::chrono::duration<float> fs = std::chrono::high_resolution_clock::now() - m_timestamp;
        return fs.count();
    }
private:
    std::chrono::high_resolution_clock::time_point m_timestamp;
};
#endif
