/**
 * Lambert-W function.
 * 
 * Special handling in the case of infinite.
 */
#ifndef LAMBERTW_H_
#define LAMBERTW_H_

#include "pvcommutils.h"
#include <boost/math/special_functions/lambert_w.hpp>
#include <omp.h>

constexpr unsigned maxIter = 1000;
constexpr int NEPS = 3;
double lambertW(double x);

#endif