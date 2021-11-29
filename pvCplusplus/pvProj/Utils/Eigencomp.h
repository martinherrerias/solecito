#ifndef EIGENCOMP_H_
#define EIGENCOMP_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "pvcommutils.h"
#include <vector>
#include <queue>
using namespace Eigen;
//typedef std::pair<double, std::pair<int, int>> ppi;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXRd;
typedef Array<bool, Dynamic, 1> ArrayXb;

//TODO: add a specific namespace
void mask(const Ref<MatrixXd> input, Ref<MatrixXd> output, Ref<ArrayXb> maskBool)
{}





#endif
