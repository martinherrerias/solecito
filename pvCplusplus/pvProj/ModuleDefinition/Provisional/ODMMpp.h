#ifndef ODMMPP_H
#define ODMMPP_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include "../Utils/xmlutils.h"
#include "../Utils/pvcommutils.h"
#include "../Utils/Rootfind.h"
#include "OneDiodeModel.h"

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/tools/roots.hpp>


using namespace Eigen;
class GetSimOptFile;
class OneDiodeModel;

typedef Eigen::Array<double, Dynamic, 1> ArrayXd;
typedef Eigen::Array<bool, Dynamic, 1> ArrayXb;


class ODMMpp
{
private:
	GetSimOptFile* filePtr;
	struct Params inputParams;
	double neps;
	double TolRefFac = { 16.0 }; // tolerance refinement factor
 
	std::vector<double> tol;

	double seedVmp;

	unsigned maxIter;
	std::shared_ptr<OneDiodeModel> objODMHandle;
//	OneDiodeModel objODMHandle;

public:
	double Vmp;
	double Imp;
	double Pmp;
	double Voc;
	double Isc;
	double c;

	ODMMpp(struct Params params_, const std::vector<double>& tol_ = {0.0}, 
		const double& seedVmp_ = {0.0}, const unsigned long& maxIter_ = {(unsigned long)INITVAL});
 
	std::shared_ptr<OneDiodeModel> getODM() {return objODMHandle;}
	void getODMMpp();
	void dispRes() {std::cout<<"Pmp is "<<Pmp<<" Vmp is "<<Vmp<<" Imp is "<<Imp<<" Voc is "<<Voc<<" Isc is "<<Isc<<" c is "<<c<<"\n";}
};

#endif
