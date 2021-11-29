/**
 * Modelling of a single PV cell.
 * 
 * Solve the basic equation (see Matlab) for getting (approximating) the I/V curve according to the given input parameters.
 * 
 * Input:  inputParams
 * Output: approxPointSet (as an input for the class MonoDescPWL)
 */
#ifndef ONEDIODEMODEL_H
#define ONEDIODEMODEL_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include "../Utils/xmlutils.h"
#include "../Utils/pvcommutils.h"
#include "../Utils/Rootfind.h"
#include "../Utils/Linearapprox.h"

#include <include/igl/slice.h>
#include <include/igl/slice_into.h>
#include <include/igl/find.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include <boost/math/tools/roots.hpp>


using namespace Eigen;
class GetSimOptFile;

typedef Eigen::Array<double, Dynamic, 1> ArrayXd;
typedef Eigen::Array<bool, Dynamic, 1> ArrayXb;

constexpr double biasFactor = 1.0/3.0; /// Used for bias correction.

class OneDiodeModel
{
private:
	GetSimOptFile* filePtr;
	struct Params inputParams; /// Its structure is defined in pvcommutils.h.
 
	std::vector<double> tol{ 0.0 };
	std::vector<double> tolBS{ 0.0 };
	
	ArrayXd iSet; /// Temporal current set.
	ArrayXd vSet; /// Temporal voltage set.

    double prod; /// Iph*di2Mutau.
    double commonClause1; /// Io*Rsh/nVth.
    double commonClause2; /// Iph*Rsh*di2mutau.
    double i0; /// Iph-prod/Vbi.
    double logIo; /// log(Io).
    double m0;

	ArrayXb isRev; /// judge if it is reverse-baised.
	ArrayXi indexIsRev;
    ArrayXd iRev; /// Reverse-baised voltage current.
    ArrayXd iRec; /// Recombination voltage current.

    size_t odmSize{ 0 }; /// Record the size of the above iSet and vSet.

    /// Obtain from SimOption.xml.
	unsigned maxIter;
	unsigned maxRecur;
	double   neps;

	ArrayXd d1Iset; /// First derivative.
	ArrayXd d2Iset; /// Second derivative.

	ArrayXd factor;
	
	GeneralRootFind* m_rf; /// A reference to the root-finding strategy.
	GeneralLinearApprox* m_la; /// A reference to the linear approximation strategy.

	int pos { 0 };

	/// Set a specific root-finding object at runtime.
	void setRootFind(GeneralRootFind* ptrRF) { m_rf = ptrRF; }
	void findRoot() const { m_rf->findRoot(); }
 	/// Set a specific linear approximation object at runtime.
	void setLinearApprox(GeneralLinearApprox* ptrLA) { m_la = ptrLA; } 
	void approxLine() const { m_la->approx(); }


	void fiveParameterODM(); /// Serve for the conversion from voltage to current.
	void approxIrev(); /// Approximate the reverse-baised voltage current.
	void fiveParameterODMInv(const ArrayXi& index); /// Serve for the conversion from current to voltage.
	void VocIsc(); /// Obtain the Voc and Isc

//	void BiApprox(double lx, double rx, double ly, double ry, MatrixXd* out, unsigned recurCount); //TODO: should be independent?? Only for the ODM? ask Martin
//	void BiApprox(double lx, double rx, double ly, double ry, unsigned recurCount); //TODO: should be independent?? Only for the ODM? ask Martin
public:
	/// Default constructor: create an empty object.
	OneDiodeModel() = default;

	/// Create an OneDiodeModel object with the given input parameters.
	OneDiodeModel(struct Params params_, const std::vector<double>& tol_ = {0.0}, 
		const unsigned long& maxIter_ = {(unsigned long)INITVAL});

    /// Characterize the MPP.
	double Vmp;
	double Imp;
	double Pmp;

	double Voc;
	double Isc;
	double c; /// TODO: what does the c for?

	MatrixXd approxPointSet; /// Ultimate output, which can be input for the class MonoDescPWL.

	/**
	 * Explicit 1st and 2nd derivatives.
	 *
	 * @param dim only 1st derivative or both of them.
	 */
	void getDerivatives(enum Dim dim);

	inline ArrayXd getIset() const { return iSet; }
	inline ArrayXd getVset() const { return vSet; }
	inline ArrayXd getD1Iset() const { return d1Iset; }
	inline ArrayXd getD2Iset() const { return d2Iset; }

	/// Calculate the corresponding iSet according to the given vSet_.
	void    VoltoCurr(ArrayXd vSet_);
	/// Calculate the corresponding vSet according to the given iSet_.
	void    CurrtoVol(ArrayXd iSet_);

    /**
     * Get the Maximum Power Point of this IV curve via the given seed Vmp.
     */
	void    getODMMpp(const double& seedVmp = {0.0}, const int& TolRefFac = {16});
    /**
     * Approximate the MDPWL curve, which is determined by the inputParams.
     * This function depends on all the other member functions. Its output is the approxPointSet.
     *
     * @param tol_ tolerance. Decide the appriximation precison. 
     * @param lim boundary. Restriction on the direction of X (voltage) and Y (current).
     * @param biasFactor bias. Used for the final bias correction.
     */
	void    pwlApprox( const std::vector<double>& tol_ = {}, const std::vector<double>& lim = {}, const double& biasK_= biasFactor);
	void 	dispRes() {std::cout<<"Pmp is "<<Pmp<<" Vmp is "<<Vmp<<" Imp is "<<Imp<<" Voc is "<<Voc<<" Isc is "<<Isc<<" c is "<<c<<"\n";}
};


/*
class OneDiodeModelCurrent : public OneDiodeModel
{
private:
	
	
public:
	OneDiodeModelCurrent(struct Params params_, ArrayXd vSet_, const std::vector<double>& tol_ = {0.0}, 
		const unsigned long& maxIter_ = {(unsigned long)INITVAL});
};

class OneDiodeModelVoltage : public OneDiodeModel
{
private:
	ArrayXd factor;
	

public:
	OneDiodeModelVoltage(struct Params params_, ArrayXd iSet_, const std::vector<double>& tol_ = {0.0}, 
		const unsigned long& maxIter_ = {(unsigned long)INITVAL});
};
*/

#endif
