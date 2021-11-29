/**
 * Implementation of a facility for monotonically descending Piece-Wise Line (MDPWL/IV curve)
 * 
 * This class consists of three constructors, which can be considered as three main components. 
 * The other member functions serve for the generation of the resulting MDPWL.
 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include "../Utils/xmlutils.h"
#include "../Utils/pvcommutils.h"
#include "../Utils/Sortbehavior.h"
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions> 

using namespace Eigen;
#ifndef MONODESCPWL_H_
#define MONODESCPWL_H_


class GetSimOptFile;
using Point = std::pair<double, double>;

typedef Eigen::Array<bool, 2, 2> Array2b; /// A square (2X2) bool array.
typedef Eigen::Array<bool, Dynamic, 1> ArrayXb; /// A bool array with one column and dynamic number of rows.
typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXRd;

class MonoDescPWL
{
private:

	std::vector<double> tol; /// Tolerance for determing the precision.
	std::vector<double> err;
	size_t psSize{ 0 }; /// The number of points in the generated IV curve (see pointSet below).
	Eigen::ArrayXd deriv; /// Store the first-order derivatives.
	std::vector<int> index; /// Used for the douglaspeucker algorithm.
	Eigen::Matrix2d lim; /// Boundary.
	Dim dim{ Dim::vert }; /// Differentiation between Curr and Vol.
	GetSimOptFile* filePtr; /// Pointer to the .xml file.
	double neps;

	GeneralSort* m_sort; /// A reference to a sorting strategy.

    /**
     * Locate the point that suits for evaluating the given xx.
     *
     * @return The postion of the suitable point.
     */
	int binScalIndex(const double& xx, const Dim& d);

	/**
	 * Linear interpolation according to the known points.
	 *
	 * @param begin the return value from binScalIndex.
	 * @return The interpolated point.
	 */
	double interpolate(int begin, double value, Dim d);

	/// The following three functions serve for douglaspeucker algorithm.
	void recursiveDP(int begin, int end);
	/// Initialization of the vector of Index.
	inline void preDP();
	/// Deal with Index vector and delete the 'tolerance-space' pointers.
	inline void postDP();

	/// First-order derivative at each point and write to deriv.
	void currDeriv();

    /// Clip the IV curve in both X(Vol) and Y(Curr) direction in terms of the specified boundaries.
	inline void lowerBdSubXClip(int& lowerBd);
	inline void upperBdSubXClip(int& upperbd);
	inline void lowerBdSubYClip(int& lowerBd);
	inline void upperBdSubYClip(int& upperBd);
	inline void XYClip(const int& lowerBd, const int& upperBd);

	/// Gaurantee that lim(0,1) is the upper-bound and lim(1,1) is the lower-bound.
	inline void swapLimI()
	{
		if (lim(0,1) < lim(1,1))
		{
			double tmp = lim(0,1);
			lim(0,1) = lim(1,1);
			lim(1,1) = tmp;
		}
	}

public:

	struct
	{
		double pmp;
		double vmp;
		double pmpBnd;
		double vmpBnd;
	} mppInfo; /// Store the MPP info

	MatrixXd pointSet; /// Store the generated MDPWL curve (in column major)
	MatrixXd backPS; /// A backup of pointSet, when necessary. This is used for the conversion of tolerance space.
	
	/// Default constructor: create an empty object.
	MonoDescPWL() = default;

	/** 
	 * Create a MDPWL curve from scratch according to the point set.
	 * The point set is of type "MatrixXd".
	 */
	MonoDescPWL(const Eigen::Ref<MatrixXd> pSet_, const std::vector<double>& tol_ = {EPSILON0},
		std::vector<double> lim_ = {-Inf, Inf, Inf, -Inf});

	/**
	 * Create a MDPWL curve according to the given two MDPWL curves and the related fraction.
	 * This equals to the interpolate(A, B, f) in the matlab.
	 */
	MonoDescPWL(MonoDescPWL c1, MonoDescPWL c2, double fraction);

	/**
	 * Create a MDPWL curve by combining a vector of curves in parallel or in series.
	 * This equals to addseries and addparallel in matlab.
	 *
	 * @param dim  Combine in parallel or in series. Must be either Dim::hori (in series) or Dim::vert (in parallel).
	 * @param tol_ Tolerance. Zero by default.
	 * @param lim_ Boundary. No boundary by default.
	 */
	MonoDescPWL(std::vector<MonoDescPWL*>& vecPS,
		const Dim& dim,
		const std::vector<double>& tol_ = {0.0},
		std::vector<double> lim_ = {-Inf, Inf, Inf, -Inf});

    /// Get the value of psSize
	size_t  getSize() const;

	/// Set a specific sorting object at runtime. see Sortbehavior.h
	void set_sort(GeneralSort* ptrS) { m_sort = ptrS; }
	void sort() const { m_sort->sort();  }

    /**
     * This functionality (sort) is not generic, but affiliated with this class.
     * Equal to unique sorting in the matlab. 
     * Partly sort the points with the identical Voltage according to the Current.
     */
	void uniqueSort();

	/// Convert xy coordinates to 'tolerance-space'
	void xy2TolSpace(); 
	inline void relTolSpace(Eigen::Ref<Eigen::ArrayXd> x, const double& threshold);
	
	/**
	 * Simplify the IV curve using Douglas-Peucker algorithm
	 *
	 * @see recrusiveDP
	 * @see preDP
	 * @see postDP
	 */
	void douglaspeucker();
	
	/// Guarantee the week monotonicity according to the voltage.
	void fixmonotonicity();

	/**
	 * Corresponds to fixverticals and fixhorizontals.
	 * Gaurantee the strict monotonicity vertically or horizontally.
	 * 
	 * @param d determine if it is fixverticals or fixhorzontals.
	 */
	void fixstrictMono(const Dim& d);

	/**
	 * Evaluate the corresponding voltage or current in the generated MDPWL curve, according to 
	 * the given xx.
	 * 
	 * @see binScalIndex.
	 * @see interpolate.
	 */
	double scalEval(double xx, const Dim& d);

	/**
	 * Obtain the maximum power point (MPP) of the generated MDPWL curve and update the mppInfo.
	 * 
	 * @param bound_ x and y boundary.
	 */
	void getMpp(std::vector<double> bound_);
    
    /**
     * Clip the curve according to the given lim/boundary.
     * 
     * @see lowerBdSubXClip.
	 * @see upperBdSubXClip.
	 * @see lowerBdSubYClip.
	 * @see upperBdSubYClip.
	 * @see XYClip.
     */
	void clip();

	//TODO: we need the destructor function??!!

#if 0

#ifdef VECVERSION
    std::vector<int> binIndex(std::vector<double> xx, Dim d); //break-point intervals
	std::vector<double> eval(std::vector<double> xx, Dim d);
#endif

#endif
	/*
	inline double& refPointSet(const int index, const int dim)
	{		return (dim == 0) ? (pointSet.at(index).first) : (pointSet.at(index).second);
	}*/
};

#endif // MONODESCPWL_H

