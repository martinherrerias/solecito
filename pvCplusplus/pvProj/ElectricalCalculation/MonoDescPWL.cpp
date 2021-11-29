#include "MonoDescPWL.h"
#include "../Utils/Compare.h"
#include <cstdlib>
#include <iomanip>
#include <random>
/// MSXML indicates the environment of visual studio.
#ifdef MSXML
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <vector>
#include <chrono>
#include <omp.h>

#include "../testUtils.h"

constexpr int ITRNUM = 10; /// the number of iterations for testing purpose
/// The below two constants serve for omp parallelization
/// Both of them are decided empirically and subject to change
constexpr int TASKTHLD = 256; /// The number of tasks
constexpr int SIZETHLD = 4000; /// The number of points
using namespace Eigen;

class MonoDescPWL;

/// Combination of curves in parallel or in series
MonoDescPWL::MonoDescPWL(std::vector<MonoDescPWL*>& vecPS,
	const Dim& dim,
	const std::vector<double>& tol_,
	std::vector<double> lim_
)
	: tol(tol_), lim(Eigen::Map<Eigen::Matrix2d, Eigen::Unaligned>(std::move(lim_.data())))
{
	ASSERT((vecPS.size() >= 2), "Expecting 2+ vector of MDPWL objects");
	filePtr = GetSimOptFile::Instance();
	neps    = (filePtr->getElement<double>)("NEPS"); 

//	std::vector<double> tol_backup = tol;
	filePtr->parseTolerance(tol);

	/// Extracting the unique PWL curves and the corresponding weight value.
	std::vector<MonoDescPWL*> uniqueVecPS;
	std::vector<int> w; /// Accumulated weight
	int counterDim = (dim + 1) % dimNum;

	uniqueVecPS.clear();
	uniqueVecPS.push_back(vecPS.front());
	w.push_back(ONENUM);

    /// Fill uniqueVecPS and w.
	for (auto it = vecPS.begin() + 1; it != vecPS.end(); it++)
	{
		int dist = std::distance(vecPS.begin(), it);
		bool flag = true;
		for (auto innerIt = uniqueVecPS.begin(); innerIt != uniqueVecPS.end(); innerIt++)
			if ((*innerIt) == (*it))
			{
				++w.at(std::distance(uniqueVecPS.begin(), innerIt));
				flag = false;
				break;
			}
		if (flag)
		{
			uniqueVecPS.push_back(*it);
			w.push_back(ONENUM);
		}
	}

	Array2b clipped;
	clipped = Array2b::Zero(2, 2); /// By default no clip for combined curve
	bool noclipped = lim.array().isInf().all();
	swapLimI();

    /**
     * Add parallel, multiply the current (Y), Add series, multiply the voltage (X).
     */
	if (w.size() == 1) 
	{
		pointSet = uniqueVecPS.front()->pointSet;
		psSize = uniqueVecPS.front()->getSize(); //BUGFIX: update set the value of psSize, otherwise it is invalid
		pointSet.col(counterDim) *= w.front();

		if (!noclipped) clip(); // TODO: recalculate the flag of noclipped

		if (std::any_of(tol.cbegin(), tol.cend(), [](double x) {return x != 0.0; }))
		{
			backPS = pointSet;
			xy2TolSpace();
			douglaspeucker();
		}
		//TODO: add the fixvertical or fixhorizontal
		return;
	}

	psSize = std::accumulate(uniqueVecPS.begin(), uniqueVecPS.end(), 0,
		[](size_t a, MonoDescPWL* b) { return a + b->getSize(); });

	pointSet.resize(psSize + 2, 3);
	int rowIndex = 0;
	unsigned uniquePWLSize = w.size();

	/*
	lim(0, 0) = 0.0;
	lim(1, 0) = 5.0;
	lim(0, 1) = 6.5;
	lim(1, 1) = 1.0;*/

	std::vector<int> rowIndexSet;
	Matrix2d minMax;
	minMax << Inf, -Inf,
			  -Inf, Inf;

	/// Combine the unique curves end to end and record the length of each curve into rowIndexSet.
	for (unsigned i = 0; i < uniquePWLSize; i++)
	{
		int lb = 0; ///lower bound
		int ub = uniqueVecPS.at(i)->getSize(); ///upper bound and the size

		double* startAddr = uniqueVecPS.at(i)->pointSet.col(dim).data();

		if (!noclipped)
		{
			if (uniqueVecPS.at(i)->pointSet(0, 1) > minMax(0, 1))
				minMax(0, 1) = uniqueVecPS.at(i)->pointSet(0, 1);
			if (uniqueVecPS.at(i)->pointSet(ub - 1, 1) < minMax(1, 1))
				minMax(1, 1) = uniqueVecPS.at(i)->pointSet(ub - 1, 1);
			if (uniqueVecPS.at(i)->pointSet(0, 0) < minMax(0, 0))
				minMax(0, 0) = uniqueVecPS.at(i)->pointSet(0, 0);
			if (uniqueVecPS.at(i)->pointSet(ub - 1, 0) > minMax(1, 0))
				minMax(1, 0) = uniqueVecPS.at(i)->pointSet(ub - 1, 0);
		}

		if (dim)
		{
		/*	if (!noclipped)
			{
				if (uniqueVecPS.at(i)->pointSet(0, dim) > minMax(0, dim))
					minMax(0, dim) = uniqueVecPS.at(i)->pointSet(0, dim);
				if (uniqueVecPS.at(i)->pointSet(ub - 1, dim) < minMax(1, dim))
					minMax(1, dim) = uniqueVecPS.at(i)->pointSet(ub - 1, dim);
			}*/
			if (uniqueVecPS.at(i)->pointSet(0, dim) >= lim(0, dim))
				lb = std::lower_bound(startAddr, startAddr + ub, lim(0, dim), customDoubleUpper{}) - startAddr;

			if (uniqueVecPS.at(i)->pointSet(ub - 1, dim) <= lim(1, dim))
				ub = std::lower_bound(startAddr, startAddr + ub, lim(1, dim), customDoubleDesc{}) - startAddr;
		}
		else
		{
		/*	if (!noclipped)
			{
				if (uniqueVecPS.at(i)->pointSet(0, dim) < minMax(0, dim))
					minMax(0, dim) = uniqueVecPS.at(i)->pointSet(0, dim);
				if (uniqueVecPS.at(i)->pointSet(ub - 1, dim) > minMax(1, dim))
					minMax(1, dim) = uniqueVecPS.at(i)->pointSet(ub - 1, dim);
			}*/
			/// The first element greater than lim(0, dim)
			if (uniqueVecPS.at(i)->pointSet(0, dim) <= lim(0, dim))
				lb = std::upper_bound(startAddr, startAddr + ub, lim(0, dim)) - startAddr;

			/// The first element greater or equal to lim(1, dim)
			if (uniqueVecPS.at(i)->pointSet(ub - 1, dim) >= lim(1, dim)) 
				ub = std::lower_bound(startAddr, startAddr + ub, lim(1, dim)) - startAddr;
		}

		pointSet.block(rowIndex, 0, ub - lb, 2) = uniqueVecPS.at(i)->pointSet.block(lb, 0, ub - lb, 2);
    	pointSet.block(rowIndex, 0, ub - lb, 2).col(counterDim).array() *= w.at(i);
		pointSet.block(rowIndex, 2, ub - lb, 1) = MatrixXd::Constant(ub - lb, 1, i);

		rowIndexSet.push_back(rowIndex);
		rowIndex += (ub - lb);
	}

	rowIndexSet.push_back(rowIndex);

//	if (dim)
	{
		clipped(0, 1) = minMax(0, 1) > lim(0, 1);
		clipped(1, 1) = minMax(1, 1) < lim(1, 1);
	}
//	else
	{
		clipped(0, 0) = minMax(0, 0) < lim(0, 0);
		clipped(1, 0) = minMax(1, 0) > lim(1, 0);
	}

    if (clipped(0, dim))
	{
		pointSet(rowIndex, dim) = lim(0, dim);
		pointSet(rowIndex, counterDim) = INITVAL;
		pointSet(rowIndex, 2) = UNDEFINED;
		++rowIndex;
		rowIndexSet.push_back(rowIndex);
	}
	if (clipped(1, dim))
	{
		pointSet(rowIndex, dim) = lim(1, dim);
		pointSet(rowIndex, counterDim) = INITVAL;
		pointSet(rowIndex, 2) = UNDEFINED;
		++rowIndex;
		rowIndexSet.push_back(rowIndex);
	}

    /// Note: here only decrease psSize to rowIndex instead of using conservativeSize on pointSet 
	psSize = rowIndex;
	if (dim) pointSet.col(dim).array() *= -1;

	/// use merge sorting strategy to sort the above combined curves.
	MergeSort merge(&pointSet, 0, psSize - 1, dim, &rowIndexSet);
	set_sort(&merge); /// Sort direction (current or voltage) is determined by dim

	if (dim) /// Dim::hori: add in series
	{
		if (!std::is_sorted(pointSet.col(dim).data(), pointSet.col(dim).data() + psSize, customDoubleDesc{}))
			sort();
	}
	else /// Dim::vert: add in parallel
	{
		if (!std::is_sorted(pointSet.col(dim).data(), pointSet.col(dim).data() + psSize))
			sort();
	}

	ArrayXd diffPointSet(psSize - 1), tolPS(psSize - 1);
	for (unsigned i = 0; i < psSize - 1; i++)
	{
		if (dim)
			diffPointSet(i) = DY(i, i + 1);
		else
			diffPointSet(i) = DX(i, i + 1);

		tolPS(i) = (filePtr->tolFun<double>)(tol, pointSet(i + 1, dim), neps, dim);
	}

	if (dim) pointSet.col(dim).array() *= -1; /// Only for the current (Y) direction	 
	ArrayXb tooclose = (diffPointSet < tolPS);
	ArrayXb used = ArrayXb::Ones(psSize);
	ArrayXb traversed = ArrayXb::Ones(psSize);

	bool comparedVal[] = { true };
	int resCount = 0;

	auto it = std::find_first_of(tooclose.data(), tooclose.data() + psSize - 1, comparedVal, comparedVal + 1);

	std::vector<std::vector<int>> closePoints(psSize);
	int dist1, initDist;
    
	while (it != tooclose.data() + psSize - 1)
	{
		auto innerIt = it;
		dist1 = std::distance(tooclose.data(), innerIt);
	    initDist = dist1;
	    closePoints[initDist].push_back(pointSet(initDist, 2));
      	++dist1;
	    traversed(initDist) = false;

		while ((dist1 != psSize) && (*innerIt))
		{/** a. delete the repeated point in the same curve (this will not happen)
		  *	 b. do not delete two vertical/horizontal points in the same curve (should be also solved in fixvertical/fixhorizontal)
	      *  c. remove one of the two close points in different curve and then adapt the value of V and I
          *  TODO: ask Martin if we should measure the different of pointSet(dist1,dim) and pointSet(initDist,dim) is with the tolerance tolPS? the corresponding matlab code is something wrong?
	      */		
			if (pointSet(dist1, 2) != UNDEFINED && std::find(closePoints[initDist].begin(), closePoints[initDist].end(), pointSet(dist1, 2)) == closePoints[initDist].end())
			{
				pointSet(initDist, dim) = (pointSet(initDist, dim) + pointSet(dist1, dim)) / 2;
				pointSet(initDist, counterDim) += pointSet(dist1, counterDim);
				closePoints[initDist].push_back(pointSet(dist1, 2));
				//	pointSet(dist1 + 1, 2) = UNDEFINED;
				used(dist1) = false;
				traversed(dist1) = false;
				++resCount;
			}
      		else
      		{
        		break;
      		}
			++innerIt;
      		++dist1;
		}

/*
		for (int d = initDist + 1; d < dist1; d++)
		{
			if (traversed(d))
			{
		i		closePoints[d].push_back(pointSet(d, 2));
				traversed(d) = false;
			    int interd = d;

				while (++interd < dist1)
				{
					if (traversed(interd))
					{
						if ((pointSet(d, dim) - pointSet(interd, dim) < tolPS(d)) && std::find(closePoints[d].begin(), closePoints[d].end(), pointSet(interd, 2)) == closePoints[d].end())
						{
							pointSet(d, dim) = (pointSet(d, dim) + pointSet(interd, dim)) / 2;
							pointSet(d, counterDim) += pointSet(interd, counterDim);
							closePoints[d].push_back(pointSet(interd, 2));
							//  TODO: if here we need to multiply w[i]?
							//	pointSet(dist1 + 1, 2) = UNDEFINED;
							used(interd) = false;
							traversed(interd) = false;
							++resCount;
						}
					}
				}
			}
		}*/

		if (dist1 != psSize)
    	{
			it = std::find_first_of(++innerIt, tooclose.data() + psSize - 1, comparedVal, comparedVal + 1);
    	} else { it = innerIt; }
	}

	/// interpolate to complete the combined curve
#pragma omp parallel if (psSize>SIZETHLD)
	{
	unsigned i, j;
#pragma omp for private(j)
	for (j = 0; j < psSize; j++)	
	{ /// the third column in pointSet is used to map the certain point to the related curve
		for (i = 0; i < uniquePWLSize; i++)
		{
			if (used(j) && pointSet(j, 2) != i)
      		{
				if ((closePoints[j].empty() || (!closePoints[j].empty()) && std::find(closePoints[j].begin(), closePoints[j].end(), i) == closePoints[j].end()))
					pointSet(j, counterDim) += uniqueVecPS.at(i)->scalEval(pointSet(j, dim), dim)*w[i];
      		}
      		/// Interpolate the jth point on the ith curve and add to the original point
			/// For UNDEFINED, then interpolate at each time
		}
	}
	}

	pointSet.conservativeResize(psSize, 2);
  	MatrixXd tmpPS(psSize - resCount, 2);
	unsigned j = 0;
	for (unsigned i = 0; i < psSize; i++)/// i is the row number that is to be removed
		if (used(i))
			tmpPS.row(j++) = pointSet.row(i);

	pointSet.resize(0, 0);
	pointSet = std::move(tmpPS);
	psSize -= resCount;

	std::string exceedMsg = std::string(", exceed the range of lim when adding curves");

	if (dim)
	{
		if (clipped(0, counterDim))
		{
			if (pointSet(psSize - 1, counterDim) < lim(0, counterDim))
				CATCHRTERROR(exceedMsg);
		}

		if (clipped(1, counterDim))
		{
			if (pointSet(0, counterDim) > lim(1, counterDim))
				CATCHRTERROR(exceedMsg);
		}
	}
	else
	{
		if (clipped(0, counterDim))
		{
			if (pointSet(psSize - 1, counterDim) > lim(0, counterDim))
				CATCHRTERROR(exceedMsg);	
		}

		if (clipped(1, counterDim))
		{
			if (pointSet(0, counterDim) < lim(1, counterDim))
				CATCHRTERROR(exceedMsg);
		}
	}

	clip();

	if (std::any_of(tol.cbegin(), tol.cend(), [](double x) { return x != 0.0; }))
	{
		backPS = pointSet;
		xy2TolSpace();
		douglaspeucker();
	} //BUGFIX: add simplification again and check with the tol if it is correct?

#ifndef NDEBUG
	std::cout << "Combined pointSet is :\n" << pointSet << "\n";
#endif
}

MonoDescPWL::MonoDescPWL(MonoDescPWL c1, MonoDescPWL c2, double fraction)
{
	ASSERT((fraction >= 0) && (fraction <= 1), "fraction must be within [0,1]");
	MonoDescPWL a, b;
	double diff;

	if (fraction <= 0.5)
	{
		a = std::move(c1);
		b = std::move(c2);
	}
	else
	{
		a = std::move(c2);
		b = std::move(c1);
		fraction = 1 - fraction;
	}

	psSize = a.getSize();
	pointSet.resize(psSize, 2);

	diff = b.scalEval(0, Dim::vert) - a.scalEval(0, Dim::vert);
	pointSet.col(1).array() = a.pointSet.col(1).array() + diff;
	for (unsigned i = 0; i < psSize; i++)
		pointSet(i, 0) = b.scalEval(pointSet(i, 1), Dim::hori);

	pointSet.col(1).array() = a.pointSet.col(1).array() + fraction * diff;
	pointSet.col(0).array() = (1 - fraction) * a.pointSet.col(0).array() +
		fraction * pointSet.col(0).array();

#ifndef NDEBUG
	std::cout << "P3: \n"<<pointSet << std::endl;
#endif
	tol = { Inf };
	std::vector<double> vecLim = { -Inf, Inf, Inf, -Inf };
	lim = Eigen::Map<Eigen::Matrix2d, Eigen::Unaligned>(std::move(vecLim.data()));
}


/** Standard constructor.
 * dim is by default Dim::vert.
 * Here the property dim in class is a feature related to the specific type of object.
 */
MonoDescPWL::MonoDescPWL(const Ref<MatrixXd> pSet_,
	const std::vector<double>& tol_, std::vector<double> lim_)
	: pointSet(pSet_), tol(tol_), lim(Eigen::Map<Eigen::Matrix2d, Eigen::Unaligned>(std::move(lim_.data())))
{
	// Fast version: skip all checks if tol = Inf
	
	bool noclipped = lim.array().isInf().all();
	swapLimI();
	if ((tol.size() != 1) || (tol.front() != Inf) || !noclipped)
	{
		std::vector<double> tol_backup = tol;
		filePtr = GetSimOptFile::Instance();
		filePtr->parseTolerance(tol); //	std::cout.precision(19);
		bool isNanInf = false;

		if ((pointSet.col(0).array().isNaN()).any()
			|| (pointSet.col(1).array().isNaN()).any()
			|| (pointSet.col(0).array().isInf()).any()
			|| (pointSet.col(1).array().isInf()).any())
			isNanInf = true;

		psSize = pointSet.rows();

		ASSERT((psSize > 1), "mdpwl:npts, At least 2-non-zero-points are required");
		ASSERT((!isNanInf), "mdpwl:NaN, NaN/Inf values in x or y");
		ASSERT((lim.size() == 4), "mdpwl:lim, Limits must be a 4-Vector");

		/// Quick sorting strategy is selected.
		QuickSort qs(&pointSet, 0, psSize - 1, dim);
		set_sort(&qs);
		if (!std::is_sorted(pointSet.col(0).data(), pointSet.col(0).data() + psSize))
			sort();
		uniqueSort();

		if (!noclipped) clip();

 		if (std::any_of(tol_backup.cbegin(), tol_backup.cend(), [](double x) {return x != 0.0; }))
		{
			backPS = pointSet;
			xy2TolSpace();
			douglaspeucker();
		}

    	fixmonotonicity();

		fixstrictMono(Dim::vert); //fixverticals
		fixstrictMono(Dim::hori); //fixhorizontals
	}
}

inline void MonoDescPWL::lowerBdSubYClip(int& lowerBd)
{
	lowerBd = std::lower_bound(pointSetData(1), pointSetData(1) + psSize, lim(0, 1), customDoubleUpper{})
		- pointSetData(1) - 1;

	if (pointSet(lowerBd, 1) != lim(0, 1))
	{
		pointSet(lowerBd, 0) = interpolate(lowerBd, lim(0, 1), Dim::hori);
		pointSet(lowerBd, 1) = lim(0, 1);
	}
}

inline void MonoDescPWL::lowerBdSubXClip(int& lowerBd)
{
	lowerBd = std::upper_bound(pointSetData(0), pointSetData(0) + psSize, lim(0, 0))
		- pointSetData(0) - 1; // /larger than lim.at(0)

	if (pointSet(lowerBd, 0) != lim(0, 0))
	{
		pointSet(lowerBd, 1) = interpolate(lowerBd, lim(0, 0), Dim::vert);
		/// interpolate the y(a) according to the point a and a+1
		pointSet(lowerBd, 0) = lim(0, 0);
	}
}

inline void MonoDescPWL::upperBdSubYClip(int& upperBd)
{
	upperBd = std::lower_bound(pointSetData(1), pointSetData(1) + psSize, lim(1, 1), customDoubleDesc{})
		- pointSetData(1);

	if (pointSet(upperBd, 1) != lim(1, 1))
	{
		pointSet(upperBd, 0) = interpolate(upperBd - 1, lim(1, 1), Dim::hori);
		pointSet(upperBd, 1) = lim(1, 1);
	}
}

inline void MonoDescPWL::upperBdSubXClip(int& upperBd)
{
	upperBd = std::lower_bound(pointSetData(0), pointSetData(0) + psSize, lim(1, 0))
		- pointSetData(0);

	if (pointSet(upperBd, 0) != lim(1, 0))
	{
		pointSet(upperBd, 1) = interpolate(upperBd - 1, lim(1, 0), Dim::vert);
		pointSet(upperBd, 0) = lim(1, 0);
	}
}


/// clip the pointSet[0:end] to pointSet[lowerBd:upperBd]
inline void MonoDescPWL::XYClip(const int& lowerBd, const int& upperBd)
{
	psSize = upperBd - lowerBd + 1;
	pointSet.block(0, 0, psSize, 2) = pointSet.block(lowerBd, 0, psSize, 2); //TODO: give explanation to the usage of block
}

void MonoDescPWL::currDeriv()
{
	ArrayXd derivY(psSize);
	deriv.resize(psSize);
	deriv(psSize - 1) = 0.0;
	derivY(psSize - 1) = 0.0;
	for (unsigned r = 0; r < psSize - 1; r++)
	{
		deriv(r) = pointSet(r + 1, 0) - pointSet(r, 0);
		derivY(r) = pointSet(r + 1, 1) - pointSet(r, 1);
	}

	deriv = derivY / deriv;
}

size_t MonoDescPWL::getSize() const
{
	return psSize;
}

void MonoDescPWL::uniqueSort()
{
	int begin, end;
	for (unsigned i = 0; i < psSize - 1; i++)
	{
		begin = -1;
		while ((i < psSize - 1)
			&& ((pointSet(i+1, 0) - pointSet(i, 0)) == 0.0))
		{
			if (begin < 0)
				begin = i;
			i++;
		}
		if ((begin > 0) && (i <= psSize - 1))
		{
			std::sort(pointSet.col(1).data() + begin, pointSet.col(1).data() + i + 1, customDoubleDesc{});
		}
	}

#ifndef NDEBUG
	for (unsigned i = 0; i < psSize; i++)
	{
		std::cout << pointSet(i, 0) << ", " << pointSet(i, 1) << "\n";
	}
	std::cout << "\n";
#endif
}

/**
 * range of x (first):  [lim(0, 0)/Xmin, lim(1, 0)/Xmax]
 * range of y (second): [lim(0, 1)/Ymax, lim(1, 1)/Ymin]
 *
 * clipped and extended correspond to lim
 */
void MonoDescPWL::clip()
{
	//TODO: if ~isscalar arrayfun....
	if (pointSet.rows() == 0) return;

	//TODO: this part could be repeated
	Array2b clipped;
	clipped(0, 0) = pointSet(0, 0) < lim(0, 0);
	clipped(1, 0) = pointSet(psSize - 1, 0) > lim(1, 0);
	clipped(0, 1) = pointSet(0, 1) > lim(0, 1);
	clipped(1, 1) = pointSet(psSize - 1, 1) < lim(1, 1);


	Array2b extended;
	extended(0, 0) = pointSet(0, 0) > lim(0, 0);
	extended(1, 0) = pointSet(psSize - 1, 0) < lim(1, 0);
	extended(0, 1) = pointSet(0, 1) < lim(0, 1);
	extended(1, 1) = pointSet(psSize - 1, 1) > lim(1, 1);

	extended = extended && lim.array().isFinite();
	//(clipped.unaryExpr([](bool x) { return !x; }))

	if (!((clipped || extended).count())) return;

	if (extended(0, 0))
	{
		pointSet(0, 1) = interpolate(0, lim(0, 0), Dim::vert);
		pointSet(0, 0) = lim(0, 0);
		clipped(0, 1) = pointSet(0, 1) > lim(0, 1);
		if (extended(0, 1))
			extended(0, 1) = false;
	}

	if (extended(1, 0))
	{
		pointSet(psSize - 1, 1) = interpolate(psSize - 2, lim(1, 0), Dim::vert);
		pointSet(psSize - 1, 0) = lim(1, 0);
		clipped(1, 1) = pointSet(psSize - 1, 1) < lim(1, 1);
		if (extended(1, 1))
			extended(1, 1) = false;
	}

	if (extended(0, 1))
	{
		pointSet(0, 0) = interpolate(0, lim(0, 1), Dim::hori);
		pointSet(0, 1) = lim(0, 1);
		clipped(0, 0) = pointSet(0, 0) < lim(0, 0);
	}

	if (extended(1, 1))
	{
		pointSet(psSize - 1, 0) = interpolate(psSize - 2, lim(1, 1), Dim::hori);
		pointSet(psSize - 1, 1) = lim(1, 1);
		clipped(1, 0) = pointSet(psSize - 1, 0) > lim(1, 0);
	}

	int a = 0;
	int b = psSize - 1;

	if (clipped(0, 0))
	{
		lowerBdSubXClip(a);
		clipped(0, 1) = pointSet(a, 1) > lim(0, 1);

		/*
		a = std::upper_bound(pointSetData(0), pointSetData(0) + n, lim(0, 0))
			- pointSetData(0) - 1; // larger than lim.at(0)

		if (pointSet(a, 0) != lim(0, 0))
		{
			pointSet(a, 1) = interpolate(a, a + 1, lim(0, 0), Dim::vert);
			// interpolate the y(a) according to the point a and a+1
			pointSet(a, 0) = lim(0, 0);
			clipped(0, 1) = pointSet(a, 1) > lim(0, 1);
		}*/
	}

	if (clipped(1, 0))
	{
		upperBdSubXClip(b);
		clipped(1, 1) = pointSet(b, 1) < lim(1, 1);

		/*
		b = std::lower_bound(pointSetData(0), pointSetData(0) + n, lim(1, 0))
			- pointSetData(0);

		if (pointSet(b, 0) != lim(1, 0))
		{ //TODO in no case there are equal
			pointSet(b, 1) = interpolate(b - 1, b, lim(1, 0), Dim::vert);
			pointSet(b, 0) = lim(1, 0);
			clipped(1, 1) = pointSet(b, 1) < lim(1, 1);
		}*/
	}

	if (clipped(0, 1))
	{
		lowerBdSubYClip(a);

		/*
		a = std::lower_bound(pointSetData(1), pointSetData(1) + n, lim(0, 1), customDoubleUpper{})
			- pointSetData(1) - 1;

		pointSet(a, 0) = interpolate(a, a + 1, lim(0, 1), Dim::hori);
		pointSet(a, 1) = lim(0, 1);*/
	}

	//TODO map from 2 to Dim::HORI
	if (clipped(1, 1))
	{
		upperBdSubYClip(b);

		/*
		b = std::lower_bound(pointSetData(1), pointSetData(1) + n, lim(1, 1), customDoubleDesc{})
			- pointSetData(1);

		if (pointSet(b, 1) != lim(1, 1))
		{//TODO: in no case there are equal
			pointSet(b, 0) = interpolate(b - 1, b, lim(1, 1), Dim::hori);
			pointSet(b, 1) = lim(1, 1);
		}*/
	}

	XYClip(a, b);
    pointSet.conservativeResize(psSize, 2);

#ifndef NDEBUG
	std::cout << "\nAfter clip:\n" << pointSet << std::endl;
#endif
}


inline void MonoDescPWL::relTolSpace(Eigen::Ref<Eigen::ArrayXd> x, const double& threshold)
{
	x /= threshold;

	// coefficient-wise operations
	x = (x.abs() > 1).select(SIGN(x)*(x.abs().log() + 1), x);

//	if (ABS(x) > 1)
//		x = SIGN(x)*(std::log(ABS(x)) + 1); //base-e log
}

void MonoDescPWL::xy2TolSpace()
{
	/// Change the value of backPS and later the pointSet is used without reverting the tolerance space.
	if (tol.at(2) <= 1)
	{
		double tolDiv1, tolDiv2;
		tolDiv1 = tol.at(0) / tol.at(2);
		tolDiv2 = tol.at(1) / tol.at(2);

		relTolSpace(backPS.col(0).array(), tolDiv1);
		relTolSpace(backPS.col(1).array(), tolDiv2);
	}
	else
	{
		backPS.col(0) /= tol.at(0);
		backPS.col(1) /= tol.at(1);
	}

/*
	if (tol.at(2) <= 1)
		for (unsigned i = 0; i < n; i++)
		{
			refPointSet(i, 0) = relTolSpace(refPointSet(i, 0), tol.at(0) / tol.at(2));
			refPointSet(i, 1) = relTolSpace(refPointSet(i, 1), tol.at(1) / tol.at(2));

		}
	else
		for (unsigned i = 0; i < n; i++)
		{
			refPointSet(i, 0) = refPointSet(i, 0) / tol.at(0);
			refPointSet(i, 1) = refPointSet(i, 1) / tol.at(1);
		}*/
}

#if 0
inline void MonoDescPWL::absTolSpace(Eigen::Ref<Eigen::ArrayXd> x, double threshold)
{
	//TODO: store the value of sign and isAbs of the pointSet , maintian it or put them into a macro, calculate it when use it??
//	if (ABS(x) > 1)
//		x = SIGN(x)*std::exp(ABS(x) - 1); // '((x>0)-(x<0))' indicates the sign of the number

	ArrayXd y = (x.abs() - 1).exp();

	x = (x.abs() > 1).select(y * ((x > 0).cast<double>() - (x < 0).cast<double>()), x);
	x = x * threshold;
}



void MonoDescPWL::xy2AbsSpace() //Inverse of xy2TolSpace
{
	if (tol.at(2) <= 1)
	{
		//TODO: Store this value into the class
		double tolDiv1, tolDiv2;
		tolDiv1 = tol.at(0) / tol.at(2);
		tolDiv2 = tol.at(1) / tol.at(2);
		absTolSpace(pointSet.col(0).array(), tolDiv1);
		absTolSpace(pointSet.col(1).array(), tolDiv2);
	}
	else
	{
		pointSet.col(0).array() *= tol.at(0);
		pointSet.col(1).array() *= tol.at(1);
	}


#ifndef NDEBUG
	for (unsigned i = 0; i < n; i++)
	{
		std::cout << pointSet(i, 0) << ", " << pointSet(i, 1) << "\n";
	}
	std::cout << "\n";
#endif

	/*

	if (tol.at(2) <= 1)
		for (unsigned i = 0; i < n; i++)
		{
			refPointSet(i, 0) = absTolSpace(refPointSet(i, 0), tol.at(0) / tol.at(2));
			refPointSet(i, 1) = absTolSpace(refPointSet(i, 1), tol.at(1) / tol.at(2));
		}
	else
		for (unsigned i = 0; i < n; i++)
		{
			refPointSet(i, 0) = refPointSet(i, 0) * tol.at(0);
			refPointSet(i, 1) = refPointSet(i, 1) * tol.at(1);
		}*/
}

double MonoDescPWL::perpendicularDistance(int i, int begin, int end)
{
	double dx = DX(begin, end);
	double dy = DY(begin, end);
	double mag = SQRT2POW(dx, dy);

	if (mag > 0.0)
	{
		dx /= mag;
		dy /= mag;
	}

	double pvx = DX(begin, i);
	double pvy = DY(begin, i);
	double pvdot = dx * pvx + dy * pvy;

	double dsx = pvdot * dx;
	double dsy = pvdot * dy;

	double ax = pvx - dsx;
	double ay = pvy - dsy;

	return SQRT2POW(ax, ay);
}
#endif

void MonoDescPWL::recursiveDP(int begin, int end)
{
	if (end - begin > 1)
	{
		double dmax = 0.0;
		int position = 0;

		ParametrizedLine<double, 2> line1;
		/// Form a line
		line1 = ParametrizedLine<double, 2>::Through(backPS.row(begin), backPS.row(end));

		for (int i = begin + 1; i < end; i++)
		{
			double dist2 = line1.distance(backPS.row(i));
			if (dist2 > dmax)
			{
				dmax = dist2;
				position = i;
			}
		}

		if (dmax > tol.at(2))
		{
			if (psSize < SIZETHLD)
			{
				recursiveDP(begin, position);
				recursiveDP(position, end);
			}
			else {
#pragma omp task firstprivate(begin, position) if (position - begin > TASKTHLD)
			{
			recursiveDP(begin, position);
			}
#pragma omp task firstprivate(position, end) if (end - position > TASKTHLD)
			{
			recursiveDP(position, end);
			}
#pragma omp taskwait
			}
		}
		else
		{
			/// Indicates the points that should be removed.
			for (int j = begin + 1; j < end; j++)
				index.at(j) = -1;
		}

/// Original method calculating the maximum distance
/*
		double dx = ps(end, 0) - ps(begin, 0);
		double dy = ps(end, 1) - ps(begin, 1);
		double mag = SQRT2POW(dx, dy);

		if (mag > 0.0)
		{
			dx /= mag;
			dy /= mag;
		}

		ArrayXd pvx = ps.col(0).block(begin + 1, 0, end - begin - 1, 1).array() - ps(begin, 0);

		ArrayXd pvy = ps.col(1).block(begin + 1, 0, end - begin - 1, 1).array() - ps(begin, 1);

	    ArrayXd pvdot = dx * pvx + dy * pvy;

		ArrayXd dsx = pvdot * dx;
		ArrayXd dsy = pvdot * dy;

		ArrayXd ax = pvx - dsx;
		ArrayXd ay = pvy - dsy;

		ArrayXd resl = (ax.pow(2) + ay.pow(2)).sqrt();
		int indexX;
		dmax = resl.maxCoeff(&indexX);
		position = indexX;*/
	}
}

void inline MonoDescPWL::preDP()
{
	index.assign(psSize, 0);
	for (unsigned i = 0; i < psSize; i++)
		index.at(i) = i;
}

void MonoDescPWL::postDP()
{
	int tail;
	tail = -1;
	for (unsigned i = 0; i < psSize; i++)
	{
		if (index.at(i) == -1)
		{
			if (tail == -1)
				tail = i;
		}
		else
		{
			if (tail != -1)
			{
				pointSet.row(tail).swap(pointSet.row(i));
				tail++;
			}
		}
	}
	psSize = (tail == -1) ? psSize : tail;
	pointSet.conservativeResize(psSize, 2);
}

void MonoDescPWL::douglaspeucker()
{
	preDP();

#ifndef NDEBUG
	for (unsigned i = 0; i < psSize; i++)
	{
		std::cout << backPS.row(i) << "\n";
	}
	std::cout << "\n";
#endif

if (psSize < SIZETHLD)
    recursiveDP(0, psSize - 1);
else{
#pragma omp parallel
	{
#pragma omp single
	{  
	recursiveDP(0, psSize - 1);
	}
	}
	}	

	postDP();


#ifndef NDEBUG
	std::cout << "after douglas-peucker: \n";
	std::cout << pointSet;
	std::cout << "\n" << psSize << "\n";
#endif
}

void MonoDescPWL::fixmonotonicity()
{
	if (std::is_sorted(pointSetData(1), pointSetData(1) + psSize, customDoubleDesc{}))
	{
 		err.assign(1, 0.0);
		return;
	}


 	backPS = pointSet;

	std::sort(pointSetData(1), pointSetData(1) + psSize, customDoubleDesc{}); //Only sort the y

	ArrayXd dx = ArrayXd::Zero(psSize);
	ArrayXd dy = ArrayXd::Zero(psSize);
	ArrayXd dr = ArrayXd::Zero(psSize);
	ArrayXd e  = ArrayXd::Zero(psSize);

	dy = (pointSet.col(1).array() - backPS.col(1).array()).abs();
	dr = (dy / pointSet.col(1).array()).abs();
	dx = (dr * pointSet.col(0).array()).abs();
	e = ((dx / tol.at(0)).min(dy / tol.at(1))).min(dr / tol.at(2));

	if ((e > 1.0).any())
	{
#ifndef NDEBUG
		std::cerr << "Warning: mdpwl:mono, Points not strictly sorted, RMS = "
			<< STDDEV(e) << " RelTol\n";
#endif
		err.clear();
		err.push_back(dx.maxCoeff());
		err.push_back(dy.maxCoeff());
		err.push_back(dr.maxCoeff());
	}
}

/// Make sure there is no equal values in the horizontal and vertical directions
void MonoDescPWL::fixstrictMono(const Dim& d)
{
	/**
	 * d: Dim::VERT
	 *       Dim::HORI
	 */
	double stride;
	double* beginPtr = pointSetData(d);
	auto it_curr = std::adjacent_find(beginPtr, beginPtr + psSize);

  	while (it_curr != (beginPtr + psSize))
	{
		stride = (d == 0) ? 1.0 : -1.0;
		auto it_begin = it_curr++;
		do { ++it_curr; } while (it_curr != (beginPtr + psSize) && *(it_curr) == *(it_begin));
		//TODO if we should use upper_bound to replace

		if (it_curr != (beginPtr + psSize))
		{
			int dist = std::distance(it_begin, it_curr);
			stride = (*(it_curr)-*(it_begin)) / dist;
		} //TODO: if this is necessary?
	
		stride *= EPSILONScal<double>(*it_begin);
        //stride *= EPSILON;
		for (auto i = it_begin + 1; i != it_curr; ++i)
			*i = *(i-1) + stride;

		if (it_curr != (beginPtr + psSize))
			it_curr = std::adjacent_find(it_curr, beginPtr + psSize);
	}
}


double MonoDescPWL::scalEval(double xx, const Dim& d)
{
	double yy;
	int index;

	switch (d)
	{
	case vert:
	case hori:
		index = binScalIndex(xx, d);

		if (index == -1)
			index = 0;
		else if (index == psSize - 1)
			index = psSize - 2;

	//	if (std::isnan(double(index)))
	//		index = psSize - 2;
		yy = interpolate(index, xx, d);
		return yy;
		break;
	default:
		std::string errorMsg = std::string(", could not support the given Dim ") +
			std::to_string(d);
		CATCHRTERROR(errorMsg);
		break;
	}
}

int MonoDescPWL::binScalIndex(const double& xx, const Dim& d)
{
	/**
	 * Dim::VERT //work on the 1st dimension (x axis)
	 * Dim::HORI //work on the 2nd dimension (y axis)
	 */
	int    idX;
//	double head, tail;
//	head = pointSet(0, d);
//	tail = pointSet(psSize - 1, d);

	if (d == hori)
	{
//		pointSet(0, d) = Inf;
//		pointSet(psSize - 1, d) = -Inf;
		idX = std::lower_bound(pointSetData(d), pointSetData(d) + psSize, xx, customDoubleUpper{}) -
			pointSetData(d) - 1;
	}
	else
	{
//		pointSet(0, d) = -Inf;
//		pointSet(psSize - 1, d) = Inf;
		idX = std::upper_bound(pointSetData(d), pointSetData(d) + psSize, xx) -
			pointSetData(d) - 1;
	}

//	pointSet(0, d) = head;
//	pointSet(psSize - 1, d) = tail;
	return idX;
}

double MonoDescPWL::interpolate(int begin, double value, Dim d)
{
	/*
		switch (d)
		{
		case Dim::VERT: //TODO: general handling, to be optimized
			if (pointSet.at(begin).first == value)
				return pointSet.at(begin).second;
			else
				return (value - pointSet.at(begin).first) * div
					+ pointSet.at(begin).second;
			break;

		case Dim::HORI:
			if (pointSet.at(begin).second == value)
				return pointSet.at(begin).first;
			else
				return (value - pointSet.at(begin).second) * (1/div)
					+ pointSet.at(begin).first;
			break;

		default:
			std::string errorMsg = std::string(", could not support the given Dim ") +
				std::to_string(static_cast<std::underlying_type<Dim>::type>(d));
			CATCHRTERROR(errorMsg);
		}
	*/
	double div;
	short  valMap;

	(d == vert) ? (valMap = 1) : (valMap = 0);
	if (pointSet(begin, d) == value)
		return pointSet(begin, valMap);
	else
	{
		if (d == vert) {
			div = DIFFDIV(begin, begin + 1);
		//	div = deriv(begin);
		}
		else {
			div = DIFFINDIV(begin, begin + 1);
		//	div = 1 / deriv(begin);
		}
		return (value - pointSet(begin, d)) * div +
			pointSet(begin, valMap);
	}
}


void MonoDescPWL::getMpp(std::vector<double> volbound_)
{
	/** volbound_ is targeted for the voltage: the range of voltage: [volbound_[0], volbound_[1]]
	 * The size of volbound_ is 2
	 * This is implemented under the condition that the vector of pointSet remains unchanged.
	 */
	std::vector<Point> mpPointSet(psSize);
	if (volbound_.empty())
	{
		volbound_.push_back(0.0);
		volbound_.push_back(Inf);
	}
	else
	{
		ASSERT((volbound_.size() == 2), "Bounds must be a 2-vector");
	}
	for (unsigned i = 0; i < volbound_.size(); i++)
		if (volbound_.at(i) < 0.0)
		{
			std::cerr << "WARNING getMpp: " <<
				"Pmp search is limited to 1Q (first quadrant), negative limits set to zero"
				<< std::endl;
			volbound_.at(i) = pvMAX(volbound_.at(i), 0.0);
		}
	//	nolimits = (bound_.at(0) == 0.0) && (bound_.at(1) == Inf);
	// Start with the set of break-points in the first quadrant

  

	unsigned j = 0;
	for (unsigned i = 0; i < psSize; i++)
	{
		if ((pointSet.row(i).array() > 0.0).all())
		{
			mpPointSet.at(j).first  = pointSet(i, 0);
			mpPointSet.at(j).second = pointSet(i, 1);
			++j;
		}
	}
	mpPointSet.resize(j);

	double xmp;

    currDeriv();// Added to calculate first order deriv.
	auto con = (pointSet.col(1).array() > pointSet.col(0).array() * deriv);
     
	for (unsigned i = 0; i < psSize - 1; i++)
	{
		if (con(i))
		{/// the Point outside the quarant 1
			xmp = (pointSet(i, 0) - pointSet(i, 1) / deriv(i)) / 2.0;

			/// First and last bounds are extended, in case MPP must be extrapolated.
			//TODO: alleviate the reallocation brought by push_back
			if (i == 0 || i == psSize - 2)
			{
				mpPointSet.push_back(Point(xmp, scalEval(xmp, Dim::vert)));
				continue;
			}
			if ((xmp > pointSet(i, 0)) && (xmp < pointSet(i + 1, 0)))
			{
				mpPointSet.push_back(Point(xmp, scalEval(xmp, Dim::vert)));
			}
		}
	}
     
	/// Make sure to include Isc and Voc.
	mpPointSet.push_back(Point(0.0, scalEval(0.0, Dim::vert))); /// add Ioc
	mpPointSet.push_back(Point(scalEval(0.0, Dim::hori), 0.0)); /// add Isc

	/// Make sure to include lower and upper bound points
	for (unsigned i = 0; i < volbound_.size(); i++)
	{
		if (std::abs(volbound_.at(i)) == Inf || !volbound_.at(i))
			continue;
		auto yy = scalEval(volbound_.at(i), Dim::vert);
		mpPointSet.push_back(Point(volbound_.at(i), yy)); //TODO: check not to put the repeated point to the pointset!!
	}
	//TODO: consider to optimize the following two max_element operations, to check if we can combine these two operations to one!!
	auto largest = std::max_element(mpPointSet.begin(), mpPointSet.end(), customDotProduct());

	mppInfo.pmp = ((*largest).first)*((*largest).second);
	mppInfo.vmp = (*largest).first;

	if (mppInfo.vmp >= volbound_.at(0) && mppInfo.vmp <= volbound_.at(1))
	{
		mppInfo.pmpBnd = mppInfo.pmp; /// maximum power point without bound
		mppInfo.vmpBnd = mppInfo.vmp;
	}
	else
	{
		if (!std::is_sorted(mpPointSet.begin(), mpPointSet.end(), customFirstLess()))
			std::sort(mpPointSet.begin(), mpPointSet.end(), customFirstLess());
		auto lower = mpPointSet.begin();
		auto upper = mpPointSet.end() - 1;
		if (mppInfo.vmp < volbound_.at(0))
			lower = std::lower_bound(mpPointSet.begin(), mpPointSet.end(),
				volbound_.at(0), customLowerFirst());
		if (mppInfo.vmp > volbound_.at(1))
			upper = std::lower_bound(mpPointSet.begin(), mpPointSet.end(),
				volbound_.at(1), customLowerFirst());
		auto largest = std::max_element(lower, upper + 1, customDotProduct());
		mppInfo.pmpBnd = ((*largest).first)*((*largest).second); /// maximum power point with bound
		mppInfo.vmpBnd = (*largest).first;
	}

#ifndef NDEBUG
	std::cout << "pmp, vmp, pmpBnd, vmpBnd are " << mppInfo.pmp << " , " << mppInfo.vmp <<
		" , " << mppInfo.pmpBnd << " , " << mppInfo.vmpBnd << std::endl;

	for (unsigned i = 0; i < mpPointSet.size(); i++)
	{
		std::cout << mpPointSet.at(i).first << ", " << mpPointSet.at(i).second << "\n";
	}
	std::cout << "\n";
#endif
}



#if 0

#ifdef VECVERSION
std::vector<int> MonoDescPWL::binIndex(std::vector<double> xx, Dim d)
{//TODO: put all those std::upper_bound and std::lower_bound functions into the mathfunctions!!
	//TODO: the parameter xx needs to be changed
	// xx could be any size and any set of points with one dimention
	std::vector<double> dimValue(n, 0);
	std::vector<int>    idX;

	switch (d)
	{
	case Dim::VERT: //work on the 1st dimention (x axis)
		dimValue.front() = -Inf;
		dimValue.back()  =  Inf;
		for (unsigned i = 1; i < n - 1; i++)
			dimValue.at(i) = refPointSet(i, 0);
		for (unsigned i = 0; i < xx.size(); i++)
			idX.push_back(std::upper_bound(dimValue.begin(), dimValue.end(),
				xx.at(i)) - dimValue.begin() - 1);
		break;
	case Dim::HORI: //work on the 2nd dimention (y axis)
		dimValue.front() = Inf;
		dimValue.back() = -Inf;
		for (unsigned i = 1; i < n - 1; i++)
			dimValue.at(i) = refPointSet(i, 1);
		for (unsigned i = 0; i < xx.size(); i++)
			idX.push_back(dimValue.rend() -
				std::lower_bound(dimValue.rbegin(), dimValue.rend(), xx.at(i)) - 1);
#ifndef NDEBUG
		std::cout << "the value is " << idX.front() << std::endl;
#endif
		break;
	default:
		std::string errorMsg = std::string(", could not support the given Dim ") +
			std::to_string(static_cast<std::underlying_type<Dim>::type>(d));
		CATCHRTERROR(errorMsg);
	}
	return idX;
}

std::vector<double> MonoDescPWL::eval(std::vector<double> xx, Dim d)
{// evaluate function value at points X(xx).
	std::vector<double> yy;

	yy.clear();
	if (xx.empty()) return yy;

	//TODO: check strict monotonicity??
	std::vector<int> indexSet = binIndex(xx, d);

	double elem;
	for (unsigned i = 0; i < indexSet.size(); i++)
	{
		if (std::isnan(double(indexSet.at(i))))
			indexSet.at(i) = n - 2;
		//TODO: if check with the -Inf or Inf, then issue errors

		elem = interpolate(indexSet.at(i), indexSet.at(i) + 1, xx.at(i), d);
		yy.push_back(elem);
	}
	return yy;
}
#endif
#endif

#ifdef MSXML
int main()
{
	GetSimOptFile* filePtr = GetSimOptFile::Instance();
#ifdef MSXML
	filePtr->openOptionFile("SimOption.xml");
#else
//	filePtr->openOptionFile("/zhome/academic/HLRS/hlrs/hpchuzho/PV/linuxPVProj/SimOption.xml");
	filePtr->openOptionFile("Resources/SimOption.xml");
#endif
	filePtr->readOptionFile();

	MatrixXd PS(ITRNUM,2);
//	MatrixXd P2(ITRNUM,2);
//	MatrixXd P3(ITRNUM, 2);
	// Dataset for douglas
	PS << 0.0, 0.0,
		1.0, 0.1,
		2.0, -0.1,
		3.0, 5.0,
		4.0, 6.0,
		5.0, 7.0,
		6.0, 8.1,
		7.0, 9.0,
		8.0, 9.0,
		9.0, 9.0;
  /*
	PS << 1.05, 1.8,
		1.08, 1.4,
		1.1, 1.25,
		3.1, 1.1;


	P2 << 1, 2.15,
		1.7, 1.95,
		2, 0.5,
		3, 0;*/


/*
	PS.push_back(Point(-1.2, 5.4));
	PS.push_back(Point(2.3, 4.3));
	PS.push_back(Point(2.5, 4.1));
	PS.push_back(Point(2.8, 3.2));
	PS.push_back(Point(3.2, 2.1));
	PS.push_back(Point(3.6, 1.6));
	PS.push_back(Point(4.3, 0.2));
	PS.push_back(Point(4.3, -1.2));


	PS.push_back(Point(-1.2, 5.4));
	PS.push_back(Point(2.3, 4.3));
	PS.push_back(Point(2.5, 4.1));
	PS.push_back(Point(3.2, 2.1));
	PS.push_back(Point(2.8, 3.2));
	PS.push_back(Point(3.6, 1.6));
	PS.push_back(Point(4.3, -1.2));
	PS.push_back(Point(4.3, 0.2));*/
#if 0
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);

	std::vector<double> PSFirst, PSSec;

	for (int iter = 0; iter < ITRNUM; iter++)
	{
		PSFirst.push_back(dis(gen));
		PSSec.push_back(dis(gen));
	}

	if (!std::is_sorted(PSFirst.begin(), PSFirst.end()))
		std::sort(PSFirst.begin(), PSFirst.end());

	if (!std::is_sorted(PSSec.begin(), PSSec.end(), customDoubleDesc()))
		std::sort(PSSec.begin(), PSSec.end(), customDoubleDesc());


	for (int iter = 0; iter < ITRNUM; iter++)
	{
		PS(iter, 0) = PSFirst.at(iter);
		PS(iter, 1) = PSSec.at(iter);
	}
#endif
#if 0
  PSFirst.clear();
	PSSec.clear();

	for (int iter = 0; iter < ITRNUM; iter++)
	{
		PSFirst.push_back(dis(gen));
		PSSec.push_back(dis(gen));
	}

	if (!std::is_sorted(PSFirst.begin(), PSFirst.end()))
		std::sort(PSFirst.begin(), PSFirst.end());

	if (!std::is_sorted(PSSec.begin(), PSSec.end(), customDoubleDesc()))
		std::sort(PSSec.begin(), PSSec.end(), customDoubleDesc());


	for (int iter = 0; iter < ITRNUM; iter++)
	{
		P2(iter, 0) = PSFirst.at(iter);
		P2(iter, 1) = PSSec.at(iter);
	}
#endif
#ifdef MSXML
	clock_t start = clock();
#else
    struct timeval sTime, eTime;
    gettimeofday(&sTime, NULL);
#endif

	MonoDescPWL pwl1(PS, {1.0});
	std::cout<<"pointset is \n"<<pwl1.pointSet<<"\n";
// 	MonoDescPWL pwl2(P2);

#ifndef NDEBUG
	std::cout << "pwl1: \n" << PS << "\n";
	std::cout << "pwl2: \n" << P2 << "\n";
#endif


//	MonoDescPWL pwl3(P3);
//	MonoDescPWL pwl3(pwl1, pwl2, 0.4);

//	std::vector<MonoDescPWL*> set;

//	set.push_back(&pwl1);
//	set.push_back(&pwl2);

//	MonoDescPWL pwl3 = MonoDescPWL(set, Dim::hori);

#ifndef NDEBUG
	std::cout << "pwl3 : \n" << pw3.pointSet << "\n";
#endif

#ifdef MSXML
	clock_t end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
#else
  	gettimeofday(&eTime, NULL);
  	long seconds = (eTime.tv_sec - sTime.tv_sec);
  	long time = ((seconds * 1000000) + eTime.tv_usec) - (sTime.tv_usec);
#endif

	std::cout << "\nThe time is " << time << " us\n";


#ifdef MSXML
	std::cin.get();
#endif
}
#endif
