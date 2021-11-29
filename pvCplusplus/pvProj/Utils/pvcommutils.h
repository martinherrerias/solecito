/**
 * Commonly-used utilities
 */
#ifndef PV_COMMUTILS_H_
#define PV_COMMUTILS_H_

#include <cmath>
#include <stdexcept>  
#include <cassert>
#include <numeric>
#include <algorithm>

#ifndef NMEX
#include "mex.h"
#endif

/// The structure of the input parameters for the OneDiodeModel.h.
struct Params
{
	double Iph;
	double Io;
	double Rs;
	double nVth;
	double Rsh;
	double di2Mutau;
	double Vbi;
	double Brev;
};

enum Dim
{
	vert = 0, /// fixvertical (voltage), addparallel, evaluate the current
	hori, /// fixhorizontal (current), addseries, evaluate the voltage
	dimNum
};

/// Constants
constexpr double INITVAL = 0.0;
constexpr int ONENUM = 1;
constexpr double UNDEFINED = -1.0; /// Beyond the given curves.
//constexpr auto TOLERANCE = 10.0;
//constexpr auto EPSILON = TOLERANCE * std::numeric_limits<double>::epsilon();
constexpr auto Inf = std::numeric_limits<double>::infinity();
constexpr auto Nan = std::numeric_limits<double>::quiet_NaN();

const auto EPSILON1 = std::nextafter(1.0, 2.0) - 1.0; /// Equals to eps(1) in Matlab.
const auto EPSILON0 = std::nextafter(0.0, 1.0); /// Equals to eps(0) in Matlab.

template<class Number>
inline Number EPSILONVec(const Number &x, int size) 
{ Number rev(size);
for (int i = 0; i < size; i++)
	rev(i) = std::nextafter(x(i), x(i)+1.0) - x(i);
return std::move(rev);}

template<class Number>
inline Number EPSILONScal(const Number x) 
{ return std::nextafter(x, x+1.0) - x; }

/// \brief Template class and macro for the math expression.
template<class Number>
inline Number ABS(const Number x) { return std::abs(x);  }
template<class Number> 
inline Number pvMIN(const Number x, const Number y) { return std::min(x, y); }
template<class Number>
inline Number pvMAX(const Number x, const Number y) { return std::max(x, y); }
template<class Number>
inline Number MIN3NUM(const Number x, const Number y, const Number z) { return std::min(std::min(x, y), z); }
template<class Number>
inline Number MAX3NUM(const Number x, const Number y, const Number z) { return std::max(std::max(x, x), z); }

template < class Real >
inline Real POW2(const Real x) { return (x*x); }
template <class Real>
inline Real SQRT2POW(const Real x, const Real y) { return std::sqrt(POW2(x) + POW2(y));  }

/// \brief Returns the sign of the number x,
/// i.e., negative: -1, positive: 1, and zero: 0.
template <class Number>
inline Number SIGN(const Number x) {
	return (x > 0) - (x < 0);
}

#define DX(begin, end) (pointSet(end, 0) - pointSet(begin, 0))
#define DY(begin, end) (pointSet(end, 1) - pointSet(begin, 1))
#define DIFFDIV(begin, end) (DY(begin, end)/DX(begin, end))
#define DIFFINDIV(begin, end) (DX(begin, end)/DY(begin, end))
#define SIGN(x) ((x > 0).cast<double>() - (x < 0).cast<double>())

#define STDDEV(arr) (std::sqrt((((arr - arr.mean()).square()).sum())/(psSize - 1)))
#define pointSetData(dim) (pointSet.col(dim).data())


/// \brief Assertion expression, which works only when the NDEBUG is disabled.
#ifndef NMEX
#define ASSERT(condition, message) do { mxAssert(condition, message); } while (false)
#else
#define ASSERT(condition, message) do { assert((condition) && (message)); } while (false)
#endif
#define TOLSIZEASSERT do { ASSERT((t.size() == 3), "The size of the tolerance vector should be three"); } while(false)

#define CATCHRTERROR(message) \
do { \
	try{ \
		throw std::runtime_error(std::string("Error: on function ") + __FUNCTION__ + message); \
    }catch(std::runtime_error& e) {\
		std::cerr << e.what() << std::endl; \
		abort();  } } while(false)

#define ASSIGN(val) \
do { \
	if (j >= length){ \
		mpPointSet.push_back(Point(val, scalEval(val, Dim::vert))); \
	} else { \
		mpPointSet.at(j).first = val; \
		mpPointSet.at(j).second = scalEval(val, Dim::vert); \
		++j; \
	} } while(false)

#define ASSIGNSC(val) \
do { \
	if (j >= length){ \
		mpPointSet.push_back(Point(scalEval(val, Dim::hori), val)); \
	} else { \
		mpPointSet.at(j).first = scalEval(val, Dim::hori); \
		mpPointSet.at(j).second = val; \
		++j; \
	} } while(false)
#endif
