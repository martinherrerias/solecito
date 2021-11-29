/**
 * Collection of the root-finding algorithms.
 *
 * Bisection: use the boost bisection method with a given interval.
 * 
 * BoostNewtonRaphson: use the boost boost Newton Raphson method (scalar processing).
 * NewtonRaphson: vector processing.
 */

#ifndef ROOTFIND_H_
#define ROOTFIND_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <random>
#include "xmlutils.h"
#include "pvcommutils.h"
#include <boost/math/tools/roots.hpp>
#include <omp.h>

using namespace Eigen;

using namespace boost::math::tools;

/** The usage of policy in boost
*using boost::math::policies::error_policy_type;
using boost::math::policies::policy;
using boost::math::policies::evaluation_error;

typedef policy<
	evaluation_error<error_policy_type::ignore_error>
> my_policy;*/

class GeneralRootFind
{
protected:
	GetSimOptFile* filePtr;
	std::vector<double> tol;
	double neps;
public:
	GeneralRootFind(GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_)
	: filePtr(filePtr_), tol(tol_), neps(neps_) {};
	virtual void findRoot() = 0;
};

class BiSection : public GeneralRootFind
{
private:
	std::function<double(double)> f;
	double leftInterval, rightInterval;
	double& root;
public:

	class tolerance {
	public:
		tolerance(double eps_) :
			eps(eps_) {}
		bool operator()(double first, double second) {
			return (fabs(second - first) <= eps);
		}
	private:
		double eps;
	};

	BiSection(const std::function<double(double)> f_, GetSimOptFile* filePtr_,
		std::vector<double> tol_, double neps_, double li, double ri, double& root_) 
	: GeneralRootFind(filePtr_, tol_, neps_), f(f_), leftInterval(li), rightInterval(ri), root(root_) {
	} //TODO: should change the constructor!! delete the std::function<double(double)> f_

	virtual void findRoot()
	{ 	
		std::pair<double, double> found;
		double midPoint; 
	    
		if (leftInterval == rightInterval)
		{
			std::cerr<<"WARNING BiSection: "
						<<"Closed (zero-distance) intervals included\n";
			root = leftInterval;
			return;
		}
		
		midPoint = (leftInterval + rightInterval)/2.0;
		root = leftInterval;
		if (f(root) > f(midPoint))
			root = midPoint;
		if (f(root) > f(rightInterval))
			root = rightInterval;

		if (f(leftInterval)*f(rightInterval) > 0.0)
		{/// First revision
			if (f(midPoint)*f(leftInterval) <= 0.0)
				{rightInterval = midPoint;}
			else if (f(midPoint)*f(rightInterval) <= 0.0)
				{leftInterval = midPoint;}
			else {
				std::cerr<<"WARNING BiSection: "
					<<"The root is not included in provided interval\n";
				return;
			}
		}
	    
	    /// The returned 'found' consists of two consecutive roots.
	    /// The terminiation condition is |found.first1 - found.second| < toli.
	    tolerance toli(filePtr->tolFun(tol, root, neps, Dim::vert));
	    /// Note: rightInterval > leftInterval, otherwise raise error.
	    ASSERT((rightInterval >= leftInterval), "The rightInterval should not be smaller than leftInterval");
		found =
			bisect(f, leftInterval, rightInterval, toli);
		root = (found.first + found.second)/2.0;
	}
};

class BoostNewtonRaphson : public GeneralRootFind
{
	template <class T>
	struct curr_functor
	{ 
  		curr_functor(T const& iPh_, T const& vP_, struct Params input_) : iPh(iPh_),vP(vP_),inputParams(input_)
  		{ 
  		}
  		std::pair<T, T> operator()(T const& x)
  		{
    		T fx = iPh - std::exp(std::log(inputParams.Io) + 
				(vP + x*inputParams.Rs)/inputParams.nVth) + inputParams.Io 
				- (vP + x*inputParams.Rs)/inputParams.Rsh - x;
			T dx = -1.0*inputParams.Rs/inputParams.Rsh - inputParams.Rs/inputParams.nVth*
					std::exp(std::log(inputParams.Io)+(vP + x*inputParams.Rs)/inputParams.nVth) - 1.0;
    		return std::make_pair(fx, dx); /// 'return' both fx and dx.
  		}
		private:
  		T iPh, vP;
  		struct Params inputParams;
  	}; /// curr_functor

  	template<class T>
  	struct vol_functor
	{ 
  		vol_functor(T const& iP_, struct Params input_) : iP(iP_),inputParams(input_)
  		{ 
  		}
  		std::pair<T, T> operator()(T const& x)
  		{
    		T fx = inputParams.Rsh*(inputParams.Iph - std::exp(std::log(inputParams.Io) + 
			(x + iP*inputParams.Rs)/inputParams.nVth) + inputParams.Io)
			- iP*(inputParams.Rs + inputParams.Rsh) - x;;
			T dx = -1.0*inputParams.Rsh/inputParams.nVth*std::exp((std::log(inputParams.Io) + 
			(x + iP*inputParams.Rs)/inputParams.nVth)) - 1.0;
    		return std::make_pair(fx, dx); /// 'return' both fx and dx.
  		}
		private:
  		T iP;
  		struct Params inputParams;
  	}; // vol_functor

private:
	unsigned maxIter;
	Ref<ArrayXd> seed;
	ArrayXd givenSet;
	ArrayXd newIph;
	struct Params inputP;
	enum Dim dim;
public:
	BoostNewtonRaphson(GetSimOptFile* filePtr_,
		std::vector<double> tol_, double neps_, const struct Params inputP_,
		ArrayXd newIph_,
		unsigned maxIter_, Ref<ArrayXd> seed_, ArrayXd givenSet_, enum Dim dim_)
	: GeneralRootFind(filePtr_, tol_, neps_), maxIter(maxIter_), seed(seed_), givenSet(givenSet_), inputP(inputP_), newIph(newIph_), dim(dim_) {}

	BoostNewtonRaphson(GetSimOptFile* filePtr_,
		std::vector<double> tol_, double neps_, const struct Params inputP_,
		unsigned maxIter_, Ref<ArrayXd> seed_, ArrayXd givenSet_, enum Dim dim_)
	: GeneralRootFind(filePtr_, tol_, neps_), maxIter(maxIter_), seed(seed_), givenSet(givenSet_), inputP(inputP_), dim(dim_) {}

	virtual void findRoot()
	{
		boost::uintmax_t mIter;
		int digits = std::numeric_limits<double>::digits; /// Maximum possible binary digits accuracy for type double.
  		int get_digits = (digits * 3)/4; /// Near maximum (3/4) possible accuracy.
  		int seedSize = seed.matrix().rows();
  		
  		switch (dim)
  		{
  			case Dim::vert:
#pragma omp parallel for default(shared) private(mIter) if (seedSize > 3000)
  			for(int i = 0; i < seedSize; i++)
    		{
    			mIter = maxIter;
    			seed(i) = 
    				newton_raphson_iterate(curr_functor<double>(newIph(i), givenSet(i), inputP), 
    									seed(i), -Inf, Inf, get_digits, mIter);
    		}
    		break;
    		case Dim::hori:
#pragma omp parallel for default(shared) private(mIter) if (seedSize > 3000)  		
    		for(int i = 0; i < seedSize; i++)
    		{
    			mIter = maxIter;
    			seed(i) = 
    				newton_raphson_iterate(vol_functor<double>(givenSet(i), inputP), 
    									seed(i), -Inf, Inf, get_digits, mIter);
    		}
    		break;
    		default: break;
  		}
  		
    }
};

class NewtonRaphson : public GeneralRootFind
{
private:
	unsigned long maxIter;
	std::function<ArrayXd(ArrayXd)> f;
	std::function<ArrayXd(ArrayXd)> df;
	Ref<ArrayXd> seed;
public:
	NewtonRaphson(const std::function<ArrayXd(ArrayXd)> f_, GetSimOptFile* filePtr_,
		std::vector<double> tol_, double neps_, const std::function<ArrayXd(ArrayXd)> df_,
		unsigned maxIter_, Ref<ArrayXd> seed_)
	: GeneralRootFind(filePtr_, tol_, neps_), f(f_), df(df_), maxIter(maxIter_), seed(seed_) {}

	virtual void findRoot()
	{
		/// Set random device
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0.0, 1.0);

		size_t seedSize;
		ArrayXd xp, fXp, dXp, fx, dx;
		ArrayXd comTmp;
		unsigned iter = 0;
		std::string convergeErrMsg = std::string(", NewtonRaphson could not converge to a solution");
		// TODO: generize this error message

		seedSize = seed.matrix().rows();
		ArrayXd toli(seedSize);
	    Ref<ArrayXd> x = seed; /// iSet is the initial seed root
		do
		{/// xp: previous x; x: update x
			if (iter > maxIter)
			{
				CATCHRTERROR(convergeErrMsg);
			}
			xp = x;
			if (iter > 0)
			{
				fXp = fx;
				dx = -1.0*fx/df(x);
				if ((dx + dXp <= toli).any())
					dx *= dis(gen);
				dXp = dx;
			}
			else{
				fXp = f(xp);
				dXp = -1.0*fXp/df(xp);
			}
			x = xp + dXp;
			fx = f(x);

			comTmp = ((fx.abs()).min((x-xp).abs())).min((fx-fXp).abs());

			toli = filePtr->tolFunVec(tol, x, neps, Dim::hori); //TODO: generize it
			++iter;
		}while((comTmp > toli).any());
	}
};

#endif