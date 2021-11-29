#include "OneDiodeModel.h"
#include "../Utils/Compare.h"
#include "../Utils/LambertW.h"

#include <cstdlib>
#include <iomanip>
#include <random>
#ifdef MSXML
#include <time.h>
#else
#include <sys/time.h>
#endif
#include <vector>
#include <chrono>
#include <omp.h>
#include "../testUtils.h"

class OneDiodeModel;

OneDiodeModel::OneDiodeModel(struct Params params_,
	const std::vector<double>& tol_,
	const unsigned long& maxIter_)
	: inputParams(params_), tol(tol_), maxIter(maxIter_)
{
	filePtr = GetSimOptFile::Instance();
	if (maxIter == (unsigned)INITVAL)
		maxIter = (filePtr->getElement<unsigned>)("MaxIter");
	neps = (filePtr->getElement<double>)("NEPS");

	filePtr->parseTolerance(tol);
	/// Replace the "minRelTol" and "minAbsTol".
	filePtr->setElement<double>("minRelTol", EPSILON1);
	filePtr->setElement<double>("minAbsTol", EPSILON0);

	filePtr->parseTolerance(tolBS);

    /// Wrap up the repeatedly-calculated variables.
	prod = (inputParams.Iph*inputParams.di2Mutau);
	commonClause1 = inputParams.Io*inputParams.Rsh/inputParams.nVth;
	commonClause2 = inputParams.Iph*inputParams.Rsh*inputParams.di2Mutau;
	m0 = commonClause1 +
		commonClause2/(POW2(inputParams.Vbi)) + 1.0;
	m0 = -1.0*(m0 / (inputParams.Rsh + inputParams.Rs*m0));

	/// The point from which the module starts to be reverser-biased (V < -IRs).
	i0 = inputParams.Iph-prod/inputParams.Vbi;
	logIo = std::log(inputParams.Io);
}

void OneDiodeModel::fiveParameterODM()
{
//	std::cout.precision(std::numeric_limits<double>::max_digits10);
	ArrayXd newIph = inputParams.Iph-iRec+iRev;
	if (inputParams.Rs == 0.0)
		iSet = newIph - inputParams.Io*((vSet/inputParams.nVth).exp() - 1.0) 
				- vSet/inputParams.Rsh;
	else
	{
		double kr = inputParams.Rsh/(inputParams.Rs + inputParams.Rsh);
		double ki = inputParams.Rs/inputParams.nVth;

		ArrayXd logPar = std::log(kr*ki*inputParams.Io) 
				+ ki*kr*(newIph + inputParams.Io + vSet/inputParams.Rs);

		ArrayXd w(odmSize);
		/// Tolerance on the current direction.
		ArrayXd toli(odmSize);
		for (int i = 0; i < odmSize; i++)
		{
			/// Use the boost lambert function as the default policy. (see LambertW.h)
			w(i) = lambertW(logPar(i));
			//w(i) = boost::math::lambert_w0(std::exp(logPar(i)));
		}

		iSet = kr*(newIph + inputParams.Io - vSet/inputParams.Rsh) - w/ki;

		toli = filePtr->tolFunVec(tol, iSet, neps, Dim::hori);
		// TODO: probably Dim::hori can be replaced
		auto fLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd iS) -> ArrayXd 
			{return newIph - (logIo + 
				(vSet + iS*inputParams.Rs)/inputParams.nVth).exp() + inputParams.Io 
				- (vSet + iS*inputParams.Rs)/inputParams.Rsh - iS;}); // solve fLambda(x) = 0
		/// For NewtonRaphson.
	/*	auto dfLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd iS) -> ArrayXd
			{return -1*(inputParams.Rs/inputParams.Rsh + 
				inputParams.Rs/inputParams.nVth*(logIo + 
				(vSet + iS*inputParams.Rs)/inputParams.nVth).exp()) - 1.0;}); */ 
        /// Use the boost Newton-Rapson algorithm.
        if ((fLambda(iSet).abs() > toli).any())
        {
        //	NewtonRaphson nrMethod(fLambda, filePtr,
        //		tol, neps, dfLambda,
        //		maxIter, iSet);
        	BoostNewtonRaphson nrMethod(filePtr,
        		tol, neps, inputParams, newIph,
        		maxIter, iSet, vSet, Dim::vert);

        	setRootFind(&nrMethod);
        	findRoot();
        }

#ifndef NDEBUG
        std::cout<<"The output of fiveParameterODM is \n"<<iSet<<"\n";
#endif
	}
}

void OneDiodeModel::getDerivatives(enum Dim dim)
{
	// for Dim::hori -> Ix
	// TODO: ask Martin to change the naming
	ArrayXd kc = commonClause1*
		((vSet + inputParams.Rs*iSet)/inputParams.nVth).exp();
	ArrayXd kd = commonClause2 / 
		(vSet - inputParams.Vbi + inputParams.Rs*iSet).pow(2);
	ArrayXd ke = kc + kd + 1.0;

	ArrayXd clause = (inputParams.Rsh + inputParams.Rs*ke);
	d1Iset = -1.0*ke/clause;


	if(dim == Dim::vert)
		d2Iset = (1.0 + inputParams.Rs*d1Iset).pow(2)*
				(2.0*kd/(vSet - inputParams.Vbi + inputParams.Rs*iSet) - kc/inputParams.nVth)
				/ clause;

    int subS;
	for (int i = 0; i < indexIsRev.matrix().rows(); i++)
	{
		subS = indexIsRev(i);
		
		double kr = 2.0*inputParams.Brev*(vSet(subS) + inputParams.Rs*iSet(subS));
		d1Iset(subS) = (d1Iset(subS) + kr)/(1.0 - inputParams.Rs*kr);

		if (dim == Dim::vert)
			d2Iset(subS) = (d2Iset(subS) + 2.0*inputParams.Brev
						*POW2(1.0 - inputParams.Rs*d1Iset(subS)))/(1.0 - inputParams.Rs*kr);
	}
}

void OneDiodeModel::approxIrev()
{
	if (inputParams.Rs == 0.0)
		iRev = isRev.select(inputParams.Brev*(vSet.pow(2)), 0.0);
	else
	{
		iSet = i0 + m0*(vSet + inputParams.Rs*i0);

		ArrayXd term = 1.0 - 4.0*inputParams.Brev*inputParams.Rs*(vSet + inputParams.Rs*iSet);

		iRev = 1.0/(2.0*inputParams.Brev*POW2(inputParams.Rs))*(1.0 - term.sqrt()) - vSet/inputParams.Rs;
		iRev = isRev.select(iRev - iSet, 0.0);
	}
}


void OneDiodeModel::fiveParameterODMInv(const ArrayXi& indexS)
{
	int subSize = indexS.matrix().rows();
	ArrayXd tmpVSet = ArrayXd::Zero(subSize);
	ArrayXd iSubSet;
	igl::slice<ArrayXd, ArrayXd>(iSet, indexS, iSubSet);
	ArrayXd Idf = inputParams.Iph - iSubSet;
	
#pragma omp parallel for if (subSize > 1000)
	for (int i = 0; i < subSize; i++)
	{
		if (Idf(i) <= 0)
			tmpVSet(i) = inputParams.Iph - Idf(i)*inputParams.Rs;
		else
		{
			double logPar = std::log(commonClause1) + (commonClause1*Idf(i)/inputParams.Io);
			
			double w = lambertW(logPar);

		//	double w = boost::math::lambert_w0(std::exp(logPar));
			tmpVSet(i) = inputParams.Rsh*Idf(i) - inputParams.Rs*iSubSet(i) - inputParams.nVth*w;
		}
	}

	ArrayXd tolv(subSize); /// Tolerance for voltage direction.
	tolv = filePtr->tolFunVec(tol, tmpVSet, neps, Dim::vert);
	
	auto fLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd vS) -> ArrayXd 
		{return inputParams.Rsh*(inputParams.Iph - (logIo + 
			(vS + iSubSet*inputParams.Rs)/inputParams.nVth).exp() + inputParams.Io)
			- iSubSet*(inputParams.Rs + inputParams.Rsh) - vS;});
	/// Reserved for NewtonRaphson.
	//auto dfLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd vS) -> ArrayXd
	//	{return -1*(inputParams.Rsh/inputParams.nVth*logIo + 
	//		(vS + iSubSet*inputParams.Rs)/inputParams.nVth)).exp() - 1;});
    if ((fLambda(tmpVSet).abs() > tolv).any())
    {
    //	NewtonRaphson nrMethod(fLambda, filePtr,
    //		tol, neps, dfLambda,
    //		maxIter, tmpVSet);
    	/// Scalar processing using boost NewtonRaphson. (see RootFind.h)
    	BoostNewtonRaphson nrMethod(filePtr,
    		tol, neps, inputParams,
    		maxIter, tmpVSet, iSubSet, Dim::hori);
    	setRootFind(&nrMethod);
    	findRoot();
    }
    igl::slice_into(tmpVSet,indexS,vSet);
}

void OneDiodeModel::VoltoCurr(ArrayXd vSet_)
{
	vSet = vSet_;
	odmSize = vSet.matrix().rows();
	iSet.resize(odmSize);

	isRev = vSet < -1.0*(inputParams.Rs*i0);
	igl::find<ArrayXb, ArrayXi>(isRev, indexIsRev);

	/// Approximate the recombination loss current.
	iRec = prod/(inputParams.Vbi - vSet);

	iRev = ArrayXd::Zero(odmSize);
	if (inputParams.Brev > 0.0 && isRev.any())
		approxIrev();

//	double test_start, test_end;
//	TIMESTAMP(test_start);
	fiveParameterODM();

//    TIMESTAMP(test_end);
//    std::cout<<"test time is "<<test_end - test_start<<"\n";
	ArrayXb triky; //TODO: change the naming of triky
	triky = (inputParams.di2Mutau > 0.0) * (vSet > inputParams.Vbi);
    ArrayXi indexTricky;
    igl::find<ArrayXb, ArrayXi>(triky, indexTricky);
		
//  double bs_start, bs_end;
//  TIMESTAMP(bs_start);
	int trickySize = indexTricky.matrix().rows();

#pragma omp parallel for if (trickySize > 1000)
    for (int i = 0; i < trickySize; i++)
    {
		double root;
		int    subS = indexTricky(i);
		double vPoint = vSet(subS);

		auto errLambda = std::function<double(double)>([=](double iS) -> double
		{return inputParams.Iph - std::exp(logIo + (vPoint + iS*inputParams.Rs)/inputParams.nVth)
			+ inputParams.Io - (vPoint + iS*inputParams.Rs)/inputParams.Rsh
			- prod/pvMAX<double>(EPSILON0, inputParams.Vbi - vPoint - inputParams.Rs*iS) - iS;});
		/// Get the root by adopting the bisection method. (see RootFind.h)
		BiSection bsMethod(std::move(errLambda), filePtr, tolBS, neps, -vPoint/inputParams.Rs, i0, root);
			
		setRootFind(&bsMethod);
		findRoot();
#ifndef NDEBUG
		std::cout<<"the root of BiSection is "<<root<<"\n";
#endif
		iSet(subS) = root;
#ifndef NDEBUG
		std::cout<<"the iSet is \n"<<iSet <<"\n";
#endif
	}
//	TIMESTAMP(bs_end);
//	std::cout<<"bs time is "<<bs_end - bs_start<<"\n";

//  double iter_start, iter_end;
//  TIMESTAMP(iter_start);
	if ((iRev != 0.0).any() || (iRec != 0.0).any())
	{
		ArrayXd preIset;
		int i;
		for (i = 0; i < maxIter + 1; i++)
		{
			preIset = iSet;
			iRec = prod/((inputParams.Vbi - vSet - inputParams.Rs*iSet).max(EPSILON0));
			iRev = isRev.select(inputParams.Brev*((vSet + inputParams.Rs*iSet).pow(2)), 0.0);
			fiveParameterODM();

			if (((iSet - preIset).abs() <= filePtr->tolFunVec(tol, iSet, neps, Dim::hori)).all())
				break;
		}
#ifndef NDEBUG
		std::cout<<"the i is "<<i<<"\n";
#endif
		if (i > maxIter)
			CATCHRTERROR(", ODM could not converge to a solution");
	}

//	TIMESTAMP(iter_end);
//	std::cout<<"iter time is "<<iter_end - iter_start<<"\n";
//	getDerivatives(Dim::vert);
#ifndef NDEBUG
	std::cout<<"The output iSet is \n"<< iSet <<"\n";
#endif
	//TODO: ask Martin about the reshape??
}

void OneDiodeModel::CurrtoVol(ArrayXd iSet_)
{
	iSet = iSet_;
	odmSize = iSet.matrix().rows();

	isRev = (iSet >= i0);
	ArrayXb notIsRev = (iSet < i0);

	vSet = ArrayXd::Zero(odmSize);
	factor = ArrayXd::Ones(odmSize);

	igl::find<ArrayXb, ArrayXi>(isRev, indexIsRev);
	ArrayXi indexNotIsRev;
    igl::find<ArrayXb, ArrayXi>(notIsRev, indexNotIsRev);

	/// Get approx. reverse-bias voltages.
	if (inputParams.Brev > 0.0 && isRev.any())
	{
		// Solution to I=i0+m0(V+Rsi)+Brev(V+Rsi)^2.
		for (int i = 0; i < indexIsRev.matrix().rows(); i++)
		{
			int subIsRev = indexIsRev(i);
			
			vSet(subIsRev) = -inputParams.Rs*iSet(subIsRev) - 
				m0/(2.0*inputParams.Brev)*(1.0-std::sqrt(1.0 + 4.0*inputParams.Brev/m0*(iSet(subIsRev)-i0)*(1.0/m0+inputParams.Rs)));
		}
	}

//  double bs_start, bs_end;
//  TIMESTAMP(bs_start);
    if (inputParams.di2Mutau > 0.0 && notIsRev.any())
    {
    	int notIsRevSize = indexNotIsRev.matrix().rows();

#pragma omp parallel for if (notIsRevSize > 1000)
    	for (int i = 0; i < notIsRevSize; i++)
    	{
			double root;
			int    subS = indexNotIsRev(i);
			double iPoint = iSet(subS);

			auto errLambda = std::function<double(double)>([=](double vS) -> double
			{return inputParams.Iph - std::exp(logIo + (vS + iPoint*inputParams.Rs)/inputParams.nVth)
				+ inputParams.Io - (vS + iPoint*inputParams.Rs)/inputParams.Rsh
				- prod/pvMAX<double>(EPSILON0, inputParams.Vbi - vS - inputParams.Rs*iPoint) - iPoint;});
			BiSection bsMethod(std::move(errLambda), filePtr, tol, neps, -i0*inputParams.Rs, inputParams.Vbi-inputParams.Rs*iPoint, root);
				
			setRootFind(&bsMethod);
			findRoot();
	#ifndef NDEBUG
			std::cout<<"the root of BiSection is "<<root<<"\n";
	#endif
			vSet(subS) = root;
			factor(subS) = 0.2;
	#ifndef NDEBUG
			std::cout<<"the vSet is \n"<<vSet <<"\n";
	#endif
		}
    }
    else
    {
    	fiveParameterODMInv(indexNotIsRev);
    }

//	TIMESTAMP(bs_end);
//	std::cout<<"bs time is "<<bs_end - bs_start<<"\n";

//  double iter_start, iter_end;
//  TIMESTAMP(iter_start);

	if (isRev.any() || (inputParams.di2Mutau > 0.0))
	{
		ArrayXd preVset;
		ArrayXd preIset;
		ArrayXd err;
		int i;
        preIset = iSet;
		for (i = 0; i < maxIter + 1; i++)
		{
			preVset = vSet;
			
			iRec = prod/(inputParams.Vbi - vSet - inputParams.Rs*preIset);
			iRev = isRev.select(inputParams.Brev*((vSet + inputParams.Rs*preIset).pow(2)), 0.0);
			fiveParameterODM();
			err = iSet - preIset;
			getDerivatives(Dim::hori); /// Get 1st derivative based on the updated iSet.
			vSet = vSet - factor*err/d1Iset;

			if ((((err.abs() <= filePtr->tolFunVec(tol, iSet, neps, Dim::hori)).cast<int>() + 
				((vSet-preVset).abs() <= filePtr->tolFunVec(tol, vSet, neps, Dim::vert)).cast<int>()) > 0).all())
				break;
		}
		iSet = preIset; /// Restore the value of iSet.
#ifndef NDEBUG
		std::cout<<"the i is "<<i<<"\n";
#endif
		if (i > maxIter)
			CATCHRTERROR(", ODMInv could not converge to a solution");
	}

//	TIMESTAMP(iter_end);
//	std::cout<<"iter time is "<<iter_end - iter_start<<"\n";
//	getDerivatives();
#ifndef NDEBUG
	std::cout<<"The output vSet is \n"<< vSet <<"\n";
#endif
}

void OneDiodeModel::VocIsc()
{
	ArrayXd givenSet(1);
	givenSet << 0.0;

	VoltoCurr(givenSet);
	Isc = iSet(0);
	CurrtoVol(givenSet);
	Voc = vSet(0);
}

void OneDiodeModel::getODMMpp(const double& seedVmp, const int& TolRefFac)
{
	// TODO: if need to restore the value of tol??
	unsigned tolSize = tol.size();
    for (int i = 0; i < tolSize; i++)
    	tol[i] /= TolRefFac;
	filePtr->parseTolerance(tol);
	tolSize = tol.size();
    for (int i = 0; i < tolSize; i++)
		tol[i] *= TolRefFac;

	VocIsc();

	if (seedVmp > 0.0 && seedVmp < Voc)
		Vmp = seedVmp;
	else
		Vmp = Voc*0.84;

	ArrayXd vmpSet(1);
	vmpSet << Vmp;

	VoltoCurr(vmpSet);
	getDerivatives(Dim::vert);

	Imp = iSet(0);

	double i, k;

	for (i = 0; i < maxIter; i++)
	{
		double oldVmp = Vmp;
		
		if (i < 2)
		{
			c = Vmp*d2Iset(0)/d1Iset(0);
			k = (1.0 - d2Iset(0)*Imp/POW2(d1Iset(0)));
		//	Vmp = (boost::math::lambert_w0(std::exp(std::log(k) + c + 1.0)) - 1.0)*Vmp/c;
			Vmp = (lambertW(std::log(k) + c + 1.0) - 1.0)*Vmp/c;
		}
		else
			Vmp = Vmp - (Vmp*d1Iset(0) + Imp)/(Vmp*d2Iset(0) + 2.0*d1Iset(0));

		vmpSet << Vmp;
		VoltoCurr(vmpSet);
		getDerivatives(Dim::vert);
		Imp = iSet(0);

		if ((ABS(oldVmp - Vmp) <= filePtr->tolFun(tol, Vmp, neps, Dim::vert)) && (ABS((oldVmp - Vmp)*d1Iset(0)) <= filePtr -> tolFun(tol, Imp, neps, Dim::hori)))
			break;
	}
	if (i >= maxIter)
		CATCHRTERROR(", ODMMPP could not converge to a solution");

    Pmp = Vmp * Imp;
    c = Vmp*d2Iset(0)/d1Iset(0);
}

void  OneDiodeModel::pwlApprox(const std::vector<double>& tol_, const std::vector<double>& lim, const double& biasK)
{
	ArrayXd limV(2);
	ArrayXd limI(2);

	tol = tol_;
	/// Adjust the tol.
	if (tol.empty())
	{
		tol = { 0.0 };
		getODMMpp();
		tol = tol_;
		filePtr->parseTolerance(tol); /// tol is passed as reference.
		tol[0] = Vmp*tol[2];
		tol[1] = Imp*tol[2];
	}
	else if (lim.empty())
	{
		VocIsc();
	}
	
	filePtr->parseTolerance(tol);

	int limSize = lim.size();
	switch (limSize)
	{
		case 0:
			limV[0] = 0.0;
			limV[1] = Voc;
			limI[0] = Isc;
			limI[1] = 0.0;
			break;
		case 2:
			if (lim[0] <= lim[1])
			{
				limV(0) = lim[0];
				limV(1) = lim[1];
			}else
			{
				limV(0) = lim[1];
				limV(1) = lim[0];
			}
			if (limV.isFinite().all())
				VoltoCurr(limV);
			else
				CATCHRTERROR(", Bad/empty limits");
			limI = iSet;
			break;
		case 4:
		{
			/// Guarantee that limV is in an ascending order.
			if (lim[0] <= lim[1])
			{
				limV(0) = lim[0];
				limV(1) = lim[1];
			}else
			{
				limV(0) = lim[1];
				limV(1) = lim[0];
			}
			/// Guarantee that limI is in a descending order.
			if (lim[2] <= lim[3])
			{
				limI(0) = lim[3];
				limI(1) = lim[2];
			}else
			{
				limI(0) = lim[2];
				limI(1) = lim[3];
			}

			CurrtoVol(limI);

			/// Special handling in the case of Nan or infinite.
			ArrayXb vSetisNan = vSet.isFinite();

			if (vSetisNan(0))
			{
				limV(0) = pvMAX(limV(0), vSet(0));
			}
		    if (vSetisNan(1))
		    {	
				limV(1) = pvMIN(limV(1), vSet(1));
			}

			if (limV.isFinite().all())
				VoltoCurr(limV);
			else
				CATCHRTERROR(", Bad/empty limits");
			limI = iSet;
		}	
		break;
		default:
			limV(0) = 0.0;
			limV(1) = 0.0;
	}

	if (limV(1) < limV(0))
		CATCHRTERROR(", Bad/empty limits");

	/// Differentiation between forward-biased (>0) and reverse-biased (<0).
	ArrayXd limZ = limV + inputParams.Rs*limI;

	/// Forward-biased.
	auto fbLambda = std::function<double(double)>([=](double z) -> double
		{return inputParams.Iph - std::exp(logIo + z/inputParams.nVth)
			+ inputParams.Io - z/inputParams.Rsh
			- prod/pvMAX<double>(EPSILON0, inputParams.Vbi - z);});

	/// Reverse-biased.
	auto rbLambda = std::function<double(double)>([=](double z) -> double
		{
			double tmp = -(inputParams.Io/inputParams.nVth + inputParams.Iph*inputParams.di2Mutau/POW2(inputParams.Vbi) + 1.0/inputParams.Rsh);
			return inputParams.Iph*(1.0 - inputParams.di2Mutau/inputParams.Vbi)
			+ tmp*z + inputParams.Brev*POW2(z);});

	if (limZ(0) < 0.0 && limZ(1) > 0.0)
	{
		// Interval contains Z = 0.0
		MatrixXd ps0, ps1;
		ArrayXd  limZ0(2), limI0(2), limZ1(2), limI1(2);
		limZ0(0) = limZ(0);
		limZ0(1) = 0.0;
		limZ1(0) = 0.0;
		limZ1(1) = limZ(1);

        /// Use the linear approximation algorithm defined in class StaticLA. (see Linearapprox.h)
		StaticLA sLA0(&ps0, filePtr, tol, neps, limZ0, biasK, rbLambda);
		setLinearApprox(&sLA0);
		approxLine();

		StaticLA sLA1(&ps1, filePtr, tol, neps, limZ1, biasK, fbLambda);
		setLinearApprox(&sLA1);
		approxLine();

		int size0 = ps0.rows();
		int size1 = ps1.rows();
		
		approxPointSet.resize(size0 + size1 - 1, 2);
		approxPointSet.topRows(size0 - 1) = ps0.topRows(size0 - 1);
		approxPointSet.bottomRows(size1) = ps1;
	}
	else
	{
		StaticLA* psLA;
		if (limZ(0) >= 0.0)
		{ 
			psLA = new StaticLA(&approxPointSet, filePtr, tol, neps, limZ, biasK, fbLambda); 
		}
		else if (limZ(1) <= 0.0)
		{ 
			psLA = new StaticLA(&approxPointSet, filePtr, tol, neps, limZ, biasK, rbLambda); 
		}
		setLinearApprox(psLA);
		approxLine();
		delete psLA;
	}
	approxPointSet.col(0) -= approxPointSet.col(1)*inputParams.Rs;
	
#ifndef NDEBUG
	std::cout<<"the size is "<<approxPointSet.rows()<<"\n";
	std::cout<<"approxPointSet is \n"<<approxPointSet<<"\n";
#endif
}


#if 0
//void OneDiodeModel::BiApprox(MatrixXd* out, ArrayXd limZ, ArrayXd limI)
void OneDiodeModel::BiApprox(double lx, double rx, double ly, double ry, unsigned recurCount)
{
	// TODO: temporarily only consider the forward biased condition and g is 0.5

    double g = 0.5;
	double mx = (rx - lx)*g + lx;
	double my = fbLambda(mx);

	double dy = my - ((mx - lx)*(ry - ly)/(rx - lx) + ly);
	double dx = mx - ((my - ly)*(rx - lx)/(ry - ly) + lx);

	if ((ABS(rx - lx) < 2*EPSILONScal(mx)) || (ABS(ry - ly) < 2*EPSILONScal(my)) 
		|| ((ABS(dy) < filePtr->tolFun(tol, my, neps, Dim::hori)) 
		&& (ABS(dx) < filePtr->tolFun(tol, mx, neps, Dim::vert))))
	{
		approxPointSet(pos, 0) = rx;
		approxPointSet(pos, 1) = ry;
        ++pos;
	}
	else
	{

		if (recurCount > maxRecur)
			std::cerr << "Warning: OneDiodeModel, exceed the maximum recursive depth = "
			<< recurCount << "\n";
		BiApprox(lx, mx, ly, my, recurCount + 1);
		BiApprox(mx, rx, my, ry, recurCount + 1);
	}
}

//void OneDiodeModel::BiApprox(MatrixXd* out, ArrayXd limZ, ArrayXd limI)
void OneDiodeModel::BiApprox(double lx, double rx, double ly, double ry, MatrixXd* out, unsigned recurCount)
{
	// TODO: temporarily only consider the forward biased condition and g is 0.5

    double g = 0.5;
	double mx = (rx - lx)*g + lx;
	double my = fbLambda(mx);

	double dy = my - ((mx - lx)*(ry - ly)/(rx - lx) + ly);
	double dx = mx - ((my - ly)*(rx - lx)/(ry - ly) + lx);

	if ((ABS(rx - lx) < 2*EPSILONScal(mx)) || (ABS(ry - ly) < 2*EPSILONScal(my)) 
		|| ((ABS(dy) < filePtr->tolFun(tol, my, neps, Dim::hori)) 
		&& (ABS(dx) < filePtr->tolFun(tol, mx, neps, Dim::vert))))
	{
		(*out).resize(2,2);
		(*out)<<lx, ly,
		        rx, ry;
	}
	else
	{
		MatrixXd recRes1, recRes2;

		if (recurCount > maxRecur)
			std::cerr << "Warning: OneDiodeModel, exceed the maximum recursive depth = "
			<< recurCount << "\n";
		BiApprox(lx, mx, ly, my, &recRes1, recurCount + 1);
		BiApprox(mx, rx, my, ry, &recRes2, recurCount + 1);
		int size1 = recRes1.rows();
		int size2 = recRes2.rows();
		(*out).resize(size1+size2-1, 2);
		(*out).topRows(size1) = recRes1;
		(*out).bottomRows(size2-1) = recRes2.bottomRows(size2-1);
	}
}
#endif