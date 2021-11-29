#include "ODMMpp.h"
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
#include "../Utils/LambertW.h"

class ODMMpp;

ODMMpp::ODMMpp(struct Params params_,
	const std::vector<double>& tol_,
	const double& seedVmp_,
	const unsigned long& maxIter_)
	: inputParams(params_), tol(tol_), seedVmp(seedVmp_), maxIter(maxIter_)
{
	filePtr = GetSimOptFile::Instance();
	if (maxIter == (unsigned)INITVAL)
		maxIter = (filePtr->getElement<unsigned>)("MaxIter");

	neps    = (filePtr->getElement<double>)("NEPS");
    
    unsigned tolSize = tol.size();
    for (int i = 0; i < tolSize; i++)
    	tol[i] /= TolRefFac;
	filePtr->parseTolerance(tol);
	tolSize = tol.size();
    for (int i = 0; i < tolSize; i++)
		tol[i] *= TolRefFac;

	objODMHandle = std::make_shared<OneDiodeModel>(inputParams, tol, maxIter);
    
//	objODMHandle = std::move(OneDiodeModel(inputParams, tol, maxIter)); 
	ArrayXd givenSet(1);
	givenSet << 0.0;

	objODMHandle->VoltoCurr(givenSet);
	Isc = (objODMHandle->getIset())(0);
	objODMHandle->CurrtoVol(givenSet);
	Voc = (objODMHandle->getVset())(0);

	if (seedVmp > 0.0 && seedVmp < Voc)
		Vmp = seedVmp;
	else
		Vmp = Voc*0.84;
//	std::cout<<"Vmp is "<<Vmp<<" Isc is "<<Isc<<" Voc is "<<Voc<<"\n";
}

void ODMMpp::getODMMpp()
{
	ArrayXd di, d2i;
	ArrayXd vmpSet(1);
	vmpSet << Vmp;
	
	objODMHandle->VoltoCurr(vmpSet);
	objODMHandle->getDerivatives(Dim::vert);

	Imp = (objODMHandle->getIset())(0);
	di = objODMHandle->getD1Iset();
	d2i = objODMHandle->getD2Iset();
	double i, k;

	for (i = 0; i < maxIter; i++)
	{
		double oldVmp = Vmp;
		
		if (i < 2)
		{
			c = Vmp*d2i(0)/di(0);
			k = (1.0 - d2i(0)*Imp/POW2(di(0)));
		//	Vmp = (boost::math::lambert_w0(std::exp(std::log(k) + c + 1.0)) - 1.0)*Vmp/c;
			Vmp = (lambertW(std::log(k) + c + 1.0) - 1.0)*Vmp/c;
		}
		else
			Vmp = Vmp - (Vmp*di(0) + Imp)/(Vmp*d2i(0) + 2.0*di(0));

		vmpSet << Vmp;
		objODMHandle->VoltoCurr(vmpSet);
		objODMHandle->getDerivatives(Dim::vert);
		Imp = (objODMHandle->getIset())(0);
		di = objODMHandle->getD1Iset();
		d2i = objODMHandle->getD2Iset();

		if ((ABS(oldVmp - Vmp) <= filePtr->tolFun(tol, Vmp, neps, Dim::vert)) && (ABS((oldVmp - Vmp)*di(0)) <= filePtr -> tolFun(tol, Imp, neps, Dim::hori)))
			break;
	}
	if (i >= maxIter)
		CATCHRTERROR(", ODMMPP could not converge to a solution");

    Pmp = Vmp * Imp;
    c = Vmp*d2i(0)/di(0);
}
