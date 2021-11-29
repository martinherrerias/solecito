/**
 * Benchmark on adding two PWL curves in parallel.
 *
 * @see benchOnePWLCurves.cpp
 */
#include "../../testUtils.h"
#include "../pvProj/ElectricalCalculation/MonoDescPWL.h"
#include "../pvProj/Utils/Compare.h"
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>

#undef NDEBUG

class MonoDescPWL;

/// The number of input points
constexpr int ITRNUM = 1000000;

int main()
{
	GetSimOptFile* filePtr = GetSimOptFile::Instance();
#ifdef MSXML
	filePtr->o

#undef NDEBUGpenOptionFile("SimOption.xml");
#else
	filePtr->openOptionFile("../../../Resources/SimOption.xml");
#endif
	filePtr->readOptionFile();

	MatrixXd PS(ITRNUM,2);
    MatrixXd P2(ITRNUM,2);

    std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0.0, 1.0);

	std::vector<double> PSFirst, PSSec, P2First, P2Sec;

	for (int iter = 0; iter < ITRNUM; iter++)
	{
		PSFirst.push_back(dis(gen));
		PSSec.push_back(dis(gen));

    P2First.push_back(dis(gen));
    P2Sec.push_back(dis(gen));
	}

	if (!std::is_sorted(PSFirst.begin(), PSFirst.end()))
		std::sort(PSFirst.begin(), PSFirst.end());

    if (!std::is_sorted(P2First.begin(), P2First.end()))
    	std::sort(P2First.begin(), P2First.end());

	if (!std::is_sorted(PSSec.begin(), PSSec.end(), customDoubleDesc()))
		std::sort(PSSec.begin(), PSSec.end(), customDoubleDesc());

  	if (!std::is_sorted(P2Sec.begin(), PSSec.end(), customDoubleDesc()))
    	std::sort(P2Sec.begin(), P2Sec.end(), customDoubleDesc());


	for (int iter = 0; iter < ITRNUM; iter++)
	{
		PS(iter, 0) = PSFirst.at(iter);
		PS(iter, 1) = PSSec.at(iter);

    	P2(iter, 0) = P2First.at(iter);
    	P2(iter, 1) = P2Sec.at(iter);
	}


  	double sTime, eTime;

  	std::cout << "Adding two PWL IV-curves with "<< ITRNUM <<" points in parallel...\n";
  	TIMESTAMP( sTime );
	MonoDescPWL pwl1(PS);
  	MonoDescPWL pwl2(P2);

  	std::vector<MonoDescPWL*> set;
  	set.push_back(&pwl1);
  	set.push_back(&pwl2);

  	MonoDescPWL pwl3 = MonoDescPWL(set, Dim::vert);
  	TIMESTAMP( eTime );

  	std::cout << "The elasped time is " << eTime-sTime << " s\n";
}
