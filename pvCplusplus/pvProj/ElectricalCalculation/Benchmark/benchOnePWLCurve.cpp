/**
 * Benchmark on generating a PWL curve in terms of the input raw points.
 *
 * @see benchTwoParPWLCurves.cpp
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


    double sTime, eTime;

    std::cout << "Generating a PWL IV-curve with "<< ITRNUM <<" points...\n";
    TIMESTAMP( sTime );
	MonoDescPWL pwl(PS);
    TIMESTAMP( eTime );

    std::cout << "The elasped time is " << eTime-sTime << " s\n";
}
