/**
 * Benchmark on the conversion from voltage to current in terms of the input parameter.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../OneDiodeModel.h"
#include <iostream>

constexpr int ITRNUM = 1000000;

class OneDiodeModel;

int main()
{
	GetSimOptFile* filePtr = GetSimOptFile::Instance();
#ifdef MSXML
	filePtr->openOptionFile("SimOption.xml");
#else
	filePtr->openOptionFile("../../../Resources/SimOption.xml");
	//TODO: to change the path
#endif
	filePtr->readOptionFile();

	ArrayXd Vol(ITRNUM);
  ArrayXd Curr(ITRNUM);
  struct Params input;
  double start, end;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(-50.0, 50.0);

  std::vector<double> volVec;

  for (int iter = 0; iter < ITRNUM; iter++)
  {
    volVec.push_back(dis(gen));
  }

  if (!std::is_sorted(volVec.begin(), volVec.end()))
    std::sort(volVec.begin(), volVec.end());


  for (int iter = 0; iter < ITRNUM; iter++)
  {
    Vol(iter) = volVec.at(iter);
  }

  std::cout<<"\nBenchmarking: OneDiodeModelCurrent (Vol -> Curr) \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0.0, 64.8, 2.102623456790123e-06};


  double sTime, eTime;

  TIMESTAMP( sTime );
  OneDiodeModel odmObj(input);
  odmObj.VoltoCurr(Vol);
  TIMESTAMP( eTime );

  std::cout << "The elasped time is " << eTime-sTime << " s\n";
}
