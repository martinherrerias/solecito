/**
 * Benchmark on approximating the IV curve in terms of the input parameters.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../OneDiodeModel.h"
#include <iostream>

constexpr int ITRNUM = 5;

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


  struct Params input;
  double start, end;
  bool   res;

  double refVmp;
  double refImp;
  double refPmp;
  double refVoc;
  double refIsc;
  double refc;

  std::cout<<"\nBenchmark: 1) PWLApprox forward-biased&reverse-biased (BB) \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0, 64.8, 2.102623456790123e-06};

  double startTime, endTime;

  TIMESTAMP( startTime ); /// see testUtils.h
  OneDiodeModel odmObj(input);

  odmObj.pwlApprox({0.0});
//  odmObj.pwlApprox({0.0}, {-20.2, 18.2});
  TIMESTAMP( endTime );
  
  std::cout<<"the elapsed time is "<<endTime - startTime<<"\n";
}
