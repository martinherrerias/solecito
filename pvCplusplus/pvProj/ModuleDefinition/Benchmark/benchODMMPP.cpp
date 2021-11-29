/**
 * Benchmark on obtaining the MPP of the IV curve characterized by the input parameter.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../OneDiodeModel.h"
#include <omp.h>
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

  ArrayXd VmpSet(ITRNUM);
  ArrayXd ImpSet(ITRNUM);
  ArrayXd PmpSet(ITRNUM);
  ArrayXd VocSet(ITRNUM);
  ArrayXd IscSet(ITRNUM);
  ArrayXd cSet(ITRNUM);

  
  std::cout<<"\nbenchmarking: ODMMPP \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0, 64.8, 2.102623456790123e-06};

  TIMESTAMP( start );
#pragma omp parallel for if (ITRNUM >= 100)
  for (int i = 0; i < ITRNUM; i++)
  {
    OneDiodeModel odmMppObj(input);
    odmMppObj.getODMMpp();
    VmpSet(i) = odmMppObj.Vmp;
    ImpSet(i) = odmMppObj.Imp;
    PmpSet(i) = odmMppObj.Pmp;
    VocSet(i) = odmMppObj.Voc;
    IscSet(i) = odmMppObj.Isc;
    cSet(i)   = odmMppObj.c;
  }
  TIMESTAMP( end );
  std::cout<<"The elapsed time is "<<end - start<<"\n";
 
}
