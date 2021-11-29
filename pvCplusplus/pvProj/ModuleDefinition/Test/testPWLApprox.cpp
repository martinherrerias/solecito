/**
 * Test on the linear approximation of the IV curve.
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
  

  std::cout<<"\nTesting: 1) PWLApprox forward-biased/reverse-biased (BB) \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0, 64.8, 2.102623456790123e-06};

  OneDiodeModel odmObj(input);
  odmObj.pwlApprox({0.0});
  
  res = EXPECTNUM_EQ((odmObj.approxPointSet.rows()), 2115);
//  std::cout<<odmObj.approxPointSet<<"\n";
  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";


  std::cout<<"\nTesting: 1) PWLApprox forward-biased/reverse-biased (BB) with limitation \n";

  odmObj.pwlApprox({0.0}, {0.3, 44.7, 4.3, -5.0});
  
  res = EXPECTNUM_EQ((odmObj.approxPointSet.rows()), 764);

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}
