/**
 * Test on the aproximation of curve in terms of two given curves.
 */
#include "../../test.h"
#include "../../Utils/pvcommutils.h"
#include "../MonoDescPWL.h"
#include <iostream>

class MonoDescPWL;

constexpr int ITRNUM1 = 8;
constexpr int ITRNUM2 = 4;

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

	MatrixXd PS1(ITRNUM1,2);
  MatrixXd PS2(ITRNUM2,2);
	MatrixXd refPS(ITRNUM1,2);

	// Dataset for douglas
	PS1 << -1.2, 5.4,
      2.3, 4.3,
      2.5, 4.1,
      2.8, 3.2,
      3.2, 2.1,
      3.6, 1.6,
      4.3, 0.2,
      4.3, -1.2;

  PS2 << 1.05, 1.8,
      1.08, 1.4,
      1.1, 1.25,
      3.1, 1.1;


	refPS << -0.731314, 9.71086,
        1.40169, 8.61086,
        1.52769, 8.41086,
        1.73469, 7.51086,
        2.00769, 6.41086,
        2.26269, 5.91086,
        2.72469, 4.51086,
        2.76669, 3.11086;


	MonoDescPWL pwl1(PS1);
  MonoDescPWL pwl2(PS2);

  std::cout<<"\nTesting: 1) combine two curves with fraction smaller than 0.5.\n";
  MonoDescPWL pwl3 = MonoDescPWL(pwl1, pwl2, 0.4);

  bool res = EXPECTMATRIX_EQ(refPS, pwl3.pointSet);

  if(res)
	  std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 2) combine two curves with fraction larger than 0.5.\n";

  MatrixXd refPS2(4,2);

  MonoDescPWL pwl4 = MonoDescPWL(pwl1, pwl2, 0.9);
  refPS2 << 1.375, 0.722286,
         1.402, 0.322286,
         1.42, 0.172286,
         3.22, 0.0222857;

  res = EXPECTMATRIX_EQ(refPS2, pwl4.pointSet);

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}