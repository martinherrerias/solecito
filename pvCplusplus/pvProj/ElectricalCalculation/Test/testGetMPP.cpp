/**
 * Test on the calculation of Maximum Power Point (MPP) of the given curve.
 */
#include "../../test.h"
#include "../../Utils/pvcommutils.h"
#include "../MonoDescPWL.h"
#include <iostream>

class MonoDescPWL;

constexpr int ITRNUM = 4;

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

	MatrixXd PS1(ITRNUM,2);
  MatrixXd PS2(ITRNUM,2);

	// Dataset for douglas
	PS1 << 1.05, 1.8,
      1.08, 1.4,
      1.1, 1.25,
      3.1, 1.1;

  PS2 << 1, 2.15,
      1.7, 1.95,
      2, 0.5,
      3, 0;

  double refPmp = 7.28979;
  double refVmp = 9.93506;
  double refPmpl = 7.2252;
  double refVmpl = 9.0;

  std::cout<<"\nTesting: 1) getMPP function with infinite boundary.\n";
	MonoDescPWL pwl1(PS1);
  MonoDescPWL pwl2(PS2);

  std::vector<MonoDescPWL*> set;
  set.push_back(&pwl1);
  set.push_back(&pwl2);
  MonoDescPWL pwl3 = MonoDescPWL(set, Dim::hori);
  pwl3.getMpp({});

  bool res1 = EXPECTNUM_EQ(pwl3.mppInfo.pmp, refPmp);
  bool res2 = EXPECTNUM_EQ(pwl3.mppInfo.vmp, refVmp);

  if(res1 && res2)
	  std::cout<<"Testing: success finished\n";
  else
  {
    std::cout<<"Testing: failed finished\n";
  }

  std::cout<<"\nTesting: 2) getMPP function with finite boundary.\n";

  pwl3.getMpp({2.0, 9.0});

  res1 = EXPECTNUM_EQ(pwl3.mppInfo.pmpBnd, refPmpl);
  res2 = EXPECTNUM_EQ(pwl3.mppInfo.vmpBnd, refVmpl);

  if(res1 && res2)
	  std::cout<<"Testing: success finished\n";
  else
  {
    std::cout<<"Testing: failed finished\n";
  }
}
