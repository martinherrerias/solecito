/**
 * Test on adding more than two curves in parallel.
 * 
 * These curves could have repeated curr/volt or not.
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
	MatrixXd refPS(8,2);

	std::cout<<"\nTesting: 1) addparallel module without repeated/close curr or volt.\n";
	PS1 << 1.05, 1.8,
      1.08, 1.4,
      1.1, 1.25,
      3.1, 1.1;

  PS2 << 1, 2.15,
      1.7, 1.95,
      2, 0.5,
      3, 0;

	refPS << 1.0, 4.61667,
        1.05, 3.93571,
        1.08, 3.52714,
        1.1, 3.37143,
        1.7, 3.155,
        2.0, 1.6825,
        3.0, 1.1075,
        3.1, 1.05;

	MonoDescPWL pwl1(PS1);
  MonoDescPWL pwl2(PS2);

  std::vector<MonoDescPWL*> set;
  set.push_back(&pwl1);
  set.push_back(&pwl2);
  MonoDescPWL pwl3 = MonoDescPWL(set, Dim::vert);

  bool res = EXPECTMATRIX_EQ(refPS, pwl3.pointSet);

  if(res)
	  std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 2) addparallel module with repeated/close curr or volt.\n";

  MatrixXd refPS2(7,2);
  PS1 << 1, 1.8,
      1.2, 1.4,
      2, 1.25,
      3.1, 1.1;

  PS2 << 1, 2.15,
      1.1, 1.95,
      1.1, 0.5,
      3, 0;

  refPS2 << 1, 3.95,
        1.1, 3.55,
        1.1, 2.1,
        1.2, 1.87368,
        2, 1.51316,
        3, 1.11364,
        3.1, 1.07368;

  MonoDescPWL pwl4(PS1);
  MonoDescPWL pwl5(PS2);

  std::vector<MonoDescPWL*> set2;
  set2.push_back(&pwl4);
  set2.push_back(&pwl5);
  MonoDescPWL pwl6 = MonoDescPWL(set2, Dim::vert);

  res = EXPECTMATRIX_EQ(refPS2, pwl6.pointSet);

  if(res)
	  std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 3) addparallel module with three repeated/close curr or volt.\n";
  MatrixXd refPS3(9,2);
  MatrixXd PS3(ITRNUM+1,2);
  PS1 << 0.6, 1.8,
      1.0, 1.4,
      1.0, 1.25,
      3, 0.3;

  PS2 << 1.0, 1.7,
      1.7, 1.65,
      2, 0.5,
      3, 0.3;

  PS3 << 1.0, 1.8,
      1.2, 1.6,
      2.4, 0.8,
      3, 0.3,
      3, 0.2;

  refPS3 << 0.6, 5.72857,
      1, 4.9,
      1, 4.75,
      1.2, 4.44071,
      1.7, 3.83417,
      2, 2.34167,
      2.4, 1.805,
      3, 0.9,
      3, 0.8;

  MonoDescPWL pwl7(PS1);
  MonoDescPWL pwl8(PS2);
  MonoDescPWL pwl9(PS3);

  std::vector<MonoDescPWL*> set3;
  set3.push_back(&pwl7);
  set3.push_back(&pwl8);
  set3.push_back(&pwl9);
  MonoDescPWL pwla = MonoDescPWL(set3, Dim::vert);

  res = EXPECTMATRIX_EQ(refPS3, pwla.pointSet);
  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 4) addparallel module with three repeated/close curr or volt using clip.\n";
  MatrixXd refPS4(9,2);

  refPS4 << 0.5, 5.93571,
        1, 4.9,
        1, 4.75,
        1.2, 4.44071,
        1.7, 3.83417,
        2, 2.34167,
        2.4, 1.805,
        3, 0.9,
        3, 0.7;

  MonoDescPWL pwlb = MonoDescPWL(set3, Dim::vert, {EPSILON0}, {0.5, 4, 6, 0.7});

  res = EXPECTMATRIX_EQ(refPS4, pwlb.pointSet);
  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

}
