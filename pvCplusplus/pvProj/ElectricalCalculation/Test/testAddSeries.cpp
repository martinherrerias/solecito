/**
 * Test on adding more than two curves in series.
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

	// Dataset for douglas
	PS1 << 1.05, 1.8,
      1.08, 1.4,
      1.1, 1.25,
      3.1, 1.1;

  PS2 << 1, 2.15,
      1.7, 1.95,
      2, 0.5,
      3, 0;

	refPS << 2.0238, 2.15,
        2.7388, 1.95,
        2.781, 1.8,
        2.8938, 1.4,
        2.9448, 1.25,
        4.9759, 1.1,
        13.1, 0.5,
        20.7667, 0;


  std::cout<<"\nTesting: 1) addseries module without repeated/close curr or volt.\n";
	MonoDescPWL pwl1(PS1);
  MonoDescPWL pwl2(PS2);

  std::vector<MonoDescPWL*> set;
  set.push_back(&pwl1);
  set.push_back(&pwl2);
  MonoDescPWL pwl3 = MonoDescPWL(set, Dim::hori);

  bool res = EXPECTMATRIX_EQ(refPS, pwl3.pointSet);

  if(res)
	  std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 2) addseries module with repeated/close curr or volt.\n";

  MatrixXd refPS2(7,2);
  PS1 << 1.05, 1.95,
      1.2, 1.8,
      1.4, 1.25,
      3.1, 1.1;

  PS2 << 1, 2.15,
      1.7, 1.95,
      2, 0.5,
      3, 0.5;

  refPS2 << 1.85, 2.15,
        2.75, 1.95,
        2.93103, 1.8,
        3.24483, 1.25,
        4.97586, 1.1,
        11.9, 0.5,
        12.9, 0.5;

  MonoDescPWL pwl4(PS1);
  MonoDescPWL pwl5(PS2);

  std::vector<MonoDescPWL*> set2;
  set2.push_back(&pwl4);
  set2.push_back(&pwl5);
  MonoDescPWL pwl6 = MonoDescPWL(set2, Dim::hori);

  res = EXPECTMATRIX_EQ(refPS2, pwl6.pointSet);
    if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  
  std::cout<<"\nTesting: 3) addseries module with three repeated/close curr or volt.\n";
  MatrixXd refPS3(10,2);
  MatrixXd PS3(ITRNUM+1,2);
  PS1 << 1.05, 1.8,
      1.08, 1.4,
      1.1, 1.25,
      3.4, 0.3;

  PS2 << 1, 1.7,
      1.7, 1.65,
      2, 0.5,
      3, 0.3;

  PS3 << 1.0, 1.8,
      1.2, 1.8,
      2.4, 0.8,
      2.8, 0.4,
      2.9, 0.3;

  refPS3 << 1.65, 1.8,
      1.85, 1.8,
      3.3775, 1.7,
      4.14125, 1.65,
      4.52522, 1.4,
      4.76435, 1.25,
      6.51121, 0.8,
      7.61579, 0.5,
      8.45789, 0.4,
      9.3, 0.3;

  MonoDescPWL pwl7(PS1);
  MonoDescPWL pwl8(PS2);
  MonoDescPWL pwl9(PS3);

  std::vector<MonoDescPWL*> set3;
  set3.push_back(&pwl7);
  set3.push_back(&pwl8);
  set3.push_back(&pwl9);
  MonoDescPWL pwla = MonoDescPWL(set3, Dim::hori);

  res = EXPECTMATRIX_EQ(refPS3, pwla.pointSet);
  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 4) addseries module with three repeated/close curr or volt using clip.\n";
  MatrixXd refPS4(8,2);

  refPS4 << 2.61375, 1.75,
      3.3775, 1.7,
      4.14125, 1.65,
      4.52522, 1.4,
      4.76435, 1.25,
      6.51121, 0.8,
      7.61579, 0.5,
      9.3, 0.3;

  MonoDescPWL pwlb = MonoDescPWL(set3, Dim::hori, {EPSILON0}, {1.8, 9.4, 1.75, 0.3});

  res = EXPECTMATRIX_EQ(refPS4, pwlb.pointSet);
  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}
