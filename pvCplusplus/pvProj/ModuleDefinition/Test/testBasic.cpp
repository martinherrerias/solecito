/**
 * Test on converting from current to voltage.
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

	ArrayXd Vol(ITRNUM);
  ArrayXd Curr(ITRNUM);
  struct Params input;
  double start, end;
  bool res;

  std::cout<<"\nTesting: 1) OneDiodeModelCurrent (Vol -> Curr) \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0, 64.8, 2.102623456790123e-06};

  Vol << 0.0, 
       22.32154,
       35.48047853789603,
       40.06177931022444,
       44.64308008255284;

  OneDiodeModel odmObj(input);
  odmObj.VoltoCurr(Vol);
  ArrayXd refIset(ITRNUM);
  refIset <<  11.833109775636084,
              11.776671149702414,
              11.059055521310057,
              7.726205852087536,
              0.0;

  res = EXPECTMATRIX_EQ(refIset.matrix(), odmObj.getIset().matrix());

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
 
  std::cout<<"\nTesting: 2) OneDiodeModelVoltage (Curr -> Vol) \n";

  Curr << 2.05, 
         1.6,
         -2.1,
         2.5,
         2.9;

  input = {1.3, 1.2, 1.0, 0.4, 3.2, -1.4, 2.3, 2.3};

  OneDiodeModel odmObj2(input);
  odmObj2.CurrtoVol(Curr);
  ArrayXd refVset(ITRNUM);
  refVset << -2.03632,
            -1.46097,
            2.70792,
            -2.64394,
            -3.18803;

  res = EXPECTMATRIX_EQ(refVset.matrix(), odmObj2.getVset().matrix()); /// see test.h

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}
