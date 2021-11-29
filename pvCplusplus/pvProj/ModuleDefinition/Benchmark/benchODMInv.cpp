/**
 * Benchmark on the conversion from current to voltage in terms of the input parameter.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../OneDiodeModel.h"
#include "../pvProj/Utils/Compare.h"
#include <iostream>

constexpr int ITRNUM = 10000;

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

//	ArrayXd Vol(ITRNUM);
  ArrayXd Curr(ITRNUM);
  struct Params input;
  double start, end;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 2.0);

  std::vector<double> currVec;

  for (int iter = 0; iter < ITRNUM; iter++)
  {
    currVec.push_back(dis(gen));
  }

  if (!std::is_sorted(currVec.begin(), currVec.end()))
    std::sort(currVec.begin(), currVec.end(), customDoubleDesc());


  for (int iter = 0; iter < ITRNUM; iter++)
  {
    Curr(iter) = currVec.at(iter);
  }

  std::cout<<"\nBenchmarking: OneDiodeModelCurrent (Curr -> Vol) \n";
  input = {0.878627715718, 0.0000000711589600036043, 2.04, 3.43435521606909, 7158.91256877715 , 0.35, 69.3, 0.00000132231404958678};


  double sTime, eTime;

  TIMESTAMP( sTime );
  OneDiodeModel odmObj(input);
  odmObj.CurrtoVol(Curr);
  TIMESTAMP( eTime );

  std::cout << "The elasped time is " << eTime-sTime << " s\n";
}
