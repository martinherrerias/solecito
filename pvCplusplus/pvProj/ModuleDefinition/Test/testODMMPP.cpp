/**
 * Test on obtaining MPP of the related IV curve.
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
  

  std::cout<<"\nTesting: 1) ODMMPP \n";
  input = {11.842363, 2.394890881115644e-09, 0.315, 2.000836411141469, 4.028255106927194e+02, 0, 64.8, 2.102623456790123e-06};

  OneDiodeModel odmMppObj(input);
  odmMppObj.getODMMpp();
  refVmp = 35.4805;
  refVoc = 44.6431;
  refIsc = 11.8331;
  refImp = 11.0591;
  refPmp = 392.381;
  refc   = 14.318;

  res = EXPECTNUM_EQ(odmMppObj.Vmp, refVmp)&&EXPECTNUM_EQ(odmMppObj.Voc, refVoc)&&EXPECTNUM_EQ(odmMppObj.Isc, refIsc)
                        &&EXPECTNUM_EQ(odmMppObj.Imp, refImp)&&EXPECTNUM_EQ(odmMppObj.Pmp, refPmp)&&EXPECTNUM_EQ(odmMppObj.c, refc);

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";

  std::cout<<"\nTesting: 2) ODMMPP \n";
  input = {0.878627715718, 0.0000000711589600036043, 2.04, 3.43435521606909, 7158.91256877715, 0.35, 69.3, 0.00000132231404958678};

  OneDiodeModel odmMppObj2(input);
  odmMppObj2.getODMMpp();

  refVmp = 45.2481;
  refVoc = 55.9562;
  refIsc = 0.873824;
  refImp = 0.798108;
  refPmp = 36.1129;
  refc   = 11.8667;

  res = EXPECTNUM_EQ(odmMppObj2.Vmp, refVmp)&&EXPECTNUM_EQ(odmMppObj2.Voc, refVoc)&&EXPECTNUM_EQ(odmMppObj2.Isc, refIsc)
                        &&EXPECTNUM_EQ(odmMppObj2.Imp, refImp)&&EXPECTNUM_EQ(odmMppObj2.Pmp, refPmp)&&EXPECTNUM_EQ(odmMppObj2.c, refc);

  if(res)
    std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}
