/**
 * Test on the Douglaspeucker function.
 */
#include "../../test.h"
#include "../MonoDescPWL.h"
#include <iostream>

class MonoDescPWL;

constexpr int ITRNUM = 10;

int main()
{
	GetSimOptFile* filePtr = GetSimOptFile::Instance();
#ifdef MSXML
	filePtr->o
    
#undef NDEBUGpenOptionFile("SimOption.xml");
#else
	filePtr->openOptionFile("../../../Resources/SimOption.xml");
	//TODO: to change the path
#endif
	filePtr->readOptionFile();

	MatrixXd PS(ITRNUM,2);
	MatrixXd refPS(5,2);

	// Dataset for douglas
	PS << 0.0, 0.0,
		1.0, 0.1,
		2.0, -0.1,
		3.0, 5.0,
		4.0, 6.0,
		5.0, 7.0,
		6.0, 8.1,
		7.0, 9.0,
		8.0, 9.0,
		9.0, 9.0;

	refPS << 0.0, 9.0,
			1.0, 5.0,
			2.0, 0.1,
			3.0, 0.0,
			9.0, -0.1;

  std::cout<<"\nTesting: douglaspeucker function.\n";
	MonoDescPWL pwl(PS, {1.0});
  bool res = EXPECTMATRIX_EQ(refPS, pwl.pointSet);
  if(res)
	  std::cout<<"Testing: success finished\n";
  else
    std::cout<<"Testing: failed finished\n";
}
