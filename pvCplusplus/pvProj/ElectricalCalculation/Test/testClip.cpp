/**
 * Test on the clip function in terms of the given boundary (lim).
 */
#include "../../test.h"
#include "../../Utils/pvcommutils.h"
#include "../MonoDescPWL.h"
#include <iostream>

class MonoDescPWL;

constexpr int ITRNUM = 9;

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

	MatrixXd PS(ITRNUM,2);
	MatrixXd refPS(4,2);

	// Dataset for douglas
	PS << 0.0, 9.0,
		2.0, 8.5,
		3.0, 8.1,
		4.0, 7.0,
		5.0, 6.0,
		6.0, 5.0,
		7.0, 0.1,
		8.0, 0.0,
		9.0, -0.1; 

	refPS << -0.8, 9.2,
   		2, 8.5,
   		3, 8.1,
 		3.5, 7.55;

    std::cout<<"\nTesting: clip function.\n";
//  MonoDescPWL pwl(PS, {EPSILON}, {-1.3, 3.5, 8.8, 0.2});

    MonoDescPWL pwl(PS, {EPSILON0}, {-1.3, 3.5, 9.2, -0.2});
    bool res = EXPECTMATRIX_EQ(refPS, pwl.pointSet);

    if(res)
		std::cout<<"Testing: success finished\n";
  	else
    	std::cout<<"Testing: failed finished\n";
}
