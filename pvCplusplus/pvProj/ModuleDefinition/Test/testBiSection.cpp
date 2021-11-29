/**
 * Test on bisection algorithm.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../../Utils/Rootfind.h"
#include <iostream>

int main()
{
    GetSimOptFile* filePtr = GetSimOptFile::Instance();
    double neps = filePtr->getElement<double>("NEPS");
#ifdef MSXML
    filePtr->openOptionFile("SimOption.xml");
#else
    filePtr->openOptionFile("../../../Resources/SimOption.xml");
  
#endif
    filePtr->readOptionFile();
    std::vector<double> tol = { 0.0 };
    filePtr->parseTolerance(tol);

    std::cout<<"\nTesting: Bisection method.\n";
    auto testLambda = std::function<double(double)>([=](double x) -> double {return x*x - 9;});
    double root;
    double start, end;
    TIMESTAMP( start );
    BiSection bsMethod(std::move(testLambda), filePtr,
                  tol, neps, 2.4, 4.5, root);
    bsMethod.findRoot();
    TIMESTAMP( end );
    std::cout<<"the elapsed time is "<< end - start <<"\n";

    double refRoot = 3.0;
    bool res = (root - refRoot < 0.0001);

    if(res)
      std::cout<<"Testing: success finished\n";
    else
      std::cout<<"Testing: failed finished\n";
}
