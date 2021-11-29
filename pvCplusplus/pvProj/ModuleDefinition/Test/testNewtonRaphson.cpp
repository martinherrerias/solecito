/**
 * Test on NewtonRaphson iteration algorithm.
 */
#include "../../test.h"
#include "../../testUtils.h"
#include "../../Utils/pvcommutils.h"
#include "../../Utils/Rootfind.h"
#include <boost/math/tools/roots.hpp>
#include <iostream>

constexpr size_t arraySize = 8;

template <class T>
struct cbrt_functor_2
{ // Functor also returning 1st derviative.
  cbrt_functor_2(T const& to_find_root_of) : value(to_find_root_of)
  { // Constructor stores value to find root of,
    // for example: calling cbrt_functor_2<T>(x) to use to get cube root of x.
  }
  std::pair<T, T> operator()(T const& x)
  { // Return both f(x) and f'(x).
    T fx = x*x - value; // Difference (estmate x^3 - value).
    T dx =  2 *x; // 1st derivative = 3x^2.
    return std::make_pair(fx, dx); // 'return' both fx and dx.
  }
private:
  T value; // to be 'cube_rooted'.
}; // cbrt_functor_2


int main()
{
    GetSimOptFile* filePtr = GetSimOptFile::Instance();
#ifdef MSXML
    filePtr->openOptionFile("SimOption.xml");
#else
    filePtr->openOptionFile("../../../Resources/SimOption.xml");
  
#endif
    filePtr->readOptionFile();
    double neps = filePtr->getElement<double>("NEPS");
    std::vector<double> tol = { 0.0 };
    filePtr->parseTolerance(tol);

    std::cout<<"\nTesting: Newton-Raphson iteration method.\n";
    auto testLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd x) -> ArrayXd {return x*x - 9;});
    auto dtestLambda = std::function<ArrayXd(ArrayXd)>([=](ArrayXd x) -> ArrayXd {return 2*x;});

    ArrayXd seedSet(arraySize);
    seedSet << 1.3, 2.3, 3.4, 2.5, 1.3, 2.3, 3.4, 2.5;
      int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits accuracy for type T.
  // digits used to control how accurate to try to make the result.
  int get_digits = (digits * 3) /4; // Near maximum (3/4) possible accuracy.
 //   unsigned maxIter = 1000;
  unsigned max = 1000;
  boost::uintmax_t maxIter = max;
    double start, end;
    TIMESTAMP( start );

    for(int i = 0; i < seedSet.matrix().rows(); i++)
    {
        seedSet(i) = boost::math::tools::newton_raphson_iterate(cbrt_functor_2<double>(9), seedSet(i), -Inf, Inf, get_digits, maxIter);
    }
/*    NewtonRaphson nrMethod(std::move(testLambda), filePtr,
                  tol, neps, std::move(dtestLambda),
                  maxIter, seedSet);
    nrMethod.findRoot(); */
    TIMESTAMP( end );
    std::cout<<"the elapsed time is "<<end - start<<"\n";

//    std::cout<<"the result is \n"<<seedSet<<"\n";
    ArrayXd refRoot(arraySize);
    refRoot << 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0;

    bool res = EXPECTMATRIX_EQ(refRoot.matrix(), seedSet.matrix());

    if(res)
      std::cout<<"Testing: success finished\n";
    else
      std::cout<<"Testing: failed finished\n";
}
