#include "LambertW.h"

double lambertW(double x)
{
	double w;
	if (boost::math::isinf(x))
	{
		if (x > 0.0)
			w = Inf;
		else
			w = 0.0;
		return w;
	}
	if (boost::math::isnan(x))
	{
		return Nan;
	}

	double expLog = std::exp(x);
	if (boost::math::isinf(expLog))
	{
		w = x - std::log(x); /// too big and then use log(th) - log(log(th)) as seed
		double err = Inf;
		int iter = 0;

		while (ABS(err) > NEPS*std::exp(w))
		{
			err = w + std::log(w) - x;
			w = w - err;
			++iter;
			if (iter >= maxIter)
				w = Nan;
		}
	}
	else
	{
		w = boost::math::lambert_w0(expLog);
	}
	
	return w;
}