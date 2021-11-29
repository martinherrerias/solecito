/**
 * Customize the operator and operand for the stl algorithms.
 */
#ifndef COMPARE_H
#define COMPARE_H

using Point = std::pair<double, double>;

struct customDoubleDesc {
	bool operator()(double x, double y) const
	{
		return (x > y);
	}
};

struct customDoubleUpper {
	bool operator()(double x, double y) const
	{
		return !(x < y);
	}
};


struct customFirstLess{
	bool operator()(Point x, Point y) const
	{
		return (x.first < y.first);
	}
};

struct customSecLess{
	bool operator()(Point x, Point y) const
	{
		return (x.second < y.second);
	}
};

struct customUpperFirst {
	bool operator() (Point it, double y) const
	{
		return !(it.first > y);
	}
};

struct customUpperSecond {
	bool operator() (Point it, double y) const
	{
		return !(it.second > y);
	}
};
struct customLowerFirst {
	bool operator() (Point it, double y) const
	{
		return (it.first < y);
	}
};
struct customLowerSecond {
	bool operator() (Point it, double y) const
	{
		return (it.second < y);
	}
};
struct customDotProduct {
	bool operator() (Point x, Point y) const
	{
		return (((x).first)*((x).second)) < (((y).first)*((y).second));
	}
};


#endif
