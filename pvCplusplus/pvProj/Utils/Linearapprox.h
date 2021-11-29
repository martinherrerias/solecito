/**
 * Collection of the linear approximation algorithms.
 * 
 * StackLA: use stack structure to replace the recursion.
 * TreeLA: use the tree structure to replace the recursion.
 * StaticLA: recursive method with minimum memory reallocation (best).
 * DynamicLA: recursive method with memory reallocation in each iteration.
 */
#ifndef LINEARAPPROX_H_
#define LINEARAPPROX_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <random>
#include "xmlutils.h"
#include "pvcommutils.h"
#include <omp.h>
#include <stack>

using namespace Eigen;
typedef Eigen::Array<double, Dynamic, 1> ArrayXd;

constexpr int TASKTHLD = 8;

class GeneralLinearApprox
{
protected:
	MatrixXd* approxPointSet;

	GetSimOptFile* filePtr;
	std::vector<double> tol;
	double neps;
	double biasK;
	double g{ 0.5 };
	ArrayXd limZ, limI;
	std::function<double(double)> fun;
	unsigned maxRecur;
	int approxPSSize;

	// Bias correction
	void biasCorr ()
	{
		ArrayXd divApproxPS(approxPSSize - 1);
		ArrayXd diffDiv(approxPSSize - 2);
		ArrayXd errMax(approxPSSize - 1);
		ArrayXd diffErrMax(approxPSSize - 2);
		ArrayXd offsetX(approxPSSize);
		ArrayXd offsetY(approxPSSize);
		for (int i = 0; i < approxPSSize - 1; i++)
		{
			divApproxPS(i) = ((*approxPointSet)(i + 1, 1) - (*approxPointSet)(i, 1))/((*approxPointSet)(i + 1, 0) - (*approxPointSet)(i, 0));
			errMax(i) = ((*approxPointSet)(i, 1) + (*approxPointSet)(i + 1, 1))/2.0
				- fun(((*approxPointSet)(i, 0) + (*approxPointSet)(i + 1, 0))/2.0);
			if (i > 0)
			{ 
				diffDiv(i - 1) = divApproxPS(i) - divApproxPS(i - 1);
				diffErrMax(i - 1) = errMax(i) - errMax(i - 1);
			}
		}

		for (int i = 0; i < approxPSSize; i++)
		{
			if (i == 0)
			{
				offsetY(i) = -1.0*biasK*errMax(i)/std::sqrt(1 + POW2<double>(divApproxPS(i)));
				offsetX(i) = -1.0*divApproxPS(i)*offsetY(i);
			}
			else if (i == (approxPSSize - 1))
			{
				offsetY(i) = -1.0*biasK*errMax(i - 1)/std::sqrt(1 + POW2(divApproxPS(i - 1)));
				offsetX(i) = -1.0*divApproxPS(i - 1)*offsetY(i);
			}
			else
			{
				offsetY(i) = biasK*(errMax(i)*divApproxPS(i - 1) - errMax(i - 1)*divApproxPS(i))/diffDiv(i - 1);
				offsetX(i) = biasK*diffErrMax(i - 1)/diffDiv(i - 1);
			}
		}

		(approxPointSet->col(0)).array() += offsetX;
		(approxPointSet->col(1)).array() += offsetY;
	}

	void initLimI ()
	{
		limI.resize(2);
		limI(0) = fun(limZ(0));
		limI(1) = fun(limZ(1));
	}

public:
	GeneralLinearApprox(MatrixXd* ps_, GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_, ArrayXd limZ_, double biasK_, std::function<double(double)> fun_)
	: approxPointSet(ps_), filePtr(filePtr_), tol(tol_), neps(neps_), limZ(limZ_), biasK(biasK_), fun(fun_)
	{ maxRecur = (filePtr->getElement<unsigned>)("MaxRecur"); };

	virtual void approx() = 0;
};

class  StackLA: public GeneralLinearApprox
{
private:
	int pos{ 0 };
	unsigned numVertices = 256; /// Empirical value.

	std::vector<std::pair<double, double>> mystack;
	int top{ -1 };

	void stackBBApprox()
	{
		while (top != -1)
		{
			std::pair<double, double> end = mystack[top];

			double mx = (end.first - (*approxPointSet)(pos, 0))*g + (*approxPointSet)(pos, 0);
			double my = fun(mx);
			double slope = (end.second - (*approxPointSet)(pos, 1))/(end.first - (*approxPointSet)(pos, 0));

			double dy = my - ((mx - (*approxPointSet)(pos, 0))*slope + (*approxPointSet)(pos, 1));

			if ((ABS(dy) <= pvMIN(filePtr->tolFun(tol, my, neps, Dim::hori), ABS(slope*(filePtr->tolFun(tol, mx, neps, Dim::vert)))))
			|| (ABS(end.first - (*approxPointSet)(pos, 0)) < 2.0*EPSILONScal(mx)) || (ABS(end.second - (*approxPointSet)(pos, 1)) < 2.0*EPSILONScal(my)))
			{
				++pos;
				if (pos == numVertices)
				{	
					numVertices *= 2;
					approxPointSet->conservativeResize(numVertices, 2);
				}
				
				(*approxPointSet)(pos, 0) = end.first;
				(*approxPointSet)(pos, 1) = end.second;
				--top;
			}
			else
			{
				++top;
				if (top == numVertices)
				{
					numVertices *= 2;
					mystack.resize(numVertices);
				}
				mystack[top] = std::make_pair(mx, my);
			}
		}
	}
public:
	StackLA(MatrixXd* ps_, GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_, ArrayXd limZ_, double biasK_, std::function<double(double)> fun_) 
	: GeneralLinearApprox(ps_, filePtr_, tol_, neps_, limZ_, biasK_, fun_) {
	} 

	virtual void approx()
	{ 	
		initLimI();
		mystack.resize(numVertices);
		approxPointSet->resize(numVertices, 2);
		(*approxPointSet)(pos, 0) = limZ(0);
		(*approxPointSet)(pos, 1) = limI(0);

	    double mx = (limZ(1) - limZ(0))*g + limZ(0);
		double my = fun(mx);
		double slope = (limI(1) - limI(0))/(limZ(1) - limZ(0));

		double dy = my - ((mx - limZ(0))*slope + limI(0));

		if ((ABS(dy) <= pvMIN(filePtr->tolFun(tol, my, neps, Dim::hori), ABS(slope*(filePtr->tolFun(tol, mx, neps, Dim::vert)))))
		|| (ABS(limZ(1) - limZ(0)) < 2.0*EPSILONScal(mx)) || (ABS(limI(1) - limI(0)) < 2.0*EPSILONScal(my)))
		{
			++pos;
			(*approxPointSet)(pos, 0) = limZ(1);
			(*approxPointSet)(pos, 1) = limI(1);
		}
		else
		{
			mystack[++top] = std::make_pair(limZ(1), limI(1));
			mystack[++top] = std::make_pair(mx, my);
			
			stackBBApprox();
			if (++pos != numVertices)
			{ approxPointSet->conservativeResize(pos, 2); }
		}

		approxPSSize = pos;
		biasCorr();
	}
};


class  TreeLA: public GeneralLinearApprox
{
	typedef struct node
	{
		double x;
		double y;
		struct node* next;
	} pwlnode;
private:
	int count { 0 };

	pwlnode* begin, *curr;

	void trBBApprox(pwlnode* curr)
	{
		while (curr->next != nullptr)
		{
			pwlnode* end = curr->next;
			double mx = (end->x - curr->x)*g + curr->x;
			double my = fun(mx);
			double slope = (end->y - curr->y)/(end->x - curr->x);

			double dy = my - ((mx - curr->x)*slope + curr->y);

			if ((ABS(dy) <= pvMIN(filePtr->tolFun(tol, my, neps, Dim::hori), ABS(slope*(filePtr->tolFun(tol, mx, neps, Dim::vert)))))
			|| (ABS(end->x - curr->x) < 2.0*EPSILONScal(mx)) || (ABS(end->y - curr->y) < 2.0*EPSILONScal(my)))
			{
				curr = curr->next;
			}
			else
			{
				++count;
				pwlnode* p = new pwlnode();
				p->x = mx;
				p->y = my;
				p->next = curr->next;
				curr->next = p;
			}
		}
	}
public:
	TreeLA(MatrixXd* ps_, GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_, ArrayXd limZ_, double biasK_, std::function<double(double)> fun_) 
	: GeneralLinearApprox(ps_, filePtr_, tol_, neps_, limZ_, biasK_, fun_) {
	} 

	virtual void approx()
	{ 
		initLimI();

		begin = new pwlnode();
		begin->x = limZ(0);
		begin->y = limI(0);
	    pwlnode* p = new pwlnode();
	    p->x = limZ(1);
	    p->y = limI(1);
	    begin->next = p;
	    p->next = nullptr;

		count += 2;

		curr = begin;
		trBBApprox(curr);

		approxPointSet->resize(count, 2);
		pwlnode* pre;
		p = begin;
		for (int i = 0; i < count; i++)
		{
			(*approxPointSet)(i,0) = p->x;
			(*approxPointSet)(i,1) = p->y;
			pre = p;
			p = p->next;
			delete pre;
		}
		approxPSSize = count;
		biasCorr();
	}
};

class  StaticLA: public GeneralLinearApprox
{
private:
	int pos {0};
	unsigned numVertices = 256; /// Empirical value.
	
	void sBBApprox(double lx, double rx, double ly, double ry, unsigned recurCount)
	{
		//TODO: move the repeated code below to a common part
		double mx = (rx - lx)*g + lx;
		double my = fun(mx);
		double slope = (ry - ly)/(rx - lx);

		double dy = my - ((mx - lx)*slope + ly);

		if ((ABS(dy) <= pvMIN(filePtr->tolFun(tol, my, neps, Dim::hori), ABS(slope*(filePtr->tolFun(tol, mx, neps, Dim::vert)))))
			|| (ABS(rx - lx) < 2.0*EPSILONScal(mx)) || (ABS(ry - ly) < 2.0*EPSILONScal(my)))
		{
			if (pos == numVertices)
			{
				numVertices *= 2;
				approxPointSet->conservativeResize(numVertices, 2);
			}
			(*approxPointSet)(pos, 0) = rx;
			(*approxPointSet)(pos, 1) = ry;
        	++pos;
		}
		else
		{
			if (recurCount > maxRecur)
				std::cerr << "Warning: OneDiodeModel, exceed the maximum recursive depth = "
				<< recurCount << "\n";

			sBBApprox(lx, mx, ly, my, recurCount + 1);			
			sBBApprox(mx, rx, my, ry, recurCount + 1);	
		}
	}
public:
	StaticLA(MatrixXd* ps_, GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_, ArrayXd limZ_, double biasK_, std::function<double(double)> fun_) 
	: GeneralLinearApprox(ps_, filePtr_, tol_, neps_, limZ_, biasK_, fun_) {
	} 

	virtual void approx()
	{ 
		initLimI();
		approxPointSet->resize(numVertices, 2);
		(*approxPointSet)(pos, 0) = limZ(0);
		(*approxPointSet)(pos, 1) = limI(0);	
		++pos;

		sBBApprox(limZ(0), limZ(1), limI(0), limI(1), 0);

		if (pos != numVertices)
			{approxPointSet->conservativeResize(pos, 2);}
		approxPSSize = pos;
		biasCorr();
	}	
};

class DynamicLA : public GeneralLinearApprox
{
private:
	void dBBApprox(double lx, double rx, double ly, double ry, MatrixXd* out, unsigned recurCount)
	{
		double mx = (rx - lx)*g + lx;
		double my = fun(mx);
		double slope = (ry - ly)/(rx - lx);

		double dy = my - ((mx - lx)*slope + ly);

		if ((ABS(dy) <= pvMIN(filePtr->tolFun(tol, my, neps, Dim::hori), ABS(slope*(filePtr->tolFun(tol, mx, neps, Dim::vert)))))
			|| (ABS(rx - lx) < 2.0*EPSILONScal(mx)) || (ABS(ry - ly) < 2.0*EPSILONScal(my)))
		{
			(*out).resize(2,2);
			(*out)<<lx, ly,
			        rx, ry;
		}
		else
		{
			MatrixXd recRes1, recRes2;

			if (recurCount > maxRecur)
				std::cerr << "Warning: OneDiodeModel, exceed the maximum recursive depth = "
				<< recurCount << "\n";

			dBBApprox(lx, mx, ly, my, &recRes1, recurCount + 1);			
			dBBApprox(mx, rx, my, ry, &recRes2, recurCount + 1);
			
			int size1 = recRes1.rows();
			int size2 = recRes2.rows();
			(*out).resize(size1+size2-1, 2);
			(*out).topRows(size1) = recRes1;
			(*out).bottomRows(size2-1) = recRes2.bottomRows(size2-1);
		}
	}
public:

	DynamicLA(MatrixXd* ps_, GetSimOptFile* filePtr_, std::vector<double> tol_, double neps_, ArrayXd limZ_, double biasK_, std::function<double(double)> fun_) 
	: GeneralLinearApprox(ps_, filePtr_, tol_, neps_, limZ_, biasK_, fun_) {
	} 

	virtual void approx()
	{ 
		initLimI();
		dBBApprox(limZ(0), limZ(1), limI(0), limI(1), approxPointSet, 0);		
		approxPSSize = approxPointSet->rows();
		biasCorr();
	}

};

#endif