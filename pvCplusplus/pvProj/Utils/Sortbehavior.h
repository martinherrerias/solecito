/**
 * Collection of sorting algorithms.
 *
 * IglSort: provided by igl library.
 * QuickSort
 * StlSort: provided by stl.
 * MergeSort: sort a series of already-ordered curves.
 */
#ifndef SORTBEHAVIOR_H_
#define SORTBEHAVIOR_H_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include "pvcommutils.h"
#include <vector>
#include <queue>
#include <include/igl/unique_rows.h>

using namespace Eigen;
typedef std::pair<double, std::pair<int, int>> ppi;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> MatrixXRd;


class GeneralSort
{
protected:
	MatrixXd* ptrPS{ nullptr };
	int low{ 0 };
	int high{ 0 };
	Dim dim{ Dim::vert };
	
public:
	/// Default constructor
	GeneralSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_)
		: ptrPS(ptrPS_), low(low_), high(high_), dim(dim_) {} /// Initialization
	virtual void sort() = 0;
};

/* 
* To be removed, too slow
class InsertSort : public GeneralSort
{
public:
	InsertSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_) : GeneralSort(ptrPS_, low_, high_, dim_) {}
	
	virtual void sort()
	{//sort rows
		for (int i = low + 1; i <= high; i++)
		{
			double temp = (*ptrPS)(i, dim);
			VectorXd vecTmp = ptrPS->row(i);
			int j = i - 1;
			while ((*ptrPS)(j, dim) > temp && j >= low)
			{
				ptrPS->row(j+1) = (ptrPS->row(j));
				j--;
			}
			ptrPS->row(j+1) = vecTmp;
		}
	}
};*/

class IglSort : public GeneralSort
{
public:
	IglSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_) : GeneralSort(ptrPS_, low_, high_, dim_) {}
	
	virtual void sort()
	{
		(ptrPS->col(1)).array() *= -1;
		Eigen::VectorXi Ia, Ic;
		MatrixXd resPS;
		igl::unique_rows(*ptrPS, resPS, Ia, Ic);
		resPS.col(1).array() *= -1;
		*ptrPS = resPS;
	}
};

class QuickSort : public GeneralSort
{
public:
	QuickSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_) : GeneralSort(ptrPS_, low_, high_, dim_) {}

	virtual void quickSort(const int& l, const int& h)
	{
		if (l < h)
		{
			//	double pivot = refPointSet(high, 1);
			double pivot = (*ptrPS)(h, dim);
			int i = (l - 1);

			for (int j = l; j <= h - 1; j++)
			{
				if ((*ptrPS)(j, dim) <= pivot)
				{
					i++;
					ptrPS->row(i).swap(ptrPS->row(j));
					//	std::swap(pointSet(i, 0), pointSet(j, 0));
					//	std::swap(pointSet(i, 1), pointSet(j, 1));
				}
			}

			ptrPS->row(i + 1).swap(ptrPS->row(h));
			//	std::swap(pointSet(i + 1, 0), pointSet(high, 0));
			//	std::swap(pointSet(i + 1, 1), pointSet(high, 1)); // TODO: replace them with swap between rows
			quickSort(l, i);
			quickSort(i + 2, h);
		}
	}

	virtual void sort()
	{
		quickSort(low, high);
	}
};

class StlSort : public GeneralSort
{
public:
	StlSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_) : GeneralSort(ptrPS_, low_, high_, dim_) {}

	virtual void sort()
	{
		std::vector<VectorXd> vec;
		for (int i = low; i <= high; ++i)
			vec.push_back(ptrPS->row(i));
		if (dim)
			std::sort(vec.begin(), vec.end(), [](const VectorXd& lhs, const VectorXd& rhs) {return lhs(1) < rhs(1); }); // when dim = hori, then according to Y
		else
			std::sort(vec.begin(), vec.end(), [](const VectorXd& lhs, const VectorXd& rhs) {return lhs(0) < rhs(0); }); // when dim = vert, then according to X

		for (int i = low; i <= high; ++i)
			ptrPS->row(i) = vec[i];
	}
};

class MergeSort : public GeneralSort
{
private:
	const std::vector<int>* rowIndexSetPtr;
public:
	MergeSort(MatrixXd* ptrPS_, const int& low_, const int& high_, const Dim& dim_, const std::vector<int>* rowIndexSetPtr_) : GeneralSort(ptrPS_, low_, high_, dim_), rowIndexSetPtr(rowIndexSetPtr_)
	{}

	virtual void sort()
	{
		int arraySize = rowIndexSetPtr->size() - 1;
		int outputIndex = 0;
		MatrixXd tempPS(ptrPS->rows(), ptrPS->cols());

		std::priority_queue<ppi, std::vector<ppi>, std::greater<ppi>> pq;

		for (int i = 0; i < arraySize; i++)
			pq.push({ (*ptrPS)((*rowIndexSetPtr)[i], dim), { i, (*rowIndexSetPtr)[i] } }); 
		/// The data structure in the pq is {double value in direction dim, {arrayIndex in the rowIndexSetPtr, row number in the pointSet}}

		while (!pq.empty())
		{
			ppi curr = pq.top();
			pq.pop();

			int setIndex = curr.second.first;
			int rowNum = curr.second.second;
			tempPS.row(outputIndex++) = ptrPS->row(rowNum);

			if ((rowNum + 1) < (*rowIndexSetPtr)[setIndex + 1])
				pq.push({ (*ptrPS)(rowNum + 1, dim), { setIndex, rowNum + 1 } });
		}

		ptrPS->resize(0, 0);
		*ptrPS = std::move(tempPS);
	}
};

#endif