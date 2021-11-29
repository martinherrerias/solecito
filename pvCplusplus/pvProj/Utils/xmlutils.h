/**
 * Operations on the SimOption.xml file. 
 *
 * Basic: open and read file, read and update its elements.
 * Further: parse tolerance (element 'tol') and tolFun to get the tolerance in the direction of current or voltage.
 */
#ifndef XMLUTILS_H_
#define XMLUTILS_H_

#include "pvcommutils.h"
#include <iostream>
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>
#include <cstdarg>

#ifdef MSXML
#define NOMINMAX
#include <comutil.h>
//#include <msxml.h>

#import "C:\Users\huanz\source\PV-PWL\PV-PWL\msxml3.dll"
#else
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <string>
using boost::property_tree::ptree;
#endif
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::Array<double, Eigen::Dynamic, 1> ArrayXd;

// Its Value is somehow dependent on user.
struct TolProp
{
	double RelTol;
	double minAbsTol;
	double minRelTol;
};

/// Singleton class.
class GetSimOptFile
{	
private:
	GetSimOptFile() {}; /// Private.
	static GetSimOptFile* m_pInstance;

	std::string fileName;

	/// Def: Default; Opt: User-defined.
	struct TolProp Def, User;
/*	double RelTol;
	double minAbsTol;
	double minRelTol; */
	/// Properties.
	double neps; //TODO: ask what the neps is? the data type?
	unsigned MaxIter;
	unsigned MaxRecur;

#ifdef MSXML
	/// Dom access handlers
	HRESULT hResult;
	IXMLDOMDocumentPtr spDocInput;
	IXMLDOMElementPtr spElemTemp;
	IXMLDOMAttributePtr spAttrTemp;
#else
	ptree pt;
#endif

#ifdef MSXML
	_bstr_t getValueByName(const char* str);
#endif

public:
	static GetSimOptFile* Instance();
	void openOptionFile(std::string optFile);
	template<class T>
  	inline	T getElement(std::string eleName_)
  	{
    	if (!strcmp(eleName_.c_str(), "RelTol"))
			return User.RelTol;
	  	else if (!strcmp(eleName_.c_str(), "minRelTol"))
			return User.minRelTol;
	  	else if (!strcmp(eleName_.c_str(), "minAbsTol"))
			return User.minAbsTol;
	  	else if (!strcmp(eleName_.c_str(), "NEPS"))
			return neps;
	  	else if (!strcmp(eleName_.c_str(), "MaxIter"))
			return MaxIter;
		else if (!strcmp(eleName_.c_str(), "MaxRecur"))
			return MaxRecur;
	  	else {
			std::string errorMsg = std::string(", could not identify the elemnet ") + eleName_;
			CATCHRTERROR(errorMsg);
	  	}
  	}

  	template<class T>
  	inline	void setElement(std::string eleName_, T val)
  	{
    	if (!strcmp(eleName_.c_str(), "RelTol"))
			User.RelTol = val;
	  	else if (!strcmp(eleName_.c_str(), "minRelTol"))
			User.minRelTol = val;
	  	else if (!strcmp(eleName_.c_str(), "minAbsTol"))
			User.minAbsTol = val;
	  	else {
			std::string errorMsg = std::string(", could not identify the elemnet ") + eleName_;
			CATCHRTERROR(errorMsg);
	  	}
  	}
  	void readOptionFile();
	void parseTolerance(std::vector<double>& tol);
	
	/// tolx/tolv: Dim::vert; toly/toli: Dim::hori.
	template<class T>
	inline T tolFun(std::vector<double> t, T x, double neps_, int dim)
	{
		TOLSIZEASSERT;
		return MAX3NUM(ABS(x)*t[2], neps_*EPSILONScal<T>(x), t[dim]);
	}

	inline ArrayXd tolFunVec(std::vector<double> t, ArrayXd x, double neps_, int dim)
	{
		TOLSIZEASSERT;
		return ((x.abs()*t[2]).max(neps_*EPSILONVec<ArrayXd>(x, x.rows()))).max(t[dim]);
	}
};
#endif
