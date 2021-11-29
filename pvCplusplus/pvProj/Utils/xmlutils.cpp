#include "xmlutils.h"

GetSimOptFile* GetSimOptFile::m_pInstance = NULL;

GetSimOptFile* GetSimOptFile::Instance()
{
	if (!m_pInstance)
		m_pInstance = new GetSimOptFile;

	return m_pInstance;
}

void GetSimOptFile::openOptionFile(std::string optFile_)
{
#ifdef MSXML
	hResult = S_OK;
	hResult = CoInitialize(NULL);

	/// Create the DOM Document object.
	hResult = spDocInput.CreateInstance(__uuidof(DOMDocument));

	_variant_t varXml(optFile_.c_str());
	VARIANT_BOOL isSucc;
	hResult = spDocInput->load(varXml, &isSucc);
	fileName = optFile_;
#else
	boost::property_tree::read_xml(optFile_, pt);
#endif
}

#ifdef MSXML
_bstr_t GetSimOptFile::getValueByName(const char* _str)
{
	IXMLDOMNodeListPtr spNLChildren;
	IXMLDOMNodePtr spNodeTemp;
	_bstr_t ele(_str);

	hResult = spDocInput->getElementsByTagName(ele, &spNLChildren);
	spNLChildren->get_item(0, &spNodeTemp);
	spNodeTemp->get_childNodes(&spNLChildren);
	spNLChildren->get_item(0, &spNodeTemp);
	BSTR val;
	spNodeTemp->get_text(&val);
	_bstr_t retVal(val);
	return retVal;
}
#endif

void GetSimOptFile::readOptionFile()
{
#ifdef MSXML /// Using  Msxml parser in Windows.
	_bstr_t v;

	v = getValueByName("minAbsTol");
	minAbsTol = wcstod(v, NULL);
	v = getValueByName("minRelTol");
	minRelTol = wcstod(v, NULL);
	v = getValueByName("RelTol");
	RelTol = wcstod(v, NULL);
	v = getValueByName("NEPS");
	neps = wcstod(v, NULL); //TODO: change it
#else 
	/// Using boost in Linux.
	BOOST_FOREACH(ptree::value_type const& node, pt.get_child(""))
	{
		if (node.first == "PrecTol")
		{
			Def.RelTol    = node.second.get<double>("RelTol");
			Def.minAbsTol = node.second.get<double>("minAbsTol");
			Def.minRelTol = node.second.get<double>("minRelTol");
			neps          = node.second.get<double>("NEPS");
			User = Def;
			MaxIter   = node.second.get<unsigned>("MaxIter");
			MaxRecur  = node.second.get<unsigned>("MaxRecur");
		}
	}
#ifndef NDEBUG
	std::cout << "RelTol, minAbsTol and minRelTol are "
		<< Def.RelTol << ", " << Def.minAbsTol << ", " << Def.minRelTol << std::endl;
#endif
#endif
}
/*
template<class T>
inline T GetSimOptFile::getElement(std::string eleName_)
{
	if (!strcmp(eleName_.c_str(), "RelTol"))
		return RelTol;
	else if (!strcmp(eleName_.c_str(), "minRelTol"))
		return minRelTol;
	else if (!strcmp(eleName_.c_str(), "minAbsTol"))
		return minAbsTol;
	else if (!strcmp(eleName_.c_str(), "NEPS"))
		return neps;
	else if (!strcmp(eleName_.c_str(), "MaxIter"))
		return MaxIter;
	else {
		std::string errorMsg = std::string(", could not identify the elemnet ") + eleName_;
		CATCHRTERROR(errorMsg);
	}
}*/

void GetSimOptFile::parseTolerance(std::vector<double>& tol_)
{
	if (tol_.empty())
		tol_.push_back(User.RelTol);

	std::vector<double>::iterator it = tol_.begin();
	switch (tol_.size())
	{
	case 1:
		it = tol_.insert(it, 0.0);
		it = tol_.insert(it, 0.0);
		break;
	case 2:
		tol_.push_back(0.0);
		break;
	case 3:
		break;
	default:
		std::cerr << "Warning: 'mkivpp.ntol','Expecting 1-3 vector tolerance'\n";
	}
	tol_.at(0) = pvMAX(tol_.at(0), User.minAbsTol);
	tol_.at(1) = pvMAX(tol_.at(1), User.minAbsTol);
	tol_.at(2) = pvMAX(tol_.at(2), User.minRelTol);

	ASSERT(((tol_.at(0) >= 0) && (tol_.at(1) >= 0) && (tol_.at(2) >= 0)), "GetSimOptFile::parseTolerance, Tolerances must be real-positive");

	User = Def; /// Reset the value of User. TODO: implicitly-defined assignment operator?
#ifndef NDEBUG
	std::cout << "the tolerance is " << tol_.at(0) << ", "
		<< tol_.at(1) << ", " << tol_.at(2) << std::endl;
#endif
} 
