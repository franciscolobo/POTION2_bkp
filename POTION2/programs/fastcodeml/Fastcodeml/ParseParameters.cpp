
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "ParseParameters.h"
#include "Exceptions.h"


ParseParameters* ParseParameters::mInstance = NULL;

ParseParameters* ParseParameters::getInstance(void)
{
	if(!mInstance) mInstance = new ParseParameters;
	return mInstance;
}


ParseParameters::ParseParameters()
{
	// Valid parameters with their default values
	mDictionary["w0"] = 0.4;
	mDictionary["k"]  = 0.2;
#ifdef USE_ORIGINAL_PROPORTIONS
	mDictionary["p0"] = 1.0;
	mDictionary["p1"] = 0.1;
#else
	mDictionary["p0"] = 0.25; // This way all the 4 proportions are equal to 0.25
	mDictionary["p1"] = 0.25;
#endif
	mDictionary["w2"] = 1.1;
}


void ParseParameters::addParameter(const char* aParamValuePair)
{
	// Parse string in the format: name=value (or name:value)
	const char* p=aParamValuePair;
	size_t i = 0;
	for(; *p != '=' && *p != ':' && *p != '\0'; ++p, ++i) {}
	if(*p == '\0') throw FastCodeMLFatal("Added pair with missing value in addParameter");
	if(i == 0) throw FastCodeMLFatal("Added pair with missing param name in addParameter");

	// Extract the name part that should be present in the dictionary
	std::string name(aParamValuePair, i);
	if(mDictionary.count(name) < 1) throw FastCodeMLFatal("Try to add a non-existing key in addParameter");

	// Change the value corresponding to the parameter name
	mDictionary[name] = atof(p+1);
}

	
double ParseParameters::getParameter(const char* aParamName) const
{
	// The name should exist in the dictionary
	std::map<std::string, double>::const_iterator ip(mDictionary.find(aParamName));
	if(ip == mDictionary.end()) throw FastCodeMLFatal("Invalid key requested in getParameter");

	// Return the corresponding value
	return ip->second;
}


std::ostream& operator<< (std::ostream& aOut, const ParseParameters* aParamsList)
{
	// Print the dictionary
	std::map<std::string, double>::const_iterator ip(aParamsList->mDictionary.begin());
	const std::map<std::string, double>::const_iterator end(aParamsList->mDictionary.end());
	for(; ip != end; ++ip)
	{
		aOut << std::setw(20) << ip->first << " = " << std::fixed << std::setprecision(6) << ip->second << std::endl;
	}

	return aOut;
}



