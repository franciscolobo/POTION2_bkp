
#ifndef PARSEPARAMETERS_H
#define PARSEPARAMETERS_H

#include <map>
#include <fstream>
#include <string>

/// Parse and accumulate initial values for parameters.
/// New values for the parameters came in the form name=value or name:value
///
///   @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///   @date 2012-02-21 (initial version)
///   @version 1.0
///
class ParseParameters
{
public:
	/// Return a pointer to the singleton instance
	///
	/// @return The pointer to the instance
	///
	static ParseParameters* getInstance(void);

	/// Change value to a given parameter.
	/// The argument is a string: name=value or name:value where the name is in (w0, k, p0, p1, w2) and the value is a double numerical value
	///
	/// @param[in] aParamValuePair The string to be parsed for name and value
	///
	/// @exception FastCodeMLFatal For invalid or malformed string
	///
	void addParameter(const char* aParamValuePair);

	/// Get the value associated to the given name
	///
	/// @param[in] aParamName The name of the parameter to be retrieved
	///
	/// @exception FastCodeMLFatal For not existent parameter name
	///
	double getParameter(const char* aParamName) const;

	/// Print the parameter list as: cout << r;
	///
	/// @param[in] aOut Output stream
	/// @param[in] aParamsList Pointer to the instance to be printed
	///
	/// @return The output stream
	///
	friend std::ostream& operator<< (std::ostream& aOut, const ParseParameters* aParamsList);

protected:
	/// Protected constructor
	///
	ParseParameters();


private:
	static ParseParameters*			mInstance;					///< Pointer to the singleton instance
	std::map<std::string, double>	mDictionary;				///< Dictionary holding the pairs (parameter name, value)
};

#endif

