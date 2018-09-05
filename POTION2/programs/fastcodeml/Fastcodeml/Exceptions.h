
#ifndef FASTCODEMLEXCEPTIONS_H
#define FASTCODEMLEXCEPTIONS_H

#include <stdexcept>
#include <sstream>

/// Fatal error in FastCodeML.
/// The message explains the reason.
///
class FastCodeMLFatal : public std::runtime_error
{
public:
	/// Exception with no message because it has already been printed.
	///
	FastCodeMLFatal(void) : runtime_error("")
	{}

	/// Exception that prints a message before termination.
	///
	/// @param[in] aMessage The message to be printed before termination
	///
	FastCodeMLFatal(const char *aMessage) : runtime_error(aMessage)
	{}

	/// Exception that prints a message (previously built into a std::ostringstream) before termination.
	///
	/// @param[in] aMessage The message to be printed before termination (has been formatted by printing into a std::ostringstream).
	///
	FastCodeMLFatal(const std::ostringstream& aMessage) : runtime_error(aMessage.str().c_str())
	{}
};


/// Early successful termination exception.
/// It does not signal a fatal error, but something like having printed the help.
///
class FastCodeMLSuccess : public std::exception
{
public:
	/// Early successful termination exception.
	///
	FastCodeMLSuccess() : exception()
	{}
};


/// Early termination due to LRT not satisfied exception.
/// It does not signal a fatal error, but that the optimization stopped because LRT cannot be satisfied.
///
class FastCodeMLEarlyStopLRT : public std::exception
{
public:
	/// Early successful termination exception.
	///
	FastCodeMLEarlyStopLRT() : exception()
	{}
};


/// Memory error in FastCodeML.
/// The message explains the reason in detail.
///
class FastCodeMLMemoryError : public std::bad_alloc
{
public:
	/// Exception that prints the standard bad_alloc message.
	///
	FastCodeMLMemoryError(void) : bad_alloc(), mMessage("std::bad_alloc")
	{}

	/// Exception that prints a message before termination.
	///
	/// @param[in] aMessage The message to be printed before termination
	///
	FastCodeMLMemoryError(const char *aMessage) : bad_alloc(), mMessage(aMessage)
	{}

	/// Override what() to display a more meaningful message.
	///
	/// @return The message text
	///
	const char* what() const throw()
	{
		return mMessage;
	}

private:
	const char* mMessage;	///< The exception message
};

#endif

