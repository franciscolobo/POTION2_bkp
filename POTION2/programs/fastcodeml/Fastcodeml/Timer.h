
#ifndef TIMER_H
#define TIMER_H

#ifdef _MSC_VER

// Define this if you want to enable the high resolution timer on Windows
#define USE_WIN_MSEC_TIMER

#ifdef USE_WIN_MSEC_TIMER
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <time.h>
#endif

/// Simple timer.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-31 (initial version)
///     @version 1.0
///
class Timer
{
public:
	/// Constructor
	///
	Timer() : mDelta(0)
	{
#ifdef USE_WIN_MSEC_TIMER
		QueryPerformanceFrequency(&mFreq);
#endif
		start();
	}

	/// Start the timer
	///
	void start(void)
	{
#ifndef USE_WIN_MSEC_TIMER
		time(&mStartTime);
#else
		QueryPerformanceCounter(&mStartTime);
#endif
	}

	/// Stop the timer
	///
	/// @return The elapsed time in milliseconds
	///
	time_t stop(void)
	{
#ifndef USE_WIN_MSEC_TIMER
		mDelta = (time(NULL) - mStartTime)*1000L;
#else
		LARGE_INTEGER end_time_w;
		QueryPerformanceCounter(&end_time_w);

		mDelta = static_cast<time_t>(((end_time_w.QuadPart - mStartTime.QuadPart) * 1000)/mFreq.QuadPart);
#endif
		return mDelta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	///
	/// @return The elapsed time in milliseconds
	///
	time_t get(void) const
	{
		return mDelta;
	}

private:
#ifndef USE_WIN_MSEC_TIMER
	time_t			mStartTime;	///< The start time
#else
    LARGE_INTEGER	mFreq;		///< The timer frequency
    LARGE_INTEGER	mStartTime;	///< The start time
#endif
	time_t			mDelta;		///< The elapsed time in milliseconds
};

#else

#include <sys/time.h> // gettimeofday

/// Simple timer
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-08-31 (initial version)
///     @version 1.0
///
class Timer
{
public:
	/// Constructor
	///
	Timer() : mDelta(0) {start();}

	/// Start the timer
	///
	void start(void)
	{
		gettimeofday(&mStartTime, NULL);
	}

	/// Stop the timer
	///
	/// @return The elapsed time in milliseconds
	///
	time_t stop(void)
	{
		struct timeval end_time;
		gettimeofday(&end_time, NULL);

		mDelta  = end_time.tv_sec*1000000L+end_time.tv_usec;
		mDelta -= mStartTime.tv_sec*1000000L+mStartTime.tv_usec;
		mDelta /= 1000L;

		return mDelta;
	}

	/// Return the elapsed time (after a start/stop cycle)
	///
	/// @return The elapsed time in milliseconds
	///
	time_t get(void) const
	{
		return mDelta;
	}

private:
	struct timeval	mStartTime;	///< The start time
	time_t			mDelta;		///< The elapsed time in milliseconds
};

#endif
#endif
