
#ifndef CMDLINE_H
#define CMDLINE_H

#include <climits>
#include "VerbosityLevels.h"

/// The default maximum number of optimization steps.
///
static const unsigned int MAX_ITERATIONS=10000;


/// Parse the command line flags.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2010-12-22 (initial version)
///     @version 1.0
///
///
class CmdLine
{
public:
	/// Constructor.
	///
	/// Here are set the default values for the command line settable parameters
	///
	CmdLine() : 
		mDeltaValueForGradient(0.0),
		mRelativeError(1e-3),
		mTreeFile(NULL),
		mGeneFile(NULL),
		mGraphFile(NULL),
		mResultsFile(NULL),
		mVerboseLevel(VERBOSE_ONLY_RESULTS),
		mSeed(0),
		mBranchStart(UINT_MAX),
		mBranchEnd(UINT_MAX),
		mExportComputedTimes(UINT_MAX),
		mComputeHypothesis(UINT_MAX),
		mOptimizationAlgo(0),
		mExtraDebug(0),
		mMaxIterations(MAX_ITERATIONS),
		mIgnoreFreq(false),
		mDoNotReduceForest(false),
		mBranchLengthsFromFile(false),
		mNoMaximization(false),
		mTrace(false),
		mForceSerial(false),
		mBranchFromFile(false),
		mInitH0fromH1(false),
		mInitFromParams(false),
		mCleanData(false),
		mStopIfNotLRT(true),
		mCmdLineImpl(NULL)
	{}

	/// Destructor.
	///
	~CmdLine();

	/// Parse the command line.
	///
	/// @param[in] aCnt Number of command line parameters (from main argc)
	/// @param[in] aVal Array of command line parameters (from main argv)
	///
	/// @exception FastCodeMLFatal For invalid command line parameters.
	/// @exception FastCodeMLSuccess After showing help to exit the application
	///
	void parseCmdLine(int aCnt, char **aVal);


public:
	double			mDeltaValueForGradient;	///< The variable increment to compute gradient (zero means use a hardcoded default value)
	double			mRelativeError;			///< Relative error to stop maximization
	const char*		mTreeFile;				///< %Newick tree file name
	const char*		mGeneFile;				///< %Genes file name
	const char*		mGraphFile;				///< If not null export the forest to this file in GML format to be visualized using R igraph package or yEd editor
	const char*		mResultsFile;			///< File to which the results should be written
	unsigned int	mVerboseLevel;			///< Verbosity level. 0: no messages; 1: basic messages; 2: messages useful for debugging; 3: really annoying
	unsigned int	mSeed;					///< Random number generator seed (0 means not set from command line)
	unsigned int	mBranchStart;			///< Compute from this branch. The numbering starts at 0 (UINT_MAX means run all branches)
	unsigned int	mBranchEnd;				///< Compute to this branch. The numbering starts at 0 (UINT_MAX means run all branches). It is >= mBranchStart
	unsigned int	mExportComputedTimes;	///< If 0 or 1 the times exported are the computed ones in H0 or H1, otherwise (UINT_MAX) the one read from file
	unsigned int	mComputeHypothesis;		///< If set to 0 compute only H0, if set to 1 compute only H1, otherwise compute both
	unsigned int	mOptimizationAlgo;		///< Select the optimization algorithm to use
	unsigned int	mExtraDebug;			///< Extra debug parameter for development tests
	unsigned int	mMaxIterations;			///< Maximum number of iterations for the maximization
	bool			mIgnoreFreq;			///< Ignore the computed codon frequencies and set them all to 1/61
	bool			mDoNotReduceForest;		///< If true do not reduce the forest merging common subtrees
	bool			mBranchLengthsFromFile;	///< The initial value of the branch lengths is taken from the phylo tree file
	bool			mNoMaximization;		///< Only the first step of the likelihood maximization is taken
	bool			mTrace;					///< Trace the optimization steps
	bool			mForceSerial;			///< Disable all parallelism
	bool			mBranchFromFile;		///< Read the foreground branch to use from the phylo tree file (it is marked as #1)
	bool			mInitH0fromH1;			///< If set starts the H0 computation from the H1 results
	bool			mInitFromParams;		///< Initialize times from phylo tree and the other from values hardcoded or entered on the command line
	bool			mCleanData;				///< Remove ambiguous or missing sites from the MSA (genes)
	bool			mStopIfNotLRT;			///< Stop H0 maximization when LRT cannot be satisfied

private:
	struct CmdLineImpl;
	CmdLineImpl* mCmdLineImpl;				///< Implementation so this structure could be used without other modules be aware of its internals
};

#endif

