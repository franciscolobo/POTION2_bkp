
#ifndef HIGHLEVELCOORDINATOR_H
#define HIGHLEVELCOORDINATOR_H

#include "Forest.h"
#include "CmdLine.h"
#include "WriteResults.h"

/// The rank of the master job
///
static const int MASTER_JOB = 0;

/// Coordinator class for high level parallelization.
/// This class encapsulates MPI usage to parallelize FastCodeML above the maximizer level.
///
/// @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
/// @date 2011-11-22 (initial version)
/// @version 1.0
///

class HighLevelCoordinator
{
public:
	/// Constructor.
	///
	/// @param[in,out] aRgc Pointer to the number of arguments
	/// @param[in,out] aRgv Pointer to the arguments' list
	///
	/// @exception FastCodeMLFatal MPI Failed to initialize
	/// @exception FastCodeMLSuccess To terminate unused worker processes
	///
	HighLevelCoordinator(int* aRgc, char*** aRgv);

	/// Destructor.
	///
	~HighLevelCoordinator();

	/// Starts the high level parallelization of the FastCodeML application
	///
	/// @param[in,out] aForest The filled forest
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	/// @return True if the execution can go parallel at this level.
	///
	bool startWork(Forest& aForest, const CmdLine& aCmdLine);

	/// Is this process the master one?
	///
	/// @return True if this is the master process
	///
	bool isMaster(void) const {return mRank == MASTER_JOB;}

	/// Return the number of MPI processes.
	///
	/// @return The number of MPI processes
	///
	int  numJobs(void) const {return mSize;}

	/// Return the current process rank.
	///
	/// @return The current process rank.
	///
	int getRank(void) const {return mRank;}

private:
	/// The master coordination job
	///
	/// @param[in] aOutputResults To collect and output results to a results file
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	/// @exception FastCodeMLFatal Invalid job request found
	///
	void doMaster(WriteResults& aOutputResults, const CmdLine& aCmdLine);

	/// The worker high level loop
	///
	/// @param[in,out] aForest The filled forest
	/// @param[in] aCmdLine The parameters from the command line of the main program
	///
	void doWorker(Forest& aForest, const CmdLine& aCmdLine);


private:
	unsigned int		mVerbose;				///< The verbose level
	int					mRank;					///< Rank of the current process (Master has rank == MASTER_JOB)
	int					mSize;					///< Number of MPI processes
	size_t				mNumInternalBranches;	///< Number of internal branches (i.e.\ the ones that can be foreground branch)

	struct WorkTable;
	WorkTable*			mWorkTable;				///< Management of the work list
};


#endif
