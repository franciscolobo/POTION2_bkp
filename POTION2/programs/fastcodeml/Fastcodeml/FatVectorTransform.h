
#ifndef FATVECTORTRANSFORM_H
#define FATVECTORTRANSFORM_H

#include <vector>
#include <utility>
#include "ForestNode.h"
#include "Types.h"

/// Manipulations on the per-branch probability vector array.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-09-31 (initial version)
///     @version 1.0
///
class FatVectorTransform
{
public:
	/// Constructor.
	///
	FatVectorTransform() : mNumSites(0), mNumBranches(0), mNoTransformations(true) {}

	/// Destructor.
	///
	~FatVectorTransform()
	{
		mNodeStatus.clear();
		mLimits.clear();
		mCopyCmds.clear();
		mReuseCmds.clear();
		mFirstForLevel.clear();
		mBranchByLevel.clear();
	}

	/// Create the dependency list between levels in the tree
	///
	/// @param[in] aNodesByLevel List of lists of pointers to nodes one for each level
	///
	void setBranchDependencies(const std::vector< std::vector<ForestNode*> >& aNodesByLevel);

	/// Initialize the class instance
	///
	/// @param[in] aNumBranches Number of branches
	/// @param[in] aNumSites Number of sites
	///
	void initNodeStatus(size_t aNumBranches, size_t aNumSites)
	{
		mNumBranches = aNumBranches;
		mNumSites = aNumSites;
		mNodeStatus.assign(aNumBranches*aNumSites, FatVectorTransform::SITE_NOT_EXISTS);
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites));
		mCopyCmds.clear();
		mReuseCmds.clear();
		mNoTransformations = false;
	}
	
	/// Initialize the class instance so can be used if no subtree pruning is present
	///
	/// @param[in] aNumBranches Number of branches
	/// @param[in] aNumSites Number of sites
	///
	void initNodeStatusMinimal(size_t aNumBranches, size_t aNumSites)
	{
		mNumBranches = aNumBranches;
		mNumSites = aNumSites;
		mLimits.assign(aNumBranches, std::make_pair(0, aNumSites));
		mCopyCmds.clear();
		mReuseCmds.clear();
		mNoTransformations = true;
	}

	/// Set the corresponding node as existing for the given site
	///
	/// @param[in] aBranch Branch for which the node has been verified as existing
	/// @param[in] aSite Site for which the node has been verified as existing
	///
	void setNodeExists(unsigned int aBranch, unsigned int aSite)
	{
		mNodeStatus[aBranch * mNumSites + aSite] = FatVectorTransform::SITE_EXISTS;
	}

	/// Set the corresponding node as taking its value from another site
	///
	/// @param[in] aBranch Branch for which the node has been checked
	/// @param[in] aSite Site for which the node has been checked
	/// @param[in] aReusedSite Site from which the node takes its value
	///
	void setNodeReuses(unsigned int aBranch, unsigned int aSite, unsigned int aReusedSite)
	{
		mNodeStatus[aBranch * mNumSites + aSite] = aReusedSite;
	}

	/// Prints (on cout) for each branch the first and last valid positions and the valid entries in this range.
	///
	/// @exception FastCodeMLFatal No SITE_EXISTS in mNodePresent at branch
	///
	void printCountGoodElements(void) const;

	/// Prints (on cout) the visit sequence of branches.
	///
	void printBranchVisitSequence(void) const;

	/// Prints (on cout) for each branch and each site if it is valid, if it is not present and if takes the value from another site
	///
	void printNodeStatus(void) const;

	/// Compute the commands needed to compact the various fat vectors (matrices, one for each branch) of probability vectors
	///
	/// @exception FastCodeMLFatal No SITE_EXISTS in mNodePresent at branch
	///
	void compactMatrix(void);

	/// Print the lists of generated commands
	///
	void printCommands(void) const;

	/// Get the first index to be used for computation
	///
	/// @param[in] aBranch Specify which branch should be returned.
	///
	/// @return The starting index in the fat vector
	///
	size_t getLowerIndex(unsigned int aBranch) const {return mLimits[aBranch].first;}

	/// Get the number of items to be used for computation
	///
	/// @param[in] aBranch Specify which branch should be returned.
	///
	/// @return The count of sites to be used
	///
	size_t getCount(unsigned int aBranch) const {return mLimits[aBranch].second;}

	/// Compact the fat vector at the leaves.
	///
	/// @param[in,out] aProbs The fat probability vector that will be changed at the leaves level
	///
	void preCompactLeaves(CacheAlignedDoubleVector& aProbs);

	/// Compact the fat vector at a certain level in the tree
	///
	void postCompact(CacheAlignedDoubleVector& aStepResults, CacheAlignedDoubleVector& aProbs, unsigned int aLevel, unsigned int aNumSets);


private:
	size_t				mNumSites;				///< The number of valid sites.
	size_t				mNumBranches;			///< The number of branches.
	std::vector<int>	mNodeStatus;			///< For each (Branch, Site) (idx = branch*NumSites+site) the values are as in BranchSitePositionStatus enum

	/// The values that mNodePresent array could take
	///
	enum BranchSitePositionStatus
	{
		SITE_EXISTS     = -2,					///< The position (Branch, Site) in mNodePresent exists
		SITE_NOT_EXISTS = -1,					///< The position (Branch, Site) in mNodePresent refers to a not existent node
		SITE_FIRST_NUM  =  0					///< if greater or equal to this value the position contains the index from which the value should be copied
	};

	/// Representation of a range to be copied and the number of items to be copied
	///
	struct Range
	{
		/// Constructor
		///
		Range(unsigned int aFrom, unsigned int aTo) {from = aFrom; to = aTo; cnt = 1;}

		unsigned int from;		///< Starting index from which to copy
		unsigned int to;		///< Starting index to which the values should be copied
		unsigned int cnt;		///< How many items (if zero, skip this entry)

		//bool operator<(Range& rhs) { return from < rhs.from; } ///< This is needed for sorting
	};

	/// Representation of the copy of one item
	///
	struct RangeNoCnt
	{
		/// Constructor
		///
		RangeNoCnt(unsigned int aFrom, unsigned int aTo) {from = aFrom; to = aTo;}

		unsigned int from;		///< Index from which to copy
		unsigned int to;		///< Index to which the value should be copied
	};

	typedef std::vector<Range> VectorOfRanges;									///< Vector of ranges (from position, to position, number of items)
	typedef std::vector<VectorOfRanges> VectorOfVectorOfRanges;					///< Vector of vectors of ranges
	typedef std::vector<RangeNoCnt> VectorOfRangesNoCnt;						///< Vector of single item copy (from position, to position)
	typedef std::vector<VectorOfRangesNoCnt> VectorOfVectorOfRangesNoCnt;		///< Vector of vectors of single item copies
	typedef std::vector< std::pair<size_t, size_t> > VectorOfPairs;				///< Vector of pairs of values

	VectorOfPairs						mLimits;				///< Lower index and total count for each branch
	VectorOfVectorOfRanges				mCopyCmds;				///< Ranges to be copied to fill the holes (one list for each branch)
	VectorOfVectorOfRangesNoCnt			mReuseCmds;				///< Ranges to be reused copying the computed value (one list for each branch)
	std::vector<bool>					mFirstForLevel;			///< One entry for branch set to true if it is the first entry for its level
	bool								mNoTransformations;		///< If set no transformation will take place (corresponds to no tree prune case)
	std::vector< std::vector<unsigned int> >
										mBranchByLevel;			///< Each level contains a list of branch numbers at this level. List start from the leaves.
	std::vector<unsigned int>			mParentNode;			///< One entry for branch set to the parent node index
};

#endif

