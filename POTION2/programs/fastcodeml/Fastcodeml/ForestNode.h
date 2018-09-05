#ifndef FORESTNODE_H
#define FORESTNODE_H

#include <vector>
#include <iostream>
#include <climits>
#include <cstring>
#include <string>
#include "MatrixSize.h"
#include "AlignedMalloc.h"
#include "Exceptions.h"

/// Max number of children for a given node (in reality they are only 2)
static const int MAX_NUM_CHILDREN = 8;
static const unsigned char ALL_CHILDREN_SAME_TREE = 0xFF;

/// Support data needed only during forest preprocessing phase. It is deleted before the computation phase.
///
///   @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///   @date 2011-11-08 (initial version)
///   @version 1.0
///
struct ForestNodeSupport
{
	/// Constructor.
	/// Do nothing.
	///
	ForestNodeSupport() {}

	std::vector<long long>			mSubtreeCodonsSignature;	///< List of codon idx for the subtree rooted at this node
};

/// One node of each tree in the forest.
///
///   @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///   @date 2011-02-23 (initial version)
///   @version 1.0
///
struct ForestNode
{
	// Field order suggested by icc
	//'mChildrenSameTreeFlags, mBranchId, mOwnTree, mChildrenCount, mParent, mProb, mInternalNodeId, mChildrenList, mOtherTreeProb'
	unsigned char				mChildrenSameTreeFlags;		///< Bit i set if child i is in the same tree
	unsigned char				mChildrenCount;				///< Number of children of this node
	short						mLeafCodon;					///< On a leaf set to the corresponding codon number, otherwise -1
	unsigned int				mBranchId;					///< A unique index to access the branch length array (starts from zero at the first non-root node)
	unsigned int				mOwnTree;					///< Per tree identifier
	ForestNode*					mParent;					///< Pointer to the node parent (null for the root)
#ifndef NEW_LIKELIHOOD
	double*						mProb[Nt];					///< Codons probability array (called g in the pseudocode) (can be computed by concurrent tree traversals)
#endif
	unsigned int				mInternalNodeId;			///< Internal node identifier to mark a branch as foreground. UINT_MAX means not an internal node
	std::vector<ForestNode *>	mChildrenList;				///< List of the node children
	ForestNodeSupport*			mPreprocessingSupport;		///< Data needed only during forest preprocessing phase
#ifndef NEW_LIKELIHOOD
	std::vector<double *>		mOtherTreeProb;				///< Pointers to other tree precomputed mProb, zero if not used, or local array if used from other tree
#endif
#ifdef NON_RECURSIVE_VISIT
	bool						mFirstChild;
	unsigned int				mChildIdx;					///< Mark the child position in the parent node
#endif

	/// Constructor.
	///
	ForestNode() : mChildrenSameTreeFlags(ALL_CHILDREN_SAME_TREE), mChildrenCount(0), mLeafCodon(-1), mBranchId(0), mOwnTree(0), mParent(NULL), mInternalNodeId(0)
#ifdef NON_RECURSIVE_VISIT
				   , mFirstChild(false), mChildIdx(0)
#endif		  
	{
#ifndef NEW_LIKELIHOOD
		memset(mProb, 0, Nt*sizeof(double*));
		mOtherTreeProb.reserve(2);
#endif
		mChildrenList.reserve(2);
		mPreprocessingSupport = new ForestNodeSupport;
	}

	/// Destructor.
	///
	~ForestNode()
	{
		// Delete children if in the same tree. Delete partial Prob arrays if not pointer to other tree partial Prob array
		const unsigned int nc = mChildrenCount;
		for(unsigned int i=0; i < nc; ++i)
		{
			if(isSameTree(i))
			{
				delete mChildrenList[i];

#ifndef NEW_LIKELIHOOD
				alignedFree(mOtherTreeProb[i]);
#endif
			}
		}

		// Clean all arrays
		mChildrenList.clear();
#ifndef NEW_LIKELIHOOD
		mOtherTreeProb.clear();
#endif
		// Delete preprocessing support data
		delete mPreprocessingSupport;
	}

	/// Copy constructor.
	///
	/// @param[in] aNode Node that has to be assigned to the current node
	///
	ForestNode(const ForestNode& aNode)
		: mChildrenSameTreeFlags(aNode.mChildrenSameTreeFlags),
		  mChildrenCount(aNode.mChildrenCount), mLeafCodon(aNode.mLeafCodon), mBranchId(aNode.mBranchId), mOwnTree(aNode.mOwnTree),
		  mParent(aNode.mParent), mInternalNodeId(aNode.mInternalNodeId), mChildrenList(aNode.mChildrenList)
#ifdef NON_RECURSIVE_VISIT
		  , mFirstChild(aNode.mFirstChild), mChildIdx(aNode.mChildIdx)
#endif
#ifndef NEW_LIKELIHOOD
		, mOtherTreeProb(aNode.mOtherTreeProb)
#endif
	{
#ifndef NEW_LIKELIHOOD
		memcpy(mProb, aNode.mProb, Nt*sizeof(double*));
#endif
		mPreprocessingSupport = new ForestNodeSupport;
		mPreprocessingSupport->mSubtreeCodonsSignature = aNode.mPreprocessingSupport->mSubtreeCodonsSignature;
	}

	/// Assignment operator.
	///
	/// @param[in] aNode Node that has to be assigned to the current node
	///
	/// @return The node itself
	///
	ForestNode& operator=(const ForestNode& aNode)
	{
		// Make sure not same object
		if(this != &aNode)
		{
			mChildrenList			= aNode.mChildrenList;
			mParent					= aNode.mParent;
			mLeafCodon				= aNode.mLeafCodon;
#ifndef NEW_LIKELIHOOD
			memcpy(mProb, aNode.mProb, Nt*sizeof(double*));
#endif
			mInternalNodeId			= aNode.mInternalNodeId;
			mBranchId				= aNode.mBranchId;
			mOwnTree				= aNode.mOwnTree;
#ifndef NEW_LIKELIHOOD
			mOtherTreeProb			= aNode.mOtherTreeProb;
#endif
			mChildrenSameTreeFlags	= aNode.mChildrenSameTreeFlags;
			mChildrenCount			= aNode.mChildrenCount;

			mPreprocessingSupport = new ForestNodeSupport;
			mPreprocessingSupport->mSubtreeCodonsSignature = aNode.mPreprocessingSupport->mSubtreeCodonsSignature;
#ifdef NON_RECURSIVE_VISIT
			mFirstChild				= aNode.mFirstChild;
			mChildIdx				= aNode.mChildIdx;
#endif
		}

		// Return ref for multiple assignment
		return *this;
	}

	/// Build a node aligned to a 64 bits boundary (cache line size).
	///
	/// @param[in] aSize The size to be allocated
	///
	/// @return Pointer to the allocated memory area
	///
	/// @exception FastCodeMLMemoryError If no memory available
	///
	void* operator new(size_t aSize)
	{
		void *m = alignedMalloc(aSize, CACHE_LINE_ALIGN);
		if(!m) throw FastCodeMLMemoryError("Error in ForestNode allocation");
		return m;
	}

	/// Release the node
	///
	/// @param[in] aPtr Pointer to the memory area to be released
	///
	void operator delete(void* aPtr)
	{
		alignedFree(aPtr);
	}

#ifndef __MTA__
	/// Placement new required by PGI compiler
	///
	/// @fn void* operator new(std::size_t aSize, ForestNode* aHere)
	///
	/// @param[in] aSize Requested size (ignored)
	/// @param[in] aHere Where the placement new should go
	///
	/// @return The placed memory
	///
	void* operator new(std::size_t /* aSize */, ForestNode* aHere)
	{
		return aHere;
	}

	/// Placement delete required by PGI compiler
	///
	/// @fn void operator delete(void* aPtr, ForestNode* aHere)
	///
	/// @param[in] aPtr Pointer to the memory area to be released (ignored)
	/// @param[in] aHere Where the placement new should go (ignored)
	///
	void operator delete(void* /* aPtr */, ForestNode* /* aHere */)
	{
		// Do nothing
	}
#else
	/// Placement new required by XMT compiler
	///
	/// @fn void* operator new(unsigned long aSize, void* aHere)
	///
	/// @param[in] aSize Requested size (ignored)
	/// @param[in] aHere Where the placement new should go
	///
	/// @return The placed memory
	///
	void* operator new(unsigned long /* aSize */, void* aHere)
	{
		return aHere;
	}

	/// Placement delete required by XMT compiler
	///
	/// @fn void operator delete(void* aPtr, void* aHere)
	///
	/// @param[in] aPtr Pointer to the memory area to be released (ignored)
	/// @param[in] aHere Where the placement new should go (ignored)
	///
	void operator delete(void* /* aPtr */, void* /* aHere */)
	{
		// Do nothing
	}
#endif

	/// Print from this node down
	///
	/// @param[in] aNodeNames The list of node names
	/// @param[in] aOut Output stream
	/// @param[in] aIndent Initial number of indent spaces
	/// @param[in] aIncrement The indent amount is incremented by this value at each level
	///
	void print(const std::vector<std::string>& aNodeNames, std::ostream& aOut=std::cout, unsigned int aIndent=0, unsigned int aIncrement=3) const
	{
		unsigned int i;

		// Indent
		for(i=0; i < aIndent; ++i) aOut << ' ';

		// Print the name
		aOut << '<' << ((mBranchId  != UINT_MAX) ? aNodeNames[mBranchId+1] : aNodeNames[0]) << "> ";
	
		// Print the ID
		if(mInternalNodeId != UINT_MAX) aOut << '(' << mInternalNodeId << '|' << mBranchId << '|' << mLeafCodon << ") ";
		else                            aOut << '('                    << '|' << mBranchId << '|' << mLeafCodon << ") ";

		// Print the indexes of the codons accumulated till this node
		if(mPreprocessingSupport)
		{
			std::vector<long long>::const_iterator ig(mPreprocessingSupport->mSubtreeCodonsSignature.begin());
			const std::vector<long long>::const_iterator end(mPreprocessingSupport->mSubtreeCodonsSignature.end());
			for(; ig != end; ++ig) aOut << *ig << ' ';
			aOut << std::endl;
		}

		// Print the subtree
		std::vector<ForestNode*>::const_iterator irn(mChildrenList.begin());
		const std::vector<ForestNode*>::const_iterator end(mChildrenList.end());
		for(i=0; irn != end; ++irn, ++i)
		{
			// If the subtree is on the same tree, then print it, otherwise print only the subtree root node name.
			if(isSameTree(i))
			{
				(*irn)->print(aNodeNames, aOut, aIndent+aIncrement, aIncrement);
			}
			else
			{
				for(i=0; i < aIndent+aIncrement; ++i) aOut << ' ';
				i = (*irn)->mBranchId+1;
				aOut << '[' << aNodeNames[i] << ']' << std::endl;
			}
		}
	}

	/// Create a list of pointers to leaves.
	///
	/// @param[out] aLeafsList Pointers to leaves are pushed to this vector
	///
	void pushLeaf(std::vector<ForestNode*>& aLeafsList)
	{
		if(mChildrenList.empty())
		{
			aLeafsList.push_back(this);
		}
		else
		{
			std::vector<ForestNode*>::const_iterator irn(mChildrenList.begin());
			const std::vector<ForestNode*>::const_iterator end(mChildrenList.end());
			for(; irn != end; ++irn) (*irn)->pushLeaf(aLeafsList);
		}
	}

	/// Fills the mPreprocessingSupport->mSubtreeCodonsSignature list with the ordered union of its children's lists.
	///
	void gatherCodons(void)
	{
		std::vector<ForestNode*>::const_iterator irn(mChildrenList.begin());
		const std::vector<ForestNode*>::const_iterator end(mChildrenList.end());
		for(; irn != end; ++irn)
		{
			(*irn)->gatherCodons();
			mPreprocessingSupport->mSubtreeCodonsSignature.insert(mPreprocessingSupport->mSubtreeCodonsSignature.end(), (*irn)->mPreprocessingSupport->mSubtreeCodonsSignature.begin(), (*irn)->mPreprocessingSupport->mSubtreeCodonsSignature.end());
		}
	}

	/// Count the total branches in the forest
	///
	/// @param[in] aAggressiveStrategy If true use the aggressive simplification strategy
	///
	/// @return The total number of branches of the forest
	///
	unsigned int countBranches(bool aAggressiveStrategy=false) const
	{
		unsigned int cnt = 0;
		unsigned int i;

		// Visit the subtrees
		std::vector<ForestNode*>::const_iterator irn(mChildrenList.begin());
		const std::vector<ForestNode*>::const_iterator end(mChildrenList.end());
		for(i=0; irn != end; ++irn, ++i)
		{
			// If the subtree is on the same tree, then print it, otherwise print only the subtree root node name.
			if(isSameTree(i))
			{
				cnt += (*irn)->countBranches(aAggressiveStrategy)+1;
			}
			else if(!aAggressiveStrategy)
			{
				cnt += 1;
			}
		}
	
		return cnt;
	}

	/// Get the computation cost from the node to the leaves.
	///
	/// @param[in] aCostAtLeaf The cost associated to a leaf
	/// @param[in] aCostIntern The cost associated to an internal node
	/// @param[in] aCostPtr The cost associated to a pointer to another node
	///
	/// @return The total cost of the tree starting from the given node
	///
	unsigned int getCost(unsigned int aCostAtLeaf, unsigned int aCostIntern, unsigned int aCostPtr) const
	{
		unsigned int cost = 0;
		unsigned int i;

		// Visit the subtrees
		std::vector<ForestNode*>::const_iterator irn(mChildrenList.begin());
		const std::vector<ForestNode*>::const_iterator end(mChildrenList.end());
		for(i=0; irn != end; ++irn, ++i)
		{
			// If the subtree is on the same tree, then print it, otherwise print only the subtree root node name.
			if(isSameTree(i))
			{

				if((*irn)->mLeafCodon >= 0)
				{
					cost += aCostAtLeaf;
				}
				else
				{
					cost += (*irn)->getCost(aCostAtLeaf, aCostIntern, aCostPtr)+aCostIntern;
				}
			}
			else
			{
				cost += aCostPtr;
			}
		}
	
		return cost;
	}

	/// Bitmask for the mChildrenSameTreeFlags bitset
	///
	static const unsigned char mMaskTable[MAX_NUM_CHILDREN];

	/// Mark child aChildIndex as not in the same tree (Reset the given flag to false)
	///
	/// @param[in] aChildIndex The index of the flag to be set to false
	///
	void markNotSameTree(unsigned int aChildIndex)
	{
		mChildrenSameTreeFlags &= static_cast<unsigned char>(~mMaskTable[aChildIndex]);
	}

	/// Test if child aChildIndex is in the same tree.
	///
	/// @param[in] aChildIndex The index of the flag to be tested
	///
	/// @return The flag status
	///
	bool isSameTree(unsigned int aChildIndex) const
	{
		return (mChildrenSameTreeFlags & mMaskTable[aChildIndex]) != 0;
	}

	/// Set all flags to true (i.e.\ all children of the node are in the same tree)
	///
	void setAllFlagsSameTree(void)
	{
		mChildrenSameTreeFlags = ALL_CHILDREN_SAME_TREE;
	}
};

#endif

