
#ifndef FOREST_H
#define FOREST_H

#include <vector>
#include "PhyloTree.h"
#include "Genes.h"
#include "ForestNode.h"
#include "TransitionMatrix.h"
#include "ProbabilityMatrixSet.h"
#include "MatrixSize.h"
#ifdef NEW_LIKELIHOOD
#include "FatVectorTransform.h"
#endif
#include "CodonFrequencies.h"
#include "Types.h"
#include "ForestExport.h"
#ifdef USE_DAG
#include "DAGScheduler.h"
#endif
#include "TreeAndSetsDependencies.h"
#include "CmdLine.h"


/// The phylogenetic tree's forest.
/// This class encapsulates the forest of phylogenetic tree that will be used for computing the tree's maximum likelihood
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2011-02-23 (initial version)
///  @version 1.0
///
class Forest
{
public:
	/// Constructor
	///
	/// @param[in] aVerbose The verbosity level
	///
	explicit Forest(unsigned int aVerbose=0)
							: mNumSites(0), mCodonFreq(NULL), mInvCodonFreq(NULL), mInv2CodonFreq(NULL),
							  mNumBranches(0), mNumInternalBranches(0), mMarkedInternalBranch(UINT_MAX), mVerbose(aVerbose) {}

	/// Destructor
	///
	~Forest()
	{
		mRoots.clear();
		mNodeNames.clear();
		mBranchLengths.clear();
		mProbs.clear();
		mSiteMultiplicity.clear();		
		mTableInternalToBranchID.clear();
#ifdef NEW_LIKELIHOOD
		mProbsOut.clear();
		mNodesByLevel.clear();
#endif
#ifdef NON_RECURSIVE_VISIT
		mVisitTree.clear();
		mVisitTreeParents.clear();
#endif
	}
	
	/// Build the forest and reduces the subtrees.
	///
	/// @param[in] aTree The phylogenetic tree
	/// @param[in] aGenes The corresponding genes
	/// @param[in] aCodonFrequencyModel Model to be used to compute the codon empirical frequencies.
	///
	/// @exception FastCodeMLFatal Invalid codon found
	///
	void loadTreeAndGenes(const PhyloTree& aTree,
						  const Genes& aGenes,
						  CodonFrequencies::CodonFrequencyModel aCodonFrequencyModel);

	/// Print the class statistics as: cout << r;.
	///
	/// @param[in] aOut Output stream
	/// @param[in] aForest The forest to be printed
	///
	/// @return The output stream
	///
	friend std::ostream& operator<< (std::ostream& aOut, const Forest& aForest);

	/// Get the first and last branches to be marked as foreground.
	///
	/// @param[in] aCmdLine The parameters from the command line of the main program
	/// @param[out] aBranchStart The first branch to be marked as foreground
	/// @param[out] aBranchEnd The last branch to be marked as foreground
	///
	/// @return True if all branches are selected
	///
	/// @exception FastCodeMLFatal Invalid range from command line
	///
	bool getBranchRange(const CmdLine& aCmdLine, size_t& aBranchStart, size_t& aBranchEnd) const;

	/// Reduce common subtrees on the whole forest.
	///
	void reduceSubtrees(void);

#ifndef NEW_LIKELIHOOD
	/// Add more aggressive subtree reduction.
	///
	/// @param[in] aNode The tree node from which the walker should start (no argument starts from the root)
	///
	/// @exception FastCodeMLMemoryError Cannot allocate mOtherTreeProb
	///
	void addAggressiveReduction(ForestNode* aNode=NULL);
#endif

	/// Remove all work data used for reduction.
	///
	/// @param[in] aNode The node from which to start. Pass zero to start from the root of all the trees in the forest.
	///
	void cleanReductionWorkingData(ForestNode* aNode=NULL);

#if !defined(NON_RECURSIVE_VISIT) && !defined(NEW_LIKELIHOOD)
	/// Compute likelihood visiting the trees in a non-recursive way
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aDependencies The dependency list between sets of trees
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, const ListDependencies& aDependencies);
#endif

#ifdef NON_RECURSIVE_VISIT
	/// Prepare the list of threading pointers for non-recursive trees visit
	///
	void prepareNonRecursiveVisit(void);

	/// Compute likelihood visiting the trees in a non-recursive way
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1)
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int aHyp);
#endif

#ifdef NEW_LIKELIHOOD
	/// Compute the log likelihood of the forest given the set of precomputed matrices.
	/// If NEW_LIKELIHOOD is defined, this routine adopts the experimental "Long Vector" approach.
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[out] aLikelihoods Values of the codon probabilities at the tree root (one set for each set of matrices)
	/// @param[in] aHyp The hypothesis to be computed (H0: 0; H1: 1) (currently ignored)
	///
	void computeLikelihoods(const ProbabilityMatrixSet& aSet, CacheAlignedDoubleVector& aLikelihoods, unsigned int /*aHyp*/);
#endif

	/// Export the forest as graph file
	///
	friend class ForestExport;

	/// Return the total number of branches
	///
	/// @return The total number of branches
	///
	size_t getNumBranches(void) const {return mNumBranches;}

	/// Return the number of internal branches (i.e.\ the ones that do not connect to leaves)
	///
	/// @return The number of internal branches
	///
	size_t getNumInternalBranches(void) const {return mNumInternalBranches;}

	/// Get the number of sites
	///
	/// @return The number of sites
	///
	size_t getNumSites(void) const {return mNumSites;}

	/// Get the marked internal branch
	///
	/// @return The internal branch index of the branch marked in the tree file. UINT_MAX otherwise.
	///
	size_t getMarkedInternalBranch(void) const {return mMarkedInternalBranch;}

	/// Get site multiplicity values.
	///
	/// @return Reference to the array of site multiplicities
	///
	const std::vector<double>& getSiteMultiplicity(void) const {return mSiteMultiplicity;}

	/// Set the times (i.e.\ the branch lengths) from the values read from the tree file
	///
	/// @param[out] aTimes The array with all the tree times
	/// @param[in] aNode The node from which to start (if zero starts from the root)
	///
	void setTimesFromLengths(std::vector<double>& aTimes, const ForestNode* aNode=NULL) const;

	/// Set the times (i.e.\ the branch lengths) on the tree from the values read from the times array
	///
	/// @param[in] aTimes The array with all the tree times
	/// @param[out] aNode The node from which to start (if zero starts from the root)
	///
	void setLengthsFromTimes(const std::vector<double>& aTimes, ForestNode* aNode=NULL);

	/// Change the internal branch identifier for the foreground branch into the corresponding internal branch index.
	///
	/// @param[in] aFgBranch Number of the foreground branch
	///
	/// @return The node index corresponding to the foreground branch
	///
	unsigned int adjustFgBranchIdx(size_t aFgBranch) const {return mTableInternalToBranchID[aFgBranch];}

	/// Access the global list of node names.
	///
	/// @return A reference to the list of node names.
	///
	const std::vector<std::string>& getNodeNames(void) const {return mNodeNames;}

	/// Get the mapping from the internal site number to the original site.
	///
	/// @return Multimap with key the internal site, and value one of the original sites.
	///
	const std::multimap<size_t, size_t>& getSitesMappingToOriginal(void) {return mSitesMappingToOriginal;}

#ifdef NEW_LIKELIHOOD
	/// All the preparatory steps needed for the Fat Vector approach.
	///
	void postLoad(void);

	/// Analyze the forest to prepare the operation to be done to restore the contiguity to the grouped vector approach.
	///
	/// @param[in] aNode The node from which to start. If null then starts from all the trees' roots.
	///
	void prepareNewReduction(ForestNode* aNode=NULL);

	/// Prepare the data for a forest that has not been reduced
	///
	void prepareNewReductionNoReuse(void);
#endif

#ifdef USE_DAG
	/// Load the forest into a DAG
	///
	/// @param[in] aMaxCopies Max number of forest copies for the various codon classes
	/// @param[in] aCopyId The current copy entered
	/// @param[in] aNode If null starts from the roots. It is used for recursive visit
	///
	void loadForestIntoDAG(unsigned int aMaxCopies, unsigned int aCopyId=0, const ForestNode* aNode=NULL);
#endif

	/// Access the dependency list.
	/// result[tj] = [t1 t2 t3] means: tj can be done after: t1 t2 t3.
	///
	/// @return List of lists of dependencies
	///
	const std::vector< std::vector<unsigned int> >& getTreeDependencies(void) const {return mTreeDependencies;}

	/// Access the reverse dependency list.
	/// result[tj] = [t1 t2 t3] means: tj should be ready before: t1 t2 t3
	///
	/// @return List of lists of reverse dependencies
	///
	const std::vector< std::vector<unsigned int> >& getTreeRevDependencies(void) const {return mTreeRevDependencies;}

	/// Get the computation cost for all sites from the root to the leaves.
	///
	/// @param[out] aEfforts The cost for each site
	/// @param[in] aCostAtLeaf The cost associated to a leaf
	/// @param[in] aCostIntern The cost associated to an internal node
	/// @param[in] aCostPtr The cost associated to a pointer to another node
	///
	void getEffortPerSite(std::vector<unsigned int>& aEfforts, unsigned int aCostAtLeaf, unsigned int aCostIntern, unsigned int aCostPtr) const;


private:
	/// Reduce the common subtree between two (sub)trees
	///
	/// @param[in] aNode The subtree to be tested (i.e.\ if it exists in both trees)
	/// @param[in] aNodeDependent The dependent tree (i.e.\ it could point to subtrees of aNode)
	///
	void reduceSubtreesWalker(ForestNode* aNode, ForestNode* aNodeDependent);

#if !defined(NON_RECURSIVE_VISIT) && !defined(NEW_LIKELIHOOD)
	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aNode The node from which the visit should start
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	///
	/// @return The vector of codons probabilities at the aNode node
	///
	double* computeLikelihoodsWalkerTC(const ForestNode* aNode, const ProbabilityMatrixSet& aSet, unsigned int aSetIdx);
#endif

#ifdef NON_RECURSIVE_VISIT
	/// Walker to prepare the non recursive visit list
	///
	/// @param[in] aNode The current node to be visited
	/// @param[in] aParentNode Parent node for aNode
	/// @param[in] aSite The current site
	/// @param[in,out] aVisitList The list of nodes to be visited in the order of visit.
	/// @param[in,out] aParentList The corresponding parent nodes
	///
	void prepareNonRecursiveVisitWalker(ForestNode* aNode,
										ForestNode* aParentNode,
										unsigned int aSite, 
										std::vector<ForestNode*>& aVisitList, 
										std::vector<ForestNode*>& aParentList);

	/// Walker for the computation of tree likelihood
	///
	/// @param[in] aSet Set of exp(Q*t) matrices
	/// @param[in] aSetIdx Identifier of the set of matrices to be used
	/// @param[in] aSiteIdx The site under computation
	///
	void computeLikelihoodsWalkerNR(const ProbabilityMatrixSet& aSet, unsigned int aSetIdx, unsigned int aSiteIdx);
#endif

	/// Walk the tree to fill the mMapInternalToBranchID map.
	///
	///	@param[in] aNode The node from which to start
	/// @param[out] aMapInternalToBranchID Maps internal branch id to branch id
	///
	void mapInternalToBranchIdWalker(const ForestNode* aNode, std::map<unsigned int, unsigned int>& aMapInternalToBranchID);


private:
	size_t					mNumSites;					///< Number of sites
	const double*			mCodonFreq;					///< Experimental codon frequencies
	const double*			mInvCodonFreq;				///< Inverse of the codon frequencies
	const double*			mInv2CodonFreq;				///< Squared inverse of the codon frequencies
	size_t					mNumBranches;				///< Total number of branches of the original tree
	std::vector<ForestNode>	mRoots;						///< The roots of the forest's trees. Its length is the number of valid sites
	std::vector<double>		mSiteMultiplicity;			///< Multiplicity of the valid sites
	size_t					mNumInternalBranches;		///< Total number of branches of the original tree
	std::vector<unsigned int>
							mTableInternalToBranchID;	///< Map from internal branch number to branch number

	/// Here are global data that will be removed from the various (site) trees
	std::vector<std::string>
							mNodeNames;					///< List of node names. Zero is the root, then its first child and so on
	std::vector<double>		mBranchLengths;				///< List of branch lengths (read from file or stored here to be exported in the tree file)
	size_t					mMarkedInternalBranch;		///< Number of the internal branch as marked in the tree file

#ifdef NEW_LIKELIHOOD

	/// New loglikelihood computation support
		
	/// The mProbs and mProbsOut layout
	///
	/// [site0][site1][site2]...  [site0][site1][site2]...               each is VECTOR_SLOT bytes long (for which only the first N are significant)
	/// [  set 0                 ][  set 1                 ]...          there are 4 (Nt) sets
	/// [    node 0                                             ]...
	///
	/// site_index = node*(Nt*NumSites*VECTOR_SLOT) + set*(NumSites*VECTOR_SLOT) + site*(VECTOR_SLOT)
	///
	CacheAlignedDoubleVector	mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
	CacheAlignedDoubleVector	mProbsOut;					///< mProbs after multiplication by exp(Qt)
	std::vector< std::vector<ForestNode*> >
								mNodesByLevel;				///< Each level contains a list of pointers to nodes at this level. List start from the root.
	FatVectorTransform			mFatVectorTransform;		///< Compute and manage the transformations to pack the "long vector" based on subtree pruning
#else
	/// Unified array for each branch probability vector
	CacheAlignedDoubleVector	mProbs;						///< The concatenation of all the probability vectors for all the nodes and all the classes
#endif
	std::vector< std::vector<unsigned int> >	mTreeDependencies;		///< mTreeDependencies[tj] = [t1 t2 t3] means: tj can be done after: t1 t2 t3
	std::vector< std::vector<unsigned int> >	mTreeRevDependencies;	///< mTreeRevDependencies[tj] = [t1 t2 t3] means: tj should be ready before: t1 t2 t3

#ifdef NON_RECURSIVE_VISIT
	std::vector< std::vector<ForestNode*> >		mVisitTree;				///< List of pointers to tree nodes (a list per site) in the non-recursive visit order
	std::vector< std::vector<ForestNode*> >		mVisitTreeParents;		///< List of parent pointers for the corresponding nodes in the mVisitTree
#endif

#ifdef USE_DAG
	DAGScheduler mDAG;		///< DAG Scheduler
#endif
	unsigned int					mVerbose;					///< If greater than zero prints more info
	std::multimap<size_t, size_t>	mSitesMappingToOriginal;	///< Map reduced site num. to list of corresponding original sites.
};

#endif

