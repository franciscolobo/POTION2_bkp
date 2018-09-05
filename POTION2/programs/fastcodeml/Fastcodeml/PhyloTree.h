
#ifndef PHYLOTREE_H
#define PHYLOTREE_H

#include <string>
#include <vector>
#include <fstream>
#include "TreeNode.h"
#include "ForestNode.h"
#include "Types.h"

/// Phylogenetic tree public interface.
///
///   @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///   @date 2010-08-30 (initial version)
///   @version 1.0
///
class PhyloTree
{
protected:
	/// Constructor.
	///
	/// @param[in] aVerboseLevel Set the verbosity level
	/// - 0: No messages
	/// - 3: Print the tree structure
	///
	explicit PhyloTree(unsigned int aVerboseLevel=0) : mVerboseLevel(aVerboseLevel) {}

	/// Destructor.
	///
	virtual ~PhyloTree();

public:
	/// Load a phylo tree definition from a Newick formatted file.
	/// This routine should be defined by a derived class.
	///
	/// @param[in] aFilename The filename
	///
	/// @exception FastCodeMLFatal For errors like cannot open the file
	///
	virtual void readFile(const char *aFilename) =0;

	/// Print the phylogenetic tree completed with all the info loaded in the same format as read in.
	///
	/// @param[in] aOut Output stream
	/// @param[in] aNode The node from which to start. If null starts from the root.
	///
	virtual void printTreeUnformatted(std::ostream& aOut, TreeNode *aNode=NULL) const =0;

	/// Print the phylogenetic tree completed with all the info loaded in the same format as read in and annotated with the internal branch number.
	///
	/// @param[in] aOut Output stream
	/// @param[in] aNode The node from which to start. If null starts from the root.
	/// @param[in] aInternalBranch Internal branch identifier to annotate the current branch.
	///
	/// @return The new internal branch id
	///
	virtual int printTreeAnnotated(std::ostream& aOut, TreeNode *aNode=NULL, int aInternalBranch=0) const =0;

	/// Return the list of species.
	///
	/// @return Reference to a vector of species names as read from the leaves of the tree
	///
	const std::vector<std::string>& getSpecies(void) const;

	/// Number of tree branches.
	///
	/// @return The number of tree branches
	///
	size_t getNumBranches(void) const {return mInternalNodes.size()+mLeavesSpecies.size();}

	/// Return the index of the first marked branch.
	///
	/// @return The index of the marked internal branch or UINT_MAX if none marked
	///
	size_t getMarkedInternalBranch(void) const;

	/// Clone the tree using ForestNode.
	/// Called without aTreeNode starts from the tree root.
	///
	/// @param[out] aForestNode The ForestNode that becomes the root of the cloned tree
	/// @param[in] aTreeId The tree running id.
	/// @param[in] aNumSites Total number of sites. Needed to assign pointers to aProbVectors.
	/// @param[in] aProbVectors Contiguous storage for the probability vectors for each branch (ignored if NEW_LIKELIHOOD defined).
	/// @param[in] aTreeNode The node from which to start the cloning in the tree. If not present starts from the root
	/// @param[in] aNodeId The node running id. For the root it is UINT_MAX.
	///
	/// @return The node id to the next node
	///
	unsigned int cloneTree(ForestNode* aForestNode, unsigned int aTreeId, size_t aNumSites, CacheAlignedDoubleVector& aProbVectors, const TreeNode* aTreeNode=NULL, unsigned int aNodeId=0) const;

	/// Extract global data (data that do not depend on the site) from the phylo tree.
	///
	/// @param[out] aNodeNames Ordered list of the node labels
	/// @param[out] aBranchLengths Ordered list of branch lists as read from the file
	/// @param[out] aMarkedIntBranch Pointer to location where the marked internal branch number is stored
	/// @param[in] aTreeNode The node from which to start the cloning in the tree. If not present starts from the root
	/// @param[in] aNodeId The node running id. For the root it is UINT_MAX.
	///
	/// @return The node id to the next node
	///
	unsigned int collectGlobalTreeData(std::vector<std::string>& aNodeNames, std::vector<double>& aBranchLengths, size_t* aMarkedIntBranch, const TreeNode* aTreeNode=NULL, unsigned int aNodeId=0) const;

	///	Check if any leaf has an associated branch length equal to zero.
	///
	/// @param[in,out] aOnLeafCnt The count of leaf branches with length zero (should be zeroed before first call)
	/// @param[in,out] aOnIntCnt The count of internal branches with length zero (should be zeroed before first call)
	/// @param[in] aTreeNode The node from which to start checking the tree. If not present starts from the root.
	///
	void countNullBranchLengths(int& aOnLeafCnt, int& aOnIntCnt, const TreeNode* aTreeNode=NULL) const;

	/// Verify if the root of the tree is well formed for the analysis.
	/// Currently it prints the number of branches emanating from the tree root and how many of them are leaves.
	///
	/// @exception FastCodeMLFatal Root has only one branch or points only to leaves.
	///
	void checkRootBranches(void) const;

protected:

	/// Fill the list of Species (the leaves of the tree).
	///
	/// @param[in,out] aNode The tree node (recursively called function starting from the tree root)
	///
	void fillSpecies(TreeNode *aNode);

	/// Fill the list of internal nodes.
	///
	/// @param[in,out] aNode The tree node (recursively called function starting from the tree root)
	///
	void fillInternalBranches(TreeNode *aNode);


protected:
	TreeNode				mTreeRoot;				///< The root of the phylogenetic tree in memory
	unsigned int			mVerboseLevel;			///< The verbosity level
	std::vector<TreeNode *>	mLeavesSpecies;			///< The list of the tree leaves
	std::vector<TreeNode *> mInternalNodes;			///< The list of the tree internal nodes
	mutable std::vector<std::string> mSpeciesList;	///< Temporary to securely return the list of species from getSpecies()
};

#endif
