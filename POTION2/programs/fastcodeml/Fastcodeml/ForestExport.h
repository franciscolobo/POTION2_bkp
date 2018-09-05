
#ifndef FORESTEXPORT_H
#define FORESTEXPORT_H

#include "ForestNode.h"
class Forest;

/// Export the phylogenetic tree's forest.
/// This class encapsulates all the routines needed to export in textual form the forest of phylogenetic trees.
/// It is a friend class of Forest.
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2012-02-13 (initial version)
///  @version 1.0
///
class ForestExport
{
public:
	/// Constructor
	///
	/// @param[in] aForest The forest to be exported
	///
	explicit ForestExport(const Forest& aForest) : mForest(aForest) {}

	/// Export the forest in GML format
	///
	/// @param[in] aFilename The filename to be written
	/// @param[in] aCounter Value to substitute \%d or \@d in filename (it is a printf format)
	///
	void exportForest(const char* aFilename, size_t aCounter=0) const;

private:
	/// Walker for the exporter
	///
	///	@param[in] aNode The node from which to start
	///	@param[in] aBranchLengths List of all branch lengths
	/// @param[out] aNodeFrom List of starting nodes
	/// @param[out] aNodeTo List of ending nodes
	/// @param[out] aLength Resulting branch lengths to label branches in exported tree
	///
	void exportForestWalker(const ForestNode* aNode,
							const std::vector<double>& aBranchLengths,
							std::vector< std::pair<int, int> >& aNodeFrom,
							std::vector< std::pair<int, int> >& aNodeTo,
							std::vector<double>& aLength) const;

	/// Disabled assignment operator to avoid warnings on Windows.
	///
	/// @fn ForestExport& operator=(const ForestExport& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	ForestExport& operator=(const ForestExport& /*aObj*/);


private:
	const Forest& mForest;	///< The Forest to be exported
};

#endif
