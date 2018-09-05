
#ifndef TREEANDSETSDEPENDENCIES_H
#define TREEANDSETSDEPENDENCIES_H

/// List (each list depends on the previous) of list (sites to be executed in parallel) of pairs (site, site class) stored as site*4+site_class 
typedef std::vector< std::vector<unsigned int> > ListDependencies;

#include "Forest.h"

/// Create the lists of trees and sets that can be computed concurrently
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-21 (initial version)
///     @version 1.0
///
///
class TreeAndSetsDependencies
{
public:
	/// Constructor.
	///
	/// @param[in] aForest The corresponding forest
	/// @param[in] aVerbose The verbosity level
	///
	explicit TreeAndSetsDependencies(const Forest& aForest, unsigned int aVerbose=0) : mForest(aForest), mNoParallel(false), mVerbose(aVerbose) {}

	/// Destructor.
	///
	~TreeAndSetsDependencies() {}

	/// Prepare the dependency list.
	///
	/// @param[in] aNumSets How many sets (that is parallel visit to the forest) should this dependency list cover. For example H0: 3 codon classes; H1: 4 codon classes
	/// @param[in] aNoParallel True if no dependencies needed to be setup
	///
	void computeDependencies(unsigned int aNumSets, bool aNoParallel);

	/// Optimize the dependency list.
	///
	void optimizeDependencies(void);

	/// Print the dependency list and the maximum effort
	///
	/// @param[in] aTitle If present a title is written on top of the table
	///
	void print(const char* aTitle=NULL) const;

	/// Access the computed dependency list
	///
	/// @return Reference to the dependency list.
	///
	const ListDependencies& getDependencies(void) {return mDependenciesClassesAndTrees;}

private:
	/// Move dependent trees from class to class trying to balance the size of the classes
	///
	/// @param[in] aGreedy If set move the trees as far as possible, otherwise move them to the first useful pot.
	///
	/// @return True if the redistribution has been successful
	///
	bool balanceDependenciesClassesAndTrees(bool aGreedy);

	/// Measure the relative effort of doTransitionAtLeaf and doTransition.
	///
	/// @return The ratio between the time spent by doTransition and doTransitionAtLeaf.
	///
	unsigned int measureRelativeEffort(void);

public:
	/// Given the combined entry in the dependency list, extract the site number.
	///
	/// @param[in] aPair The combined entry
	///
	/// @return The site number
	///
	inline static unsigned int getSiteNum(unsigned int aPair) {return aPair >> 2;}

	/// Given the combined entry in the dependency list, extract the set number.
	///
	/// @param[in] aPair The combined entry
	///
	/// @return The set number
	///
	inline static unsigned int getSetNum(unsigned int aPair) {return aPair & 03u;}

private:
	/// Create the combined entry in the dependency list from site and set numbers.
	///
	/// @param[in] aSite The site number
	/// @param[in] aSet The set number
	///
	/// @return The combined entry
	///
	inline static unsigned int makePair(unsigned int aSite, unsigned int aSet) {return aSite << 2 | aSet;}

	/// Disabled assignment operator to avoid warnings on Windows.
	///
	/// @fn TreeAndSetsDependencies& operator=(const TreeAndSetsDependencies& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	TreeAndSetsDependencies& operator=(const TreeAndSetsDependencies& /*aObj*/);

private:
	const Forest&				mForest;						///< The forest for which dependencies should be calculated
	ListDependencies			mDependenciesClassesAndTrees;	///< The groups of dependencies between trees
	bool						mNoParallel;					///< Set if the execution is sequential, so no tree dependencies needed
	unsigned int				mVerbose;						///< The verbose level
	std::vector<unsigned int>	mEffortPerSite;					///< Effort per site
};


#endif
