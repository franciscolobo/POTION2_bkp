
#ifndef DAGSCHEDULER_H
#define DAGSCHEDULER_H

#include <set>
#include <utility>
#include <vector>
#include <map>

/// The high-performance DAG scheduler. Visits the DAG obeying the dependencies.
///
///   @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///   @date 2012-05-14 (initial version)
///   @version 1.0
///
class DAGScheduler
{
public:
	/// Constructor.
	///
	DAGScheduler() {}

	/// Destructor.
	///
	~DAGScheduler() {}

	/// Clean everything
	///
	void clear(void);

	/// Add a piece of DAG
	///
	/// @param[in] aCopyId The current copy of the forest in which the dependency should be recorded
	/// @param[in] aDependsOn The object on which aDependant depends
	/// @param[in] aDependant The object to be computed only when aDependsOn is ready
	///
	void loadDependency(unsigned int aCopyId, const void* aDependsOn, const void* aDependant);
	
	/// Signals dependencies load has finished.
	///
	/// @param[in] aNumCodonSets Number of DAG copies (one for each codon classes)
	///
	void endLoadDependencies(unsigned int aNumCodonSets=1);

	/// Get next item to process.
	///
	/// @return The next object with all dependencies satisfied. Return NULL if no more objects to be executed.
	///
	void* getNext(void);
	
	/// Signals finished execution on aItem.
	///
	/// @param[in] aItem Object on which execution has ended.
	///
	void setDone(const void* aItem);
	
	/// Reset all dependencies to initial state
	///
	void resetDependencies(void);

	/// Dump DAG in textual form
	///
	void dumpDAG(std::ostream& aOut) const;
	

private:
	std::set<const void*> mNodes;								///< Load the distinct node addresses
	std::vector<std::pair<const void*, const void*> > mEdges;	///< Load the edges
	std::map<const void*, int> mRefCounter;						///< Reference counter
	std::map<const void*, int> mRefCounterSave;					///< Save the reference counter map 
};

#endif

