
#ifdef USE_DAG

#include <iostream>
#include <iomanip>
#include "DAGScheduler.h"

void DAGScheduler::clear(void)
{
	mNodes.clear();
	mEdges.clear();
	mRefCounter.clear();
	mRefCounterSave.clear();
}


void DAGScheduler::loadDependency(unsigned int aCopyId, const void* aDependsOn, const void* aDependant)
{
	// for debug load only id == 0
	if(aCopyId != 0) return;

	// For debug put inside all the nodes
	mNodes.insert(aDependsOn);
	mNodes.insert(aDependant);

	// For debug save all the edges (aDependsOn <-- aDependant)
	mEdges.push_back(std::make_pair(aDependsOn, aDependant));

	// Insert the pair in the reference counter list (increment if present the dependant; do nothing if the dependent is already here)
	mRefCounter.insert(std::pair<const void*, int>(aDependsOn, 0));
	std::pair<std::map<const void*, int>::iterator, bool> ret = mRefCounter.insert(std::pair<const void*, int>(aDependant, 1));
	if(ret.second == false) ++(ret.first->second);
}


void DAGScheduler::endLoadDependencies(unsigned int aNumCodonSets)
{
	// Save the reference counter list
	mRefCounterSave = mRefCounter;

	// For debugging dump the DAG content
	//dumpDAG(std::cout);
}


void* DAGScheduler::getNext(void)
{
	return 0;
}


void DAGScheduler::setDone(const void* aItem)
{
}


void DAGScheduler::resetDependencies(void)
{
	// Reset reference counts
	mRefCounter = mRefCounterSave;
}


void DAGScheduler::dumpDAG(std::ostream& aOut) const
{
	// Output the DAG in TGF format (can be read back by yEd)
	std::set<const void*>::const_iterator ind = mNodes.begin();
	for(; ind != mNodes.end(); ++ind)
	{
		// Output all nodes with id and label
		size_t v = reinterpret_cast<size_t>(*ind);
		aOut << std::hex << *ind << ' ' << std::setw(6) << v << std::endl;
	}
	aOut << '#' << std::endl;
	
	// Then output all edges
	std::vector<std::pair<const void*, const void*> >::const_iterator ied = mEdges.begin();
	for(; ied != mEdges.end(); ++ied)
	{
		aOut << ied->first << ' ' << ied->second << std::endl;
	}
	
	// For debug print the reference counts
	std::map<const void*, int>::const_iterator irc = mRefCounter.begin();
	for(; irc != mRefCounter.end(); ++irc)
	{
		size_t k = reinterpret_cast<size_t>(irc->first);
		int cnt = irc->second;
		std::cout << std::hex << k << ' ' << std::dec << std::setw(6) << cnt << std::endl;
	}
}

#endif
