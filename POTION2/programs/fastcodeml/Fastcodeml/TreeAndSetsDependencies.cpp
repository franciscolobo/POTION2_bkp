
#include <iostream>
#include <iomanip>
#include <set>
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wlong-long"
#endif
#include <boost/dynamic_bitset.hpp>
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#pragma GCC diagnostic pop
#endif
#include "TreeAndSetsDependencies.h"
#include "Timer.h"
#include "VerbosityLevels.h"

#ifdef USE_LAPACK
#include "blas.h"
#endif

#ifndef VTRACE
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

void TreeAndSetsDependencies::computeDependencies(unsigned int aNumSets, bool aNoParallel)
{
	size_t  j;

	// Save for the optimization phase
	mNoParallel = aNoParallel;

	// Take values from forest
	size_t num_sites = mForest.getNumSites();
	std::vector< std::vector<unsigned int> > tree_dependencies = mForest.getTreeDependencies();

	// Collect the class dependencies
	std::vector< std::vector<unsigned int> > tree_groups_dependencies;

	// If no dependencies
	if(aNoParallel)
	{
		std::vector<unsigned int> v(num_sites);

		for(size_t k=0; k < num_sites; ++k) v[k] = static_cast<unsigned int>(num_sites-k-1); // Remember: prior (could) point to subsequent

		tree_groups_dependencies.push_back(v);
	}
	else
	{
		// Prepare the search of dependencies
		boost::dynamic_bitset<> done(num_sites);	// The sites that has dependencies satisfied in the previous level
		boost::dynamic_bitset<> prev;				// Dependencies till the previous level
		std::vector<unsigned int> v;				// Temporary list of sites

		// Mark trees without dependencies
		// tree_dependencies[tj] can be done after: t1 t2 t3
		for(size_t site=0; site < num_sites; ++site)
		{
			if(tree_dependencies[site].empty())
			{
				done.set(site);
				v.push_back(static_cast<unsigned int>(site));
			}
		}

		// Prepare the dependency list
		tree_groups_dependencies.push_back(v);
		prev = done;

		// Start to find trees with one, two, ... dependencies
		for(unsigned int numdep=1;; ++numdep)
		{
			v.clear();
			bool all_done = true;
			for(size_t site=0; site < num_sites; ++site)
			{
				// If tree i has been already processed skip it
				if(prev[site]) continue;
				all_done = false;

				size_t nc = tree_dependencies[site].size();
				bool all = true;
				for(j=0; j < nc; ++j) if(!prev[tree_dependencies[site][j]]) {all = false; break;}
				if(all)
				{
					v.push_back(static_cast<unsigned int>(site));
					done.set(site);
				}
			}
			if(all_done) break;
			tree_groups_dependencies.push_back(v);
			prev = done;
		}
	}

	// Number of dependencies classes
	size_t nc = tree_groups_dependencies.size();
	
	// One dependency classes
	std::vector<unsigned int> one_class;

	// Transform the list multiplying the entries by the number of codon classes
	mDependenciesClassesAndTrees.clear();
	for(size_t dep_class=0; dep_class < nc; ++dep_class)
	{
		// Prepare the dependency classe
		one_class.clear();

		// Number of trees in the class
		const size_t nt = tree_groups_dependencies[dep_class].size();
		for(unsigned int set=0; set < aNumSets; ++set)
		{
			for(j=0; j < nt; ++j)
			{
				one_class.push_back(makePair(tree_groups_dependencies[dep_class][j], set));
			}
		}
		mDependenciesClassesAndTrees.push_back(one_class);
	}
}


unsigned int TreeAndSetsDependencies::measureRelativeEffort(void)
{
#ifdef USE_LAPACK
	// Number of measurement iterations
	static const int NR = 10000;

	// Prepare dummy data
	double dummy[N];
	double m[N*N];
	for(int i=0; i < N*N; ++i) m[i] = 0.1;
	Timer timer;

	// Measure doTransitionAtLeaf()
	timer.start();
	for(int k=0; k < NR; ++k)
	{
		for(int c=0; c < N; ++c)
		{
			int i = 0;
			for(; i < c; ++i) dummy[i] = m[c*N+i];
			for(; i < N; ++i) dummy[i] = m[i*N+c];
		}
	}
	double time_leaf = static_cast<double>(timer.stop());

	// Measure doTransition()
	timer.start();
	for(int i=0; i < NR; ++i)
	{
		for(int c=0; c < N; ++c)
		{
			dsymv_("U", &N, &D1, m, &N, dummy, &I1, &D0, dummy, &I1);
		}
	}
	double time_non_leaf = static_cast<double>(timer.stop());
	unsigned int effort_ratio = static_cast<unsigned int>(time_non_leaf/time_leaf+0.5);
#else
	unsigned int effort_ratio = 16;
#endif
	return effort_ratio;
}


void TreeAndSetsDependencies::optimizeDependencies(void)
{
#if 0
	unsigned int relative_effort = measureRelativeEffort();

	// Compute effort per site
	if(!mNoParallel)
	{
		mForest.getEffortPerSite(mEffortPerSite, 10, 10*relative_effort, 1);
	}
#endif
	balanceDependenciesClassesAndTrees(true);
}

void TreeAndSetsDependencies::print(const char* aTitle) const
{
	// Do nothing if the debug level is insufficient
	if(mVerbose < VERBOSE_MORE_DEBUG) return;

	// If present print the title
	std::cout << std::endl;
	if(aTitle) std::cout << aTitle << std::endl;

	// Print the number of executions for each class
	const size_t num_classes = mDependenciesClassesAndTrees.size();
	for(size_t tree_class=0; tree_class < num_classes; ++tree_class)
	{
		std::cout << "Trees in class " << std::setw(3) << tree_class << ": " << std::setw(5) << mDependenciesClassesAndTrees[tree_class].size() << std::endl;
	}
#if 0
	// Compute the max effort of each class
	unsigned int total_effort = 0;
	for(size_t tree_class=0; tree_class < num_classes; ++tree_class)
	{
		size_t nops = mDependenciesClassesAndTrees[tree_class].size();
		unsigned int tree_class_effort = 0;
		for(size_t op=0; op < nops; ++op)
		{
			unsigned int v = mDependenciesClassesAndTrees[tree_class][op];
			if(getSetNum(v) > 0) continue;
			unsigned int site = getSiteNum(v);

			unsigned int effort = mEffortPerSite[site];
		}
	}
#endif
}

bool TreeAndSetsDependencies::balanceDependenciesClassesAndTrees(bool aGreedy)
{
	// Do nothing if no dependencies or no parallel execution
#ifndef _OPENMP
	return false;
#else
	if(mNoParallel) return false;

	const std::vector< std::vector<unsigned int> >& tree_rev_dependencies = mForest.getTreeRevDependencies();

	// At each level collect the 'jolly' threads (trees that are not preconditions for trees in classes above)
	// This step makes sense only if run multithread and if there are more than one class
	const size_t num_threads = static_cast<size_t>(omp_get_max_threads());
	const size_t num_classes = mDependenciesClassesAndTrees.size();
	if(num_threads < 2 || num_classes < 2) return false;

	// This set will contain the site & class pairs that can be postponed without problem
	std::set<unsigned int> jolly_sites;

	// Check if the current level can be balanced
	for(size_t k=0; k < num_classes; ++k)
	{
		// Try to have num sites at this level multiple of number of threads
		const size_t num_jolly = jolly_sites.size();
		const size_t class_num_sites = mDependenciesClassesAndTrees[k].size();
		const unsigned int over = static_cast<unsigned int>(class_num_sites % num_threads);

		// Compute how many sites to add to have a multiple of num threads
		size_t needed_add = 0;
		size_t resid_jolly = num_jolly;
		if(over)
		{
			size_t min_add = num_threads - over;
			if(min_add <= num_jolly)
			{
				needed_add = min_add;
				resid_jolly -= needed_add;
			}
		}
		if(aGreedy)
		{
			needed_add += (resid_jolly/num_threads)*num_threads;
		}

		// Compute how many sites could became a jolly
		unsigned int new_possible_jolly_sites = 0;
		for(unsigned int s=0; s < class_num_sites; ++s)
		{
			// mTreeRevDependencies[tj] should be ready before: t1 t2 t3
			unsigned int site = getSiteNum(mDependenciesClassesAndTrees[k][s]);
			if(tree_rev_dependencies[site].empty()) ++new_possible_jolly_sites;
		}

		// Compute how many sites to remove to have a multiple of num threads
		size_t needed_remove = 0;
		resid_jolly = static_cast<size_t>(new_possible_jolly_sites);
		if(class_num_sites > num_threads && new_possible_jolly_sites >= over)
		{
			needed_remove = static_cast<size_t>(over);
			resid_jolly -= needed_remove;
		}
		if(aGreedy && resid_jolly >= 2*num_threads)
		{
			needed_remove += ((resid_jolly - num_threads)/num_threads)*num_threads;
		}

		// There are at least the number needed to add in the jolly list, so add them
		if(needed_add > 0)
		{
			for(unsigned int j=0; j < needed_add; ++j)
			{
				std::set<unsigned int>::iterator it(jolly_sites.begin());
				mDependenciesClassesAndTrees[k].push_back(*it);
				jolly_sites.erase(it);
			}
		}
		else if(needed_remove > 0)
		{
			// Save what remains at this level
			std::vector<unsigned int> new_content;

			for(unsigned int s=0; s < class_num_sites; ++s)
			{
				unsigned int site = getSiteNum(mDependenciesClassesAndTrees[k][s]);
				if(needed_remove && tree_rev_dependencies[site].empty())
				{
					jolly_sites.insert(mDependenciesClassesAndTrees[k][s]);
					--needed_remove;
				}
				else
				{
					new_content.push_back(mDependenciesClassesAndTrees[k][s]);
				}
			}
			mDependenciesClassesAndTrees[k].swap(new_content);
		}
	}

	// If there are still jolly sites, add them to the last class
	if(!jolly_sites.empty())
	{
		mDependenciesClassesAndTrees[num_classes-1].insert(mDependenciesClassesAndTrees[num_classes-1].end(), jolly_sites.begin(), jolly_sites.end());
	}

	return true;
#endif
}
