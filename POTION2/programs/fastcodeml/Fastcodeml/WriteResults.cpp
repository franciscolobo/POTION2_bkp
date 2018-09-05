
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include "WriteResults.h"

void WriteResults::outputResults(void)
{
	// If no file set, then do nothing
	if(!mFilename) return;

	// Range of branches to be output (for H0 and H1)
	std::map<size_t, double>::const_iterator im;
	size_t min_branch = std::numeric_limits<size_t>::max();
	size_t max_branch = 0;
	for(im = mLnL[0].begin(); im != mLnL[0].end(); ++im)
	{
		size_t v = im->first;
		if(v < min_branch) min_branch = v;
		if(v > max_branch) max_branch = v;
	}
	for(im = mLnL[1].begin(); im != mLnL[1].end(); ++im)
	{
		size_t v = im->first;
		if(v < min_branch) min_branch = v;
		if(v > max_branch) max_branch = v;
	}

	// No entries, so do nothing
	if(min_branch == std::numeric_limits<size_t>::max()) return;

	// Open the output file
	std::ofstream out(mFilename, std::ios_base::trunc | std::ios_base::out);
	if(!out.good())
	{
		std::cout << "Cannot create results file <" << mFilename << ">" << std::endl;
		return;
	}

	// Write the log-likelihood values (if a value is not present, write NA)
	for(size_t branch = min_branch; branch <= max_branch; ++branch)
	{
		out << "Branch: " << std::setw(4) << branch << "  LnL0: ";

		// Prints LnL for H0 if present
		im = mLnL[0].find(branch);
		if(im == mLnL[0].end())
		{
			out << std::setw(22) << "NA";
		}
		else
		{
			out << std::setw(22) << std::setprecision(15) << std::fixed << im->second;
		}
		out << "  LnL1: ";

		// Prints LnL for H1 if present
		im = mLnL[1].find(branch);
		if(im == mLnL[1].end())
		{
			out << std::setw(22) << "NA";
		}
		else
		{
			out << std::setw(22) << std::setprecision(15) << std::fixed << im->second;
		}
		out << std::endl;
	}

	// Write the positive selection sites (adding one to the site number because they start from 1 and not zero)
	for(size_t branch = min_branch; branch <= max_branch; ++branch)
	{
		std::map<size_t, std::pair<std::vector<unsigned int>, std::vector<double> > >::const_iterator ipss;
		ipss = mPositiveSelSites.find(branch);
		if(ipss != mPositiveSelSites.end())
		{
			const std::vector<unsigned int>& site = ipss->second.first;
			const std::vector<double>& prob       = ipss->second.second;

			const std::vector<size_t>& idx = orderSites(site);

			size_t ns = site.size();
			for(size_t s=0; s < ns; ++s)
			{
				size_t i = idx[s];
				out << "PositiveSelectionSite for branch: " << std::setw(4) << branch;
				out << "  Site: " << std::setw(6) << site[i] + 1 << "  Prob: " << std::setw(9) << std::setprecision(6) << std::fixed << prob[i] << std::endl;
			}
		}
	}

	out.close();
}

const std::vector<size_t>& WriteResults::orderSites(const std::vector<unsigned int>& aSites) const
{
	// Fill the map with pairs(site, its index)
	std::multimap<unsigned int, size_t> ordered_map;
	size_t end(aSites.size());
	for(size_t i=0; i < end; ++i) ordered_map.insert(std::pair<unsigned int, size_t>(aSites[i], i));

	// Take the list of indices after ordering the list of sites
	mSiteOrder.clear();
	std::multimap<unsigned int, size_t>::const_iterator im(ordered_map.begin());
	std::multimap<unsigned int, size_t>::const_iterator endm(ordered_map.end());
	for(; im != endm; ++im)  mSiteOrder.push_back(im->second);

	return mSiteOrder;
}

void WriteResults::saveLnL(size_t aFgBranch, double aLnL, unsigned int aHypothesis)
{
	// If no file set, then do nothing
	if(!mFilename) return;

	// Sanity check
	if(aHypothesis > 1) return;

	// Save the likelihood for later printing
	mLnL[aHypothesis][aFgBranch] = aLnL;
}

void WriteResults::savePositiveSelSites(size_t aFgBranch, const std::vector<unsigned int>& aPositiveSelSites, const std::vector<double>& aPositiveSelSitesProb)
{
	// If no file set, then do nothing
	if(!mFilename) return;

	// Save the positive selection sites and corresponding probabilities for later printing
	mPositiveSelSites[aFgBranch] = std::make_pair(aPositiveSelSites, aPositiveSelSitesProb);
}

