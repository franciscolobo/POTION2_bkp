
#ifndef BAYESTEST_H
#define BAYESTEST_H

#include <vector>
#include <cstdlib>
#include "TreeAndSetsDependencies.h"
#include "ProbabilityMatrixSet.h"
#include "VerbosityLevels.h"

/// The minimum value for class 2 sites probability to be a positive selection site.
///
const static double MIN_PROB       = 0.50;
const static double ONE_STAR_PROB  = 0.95;
const static double TWO_STARS_PROB = 0.99;

///@cond Private

/// Helper class to compute BEB_N1D^BEB_DIMS at compile time (that is Y^N)
///
template<unsigned int Y, unsigned int N>
class Pow
{
public:
	static const int value = Y * Pow<Y, N-1>::value;
};

/// Specialization of the above class to end the recursion
///
template<unsigned int Y>
class Pow<Y, 1>
{
public:
	static const int value = Y;
};

///@endcond


/// Tests to find the sites under positive selection.
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2010-12-22 (initial version)
///  @version 1.0
///
class BayesTest
{
public:
	/// Constructor.
	///
	/// @param[in] aForest The forest
	/// @param[in] aVerbose The verbosity level
	/// @param[in] aNoReduction If true the dependencies computed are for no reduced forests
	///
	explicit BayesTest(Forest& aForest, unsigned int aVerbose=0, bool aNoReduction=true)
		: mForest(aForest), mNumSites(mForest.getNumSites()), mSiteClassProb(BEB_DIMS*mNumSites),
		  mVerbose(aVerbose), mPriors(mNumSites*BEB_NUM_CAT), mDependencies(mForest, aVerbose), mBEBset(mForest.getNumBranches())
	{
		// Create the dependency list for forest likelihood computation
		mDependencies.computeDependencies(1, aNoReduction);
		if(mVerbose >= VERBOSE_ONLY_RESULTS) mDependencies.print("TEST FOR BEB (before optimization)");
		mDependencies.optimizeDependencies();
		if(mVerbose >= VERBOSE_ONLY_RESULTS) mDependencies.print("TEST FOR BEB");
	}

	/// Destructor.
	///
	~BayesTest() {}

	/// Bayes Empirical Bayes (BEB) test.
	///
	/// @param[in] aVars The variables optimized at the end of the H1 run.
	/// @param[in] aFgBranch The foreground branch under test.
	/// @param[in] aScales The two scales ([0] bg; [1] fg) to rescale the branch lengths. They are computed in H1.
	///
	void computeBEB(const std::vector<double>& aVars, size_t aFgBranch, const std::vector<double>& aScales);
	
	/// Print the sites under positive selection.
	///
	/// @param[in] aFgBranch Identifier of the branch marked as foreground branch
	///
	void printPositiveSelSites(size_t aFgBranch) const;

	/// Extract the sites under positive selection and the corresponding probabilities.
	///
	/// @param[out] aPositiveSelSites Vector of sites under positive selection
	/// @param[out] aPositiveSelSitesProb Corresponding probabilities
	///
	void extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const;

private:
	/// This sets up the grid (mPara[][]) according to the priors.  
	/// It calculates the probability of data at each site given w: f(f_h|w).  
	/// This is calculated using the branch model (NSsites = 0 model = 2), with 
	/// BayesEB=2 used to force the use of the correct scale factors in GetPMatBranch().
	///
	///@verbatim
	/// Order of site classes for iw or f(x_h|w):
	///                     back   fore     num.sets
	/// Branchsite A (121 sets)
	///   site class 0:      w0     w0        10
	///   site class 1:      w1=1   w1=1       1
	///   site class 2a:     w0     w2       100
	///   site class 2b:     w1=1   w2        10
	///@endverbatim
	///
	/// @param[in] aVars  The variables optimized at the end of the H1 run.
	/// @param[in] aSiteMultiplicity  The site multiplicity vector.
	/// @param[in] aFgBranch  The foreground branch under test.
	/// @param[in] aScales  The two scales ([0] bg; [1] fg) to rescale the branch lengths. They are computed in H1.
	///
	/// @return The computed scale.
	///
	double getGridParams(const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch, const std::vector<double>& aScales);

	/// This gives the indices (ix, iy) and the coordinates (aProbX, aProbY, 1-aProbX-aProbY) for 
	/// the aTriangleIdx-th triangle, with aTriangleIdx from 0, 1, ..., BEB_N1D*BEB_N1D-1.
	///
	/// The ternary graph (0-1 on each axis) is partitioned into BEB_N1D*BEB_N1D equal-sized triangles.  
	/// In the first row (ix=0), there is one triangle (iy=0);
	/// In the second row (ix=1), there are 3 triangles (iy=0,1,2);
	/// In the i-th row (ix=i), there are 2*i+1 triangles (iy=0,1,...,2*i).
	///
	/// aProbX rises when ix goes up, but aProbY decreases when iy increases.  (aProbX, aProbY) is the 
	/// centroid in the ij-th small triangle. aProbX and aProbY each takes on 2*BEB_N1D-1 possible values.
	///
	/// @param[out] aProbX The p0 value on the X axis of the triangular grid.
	/// @param[out] aProbY The p1 value on the Y axis of the triangular grid.
	/// @param[in] aTriangleIdx The index inside the triangular grid.
	///
	void getIndexTernary(double* aProbX, double* aProbY, unsigned int aTriangleIdx);

private:
	/// Disabled assignment operator to avoid warnings on Windows.
	///
	/// @fn BayesTest& operator=(const BayesTest& aObj)
	///
	/// @param[in] aObj The object to be assigned
	///
	/// @return The object receiving the assignment
	///
	BayesTest& operator=(const BayesTest&);

private:
	const static unsigned int BEB_N1D = 10;												///< Number of intervals for w0 and w2
	const static unsigned int BEB_DIMS = 4;												///< Number of codon classes (0, 1, 2a, 2b)
	const static unsigned int BEB_NUM_CAT = BEB_N1D + 1 + BEB_N1D*BEB_N1D + BEB_N1D;	///< Total number of categories for w0 and w2 (it is com.ncatG in codeml.c)
	const static unsigned int BEB_NGRID = Pow<BEB_N1D, BEB_DIMS>::value;				///< Number of points in the grid used to evaluate the integral. It is BEB_N1D^BEB_DIMS

private:
	Forest&					mForest;						///< The forest.
	size_t					mNumSites;						///< Number of sites.
	std::vector<double>		mSiteClassProb;					///< Probability of a site to pertain to a given class (one row per class (4 classes), one column per site).
	unsigned int			mVerbose;						///< If greater than zero prints more info
	std::vector<double>		mPriors;						///< Computed priors (each points to a list, one for each site)
	TreeAndSetsDependencies	mDependencies;					///< Dependency list for likelihood computation
	ProbabilityMatrixSetBEB	mBEBset;						///< Probability matrix set to be used for likelihood computation
};

#endif

