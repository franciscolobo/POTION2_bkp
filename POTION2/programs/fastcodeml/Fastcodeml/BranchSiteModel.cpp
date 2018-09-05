
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cmath>
#include <memory>

#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4267) // warning C4267: 'argument' : conversion from 'size_t' to 'unsigned int', possible loss of data
#endif

#include "nlopt.hpp"

#ifdef _MSC_VER
    #pragma warning(pop)
#endif

#include "BranchSiteModel.h"
#include "MathSupport.h"
#include "Exceptions.h"
#include "CodeMLoptimizer.h"
#include "ParseParameters.h"

/// Starting value for the computed maximum likelihood.
/// Beware: -HUGE_VAL is too low and on Linux it is converted to -Inf (with subsequent NLopt crash)
static const double VERY_LOW_LIKELIHOOD = -1e14;


void BranchSiteModel::printFinalVars(std::ostream& aOut) const
{
	// To nicely format num branch lengths per line
	static const unsigned int VARS_PER_LINE = 8;
	unsigned int count_per_line = 0;
	static const std::streamsize VARS_PRECISION = 7;
	static const std::streamsize VARS_WIDTH     = 11;
	
	// Write the data with uniform precision
	std::streamsize prec = aOut.precision(VARS_PRECISION);
	aOut.setf(std::ios::fixed, std::ios::floatfield);

	// Print all variables formatted to be readable
	double v0 = 0;
	std::vector<double>::const_iterator ix(mVar.begin());
	const std::vector<double>::const_iterator end(mVar.end());
	for(int k = -static_cast<int>(mNumTimes); ix != end; ++ix,++k)
	{
		switch(k)
		{
		case 0:
			if(count_per_line) aOut << std::endl;
			v0 = *ix;
			break;
		case 1:
			{
				double p[4];
				getProportions(v0, *ix, p);
				aOut << std::endl;
				aOut <<   "p0:" << std::setw(VARS_WIDTH) << p[0];
				aOut << "  p1:" << std::setw(VARS_WIDTH) << p[1];
				aOut << "  p2a:" << std::setw(VARS_WIDTH) << p[2];
				aOut << "  p2b:" << std::setw(VARS_WIDTH) << p[3];
				aOut << std::endl;
			}
			break;
		case 2:
			aOut << "w0:" << std::setw(VARS_WIDTH) << *ix;
			break;
		case 3:
			aOut << "  k: " << std::setw(VARS_WIDTH) << *ix;
			break;
		case 4:
			aOut << "  w2: " << std::setw(VARS_WIDTH) << *ix;
			break;
		default:
			aOut << std::setw(VARS_WIDTH) << *ix;
			++count_per_line;
			if(count_per_line == VARS_PER_LINE)
			{
				count_per_line = 0;
				aOut << std::endl;
			}
			break;
		}
	}
	aOut << std::endl;
	aOut.precision(prec);
}

void BranchSiteModel::printVar(const std::vector<double>& aVars, double aLnl, std::ostream& aOut) const
{
	// Write the data with uniform precision
	std::streamsize prec = aOut.precision(7);
	aOut.setf(std::ios::fixed, std::ios::floatfield);

	// Write the LnL value (if set)
	if(aLnl < DBL_MAX) aOut << std::endl << aLnl << std::endl;

	// Print all variables formatted to be readable
	double v0 = 0;
	int per_linea = 0;
	std::vector<double>::const_iterator ix(aVars.begin());
	const std::vector<double>::const_iterator end(aVars.end());
	for(int k = -static_cast<int>(mNumTimes); ix != end; ++ix,++k)
	{
		switch(k)
		{
		case 0:
			aOut << std::endl;
			aOut <<   "v0: " << *ix;
			v0 = *ix;
			break;
		case 1:
			aOut << "  v1: " << *ix;
			{
				double p[4];
				getProportions(v0, *ix, p);
				aOut << "  [" << p[0] << ", " << p[1] << ", " << p[2] << ", " << p[3] << "]";
			}
			break;
		case 2:
			aOut << "  w0: " << *ix;
			break;
		case 3:
			aOut <<  "  k: " << *ix;
			break;
		case 4:
			aOut << "  w2: " << *ix;
			break;
		default:
			aOut << *ix << ' ';
			++per_linea;
			if(per_linea == 10) {if(k != -1) aOut << std::endl; per_linea = 0;}
			break;
		}
	}
	aOut << std::endl;
	aOut.precision(prec);
}


void BranchSiteModel::setLimits(size_t aNumTimes, size_t aNumVariables)
{
	// Reserve space
	mLowerBound.reserve(aNumTimes+aNumVariables);	mUpperBound.reserve(aNumTimes+aNumVariables);
	
	// Set lower constrains							// Set upper constrains
	mLowerBound.assign(aNumTimes, 0);				mUpperBound.assign(aNumTimes, 50.0);	// T
#ifdef USE_ORIGINAL_PROPORTIONS
	mLowerBound.push_back(-99.0);					mUpperBound.push_back(99.0);			// x0 -> p0
	mLowerBound.push_back(-99.0);					mUpperBound.push_back(99.0);			// x1 -> p1
#else
	mLowerBound.push_back(0.0);						mUpperBound.push_back(1.0);				// p0+p1
	mLowerBound.push_back(0.0);						mUpperBound.push_back(1.0);				// p0/(p0+p1)
#endif
	mLowerBound.push_back(1e-6);					mUpperBound.push_back(1.0);				// w0
	mLowerBound.push_back(0.0001);					mUpperBound.push_back(20.0);			// k
	if(aNumVariables >= 5)
	{
		mLowerBound.push_back(1.0);					mUpperBound.push_back(999.0);			// w2
	}
}


void BranchSiteModel::initFromTree(void)
{
	// Initialize branch lengths from the phylo tree
	mForest.setTimesFromLengths(mVar);

	// Ask for initialization completion
	mInitStatus |= INIT_TIMES|INIT_TIMES_FROM_FILE;
}


void BranchSiteModel::initFromParams(void)
{
	// Get the parameters (the default values or the values set on the command line)
	ParseParameters* params = ParseParameters::getInstance();

	// Initialization as in CodeML (seems)
	mVar[mNumTimes+2] = params->getParameter("w0");							// w0
	mVar[mNumTimes+3] = params->getParameter("k");							// k

	double p0 = params->getParameter("p0");
	double p1 = params->getParameter("p1");
#ifdef USE_ORIGINAL_PROPORTIONS
	if(p0 <= 0 || p1 <= 0) throw FastCodeMLFatal("Invalid p0 and p1 values");
	mVar[mNumTimes+0] = log(p0);											// p0 -> x0
	mVar[mNumTimes+1] = log(p1);											// p1 -> x1
#else
	if(p0 < 0 || p1 < 0 || (p0+p1) < 1e-15) throw FastCodeMLFatal("Invalid p0 and p1 values");
	mVar[mNumTimes+0] = p0+p1;												// p0+p1
	mVar[mNumTimes+1] = p0/(p0+p1);											// p0/(p0+p1)
#endif
	if(mNumVariables == 5) mVar[mNumTimes+4] = params->getParameter("w2");	// w2

	// The parameters have been initializated
	mInitStatus |= INIT_PARAMS;
}

void BranchSiteModel::initFromResult(const std::vector<double>& aPreviousResult, unsigned int aValidLen)
{
	// Adjust the length to be copied
	if(aValidLen == 0) aValidLen = static_cast<unsigned int>(aPreviousResult.size());

	// Too long, cut. Too short, ignore. Remember H0 has 4 variables.
	if(aValidLen > mNumTimes+mNumVariables) aValidLen = mNumTimes+mNumVariables;
	else if(aValidLen < mNumTimes)
	{
		mInitStatus = INIT_NONE;
		return;
	}
	else if(aValidLen < mNumTimes+4) aValidLen = mNumTimes;

	// Copy the requested values
	mVar.assign(aPreviousResult.begin(), aPreviousResult.begin()+static_cast<size_t>(aValidLen));
	mVar.resize(static_cast<size_t>(mNumTimes+mNumVariables));

	// Ask for initialization completion
	if(aValidLen == mNumTimes)        mInitStatus = INIT_TIMES;
	else if(aValidLen == mNumTimes+4) mInitStatus = INIT_TIMES|INIT_PARAMS_H1;
	else                              mInitStatus = INIT_TIMES|INIT_PARAMS;
}


void BranchSiteModel::initVariables(void)
{
	unsigned int i;

	// Initialize times (if not already initialized)
	if((mInitStatus & INIT_TIMES) != INIT_TIMES)
	{
		for(i=0; i < mNumTimes; ++i) mVar[i] = 0.1 + 0.5 * randFrom0to1();	// T
	}

	// Initialize w0, k, v1, v2 (if not already initialized)
	if((mInitStatus & INIT_PARAMS_H1) != INIT_PARAMS_H1)
	{
		if((mInitStatus & INIT_TIMES_FROM_FILE) == INIT_TIMES_FROM_FILE)
		{
#ifdef USE_ORIGINAL_PROPORTIONS
			mVar[mNumTimes+0] = 1.0 + 0.2 * randFrom0to1();						// x0 -> p0
			mVar[mNumTimes+1] = 0.0 + 0.2 * randFrom0to1();						// x1 -> p1
#else
			double x0 =  exp(1.0 + 0.2 * randFrom0to1());
			double x1 =  exp(0.0 + 0.2 * randFrom0to1());
			double tot = x0 + x1 + 1.0;
			double p0 = x0/tot;
			double p1 = x1/tot;
			if(p0+p1 < 1e-15) {p0 = 1e-6; p1 = 1e-6;}

			mVar[mNumTimes+0] = p0+p1;											// p0+p1
			mVar[mNumTimes+1] = p0/(p0+p1);										// p0/(p0+p1)
#endif
			mVar[mNumTimes+2] = 0.2 + 0.1 * randFrom0to1();						// w0
			mVar[mNumTimes+3] = 0.4;											// k
		}
		else
		{
#ifdef USE_ORIGINAL_PROPORTIONS
			mVar[mNumTimes+0] = 0.5  +       randFrom0to1();					// x0 -> p0
			mVar[mNumTimes+1] = 0.5  +       randFrom0to1();					// x1 -> p1
#else
			double x0 =  exp(0.5 + randFrom0to1());
			double x1 =  exp(0.5 + randFrom0to1());
			double tot = x0 + x1 + 1.0;
			double p0 = x0/tot;
			double p1 = x1/tot;
			if(p0+p1 < 1e-15) {p0 = 1e-6; p1 = 1e-6;}

			mVar[mNumTimes+0] = p0+p1;											// p0+p1
			mVar[mNumTimes+1] = p0/(p0+p1);										// p0/(p0+p1)
#endif
			mVar[mNumTimes+2] = 0.5  +       randFrom0to1();					// w0
			mVar[mNumTimes+3] = 0.5  +       randFrom0to1();					// k
		}
	}

	// Initialize w2 if needed
	if(mNumVariables == 5 && (mInitStatus & INIT_PARAM_W2) != INIT_PARAM_W2)
	{
		if((mInitStatus & INIT_TIMES_FROM_FILE) == INIT_TIMES_FROM_FILE)
		{
			mVar[mNumTimes+4] = 0.5 +       randFrom0to1();						// w2
		}
		else
		{
			mVar[mNumTimes+4] = 1.0 + 0.5 * randFrom0to1();						// w2
		}
	}

	// Re-initialize the next time
	mInitStatus = INIT_NONE;

	// Check the initial values to be inside the domain (otherwise use the same clamp as in CodeML)
	// Don't clamp the results if they came from H1
	if((mInitStatus & (INIT_TIMES|INIT_PARAMS_H1)) != (INIT_TIMES|INIT_PARAMS_H1))
	{
		unsigned int nv = mNumTimes+mNumVariables;
		for(i=0; i < nv; ++i)
		{
			if(mVar[i] < mLowerBound[i] * 1.05) mVar[i] = mLowerBound[i] * 1.05;
			if(mVar[i] > mUpperBound[i] / 1.05) mVar[i] = mUpperBound[i] / 1.05;
		}
		for(i=0; i < nv; ++i)
		{
			if(mVar[i] < mLowerBound[i]) mVar[i] = mLowerBound[i] * 1.2;
			if(mVar[i] > mUpperBound[i]) mVar[i] = mUpperBound[i] * 0.8;
		}
	}
}


double BranchSiteModelNullHyp::operator()(size_t aFgBranch, bool aStopIfBigger, double aThreshold)
{
	// Initialize the variables to be optimized
	initVariables();

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;

	// Initialize the matrix set and the matrix set used for gradient computation
	mSet.initializeSet(mForest.adjustFgBranchIdx(aFgBranch));
	mSetForGradient.initializeFgBranch(mForest.adjustFgBranchIdx(aFgBranch));

	// Run the optimizer
	return maximizeLikelihood(aFgBranch, aStopIfBigger, aThreshold);
}


double BranchSiteModelAltHyp::operator()(size_t aFgBranch)
{
	// Initialize the variables to be optimized
	initVariables();

	// Initialize the variables used to avoid unneeded recomputing
	mPrevK      = DBL_MAX;
	mPrevOmega0 = DBL_MAX;
	mPrevOmega2 = DBL_MAX;

	// Initialize the matrix set and the matrix set used for gradient computation
	mSet.initializeSet(mForest.adjustFgBranchIdx(aFgBranch));
	mSetForGradient.initializeFgBranch(mForest.adjustFgBranchIdx(aFgBranch));

	// Run the optimizer
	return maximizeLikelihood(aFgBranch, false, 0.);
}

double BranchSiteModelNullHyp::computeLikelihoodForGradient(const std::vector<double>& aVar, bool aTrace, size_t aGradientVar)
{
	// One more function invocation
	++mNumEvaluations;

	// Compute the following values for gradient only if anything different from branch length has changed
	if(aGradientVar >= mNumTimes)
	{
		// Save the values to local variables to speedup access
		const double* params = &aVar[mNumTimes];
		const double  omega0 = params[2];
		const double  kappa  = params[3];

		// The values for gradient are computed in order, use this fact to reduce computations
		switch(aGradientVar - mNumTimes)
		{
		case 0:
		case 1:
			// Recompute all proportions if v0 or v1 change
			getProportions(params[0], params[1], mProportions);
			break;

		case 2:
			// Return to the original values
			getProportions(params[0], params[1], mProportions);

			mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
			mQw0.eigenQREV();
			break;

		case 3:
#ifdef _MSC_VER
			#pragma omp parallel sections default(none) shared(omega0, kappa)
#else
			#pragma omp parallel sections default(shared)
#endif
			{
				#pragma omp section
				{
					mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
					mQw0.eigenQREV();
				} 
				#pragma omp section
				{
					mScaleQ1  = mQ1.fillMatrix(kappa);
					mQ1.eigenQREV();
				}
			}

			// Initialize the variables used to refill the values next time
			mPrevK      = DBL_MAX;
			mPrevOmega0 = DBL_MAX;
			break;
		}

		// Compute the scale values
#ifdef USE_ORIGINAL_PROPORTIONS
		mFgScale = mProportions[0]*mScaleQw0 +
				   mProportions[1]*mScaleQ1  +
				   mProportions[2]*mScaleQ1  +
				   mProportions[3]*mScaleQ1;
		mBgScale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);
#else
		mFgScale = mProportions[0]*mScaleQw0 + (1.0-mProportions[0])*mScaleQ1;
		mBgScale = params[1]*mScaleQw0 + (1.0-params[1])*mScaleQ1;
#endif

		// Fill the set of Probability Matrices
		mSet.fillMatrixSet(mQw0, mQ1, mBgScale, mFgScale, aVar);

		// Compute likelihoods
		mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());
	}
	else
	{
		// Compute all the matrices for all branches. aVar contains the branch lengths varied each by its delta.
		if(aGradientVar == 0) mSetForGradient.fillMatrixSet(mQw0, mQ1, mBgScale, mFgScale, aVar);

		// Save and change one matrix
		mSet.saveMatrix(aGradientVar);
		mSet.setMatrices(aGradientVar, mSetForGradient.getChangedMatrices(aGradientVar));

		// Compute likelihoods
		mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());

		// Restore the previous value of the matrix
		mSet.restoreSavedMatrix(aGradientVar);
	}

	if(mExtraDebug > 0)
	{
		std::cout << "FG: " << std::setprecision(8) << mFgScale << " BG: " << mBgScale << std::endl;
		std::cout << "The following is the value printed by CodeML" << std::endl;
		std::cout << "FG: " << std::setprecision(8) << 1./mFgScale << " BG: " << 1./mBgScale << std::endl;
		std::cout << "Q0 " << mScaleQw0 << std::endl;
		std::cout << "Q1 " << mScaleQ1 << std::endl << std::endl;
	}

	// Combine the site likelihood into a single value
	double lnl = combineSiteLikelihoods();

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}

	return lnl;
}


double BranchSiteModelNullHyp::computeLikelihood(const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Save the values to local variables to speedup access
	const double* params = &aVar[mNumTimes];
	const double  omega0 = params[2];
	const double  kappa  = params[3];

	// Check if steps can be skipped
	const bool changed_w0 = isDifferent(omega0, mPrevOmega0);
	const bool changed_k  = isDifferent(kappa, mPrevK);
	if(changed_w0) mPrevOmega0 = omega0;
	if(changed_k)  mPrevK      = kappa;

	// Fill the matrices and compute their eigendecomposition.
	if(changed_k)
	{
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega0, kappa)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
				mQw0.eigenQREV();
			} 
			#pragma omp section
			{
				mScaleQ1  = mQ1.fillMatrix(kappa);
				mQ1.eigenQREV();
			}
		}
	}
	else if(changed_w0)
	{
		mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
		mQw0.eigenQREV();
	}

	// Compute all proportions
	getProportions(params[0], params[1], mProportions);

	// Compute the scale values
#ifdef USE_ORIGINAL_PROPORTIONS
	mFgScale = mProportions[0]*mScaleQw0 +
			   mProportions[1]*mScaleQ1  +
			   mProportions[2]*mScaleQ1  +
			   mProportions[3]*mScaleQ1;
	mBgScale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);
#else
	mFgScale = mProportions[0]*mScaleQw0 + (1.0-mProportions[0])*mScaleQ1;
	mBgScale = params[1]*mScaleQw0 + (1.0-params[1])*mScaleQ1;
#endif

	// Fill the set of Probability Matrices
	mSet.fillMatrixSet(mQw0, mQ1, mBgScale, mFgScale, aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());

	if(mExtraDebug > 0)
	{
		std::cout << "FG: " << std::setprecision(8) << mFgScale << " BG: " << mBgScale << std::endl;
		std::cout << "The following is the value printed by CodeML" << std::endl;
		std::cout << "FG: " << std::setprecision(8) << 1./mFgScale << " BG: " << 1./mBgScale << std::endl;
		std::cout << "Q0 " << mScaleQw0 << std::endl;
		std::cout << "Q1 " << mScaleQ1 << std::endl << std::endl;
	}

	// Combine the site likelihood into a single value
	double lnl = combineSiteLikelihoods();

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}

	return lnl;
}


double BranchSiteModelNullHyp::combineSiteLikelihoods(void)
{
	// Precompute the proportions to be used
	const double p0 = mProportions[0];
	const double p1_p2b = mProportions[1]+mProportions[3];
	const double p2a = mProportions[2];

	// For all (valid) sites. Don't parallelize: time increases and results are errant
	const size_t num_sites = mForest.getNumSites();
	const std::vector<double>& mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	double* likelihoods = &mLikelihoods[0];
	for(size_t site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		// double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		     (mProportions[1]+mProportions[3])*mLikelihoods[1*num_sites+site] +
		//		      mProportions[2]*mLikelihoods[2*num_sites+site];
		double p = likelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= p0;
		double x = likelihoods[1*num_sites+site];
		if(x > 0) p += p1_p2b*x;
		x = likelihoods[2*num_sites+site];
		if(x > 0) p += p2a*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];

		if(mExtraDebug > 1)
		{
			std::cout << std::setw(4) << site << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[0*num_sites+site] << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[1*num_sites+site] << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[2*num_sites+site] << " -> ";
			std::cout << std::fixed << std::setw(14) << x*mult[site] << std::endl;
		}
	}

	return lnl;
}

double BranchSiteModelAltHyp::computeLikelihoodForGradient(const std::vector<double>& aVar, bool aTrace, size_t aGradientVar)
{
	// One more function invocation
	++mNumEvaluations;

	// Compute the following values for gradient only if anything different from branch length has changed
	if(aGradientVar >= mNumTimes)
	{
		// Save the values to local variables to speedup access
		const double* params = &aVar[mNumTimes];
		const double  omega0 = params[2];
		const double  omega2 = params[4];
		const double  kappa  = params[3];

		// The values for gradient are computed in order, use this fact to reduce computations
		switch(aGradientVar - mNumTimes)
		{
		case 0:
		case 1:
			// Recompute all proportions if v0 or v1 change
			getProportions(params[0], params[1], mProportions);
			break;

		case 2:
			// Save Qw0 and Q1
#ifdef _MSC_VER
			#pragma omp parallel sections default(none) shared(params)
#else
			#pragma omp parallel sections default(shared)
#endif
			{
				#pragma omp section
				{
					mQw0.saveCheckpoint(mScaleQw0);
				} 
				#pragma omp section
				{
					mQ1.saveCheckpoint(mScaleQ1);
				}
				#pragma omp section
				{
					// Return to the original values
					getProportions(params[0], params[1], mProportions);
				}
			}
			mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
			mQw0.eigenQREV();
			break;

		case 3:
#ifdef _MSC_VER
			#pragma omp parallel sections default(none) shared(omega0, omega2, kappa)
#else
			#pragma omp parallel sections default(shared)
#endif
			{
				#pragma omp section
				{
					mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
					mQw0.eigenQREV();
				} 
				#pragma omp section
				{
					mScaleQ1  = mQ1.fillMatrix(kappa);
					mQ1.eigenQREV();
				}
				#pragma omp section
				{
					mScaleQw2 = mQw2.fillMatrix(omega2, kappa);
					mQw2.eigenQREV();
				}
			}
			break;

		case 4:
#ifdef _MSC_VER
			#pragma omp parallel sections default(none) shared(params)
#else
			#pragma omp parallel sections default(shared)
#endif
			{
				#pragma omp section
				{
					mScaleQw0 = mQw0.restoreCheckpoint();
				} 
				#pragma omp section
				{
					mScaleQ1 = mQ1.restoreCheckpoint();
				}
				#pragma omp section
				{
					mScaleQw2 = mQw2.fillMatrix(omega2, kappa);
					mQw2.eigenQREV();
				}
			}

			// Initialize the variables used to refill the values next time
			mPrevK      = DBL_MAX;
			mPrevOmega0 = DBL_MAX;
			mPrevOmega2 = DBL_MAX;
			break;
		}

		// Compute the scale values
#ifdef USE_ORIGINAL_PROPORTIONS
		mFgScale = mProportions[0]*mScaleQw0 +
				   mProportions[1]*mScaleQ1  +
				   mProportions[2]*mScaleQw2 +
				   mProportions[3]*mScaleQw2;
		mBgScale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);
#else
		mFgScale = mProportions[0]*mScaleQw0 + mProportions[1]*mScaleQ1 + (1.0-params[0])*mScaleQw2;
		mBgScale = params[1]*mScaleQw0+(1.0-params[1])*mScaleQ1;
#endif
		// Fill the set of Probability Matrices
		mSet.fillMatrixSet(mQw0, mQ1, mQw2, mBgScale, mFgScale, aVar);

		// Compute likelihoods
		mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());
	}
	else
	{
		// Compute all the matrices for all branches
		if(aGradientVar == 0) mSetForGradient.fillMatrixSet(mQw0, mQ1, mQw2, mBgScale, mFgScale, aVar);

		// Save and change one matrix
		mSet.saveMatrix(aGradientVar);
		mSet.setMatrices(aGradientVar, mSetForGradient.getChangedMatrices(aGradientVar));

		// Compute likelihoods
		mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());

		// Restore the previous value of the matrices
		mSet.restoreSavedMatrix(aGradientVar);
	}

	if(mExtraDebug > 0)
	{
		std::cout << "FG: " << std::setprecision(8) << mFgScale << " BG: " << mBgScale << std::endl;
		std::cout << "The following is the value printed by CodeML" << std::endl;
		std::cout << "FG: " << std::setprecision(8) << 1./mFgScale << " BG: " << 1./mBgScale << std::endl;
		std::cout << "Q0 " << mScaleQw0 << std::endl;
		std::cout << "Q1 " << mScaleQ1 << std::endl;
		std::cout << "Q2 " << mScaleQw2 << std::endl << std::endl;
	}

	// Combine the site likelihood into a single value
	double lnl = combineSiteLikelihoods();

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}

	return lnl;
}


double BranchSiteModelAltHyp::computeLikelihood(const std::vector<double>& aVar, bool aTrace)
{
	// One more function invocation
	++mNumEvaluations;

	// Save the values to local variables to speedup access
	const double* params = &aVar[mNumTimes];
	const double  omega0 = params[2];
	const double  omega2 = params[4];
	const double  kappa  = params[3];

	// Check if steps can be skipped
	const bool changed_w0 = isDifferent(omega0, mPrevOmega0);
	const bool changed_w2 = isDifferent(omega2, mPrevOmega2);
	const bool changed_k  = isDifferent(kappa, mPrevK);
	if(changed_w0) mPrevOmega0 = omega0;
	if(changed_w2) mPrevOmega2 = omega2;
	if(changed_k)  mPrevK      = kappa;

	// Fill the matrices and compute their eigendecomposition.
	if(changed_k)
	{
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega0, omega2, kappa)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
				mQw0.eigenQREV();
			} 
			#pragma omp section
			{
				mScaleQ1  = mQ1.fillMatrix(kappa);
				mQ1.eigenQREV();
			}
			#pragma omp section
			{
				mScaleQw2 = mQw2.fillMatrix(omega2, kappa);
				mQw2.eigenQREV();
			}
		}
	}
	else if(changed_w0 && changed_w2)
	{
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega0, omega2, kappa)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
				mQw0.eigenQREV();
			} 
			#pragma omp section
			{
				mScaleQw2 = mQw2.fillMatrix(omega2, kappa);
				mQw2.eigenQREV();
			}
		}
	}
	else
	{
		if(changed_w0)
		{
			mScaleQw0 = mQw0.fillMatrix(omega0, kappa);
			mQw0.eigenQREV();
		}
		if(changed_w2)
		{
			mScaleQw2 = mQw2.fillMatrix(omega2, kappa);
			mQw2.eigenQREV();
		}
	}

	// Compute all proportions
	getProportions(params[0], params[1], mProportions);

	// Compute the scale values
#ifdef USE_ORIGINAL_PROPORTIONS
	mFgScale = mProportions[0]*mScaleQw0 +
			   mProportions[1]*mScaleQ1  +
			   mProportions[2]*mScaleQw2 +
			   mProportions[3]*mScaleQw2;
	mBgScale = (mProportions[0]*mScaleQw0+mProportions[1]*mScaleQ1)/(mProportions[0]+mProportions[1]);
#else
	mFgScale = mProportions[0]*mScaleQw0 + mProportions[1]*mScaleQ1 + (1.0-params[0])*mScaleQw2;
	mBgScale = params[1]*mScaleQw0+(1.0-params[1])*mScaleQ1;
#endif

	// Fill the set of Probability Matrices
	mSet.fillMatrixSet(mQw0, mQ1, mQw2, mBgScale, mFgScale, aVar);

	// Compute likelihoods
	mForest.computeLikelihoods(mSet, mLikelihoods, mDependencies.getDependencies());

	if(mExtraDebug > 0)
	{
		std::cout << "FG: " << std::setprecision(8) << mFgScale << " BG: " << mBgScale << std::endl;
		std::cout << "The following is the value printed by CodeML" << std::endl;
		std::cout << "FG: " << std::setprecision(8) << 1./mFgScale << " BG: " << 1./mBgScale << std::endl;
		std::cout << "Q0 " << mScaleQw0 << std::endl;
		std::cout << "Q1 " << mScaleQ1 << std::endl;
		std::cout << "Q2 " << mScaleQw2 << std::endl << std::endl;
	}

	// Combine the site likelihood into a single value
	double lnl = combineSiteLikelihoods();

	// Output the trace message and update maxima found
	if(aTrace && lnl > mMaxLnL)
	{
		mMaxLnL = lnl;
		printVar(aVar, lnl);
	}

	return lnl;
}


double BranchSiteModelAltHyp::combineSiteLikelihoods(void)
{
	// Precompute the proportions to be used
	const double p0  = mProportions[0];
	const double p1  = mProportions[1];
	const double p2a = mProportions[2];
	const double p2b = mProportions[3];

	// For all (valid) sites. Don't parallelize: time increase and the results are errant
	const size_t num_sites = mForest.getNumSites();
	const std::vector<double>& mult = mForest.getSiteMultiplicity();
	double lnl = 0;
	double* likelihoods = &mLikelihoods[0];
	for(size_t site=0; site < num_sites; ++site)
	{
		// The following computation is split to avoid negative values
		//double p = mProportions[0]*mLikelihoods[0*num_sites+site] +
		//		     mProportions[1]*mLikelihoods[1*num_sites+site] +
		//		     mProportions[2]*mLikelihoods[2*num_sites+site] +
		//		     mProportions[3]*mLikelihoods[3*num_sites+site];
		//
		double p = likelihoods[0*num_sites+site];
		if(p < 0) p = 0;
		else      p *= p0;
		double x = likelihoods[1*num_sites+site];
		if(x > 0) p += p1*x;
		x = likelihoods[2*num_sites+site];
		if(x > 0) p += p2a*x;
		x = likelihoods[3*num_sites+site];
		if(x > 0) p += p2b*x;

		x = (p > 0) ? log(p) : mMaxLnL-100000;
		lnl += x*mult[site];

		if(mExtraDebug > 1)
		{
			std::cout << std::setw(4) << site << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[0*num_sites+site] << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[1*num_sites+site] << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[2*num_sites+site] << ' ';
			std::cout << std::scientific << std::setw(14) << likelihoods[3*num_sites+site] << " -> ";
			std::cout << std::fixed << std::setw(14) << x*mult[site] << std::endl;
		}
	}

	return lnl;
}

void BranchSiteModel::verifyOptimizerAlgo(unsigned int aOptimizationAlgo)
{
	switch(aOptimizationAlgo)
	{
	case OPTIM_LD_MING2:
	case OPTIM_LD_LBFGS:
	case OPTIM_LD_VAR1:
	case OPTIM_LD_VAR2:
	case OPTIM_LD_SLSQP:
	case OPTIM_LN_BOBYQA:
	case OPTIM_MLSL_LDS:
		return;

	default:
		throw FastCodeMLFatal("Invalid optimization algorithm identifier on the command line.");
	}
}


/// Adapter class to pass the routine to the optimizer.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-01-13 (initial version)
///     @version 1.0
///
///
class MaximizerFunction
{
public:
	/// Constructor.
	/// Saves the parameters for the routine.
	///
	/// @param[in] aModel The model to be evaluated
	/// @param[in] aTrace If set the optimization progress is traced
	/// @param[in] aUpper Upper limit for the variables (to constrain the gradient computation)
	/// @param[in] aDeltaForGradient The variable increment to compute gradient
	/// @param[in] aNumMatrixParams Number of variables besides the branch lengths
	/// @param[in] aStopIfBigger If true stop computation as soon as value is over aThreshold
	/// @param[in] aThreshold The threshold at which the maximization should be stopped
	///
	MaximizerFunction(BranchSiteModel* aModel, bool aTrace, const std::vector<double>& aUpper, double aDeltaForGradient, size_t aNumMatrixParams,
		              bool aStopIfBigger, double aThreshold)
		              : mModel(aModel), mTrace(aTrace), mUpper(aUpper), mDeltaForGradient(aDeltaForGradient),
					    mTotalNumVariables(aUpper.size()), mNumBranchLengths(aUpper.size() - aNumMatrixParams),
						mStopIfBigger(aStopIfBigger), mThreshold(aThreshold)
						{}

	/// Wrapper to be passed to the optimizer
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	/// @param[in] aData Opaque pointer containing the function to be passed to the optimizer
	///
	/// @return The evaluated function
	///
	static double wrapFunction(const std::vector<double>& aVars, std::vector<double>& aGrad, void* aData)
	{
		return (*reinterpret_cast<MaximizerFunction*>(aData))(aVars, aGrad);
	}

	/// Computes the function and the gradient if needed.
	///
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values
	///
	/// @return The evaluated function
	///
	/// @exception nlopt::forced_stop To force halt the maximization because LRT is already not satisfied
	///
	double operator()(const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		// Compute the function at the requested point
		double f0 = mModel->computeLikelihood(aVars, mTrace);

		// Stop optimization if value is greater or equal to threshold
		if(mStopIfBigger && f0 >= mThreshold) throw nlopt::forced_stop();

		// If requested compute the gradient 
		if(!aGrad.empty()) computeGradient(f0, aVars, aGrad);

		return f0;
	}


private:
	/// Compute the function gradient
	///
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void computeGradient(double aPointValue, const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		std::vector<double> vars_working_copy(aVars);
		std::vector<double> delta(mNumBranchLengths);

		// Modify all the branch lengths and compute the corresponding delta
		size_t i = 0;
		for(; i < mNumBranchLengths; ++i)
		{
			double eh = mDeltaForGradient * (vars_working_copy[i]+1.);

			// If it is going over the upper limit reverse the delta
			vars_working_copy[i] += eh;
			if(vars_working_copy[i] >= mUpper[i])
			{
				vars_working_copy[i] -= 2*eh;
				delta[i] = -eh;
			}
			else
			{
				delta[i] = eh;
			}
		}

		// Compute the partial derivative of the likelihood function for each (modified) branch length
		for(i=0; i < mNumBranchLengths; ++i)
		{
			const double f1 = mModel->computeLikelihoodForGradient(vars_working_copy, false, i);

			aGrad[i] = (f1-aPointValue)/delta[i];
		}

		// Return the branch lengths to their original values and modify the matrix parameters one at the time to compute the remaining gradient components
		vars_working_copy = aVars;
		for(; i < mTotalNumVariables; ++i)
		{
#ifdef USE_ORIGINAL_PROPORTIONS
			double eh = mDeltaForGradient * (fabs(aVars[i])+1.);
#else
			double eh = mDeltaForGradient * (aVars[i]+1.);
#endif

			vars_working_copy[i] += eh;
			if(vars_working_copy[i] >= mUpper[i]) {vars_working_copy[i] -= 2*eh; eh = -eh;}

			const double f1 = mModel->computeLikelihoodForGradient(vars_working_copy, false, i);

			aGrad[i] = (f1-aPointValue)/eh;

			vars_working_copy[i] = aVars[i];
		}
	}

#if 0
	/// Compute the function gradient (original method)
	///
	/// @param[in] aPointValue The value of the function at aVars
	/// @param[in] aVars Variables to be optimized
	/// @param[out] aGrad Gradient values computed
	///
	void computeGradient(double aPointValue, const std::vector<double>& aVars, std::vector<double>& aGrad) const
	{
		std::vector<double> vars_working_copy(aVars);
		for(size_t i=0; i < mTotalNumVariables; ++i)
		{
#ifdef USE_ORIGINAL_PROPORTIONS
			double eh = mDeltaForGradient * (fabs(aVars[i])+1.);
#else
			double eh = mDeltaForGradient * (aVars[i]+1.);
#endif
			vars_working_copy[i] += eh;
			if(vars_working_copy[i] >= mUpper[i]) {vars_working_copy[i] -= 2*eh; eh = -eh;}

			const double f1 = mModel->computeLikelihood(vars_working_copy, false);

			aGrad[i] = (f1-aPointValue)/eh;

			vars_working_copy[i] = aVars[i];
		}
	}
#endif

private:
	BranchSiteModel*	mModel;				///< Pointer to the model to be evaluated
	bool				mTrace;				///< If set traces the optimization progresses
	std::vector<double>	mUpper;				///< Upper limit of the variables to constrain the interval on which the gradient should be computed
	double				mDeltaForGradient;	///< The variable increment to compute gradient
	size_t				mTotalNumVariables;	///< Total number of variables to be used to compute gradient
	size_t				mNumBranchLengths;	///< Number of branch lengths (total number of variables to be used to compute gradient minus matrix params)
	bool				mStopIfBigger;		///< If true stop optimization as soon as function value is above mThreshold
	double				mThreshold;			///< Threshold value to stop optimization
};


double BranchSiteModel::maximizeLikelihood(size_t aFgBranch, bool aStopIfBigger, double aThreshold)
{
	// Print starting values
	if(mTrace)
	{
		std::cout << std::endl;
		std::cout << "*****************************************" << std::endl;
		std::cout << "*** Starting branch " << aFgBranch << std::endl;
		printVar(mVar);
		std::cout << "*** Upper" << std::endl;
		printVar(mUpperBound);
		std::cout << "*** Lower" << std::endl;
		printVar(mLowerBound);
		std::cout << std::endl;
	}

	// Initialize the maximum value found and the function evaluations counter
	mMaxLnL = VERY_LOW_LIKELIHOOD;
	mNumEvaluations = 0;

	// If only the initial step is requested, do it and return
	if(mOnlyInitialStep) return computeLikelihood(mVar, mTrace);

	// Special case for the CodeML optimizer
	if(mOptAlgo == OPTIM_LD_MING2)
	{
		try
		{
			// Create the optimizer (instead of mRelativeError is used the fixed value from CodeML)
			Ming2 optim(this, mTrace, mVerbose, mLowerBound, mUpperBound, mDeltaForGradient, 1e-8, aStopIfBigger, aThreshold, mMaxIterations);

			// Do the maximization
			double maxl = optim.minimizeFunction(mVar);
			
			if(mTrace)
			{
				std::cout << std::endl << "Function invocations:       " << mNumEvaluations << std::endl;
				std::cout <<              "Final log-likelihood value: " << maxl << std::endl;
				printVar(mVar);
			}
			return maxl;
		}
		catch(FastCodeMLEarlyStopLRT&)
		{
			if(mTrace) std::cout << "Optimization stopped because LRT not satisfied" << std::endl;
			return DBL_MAX;
		}
		catch(std::exception& e)
		{
			std::ostringstream o;
			o << "Exception in Ming2 computation: " << e.what() << std::endl;
			throw FastCodeMLFatal(o);
		}
	}

	// Select the maximizer algorithm (the listed ones works and are reasonably fast for FastCodeML)
	std::auto_ptr<nlopt::opt> opt;
	switch(mOptAlgo)
	{
	case OPTIM_LD_LBFGS:
		opt.reset(new nlopt::opt(nlopt::LD_LBFGS,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_VAR1:
		opt.reset(new nlopt::opt(nlopt::LD_VAR1,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_VAR2:
		opt.reset(new nlopt::opt(nlopt::LD_VAR2,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LD_SLSQP:
		opt.reset(new nlopt::opt(nlopt::LD_SLSQP,   mNumTimes+mNumVariables));
		opt->set_vector_storage(20);
		break;

	case OPTIM_LN_BOBYQA:
		opt.reset(new nlopt::opt(nlopt::LN_BOBYQA,  mNumTimes+mNumVariables));
		break;

	case OPTIM_MLSL_LDS:
		opt.reset(new nlopt::opt(nlopt::G_MLSL_LDS, mNumTimes+mNumVariables));
		{
		// For global optimization put a timeout of one day
		opt->set_maxtime(24*60*60);

		// This algorithm requires a local optimizer, add it
		nlopt::opt local_opt(nlopt::LN_BOBYQA, mNumTimes+mNumVariables);
		opt->set_local_optimizer(local_opt);
		}
		break;

	default:
		throw FastCodeMLFatal("Invalid optimization algorithm identifier on the command line.");
	}

	// Initialize bounds and termination criteria
	opt->set_lower_bounds(mLowerBound);
	opt->set_upper_bounds(mUpperBound);
    opt->set_ftol_rel(mRelativeError);
	nlopt::srand(static_cast<unsigned long>(mSeed));

	// Optimize the function
	double maxl = 0;
	try
	{
		MaximizerFunction compute(this, mTrace, mUpperBound, mDeltaForGradient, mNumVariables, aStopIfBigger, aThreshold);

		opt->set_max_objective(MaximizerFunction::wrapFunction, &compute);

		// If the user has set a maximum number of iterations set it
		if(mMaxIterations != MAX_ITERATIONS) opt->set_maxeval(mMaxIterations);

		nlopt::result result = opt->optimize(mVar, maxl);

		// Print the final optimum value
		if(mTrace)
		{
			std::cout << std::endl << "Function invocations:       " << mNumEvaluations << std::endl;
			switch(result)
			{
			case nlopt::SUCCESS:
				break;

			case nlopt::STOPVAL_REACHED:
				std::cout << "Optimization stopped because stopval was reached." << std::endl;
				break;

			case nlopt::FTOL_REACHED:
				std::cout << "Optimization stopped because ftol_rel or ftol_abs was reached." << std::endl;
				break;

			case nlopt::XTOL_REACHED:
				std::cout << "Optimization stopped because xtol_rel or xtol_abs was reached." << std::endl;
				break;

			case nlopt::MAXEVAL_REACHED:
				std::cout << "Optimization stopped because maxeval was reached." << std::endl;
				break;

			case nlopt::MAXTIME_REACHED:
				std::cout << "Optimization stopped because maxtime was reached." << std::endl;	
				break;

			default:
				std::cout << "Other reason: " << static_cast<unsigned int>(result) << std::endl;	
				break;
			}
			std::cout << "Final log-likelihood value: " << maxl << std::endl;
			printVar(mVar);
		}
	}
	catch(const nlopt::forced_stop&)
	{
		if(mTrace) std::cout << "Optimization stopped because LRT not satisfied" << std::endl;
		return DBL_MAX;
	}
	catch(const nlopt::roundoff_limited&)
	{
		throw FastCodeMLFatal("Exception in computation: Halted because roundoff errors limited progress, equivalent to NLOPT_ROUNDOFF_LIMITED.");
	}
	catch(const std::runtime_error&)
	{
		throw FastCodeMLFatal("Exception in computation: Generic failure, equivalent to NLOPT_FAILURE.");
	}
	catch(const std::invalid_argument&)
	{
		throw FastCodeMLFatal("Exception in computation: Invalid arguments (e.g. lower bounds are bigger than upper bounds, an unknown algorithm was specified, etcetera), equivalent to NLOPT_INVALID_ARGS.");
	}
	catch(const std::bad_alloc&)
	{
		throw FastCodeMLFatal("Exception in computation: Ran out of memory (a memory allocation failed), equivalent to NLOPT_OUT_OF_MEMORY.");
	}
	catch(const std::exception& e)
	{
		std::ostringstream o;
		o << "Exception in computation: " << e.what();
		throw FastCodeMLFatal(o);
	}

	return maxl;
}

/// @page vars_page Layout of free variables
/// The vector containing the independent variables has the following layout
///
/// @section blen_sect Branch lengths
/// The first `mNumTimes` positions contain the branch lengths. The index varies from `0` to `mNumTimes-1`.
///
/// @section v0_sect Combined proportions v0
/// The `v0 = (p0+p1)` value is at index `mNumTimes+0`
///
/// @section v1_sect Combined proportions v1
/// The `v1 = (p0/(p0+p1))` value is at index `mNumTimes+1`
///
/// @section w0_sect Omega 0
/// The `w0` value is at index `mNumTimes+2`
///
/// @section k_sect  Kappa
/// The `k` value is at index `mNumTimes+3`
///
/// @section w2_sect Omega 2
/// The `w2` value (if present) is at index `mNumTimes+4`
///
