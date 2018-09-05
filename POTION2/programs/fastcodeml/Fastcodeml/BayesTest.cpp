
#include <iostream>
#include <iomanip>
#include <cmath>
#include "BayesTest.h"

double BayesTest::getGridParams(const std::vector<double>& aVars, const std::vector<double>& aSiteMultiplicity, size_t aFgBranch, const std::vector<double>& aScales)
{
	// Parameters for w0, w1, w2 prior computation
	double prior_params[BEB_DIMS][BEB_N1D];

	// Searching range for w0
	const static double w0b0 = 0.;
	const static double w0b1 = 1.;

	// Searching range for w2
	const static double w2b0 = 1.;
	const static double w2b1 = 11.;

	// Initialize the w0 and w2 values to be tested
	for(unsigned int i=0; i < BEB_N1D; i++)
	{
		prior_params[0][i] = prior_params[1][i] = -1.;				// p0 & p1
		prior_params[2][i] = w0b0 + (i+0.5)*(w0b1-w0b0)/BEB_N1D;	// w0
		prior_params[3][i] = w2b0 + (i+0.5)*(w2b1-w2b0)/BEB_N1D;	// w2
	}

	// Omega for foreground and background branches (the bools are for optimization)
	double omega_fg;
	double omega_bg;
	bool   omega_fg_is_one;
	bool   omega_bg_is_one;

	// Get the optimized kappa from the last H1 computation (that is, I know the kappa is the forth value after the branch lengths)
	size_t kappa_idx = mForest.getNumBranches() - 1 + 4;
	const double kappa = aVars[kappa_idx];

	// Compute the corresponding Q matrices for foreground and background branches
	TransitionMatrix q_fg, q_bg;

	// Initialize the probability list. Only the fg branch needs to be set
	mBEBset.initializeFgBranch(mForest.adjustFgBranchIdx(aFgBranch));

	// Calculating f(x_h|w) 
	// Order of site classes for iw or f(x_h|w):
	//                     back   fore     #sets
	//   site class 0:      w0     w0        10
	//   site class 1:      w1=1   w1=1       1
	//   site class 2a:     w0     w2       100
	//   site class 2b:     w1=1   w2        10
	//
	for(unsigned int iw=0; iw < BEB_NUM_CAT; ++iw)
	{
		if(iw < BEB_N1D)										// class 0: w0 w0
		{
			omega_bg = omega_fg = prior_params[2][iw];
			omega_fg_is_one = omega_bg_is_one = false;
		}
		else if(iw == BEB_N1D)									// class 1: w1 w1
		{
			omega_bg = omega_fg = 1.;
			omega_fg_is_one = omega_bg_is_one = true;
		}
		else if(iw < (BEB_N1D+1+BEB_N1D*BEB_N1D))				// class 2a: w0 w2
		{                                  
            omega_bg = prior_params[2][(iw-BEB_N1D-1) / BEB_N1D];
            omega_fg = prior_params[3][(iw-BEB_N1D-1) % BEB_N1D];
			omega_fg_is_one = omega_bg_is_one = false;
		}
		else													// class 2b: w1 w2
		{                                                       
            omega_bg = 1.;
            omega_fg = prior_params[3][iw-BEB_N1D-1-BEB_N1D*BEB_N1D];
			omega_fg_is_one = false; omega_bg_is_one = true;
		}

		// Fill the matrices and compute their eigendecomposition.
#ifdef _MSC_VER
		#pragma omp parallel sections default(none) shared(omega_fg, omega_bg, kappa, q_fg, q_bg, omega_fg_is_one, omega_bg_is_one)
#else
		#pragma omp parallel sections default(shared)
#endif
		{
			#pragma omp section
			{
				if(omega_fg_is_one)
					q_fg.fillMatrix(kappa);
				else
					q_fg.fillMatrix(omega_fg, kappa);

				q_fg.eigenQREV();
			} 
			#pragma omp section
			{
				if(omega_bg_is_one)
					q_bg.fillMatrix(kappa);
				else
					q_bg.fillMatrix(omega_bg, kappa);

				q_bg.eigenQREV();
			}
		}

		// Fill the matrix set with the new matrices and the times computed before (and scales from the latest H1 optimization)
		mBEBset.fillMatrixSet(q_fg, q_bg, aScales[0], aScales[1], aVars);

		// Compute likelihoods
		CacheAlignedDoubleVector likelihoods(mNumSites);
		mForest.computeLikelihoods(mBEBset, likelihoods, mDependencies.getDependencies());
		for(size_t site=0; site < mNumSites; ++site)
		{
			double p = likelihoods[site];
			mPriors[iw*mNumSites+site] = (p > 0) ? log(p) : -184.2068074395237; // If p < 0 then the value in codeml.c is: log(1e-80);
 		}
	}

	double scale = 0.;
	for(size_t site=0; site < mNumSites; ++site)
	{
		// Find the maximum likelihood for the given site between all the try values (121 values)
		double fh = mPriors[site];
		for(size_t k=1; k < BEB_NUM_CAT; ++k)
		{
			if(mPriors[k*mNumSites+site] > fh) fh = mPriors[k*mNumSites+site];
		}

		// Normalize the priors so they are less or equal to 0, then exponent to remove the previous log.
		for(size_t k=0; k < BEB_NUM_CAT; ++k)
		{
			mPriors[k*mNumSites+site] = exp(mPriors[k*mNumSites+site]-fh);
		}
		scale += fh*aSiteMultiplicity[site];
	}

	return scale;
}

void BayesTest::getIndexTernary(double* aProbX, double* aProbY, unsigned int aTriangleIdx)
{
	unsigned int ix = static_cast<unsigned int>(sqrt(static_cast<double>(aTriangleIdx)));
	unsigned int iy = aTriangleIdx - ix*ix;

	*aProbX = (1. + floor(iy/2.)*3. + (iy%2))/(3.*BEB_N1D);
	*aProbY = (1. + (BEB_N1D - 1 - ix)*3. + (iy%2))/(3.*BEB_N1D);
}

void BayesTest::computeBEB(const std::vector<double>& aVars, size_t aFgBranch, const std::vector<double>& aScales)
{
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << std::endl << "Computing BEB for fg branch " << aFgBranch << std::endl;

	// Get the site multiplicity
	const std::vector<double>& site_multiplicity = mForest.getSiteMultiplicity();

	// Prepare the priors list
	double scale1 = getGridParams(aVars, site_multiplicity, aFgBranch, aScales);

	// Set up iw and codon_class_proportion, for each igrid and codon_class
	unsigned int iw[BEB_NGRID*BEB_DIMS];
	double codon_class_proportion[BEB_NGRID*BEB_DIMS];
	for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
	{
		// Get one point on the grid
		unsigned int ip[BEB_DIMS];
		unsigned int it = igrid;
		for(int j=BEB_DIMS-1; j >= 0; --j)
		{
			ip[j]= it % BEB_N1D;
			it /= BEB_N1D;
		}

		// In the original code instead of BEB_DIMS there is nclassM. Both values are 4.
		for(unsigned int codon_class=0; codon_class < BEB_DIMS; ++codon_class)
		{
			// Given the point on the grid ip[] and codon_class, this returns iw and codon_class_proportion, 
			// where iw locates the correct f(x_h|w) stored in com.fhK[], and codon_class_proportion is 
			// the proportion of the site class under the model.
			// The BEB_N1D*BEB_N1D grid for p0-p1 is mapped onto the ternary graph for p0-p1-p2.  
			//
			unsigned int idx;
			switch(codon_class)
			{
			case 0: idx = ip[2]; break;								// class 0: w0
			case 1: idx = BEB_N1D; break;							// class 1: w1
			case 2: idx = BEB_N1D+1+ip[2]*BEB_N1D+ip[3]; break;		// class 2a model A: w0 & w2
			case 3: idx = BEB_N1D+1+BEB_N1D*BEB_N1D+ip[3]; break;	// class 2b model A: w1 & w2
#if BEB_DIMS > 4
			default: throw FastCodeMLFatal("Impossible case in computeBEB");
#else
			default: idx = 0; // To stop compiler complains
#endif
			}
			iw[igrid*BEB_DIMS+codon_class] = idx;

			// Fill the volume with the probabilities for each class
			double p[3];
			getIndexTernary(&p[0], &p[1], ip[0]*BEB_N1D+ip[1]);
			p[2] = 1.0 - p[0] - p[1];

			codon_class_proportion[igrid*BEB_DIMS+codon_class] = (codon_class < 2) ? p[codon_class] : p[2]*p[codon_class-2]/(1.0-p[2]);
		}
	}

	// Calculate marginal prob of data, fX, and postpara[].  scale2 is scale.
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << std::endl << "Calculating f(X), the marginal probability of data." << std::endl;
	double fX = 1.;
	double scale2 = -1e300;
	double lnfXs[BEB_NGRID];

	// Parameters for w0, w1, w2 posterior computation (postpara[0-1] for p0p1 ignored)
	double post_params[BEB_DIMS][BEB_N1D];
	for(unsigned int j=0; j < BEB_DIMS; ++j) for(unsigned int k=0; k < BEB_N1D; ++k) post_params[j][k] = 1.;

	double postp0p1[BEB_N1D*BEB_N1D];
	for(unsigned int k=0; k < BEB_N1D*BEB_N1D; k++) postp0p1[k] = 1.;

   	for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
	{
		// Get one point on the grid
		unsigned int ip[BEB_DIMS];
		unsigned int it = igrid;
		for(int j=BEB_DIMS-1; j >= 0; --j)
		{
			ip[j]= it % BEB_N1D;
			it /= BEB_N1D;
		}

		lnfXs[igrid] = 0.;
		for(unsigned int site=0; site < mNumSites; ++site)
		{
			double fh = 0.;
			for(unsigned int k=0; k < BEB_DIMS; ++k)
				fh += codon_class_proportion[igrid*BEB_DIMS+k]*mPriors[iw[igrid*BEB_DIMS+k]*mNumSites+site];
			if(fh < 1e-300)
			{
				if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << "Strange: f[" << site << "] = " << fh << " very small." << std::endl;
				continue;
			}
			lnfXs[igrid] += log(fh)*site_multiplicity[site];
		}

		double t = lnfXs[igrid]-scale2;
		if(t > 0)
		{
			/* change scale factor scale2 */
			t = (t<200) ? exp(-t) : 0;
			fX = fX*t+1;
			for(unsigned int j=0; j < BEB_DIMS; ++j) for(unsigned int k=0; k < BEB_N1D; ++k) post_params[j][k] *= t;
			for(unsigned int k=0; k < BEB_N1D*BEB_N1D; ++k) postp0p1[k] *= t;

			for(unsigned int j=0; j < BEB_DIMS; ++j) post_params[j][ip[j]]++;
			postp0p1[ip[0]*BEB_N1D+ip[1]]++;
			scale2 = lnfXs[igrid];
		}
		else if(t > -200)
		{
			t = exp(t);
			fX += t;
			for(unsigned int j=0; j < BEB_DIMS; ++j) post_params[j][ip[j]] += t;
			postp0p1[ip[0]*BEB_N1D+ip[1]] += t;
		}
	}

	// Normalize probabilities
	for(unsigned int j=0; j<BEB_DIMS; ++j) for(unsigned int k=0; k < BEB_N1D; ++k) post_params[j][k] /= fX;
	for(unsigned int k=0; k<BEB_N1D*BEB_N1D; k++) postp0p1[k] /=fX;

	// Print
	if(mVerbose >= VERBOSE_ONLY_RESULTS)
	{
		std::cout << std::endl << "Posterior on the grid" << std::endl << std::endl;
		const char* paras[5] = {"p0","p1","w0","w2","w3"};

		for(unsigned int j=2; j < BEB_DIMS; ++j)
		{
			std::cout << paras[j] << ": ";

			for(unsigned int k=0; k < BEB_N1D; ++k) std::cout << std::fixed << std::setw(7) << std::setprecision(3) << post_params[j][k];
			std::cout << std::endl;
		}

		std::cout << std::endl << "Posterior for p0-p1 (see the ternary graph)" << std::endl << std::endl;

		double sum_postp0p1 = 0.;
		for(unsigned int k=0; k<BEB_N1D*BEB_N1D; ++k)
		{
			std::cout << std::fixed << std::setw(6) << std::setprecision(3) << postp0p1[k];

			int sq = static_cast<int>(sqrt(k+1.));

			if(fabs(sq*sq-(k+1.))<1e-5) std::cout << std::endl;

			sum_postp0p1 += postp0p1[k];
		}
		std::cout << std::endl << "Sum of density on p0-p1 = " << std::setprecision(4) << sum_postp0p1 << std::endl << std::endl;
	}

	fX = log(fX)+scale2;

	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << "log(fX) = " << (fX+scale1-BEB_DIMS*log(BEB_N1D*1.))
		                                           << "  Scales = " << scale1 << " " << scale2 << std::endl;

	// Calculate posterior probabilities for sites. scale1 is scale factor
	if(mVerbose >= VERBOSE_ONLY_RESULTS) std::cout << std::endl << "Calculating f(w|X), posterior probs of site classes." << std::endl;

	for(unsigned int site=0; site < mNumSites; ++site)
	{
		scale1 = -1e300;

		for(unsigned int j=0; j < BEB_DIMS; ++j) mSiteClassProb[j*mNumSites+site] = 1.;

   		for(unsigned int igrid=0; igrid < BEB_NGRID; ++igrid)
		{
			double fh = 0.;
			double fhk[BEB_DIMS];
			unsigned int codon_class;
			for(codon_class=0; codon_class < BEB_DIMS; ++codon_class) /* duplicated calculation */
			{
				fhk[codon_class] = codon_class_proportion[igrid*BEB_DIMS+codon_class]*mPriors[iw[igrid*BEB_DIMS+codon_class]*mNumSites+site];
				fh += fhk[codon_class];
			}

			for(codon_class=0; codon_class < BEB_DIMS; ++codon_class)
			{
				fhk[codon_class] /= fh;

				// t is log of term on grid
				double t = log(fhk[codon_class]) + lnfXs[igrid];
				if(t > scale1 + 50)
				{ 
					// Change scale factor scale1
					for(unsigned int j=0; j < BEB_DIMS; ++j) mSiteClassProb[j*mNumSites+site] *= exp(scale1-t);
					scale1 = t;
				}
				mSiteClassProb[codon_class*mNumSites+site] += exp(t-scale1);
			}
		}
		for(unsigned int j=0; j<BEB_DIMS; ++j) mSiteClassProb[j*mNumSites+site] *= exp(scale1-fX);
	}
}


void BayesTest::printPositiveSelSites(size_t aFgBranch) const
{
	bool print_title = true;

	// Access the mapping from internal to original site numbers
	const std::multimap<size_t, size_t>& internal_to_original_site = mForest.getSitesMappingToOriginal();

	// To order the sites
	std::multimap<size_t, size_t> ordered_map;
	std::vector<double> probs;
	size_t current_idx = 0;

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if is a type 2 site with prob > 50%
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
		{
			// Put a title the firts time
			if(print_title)
			{
				std::cout << "Printing positive sel sites for branch " << aFgBranch << std::endl;
				print_title = false;
			}

			// Save site number and probability after mapping the site number to the original value
			std::pair<std::multimap<size_t, size_t>::const_iterator,std::multimap<size_t, size_t>::const_iterator> ret(internal_to_original_site.equal_range(site));
			std::multimap<size_t, size_t>::const_iterator it(ret.first);
			const std::multimap<size_t, size_t>::const_iterator end(ret.second);
			for(; it != end; ++it)
			{
				ordered_map.insert(std::pair<size_t, size_t>(it->second, current_idx));
				probs.push_back(prob);
				++current_idx;
			}
		}
	}

	// Print site number and probability after mapping the site number to the original value (and changing numbering so it starts from 1 and not zero)
	std::multimap<size_t, size_t>::const_iterator im(ordered_map.begin());
	std::multimap<size_t, size_t>::const_iterator endm(ordered_map.end());
	for(; im != endm; ++im)
	{
		double prob = probs[im->second];

		// Set significance
		const char* sig;
		if(prob > TWO_STARS_PROB)     sig = "**";
		else if(prob > ONE_STAR_PROB) sig = "*";
		else                          sig = "";

		std::cout << std::setw(6) << im->first + 1 << ' ' << std::fixed << std::setprecision(6) << prob << sig << std::endl;
	}
}


void BayesTest::extractPositiveSelSites(std::vector<unsigned int>& aPositiveSelSites, std::vector<double>& aPositiveSelSitesProb) const
{
	// Prepare the output vectors
	aPositiveSelSites.clear();
	aPositiveSelSitesProb.clear();

	// Access the mapping from internal to original site numbers
	const std::multimap<size_t, size_t>& internal_to_original_site = mForest.getSitesMappingToOriginal();

	// For all sites
	for(unsigned int site=0; site < mNumSites; ++site)
	{
		// Check if it is a type 2 site with prob > minimum cutoff
		double prob = mSiteClassProb[2*mNumSites+site] + mSiteClassProb[3*mNumSites+site];
		if(prob > MIN_PROB)
		{
			std::pair<std::multimap<size_t, size_t>::const_iterator,std::multimap<size_t, size_t>::const_iterator> ret(internal_to_original_site.equal_range(site));
			std::multimap<size_t, size_t>::const_iterator it(ret.first);
			const std::multimap<size_t, size_t>::const_iterator end(ret.second);
			for(; it != end; ++it)
			{
				aPositiveSelSites.push_back(static_cast<unsigned int>(it->second));
				aPositiveSelSitesProb.push_back(prob);
			}
		}
	}
}
