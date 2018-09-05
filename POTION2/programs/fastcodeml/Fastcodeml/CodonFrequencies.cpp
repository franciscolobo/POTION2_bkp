
#include <iostream>
#include <iomanip>
#include <cstring>
#include "CodonFrequencies.h"
#include "Exceptions.h"

/// Max number of iterations to resolve codon frequencies in presence of ambiguous codons.
static const int    CODON_FREQ_MAX_ITER  = 20;

/// Max cumulative difference between an iteration and the next to stop iterations to resolve codon frequencies in presence of ambiguous codons.
static const double CODON_FREQ_MAX_ERROR = 1e-12;

CodonFrequencies* CodonFrequencies::mInstance = NULL;

CodonFrequencies* CodonFrequencies::getInstance(void)
{
	if(!mInstance) mInstance = new CodonFrequencies;
	return mInstance;
}

void CodonFrequencies::setCodonFrequencies(const std::vector<std::vector<unsigned int> >& aCodons, CodonFrequencyModel aModel, bool aShowMessages)
{
	// Compute mCodonFrequencies based on the selected model
	if(aModel == CODON_FREQ_MODEL_F3X4)
	{
		// Check if there is at least one ambiguous codon
		size_t cnt = aCodons.size();
		bool no_ambiguous_codons = true;
		for(size_t i=0; i < cnt; ++i) if(aCodons[i].size() > 2) {no_ambiguous_codons = false; break;}

		// Compute the usual way if no ambiguities (count multiplied by site multiplicity)
		if(no_ambiguous_codons)
		{
			std::vector<double> codon_count(N, 0.);
			for(size_t i=0; i < cnt; ++i) codon_count[aCodons[i][1]] += static_cast<double>(aCodons[i][0]);

			setCodonFrequenciesF3x4(codon_count, aShowMessages);
		}
		else
		{
			// Start counting ignoring ambiguous codons
			std::vector<double> codon_count(N, 0.);
			for(size_t i=0; i < cnt; ++i) if(aCodons[i].size() == 2) codon_count[aCodons[i][1]] += static_cast<double>(aCodons[i][0]);

			// Compute a first mCodonFrequencies array
			setCodonFrequenciesF3x4(codon_count, false);

			// Iterate till the values does not change
			for(int k=0; k < CODON_FREQ_MAX_ITER; ++k)
			{
				// Save the previous value
				SSEAlignedDoubleVector prev_codon_freq(mCodonFrequencies);

				// Recount the codons
				updateCodonCount(aCodons, codon_count);

				// Recompute the frequencies
				setCodonFrequenciesF3x4(codon_count, aShowMessages);

				// Compute cumulative change over previous value
				double err = 0.;
				for(int h=0; h < N; ++h) err += fabs(prev_codon_freq[h] - mCodonFrequencies[h]);
				if(aShowMessages) std::cout << "  mCodonFrequencies change after iteration " << k << " = " << std::scientific << err << std::fixed << std::endl << std::endl;

				// If the frequencies values do not change, stop
				if(err < CODON_FREQ_MAX_ERROR) break;
			}
		}
	}
	else if(aModel == CODON_FREQ_MODEL_UNIF)
	{
		mCodonFrequencies.assign(N, 1./static_cast<double>(N));
	}
	else
	{
		throw FastCodeMLFatal("Invalid codon frequency model requested.");
	}

	// Support values needed for the eigensolver
	mNumGoodCodons = 0;
	for(size_t k=0; k < static_cast<size_t>(N); ++k)
	{
		mCodonFreqSqrt[k] = sqrt(mCodonFrequencies[k]);

		// Count the number of valid codons
		if(mCodonFrequencies[k] > GOOD_CODON_THRESHOLD)
		{
			mGoodCodon.set(k);
			++mNumGoodCodons;
			mCodonFreqInv[k]  = 1./mCodonFrequencies[k];
			mCodonFreqInv2[k] = mCodonFreqInv[k]*mCodonFreqInv[k];
		}
		else
		{
			mGoodCodon.reset(k);
			mCodonFreqInv[k]  = 0.; // To have zero in non valid positions so vector norm does not diverge
			mCodonFreqInv2[k] = 0.;
		}
	}
}

void CodonFrequencies::updateCodonCount(const std::vector<std::vector<unsigned int> >& aCodons, std::vector<double>& aCodonCount) const
{
	// Zero the count array
	aCodonCount.assign(N, 0.);

	// Count again
	size_t cnt = aCodons.size();
	for(size_t i=0; i < cnt; ++i) 
	{
		//if no ambiguities count as usual (add site multiplicity)
		if(aCodons[i].size() == 2)
		{
			aCodonCount[aCodons[i][1]] += static_cast<double>(aCodons[i][0]);
		}
		else
		{
			// If ambiguous add the codons resolved from the ambiguity in proportion to their frequencies
			std::vector<double> perc;
			size_t np = aCodons[i].size()-1;
			for(size_t j=1; j <= np; ++j) perc.push_back(mCodonFrequencies[aCodons[i][j]]);
			double t = 0.;
			for(size_t j=0; j < np; ++j) t += perc[j];

			for(size_t j=0; j < np; ++j) aCodonCount[aCodons[i][j+1]] += (static_cast<double>(aCodons[i][0]) * perc[j]/t);
		}
	}
}

int CodonFrequencies::codon64to61(int aId64) const
{
	if(aId64 > 63 || aId64 == 10 || aId64 == 11 || aId64 == 14) return -1;

	if(aId64 > 14) return aId64-3;
	if(aId64 > 11) return aId64-2;
	return aId64;
}


void CodonFrequencies::setCodonFrequenciesF3x4(const std::vector<double>& aCodonCount, bool aShowMessages)
{
	int k, j;

	if(aShowMessages)
	{
		// Print the table of codon counts
		for(k=0; k < N64; ++k)
		{
			int id = codon64to61(k);
			if(id < 0) std::cout << std::setw(12) << 0;
			else       std::cout << std::setw(12) << std::setprecision(3) << std::fixed << aCodonCount[id];
			if(k % 4 == 3) std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// Compute the 3x4 table
	double fb3x4sg[12];

	memset(fb3x4sg, 0, 12*sizeof(double));

    for(k = 0; k < N64; k++)
    {
		int kk = codon64to61(k);
		if(kk < 0) continue;

        fb3x4sg[0 * 4 + k / 16]      += aCodonCount[kk];
        fb3x4sg[1 * 4 + (k / 4) % 4] += aCodonCount[kk];
        fb3x4sg[2 * 4 + k % 4]       += aCodonCount[kk];
    }

    for(j = 0; j < 3; j++)
    {
        double t = 0;
		for(k=0; k < 4; ++k) t += fb3x4sg[j*4+k];
		for(k=0; k < 4; ++k) fb3x4sg[j*4+k] /= t;
    }

	if(aShowMessages)
	{
		for(k=0; k < 12; ++k)
		{
			std::cout << std::fixed << std::setprecision(5) << std::setw(12) << fb3x4sg[k];
			if(k % 4 == 3) std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	// Compute codon frequency from the 3x4 table
	for(k=0; k < N64; ++k)
	{
		int kk = codon64to61(k);
		if(kk < 0) continue;

		mCodonFrequencies[kk] = fb3x4sg[k / 16] * fb3x4sg[4 + (k / 4) % 4] * fb3x4sg[8 + k % 4];
	}
	double t = 0;
	for(k=0; k < N; ++k) t += mCodonFrequencies[k];
	for(k=0; k < N; ++k) mCodonFrequencies[k] /= t;

	if(aShowMessages)
	{
		for(k=0; k < N64; ++k)
		{
			int kk = codon64to61(k);
			double v = (kk < 0) ? 0.0 : mCodonFrequencies[kk];
			std::cout << std::fixed << std::setprecision(8) << std::setw(12) << v;
			if(k % 4 == 3) std::cout << std::endl;
		}
	}
}

