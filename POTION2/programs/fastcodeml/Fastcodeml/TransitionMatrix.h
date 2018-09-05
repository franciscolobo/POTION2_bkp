
#ifndef TRANSITION_MATRIX_H
#define TRANSITION_MATRIX_H

#include <cstring>
#include <cmath>
#include <vector>
#include <bitset>
#include "MatrixSize.h"
#include "CompilerHints.h"
#include "CodonFrequencies.h"

#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

#ifdef USE_LAPACK
#include "blas.h"
#endif

/// The transition matrix plus its eigendecomposition.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-02-23 (initial version)
///     @version 1.0
///
class TransitionMatrix
{
public:
	/// Constructor.
	///
	TransitionMatrix()
	{
		// Initialize Q matrix to all zeroes (so only non-zero values are written)
#ifdef USE_LAPACK
		memset(mS, 0, N*N*sizeof(double));
#else
		memset(mQ, 0, N*N*sizeof(double));
#endif		
		// Initialize the codons' frequencies
		CodonFrequencies* cf = CodonFrequencies::getInstance();
		mCodonFreq = cf->getCodonFrequencies();
		mNumGoodFreq = static_cast<int>(cf->getNumGoodCodons());
		mSqrtCodonFreq = cf->getSqrtCodonFrequencies();
		cf->cloneGoodCodonIndicators(mGoodFreq);

#ifdef FORCE_IDENTITY_MATRIX
		// Fill the identity matrix (to be used when time is zero)
		memset(mIdentity, 0, N*N*sizeof(double));
		for(int i=0; i < N; ++i) mIdentity[i*(1+N)] = 1.0;
#endif
	}

	/// Fill the Q (or the S) matrix and return the matrix scale value.
	///
	/// @param[in] aOmega The omega value.
	/// @param[in] aK The k value.
	///
	/// @return The Q matrix scale value.
	///
	double fillMatrix(double aOmega, double aK);

	/// Fill the Q (or the S) matrix and return the matrix scale value. Optimized routine to be used for omega == 1
	///
	/// @param[in] aK The k value.
	///
	/// @return The Q matrix scale value.
	///
	double fillMatrix(double aK);

	/// Compute the eigendecomposition of the Q matrix.
	/// The used codon frequencies should be already loaded using setCodonFrequencies()
	/// The results are stored internally
	///
	void eigenQREV(void);

	/// Store in an external matrix the result of exp(Q*t)
	///
	/// @param[out] aOut The matrix where the result should be stored (size: N*N) under USE_LAPACK it is stored transposed
	/// @param[in] aT The time to use in the computation (it is always > 0)
	///
	void computeFullTransitionMatrix(double* RESTRICT aOut, double aT) const
	{
#ifdef FORCE_IDENTITY_MATRIX
		// if time is zero or almost zero, the transition matrix become an identity matrix
		if(aT < 1e-100)
		{
			memcpy(aOut, mIdentity, N*N*sizeof(double));
			return;
		}
#endif

#ifdef USE_LAPACK
// vdExp creates problems to one MPI implementation (even if the segfault seems harmless and happens after the program end)
#undef USE_MKL_VML
		double ALIGN64 tmp[N*N64]; // The rows are padded to 64 to increase performance
		double ALIGN64 expt[N];

		double tm = aT / 2.;
#ifndef USE_MKL_VML
		// Manual unrolling gives the best results here.
		// So it is exp(D*T/2). Remember, the eigenvalues are stored in reverse order
		for(int c=0; c < N-1; )
		{
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
			expt[c] = exp(tm*mD[N-1-c]); ++c;
		}
		expt[N-1] = exp(tm*mD[0]);

		for(int r=0; r < N; ++r)
		{
			for(int c=0; c < N; ++c)
			{
				tmp[r*N64+c] = expt[c]*mV[r*N+c];
				//tmp[r*N+c] = expt[c]*mV[r*N+c];
			}
		}
#else
		// Manual unrolling gives the best results here.
		// Remember, the eigenvalues are stored in reverse order
        for(int j=0; j < N-1; )
        {
            tmp[j] = tm*mD[N-1-j]; ++j; 
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
            tmp[j] = tm*mD[N-1-j]; ++j;
        }
		tmp[60] = tm*mD[0];

		vdExp(N, tmp, expt);
		for(int r=0; r < N; ++r)
		{
			vdMul(N, expt, &mV[r*N], &tmp[r*N64]);
			//vdMul(N, expt, &mV[r*N], &tmp[r*N]);
		}
#endif

		dsyrk_("U", "T", &N, &N, &D1, tmp, &N64, &D0, aOut, &N);
		//dsyrk_("U", "T", &N, &N, &D1, tmp, &N, &D0, aOut, &N);

#else
		// The first iteration of the loop (k == 0) is split out to initialize aOut
		double *p = aOut;
		double expt = exp(aT * mD[N-1]); // Remember, the eigenvalues are stored in reverse order
		for(int i=0; i < N; ++i)
		{
			const double uexpt = mU[i*N] * expt;

			for(int j=0; j < N; ++j)
			{
				*p++ = uexpt * mV[j];
			}
		}

		// The subsequent iterations are computed normally
		for(int k = 1; k < N; ++k)
		{
			p = aOut;
			expt = exp(aT * mD[N-1-k]); // Remember, the eigenvalues are stored in reverse order

			for(int i=0; i < N; ++i)
			{
				const double uexpt = mU[i*N + k] * expt;

				for(int j=0; j < N; ++j)
				{
					*p++ += uexpt * mV[k*N + j];
				}
			}
		}
#endif
	}


private:
	/// Compute the eigendecomposition
	///
	/// @param[in,out] aU Matrix to be decomposed on input, Eigenvectors on output
	/// @param[in] aDim The matrix dimension
	/// @param[out] aR The eigenvalues
	/// @param[out] aWork A working area used only for non lapack version
	///
	/// @exception std::range_error Error in EigenTridagQLImplicit (no lapack used)
	/// @exception FastCodeMLMemoryError Error sizing workareas
	/// @exception std::range_error No convergence in dsyevr
	///
	void eigenRealSymm(double* RESTRICT aU, int aDim, double* RESTRICT aR, double* RESTRICT aWork);


#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4324) // Padding added due to alignment request
#endif

protected:
	// Order suggested by icc to improve locality
	// 'mV, mU, mSqrtCodonFreq, mNumGoodFreq, mQ, mD, mCodonFreq, mGoodFreq'
	double ALIGN64	mV[N*N];		///< The right adjusted eigenvectors matrix (with the new method instead contains pi^1/2*R where R are the eigenvectors)
	double ALIGN64	mU[N*N];		///< The left adjusted eigenvectors matrix
#ifdef FORCE_IDENTITY_MATRIX 
	double ALIGN64	mIdentity[N*N];	///< Pre-filled identify matix
#endif
	const double*	mSqrtCodonFreq;	///< Square root of experimental codon frequencies
	int				mNumGoodFreq;	///< Number of codons whose frequency is not zero (must be int)
#ifndef USE_LAPACK
	double ALIGN64	mQ[N*N];		///< The Q matrix
#else
	double ALIGN64	mS[N*N];		///< The S matrix (remember Q = S*pi)
#endif
	double ALIGN64	mD[N];			///< The matrix eigenvalues stored in reverse order
	const double*	mCodonFreq;		///< Experimental codon frequencies
	std::bitset<N>	mGoodFreq;		///< True if the corresponding codon frequency is not small
	
#ifdef _MSC_VER
    #pragma warning(pop)
#endif
};


/// The transition matrix that can be saved ad restored afterwards.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-07 (initial version)
///     @version 1.0
///
class CheckpointableTransitionMatrix : public TransitionMatrix
{
public:
	/// Constructor.
	///
	CheckpointableTransitionMatrix() : TransitionMatrix(), mSavedScale(1.) {}
	
	/// Save a checkpoint of the matrices
	///
	/// @param[in] aScale The matrix scale to be saved
	///
	void saveCheckpoint(double aScale);

	/// Restore the status at the last checkpoint
	/// Note: no check if the saved status is valid
	///
	/// @return The original matrix scale
	///
	double restoreCheckpoint(void);

#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4324) // Padding added due to alignment request
#endif

private:
	double			mSavedScale;		///< Saved matrix scale value
	double ALIGN64	mSavedV[N*N];		///< The right adjusted eigenvectors matrix (with the new method instead contains pi^1/2*R where R are the eigenvectors)
	double ALIGN64	mSavedU[N*N];		///< The left adjusted eigenvectors matrix
#ifndef USE_LAPACK
	double ALIGN64	mSavedQ[N*N];		///< The Q matrix
#else
	double ALIGN64	mSavedS[N*N];		///< The S matrix (remember Q = S*pi)
#endif
	double ALIGN64	mSavedD[N];			///< The matrix eigenvalues stored in reverse order
	
#ifdef _MSC_VER
    #pragma warning(pop)
#endif
};


#endif

