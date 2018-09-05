
#ifndef PROBABILITYMATRIXSET_H
#define PROBABILITYMATRIXSET_H

#include <vector>
#include "MatrixSize.h"
#include "TransitionMatrix.h"
#include "AlignedMalloc.h"
#include "CodonFrequencies.h"
#include "MathSupport.h"
#include "Exceptions.h"

#ifdef USE_LAPACK
#include "blas.h"
#endif
#ifdef USE_MKL_VML
#include <mkl_vml_functions.h>
#endif

/// Set of probability matrices for all branches of a tree.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2011-04-05 (initial version)
///     @version 1.0
///
class ProbabilityMatrixSet
{
protected:
	/// Create matrix set. This is called only by subclasses.
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (is the number of branches of the tree)
	/// @param[in] aNumSets How many sets to allocate (one set is composed by the bg and fg matrices for one of the tree traversals)
	/// @param[in] aNumMatSets How many sets of matrices to allocate
	///
	/// @exception FastCodeMLMemoryError Cannot allocate mMatrixSpace or mMatrices
	///
	ProbabilityMatrixSet(size_t aNumMatrices, unsigned int aNumSets, unsigned int aNumMatSets)
		: mMatrixSpace(NULL), mMatrices(NULL), mNumMatrices(static_cast<int>(aNumMatrices)), mNumSets(aNumSets), mFgBranch(0)
	{
		mMatrixSpace  = static_cast<double*>(alignedMalloc(sizeof(double)*aNumMatSets*aNumMatrices*MATRIX_SLOT, CACHE_LINE_ALIGN));
		if(!mMatrixSpace) throw FastCodeMLMemoryError("Cannot allocate mMatrixSpace");
		mMatrices     = static_cast<double**>(alignedMalloc(sizeof(double*)*aNumSets*aNumMatrices, CACHE_LINE_ALIGN));
		if(!mMatrices) throw FastCodeMLMemoryError("Cannot allocate mMatrices");
#ifdef NEW_LIKELIHOOD
		mInvCodonFreq = CodonFrequencies::getInstance()->getInvCodonFrequencies();
#endif
	}

public:
	/// Destructor.
	///
	~ProbabilityMatrixSet()
	{
		alignedFree(mMatrixSpace);
		alignedFree(mMatrices);
	}

	/// Return the number of sets contained in this ProbabilityMatrixSet
	///
	/// @return The number of sets
	///
	unsigned int size(void) const {return mNumSets;}

	/// Initialize only the foreground branch value.
	/// It leaves the mMatrix array of pointers uninitialized. This is OK if the set is used only to store the matrices for later usage.
	///
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	///
	void initializeFgBranch(unsigned int aFgBranch) {mFgBranch = static_cast<int>(aFgBranch);} 

#ifndef NEW_LIKELIHOOD
	///	Multiply the aGin vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aGin The input vector to be multiplied by the matrix exponential
	/// @param[out] aGout The resulting vector
	///
	void doTransition(unsigned int aSetIdx, unsigned int aBranch, const double* aGin, double* aGout) const
	{
#ifdef USE_LAPACK
		dsymv_("U", &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);

		// The element wise multiplication has been moved up one level
		//#if !defined(BUNDLE_ELEMENT_WISE_MULT)
		//	elementWiseMult(aGout, mInvCodonFreq);
		//#endif
#else
		for(int r=0; r < N; ++r)
		{
			double x = 0;
			for(int c=0; c < N; ++c) x += mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+c]*aGin[c];
			aGout[r] = x;
		}
#endif
	}

	///	Multiply the aGin vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aCodon The leaf codon id (will supplant aGin)
	/// @param[out] aGout The resulting vector
	///
	void doTransitionAtLeaf(unsigned int aSetIdx, unsigned int aBranch, int aCodon, double* aGout) const
	{
#ifdef USE_LAPACK
		// Simply copy the symmetric matrix column
		// instead of: dsymv_("U", &N, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aGin, &I1, &D0, aGout, &I1);
		// The method is:
	    //	for(i=0; i < aCodon; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][aCodon*N+i];
		//	for(i=aCodon; i < N; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][i*N+aCodon];
		// The special cases below are for accellerate the routine
		switch(aCodon)
        {
        case 0:
            for(int i=0; i < N; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][i*N];
            break;
		case 1:
            aGout[0] = mMatrices[aSetIdx*mNumMatrices+aBranch][N];
            for(int i=1; i < N; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][i*N+1];
            break;
		default:
			memcpy(aGout, mMatrices[aSetIdx*mNumMatrices+aBranch]+aCodon*N, aCodon*sizeof(double));
			for(int i=aCodon; i < N; ++i) aGout[i] = mMatrices[aSetIdx*mNumMatrices+aBranch][i*N+aCodon];
			break;
		}

		// The element wise multiplication has been moved up one level
		//#if !defined(BUNDLE_ELEMENT_WISE_MULT)
		//		elementWiseMult(aGout, mInvCodonFreq);
		//#endif
#else
		for(int r=0; r < N; ++r)
		{
			aGout[r] = mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+aCodon];
		}
#endif
	}

#else

	///	Multiply the aMin fat-vector by the precomputed exp(Q*t) matrix
	///
	/// @param[in] aSetIdx Which set to use (starts from zero)
	/// @param[in] aBranch Which branch
	/// @param[in] aNumSites Number of sites composing the fat-vector
	/// @param[in] aMin The input fat-vector to be multiplied by the matrix exponential
	/// @param[out] aMout The resulting fat-vector
	///
	void doTransition(unsigned int aSetIdx, unsigned int aBranch, int aNumSites, const double* aMin, double* aMout) const
	{
#ifdef USE_LAPACK
	
	dsymm_("L", "U", &N, &aNumSites, &D1, mMatrices[aSetIdx*mNumMatrices+aBranch], &N, aMin, &VECTOR_SLOT, &D0, aMout, &VECTOR_SLOT);

#ifdef USE_MKL_VML
	for(int c=0; c < aNumSites; ++c)
	{
		vdMul(N, &aMout[c*VECTOR_SLOT], mInvCodonFreq, &aMout[c*VECTOR_SLOT]);
	}
#else
	for(int r=0; r < N; ++r)
	{
		dscal_(&aNumSites, &mInvCodonFreq[r], aMout+r, &VECTOR_SLOT);
	}
#endif

#else
		for(int r=0; r < N; ++r)
		{
			for(int c=0; c < aNumSites; ++c)
			{
				double x = 0;
				for(int k=0; k < N; ++k) x += mMatrices[aSetIdx*mNumMatrices+aBranch][r*N+k]*aMin[c*VECTOR_SLOT+k]; // aMin is transposed
				aMout[c*VECTOR_SLOT+r] = x; // also aMout is transposed
			}
		}
#endif
	}
#endif


protected:
	double*			mMatrixSpace;		///< Starts of the matrix storage area
	double**		mMatrices;			///< Access to the matrix set (contains pointers to mMatrixSpaces matrices)
#ifdef NEW_LIKELIHOOD
	const double*	mInvCodonFreq;		///< Inverse of the codon frequencies
#endif
	int				mNumMatrices;		///< Number of matrices in each set (should be int)
	unsigned int	mNumSets;			///< Number of sets
	int				mFgBranch;			///< Foreground branch number (should be int)
};


/// Set of probability matrices for all branches of a tree for the null hypothesis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-07 (initial version)
///     @version 1.0
///
class ProbabilityMatrixSetH0 : public ProbabilityMatrixSet
{
public:
	/// Create matrix set. It allocates 3 sets.
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (it is the number of branches of the tree)
	///
	explicit ProbabilityMatrixSetH0(size_t aNumMatrices) : ProbabilityMatrixSet(aNumMatrices, 3, 2) {}

	/// Initialize the set for a given foreground branch number for H0
	///
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	///
	void initializeSet(unsigned int aFgBranch);

	/// Compute the three sets of matrices for the H0 hypothesis.
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w1
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (p0+p1, p0/(p0+p1), w0, k, w2)
	///
	void fillMatrixSet(const TransitionMatrix& aQw0,
						    const TransitionMatrix& aQ1,
							double aSbg,
							double aSfg,
						    const std::vector<double>& aParams);

	/// Restore the previous value for the aBranch matrices.
	///
	/// @param[in] aBranch Branch for which the matrices should be restored
	///
	void restoreSavedMatrix(size_t aBranch);

	/// Save the previous value for the aBranch matrices.
	///
	/// @param[in] aBranch Branch for which the matrices should be saved
	///
	void saveMatrix(size_t aBranch);

	/// Access the matrices corresponding to the given branch.
	///
	/// @param[in] aBranch The branch for which the matrices are to be accessed
	///
	/// @return Pointer to an array of two pointers to the matrices to be accessed
	///
	const double** getChangedMatrices(size_t aBranch);

	/// Set the matrices for branch aBranch from the return value of a getChangedMatrices routine
	///
	/// @param[in] aBranch The branch for which the matrices are to be accessed
	/// @param[in] aMatricesPtr array of pointers as returned by the getChangedMatrices routine
	///
	void setMatrices(size_t aBranch, const double** aMatricesPtr);

private:
	const double*	mMatricesPtr[2];	///< Pointers to the changed matrices to be restored
	double			mSaveQw0[N*N];		///< Save the previous value for the Qw0 matrix
	double			mSaveQ1[N*N];		///< Save the previous value for the Q1 matrix
};


/// Set of probability matrices for all branches of a tree for the alternate hypothesis.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-07 (initial version)
///     @version 1.0
///
class ProbabilityMatrixSetH1 : public ProbabilityMatrixSet
{
public:
	/// Create matrix set
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (it is the number of branches of the tree)
	///
	explicit ProbabilityMatrixSetH1(size_t aNumMatrices) : ProbabilityMatrixSet(aNumMatrices, 4, 3) {}

	/// Initialize the set for a given foreground branch number for H1
	///
	/// @param[in] aFgBranch Number of the foreground branch (as branch number not as internal branch number!)
	///
	void initializeSet(unsigned int aFgBranch);

	/// Compute the four sets of matrices for the H1 hypothesis.
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w2
	/// - set 3: w1, w2
	///
	///	@param[in] aQw0 The mQw0 transition matrix
	///	@param[in] aQ1 The mQ1 transition matrix
	///	@param[in] aQw2 The mQw2 transition matrix
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (p0+p1, p0/(p0+p1), w0, k, w2)
	///
	void fillMatrixSet(const  TransitionMatrix& aQw0,
						    const  TransitionMatrix& aQ1,
						    const  TransitionMatrix& aQw2,
							double aSbg,
							double aSfg,
						    const std::vector<double>& aParams);
	
	/// Restore the previous value for the aBranch matrices.
	///
	/// @param[in] aBranch Branch for which the matrices should be restored
	///
	void restoreSavedMatrix(size_t aBranch);

	/// Save the previous value for the aBranch matrices.
	///
	/// @param[in] aBranch Branch for which the matrices should be saved
	///
	void saveMatrix(size_t aBranch);

	/// Access the matrices corresponding to the given branch.
	///
	/// @param[in] aBranch The branch for which the matrices are to be accessed
	///
	/// @return Pointer to an array of three pointers to the matrices to be accessed
	///
	const double** getChangedMatrices(size_t aBranch);

	/// Set the matrices for branch aBranch from the return value of a getChangedMatrices routine
	///
	/// @param[in] aBranch The branch for which the matrices are to be accessed
	/// @param[in] aMatricesPtr array of pointers as returned by the getChangedMatrices routine
	///
	void setMatrices(size_t aBranch, const double** aMatricesPtr);

private:
	const double*	mMatricesPtr[3];	///< Pointers to the changed matrices to be restored
	double			mSaveQw0[N*N];		///< Save the previous value for the Qw0 matrix
	double			mSaveQ1[N*N];		///< Save the previous value for the Q1 matrix
	double			mSaveQw2[N*N];		///< Save the previous value for the Qw2 matrix
};


/// Set of probability matrices for all branches of a tree for the BEB computation.
///
///     @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///     @date 2012-09-20 (initial version)
///     @version 1.0
///
class ProbabilityMatrixSetBEB : public ProbabilityMatrixSet
{
public:
	/// Create matrix set.
	///
	/// @param[in] aNumMatrices The number of matrices to be managed (it is the number of branches of the tree)
	///
	explicit ProbabilityMatrixSetBEB(size_t aNumMatrices) : ProbabilityMatrixSet(aNumMatrices, 1, 1)
	{
		for(int branch=0; branch < mNumMatrices; ++branch)
		{
			mMatrices[branch] = &mMatrixSpace[branch*MATRIX_SLOT];
		}
	}

	/// Compute the sets of matrices for the BEB computation.
	/// The sets are (these are the bg and fg matrices): 
	/// - set 0: w0, w0
	/// - set 1: w1, w1
	/// - set 2: w0, w2
	/// - set 3: w1, w2
	///
	///	@param[in] aQfg The transition matrix for the fg branch
	///	@param[in] aQbg The transition matrix for the bg branch
	/// @param[in] aSbg Background Q matrix scale
	/// @param[in] aSfg Foreground Q matrix scale
	/// @param[in] aParams Optimization parameters. First the branch lengths, then the variable parts (p0+p1, p0/(p0+p1), w0, k, w2)
	///
	void fillMatrixSet(const  TransitionMatrix& aQfg,
					   const  TransitionMatrix& aQbg,
					   double aSbg,
					   double aSfg,
					   const std::vector<double>& aParams);
};

#endif

