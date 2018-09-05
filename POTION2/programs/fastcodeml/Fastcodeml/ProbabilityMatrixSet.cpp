
#include "ProbabilityMatrixSet.h"

void ProbabilityMatrixSetH0::initializeSet(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		mMatrices[branch+num_matrices*0] = mMatrices[branch+num_matrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+num_matrices*1] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];
	}
	mMatrices[mFgBranch+num_matrices*2] = &mMatrixSpace[(num_matrices+mFgBranch)*MATRIX_SLOT];
}


void ProbabilityMatrixSetH1::initializeSet(unsigned int aFgBranch)
{
	mFgBranch = static_cast<int>(aFgBranch);

	int num_matrices = mNumMatrices;
	for(int branch=0; branch < num_matrices; ++branch)
	{
		mMatrices[branch+num_matrices*0] = &mMatrixSpace[branch*MATRIX_SLOT];
		mMatrices[branch+num_matrices*1] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];

		if(branch != mFgBranch)
		{
			mMatrices[branch+num_matrices*2] = &mMatrixSpace[branch*MATRIX_SLOT];
			mMatrices[branch+num_matrices*3] = &mMatrixSpace[(num_matrices+branch)*MATRIX_SLOT];
		}
		else
		{
			mMatrices[mFgBranch+num_matrices*2] = mMatrices[mFgBranch+num_matrices*3] = &mMatrixSpace[(2*num_matrices+mFgBranch)*MATRIX_SLOT];
		}
	}
}


void ProbabilityMatrixSetH0::fillMatrixSet(const TransitionMatrix& aQw0,
											 const  TransitionMatrix& aQ1,
											 double aSbg,
											 double aSfg,
											 const  std::vector<double>& aParams)
{
	const int num_matrices = mNumMatrices;
	const double* params = &aParams[0];

#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQw0, aQ1, aSbg, aSfg, params, num_matrices) schedule(guided)
#else
	#pragma omp parallel for default(shared)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
#ifndef _MSC_VER
		#pragma omp task untied
#endif
		{
		const double t = (branch == mFgBranch) ? params[branch]/aSfg : params[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		}
	}
}


void ProbabilityMatrixSetH0::restoreSavedMatrix(size_t aBranch)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQw0, N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQ1,  N*N*sizeof(double));
}

void ProbabilityMatrixSetH0::saveMatrix(size_t aBranch)
{
	memcpy(mSaveQw0, &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
	memcpy(mSaveQ1,  &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
}

const double** ProbabilityMatrixSetH0::getChangedMatrices(size_t aBranch)
{
	mMatricesPtr[0] = &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT];
	mMatricesPtr[1] = &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT];

	return mMatricesPtr;
}
	
void ProbabilityMatrixSetH0::setMatrices(size_t aBranch, const double** aMatricesPtr)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], aMatricesPtr[0], N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], aMatricesPtr[1], N*N*sizeof(double));
}

void ProbabilityMatrixSetH1::fillMatrixSet(const  TransitionMatrix& aQw0,
										   const  TransitionMatrix& aQ1,
										   const  TransitionMatrix& aQw2,
										   double aSbg,
										   double aSfg,
										   const  std::vector<double>& aParams)
{
	// To speedup access to variables
	const int num_matrices = mNumMatrices;
	const double* params = &aParams[0];

#ifdef _MSC_VER
	#pragma omp parallel default(none) shared(aQw0, aQ1, aQw2, aSbg, aSfg, params, num_matrices)
#else
	#pragma omp parallel default(shared)
#endif
	{
#ifdef _MSC_VER
	#pragma omp for schedule(guided) nowait
#else
	#pragma omp for nowait
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
#ifndef _MSC_VER
		#pragma omp task untied
#endif
		{
		const double t = (branch == mFgBranch) ? params[branch]/aSfg : params[branch]/aSbg;

		aQw0.computeFullTransitionMatrix(&mMatrixSpace[0*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		aQ1.computeFullTransitionMatrix( &mMatrixSpace[1*num_matrices*MATRIX_SLOT+branch*MATRIX_SLOT], t);
		}
	}

#pragma omp single
	{
#ifndef _MSC_VER
		#pragma omp task untied
#endif
		aQw2.computeFullTransitionMatrix(&mMatrixSpace[2*num_matrices*MATRIX_SLOT+mFgBranch*MATRIX_SLOT], params[mFgBranch]/aSfg);
	}
	}
}


void ProbabilityMatrixSetH1::restoreSavedMatrix(size_t aBranch)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQw0, N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], mSaveQ1,  N*N*sizeof(double));
	if(static_cast<int>(aBranch) == mFgBranch)
	{
		memcpy(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], mSaveQw2, N*N*sizeof(double));
	}
}

void ProbabilityMatrixSetH1::saveMatrix(size_t aBranch)
{
	memcpy(mSaveQw0, &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
	memcpy(mSaveQ1,  &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], N*N*sizeof(double));
	if(static_cast<int>(aBranch) == mFgBranch)
	{
		memcpy(mSaveQw2, &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], N*N*sizeof(double));
	}
}

const double** ProbabilityMatrixSetH1::getChangedMatrices(size_t aBranch)
{
	mMatricesPtr[0] = &mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT];
	mMatricesPtr[1] = &mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT];
	if(static_cast<int>(aBranch) == mFgBranch)
	{
		mMatricesPtr[2] = &mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT];
	}
	return mMatricesPtr;
}

void ProbabilityMatrixSetH1::setMatrices(size_t aBranch, const double** aMatricesPtr)
{
	memcpy(&mMatrixSpace[0*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], aMatricesPtr[0], N*N*sizeof(double));
	memcpy(&mMatrixSpace[1*mNumMatrices*MATRIX_SLOT+aBranch*MATRIX_SLOT], aMatricesPtr[1], N*N*sizeof(double));
	if(static_cast<int>(aBranch) == mFgBranch)
	{
		memcpy(&mMatrixSpace[(2*mNumMatrices+mFgBranch)*MATRIX_SLOT], aMatricesPtr[2], N*N*sizeof(double));
	}
}


void ProbabilityMatrixSetBEB::fillMatrixSet(const TransitionMatrix& aQfg, const TransitionMatrix& aQbg, double aSbg, double aSfg, const std::vector<double>& aParams)
{
	const int num_matrices = mNumMatrices;
	const double* params = &aParams[0];

#ifdef _MSC_VER
	#pragma omp parallel for default(none) shared(aQfg, aQbg, aSbg, aSfg, params, num_matrices) schedule(guided)
#else
	#pragma omp parallel for default(shared)
#endif
	for(int branch=0; branch < num_matrices; ++branch)
	{
#ifndef _MSC_VER
		#pragma omp task untied
#endif
		{
			if(branch == mFgBranch)
			{
				aQfg.computeFullTransitionMatrix(&mMatrixSpace[branch*MATRIX_SLOT], params[branch]/aSfg);
			}
			else
			{
				aQbg.computeFullTransitionMatrix(&mMatrixSpace[branch*MATRIX_SLOT], params[branch]/aSbg);
			}
		}
	}
}
