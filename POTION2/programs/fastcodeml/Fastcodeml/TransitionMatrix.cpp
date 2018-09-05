
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "MatrixSize.h"

#include "TransitionMatrix.h"

// Now the size of the workareas for DSYEVR are hardcoded.
// If you want to reintroduce the optimal size computation, then uncomment the following line.
//#define OPTIMAL_WORKAREAS


#ifndef USE_LAPACK

#ifdef _MSC_VER
//static inline double copysign(double a, double b) {return (b >= 0.0) ? fabs(a) : -fabs(a);}
//static inline double copysign(double x, double y) {return ((x < 0 && y > 0) || (x > 0 && y < 0)) ? -x : x;}
#define copysign _copysign
#endif

static void EigenSort(double d[], double U[], int n)
{
    // This sorts the eigenvalues d[] and rearrange the (right) eigen vectors U[]
    int k, j, i;
    double p;

    for(i = 0; i < n - 1; i++)
    {
        p = d[k = i];

        for(j = i + 1; j < n; j++)
		{
            if(d[j] >= p)
            {
                p = d[k = j];
            }
		}

        if(k != i)
        {
            d[k] = d[i];
            d[i] = p;

            for(j = 0; j < n; j++)
            {
                p = U[j * n + i];
                U[j *n + i] = U[j * n + k];
                U[j *n + k] = p;
            }
        }
    }
}


static void HouseholderRealSym(double a[], int n, double d[], double e[])
{
    /* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
       a[n*n] into a tridiagonal matrix represented by d and e.
       d[] is the diagonal (eigends), and e[] the off-diagonal.
    */
    int m, k, j, i;
    double scale, hh, h, g, f;

    for (i = n - 1; i >= 1; i--)
    {
        m = i - 1;
        h = scale = 0;

        if (m > 0)
        {
            for (k = 0; k <= m; k++)
            {
                scale += fabs(a[i * n + k]);
            }

            if (scale == 0)
            {
                e[i] = a[i * n + m];
            }
            else
            {
                for (k = 0; k <= m; k++)
                {
                    a[i *n + k] /= scale;
                    h += a[i * n + k] * a[i * n + k];
                }

                f = a[i * n + m];
                g = (f >= 0 ? -sqrt(h) : sqrt(h));
                e[i] = scale * g;
                h -= f * g;
                a[i *n + m] = f - g;
                f = 0;

                for (j = 0; j <= m; j++)
                {
                    a[j *n + i] = a[i * n + j] / h;
                    g = 0;

                    for (k = 0; k <= j; k++)
                    {
                        g += a[j * n + k] * a[i * n + k];
                    }

                    for (k = j + 1; k <= m; k++)
                   {
                        g += a[k * n + j] * a[i * n + k];
                    }

                    e[j] = g / h;
                    f += e[j] * a[i * n + j];
                }

                hh = f / (h * 2);

                for (j = 0; j <= m; j++)
                {
                    f = a[i * n + j];
                    e[j] = g = e[j] - hh * f;

                    for (k = 0; k <= j; k++)
                    {
                        a[j *n + k] -= (f * e[k] + g * a[i * n + k]);
                    }
                }
            }
        }
        else
        {
            e[i] = a[i * n + m];
        }

        d[i] = h;
    }

    d[0] = e[0] = 0;

    /* Get eigenvectors */
    for (i = 0; i < n; i++)
    {
        m = i - 1;

        if (d[i])
        {
            for (j = 0; j <= m; j++)
            {
                g = 0;

                for (k = 0; k <= m; k++)
                {
                    g += a[i * n + k] * a[k * n + j];
                }

                for (k = 0; k <= m; k++)
                {
                    a[k *n + j] -= g * a[k * n + i];
                }
            }
        }

        d[i] = a[i * n + i];
        a[i *n + i] = 1;

        for (j = 0; j <= m; j++)
        {
            a[j *n + i] = a[i *n + j] = 0;
        }
    }
}

static int EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
    /* This finds the eigensolution of a tridiagonal matrix represented by d and e.
       d[] is the diagonal (eigenvalues), e[] is the off-diagonal
       z[n*n]: as input should have the identity matrix to get the eigensolution of the
       tridiagonal matrix, or the output from HouseholderRealSym() to get the
       eigensolution to the original real symmetric matrix.
       z[n*n]: has the orthogonal matrix as output

       Adapted from routine tqli in Numerical Recipes in C, with reference to
       LAPACK fortran code.
       Ziheng Yang, May 2001
    */
    int m, j, iter, niter = 30, status = 0, i, k;
    double s, r, p, g, f, dd, c, b, aa, bb;

    for (i = 1; i < n; i++)
    {
        e[i - 1] = e[i];
    }

    e[n - 1] = 0;

    for (j = 0; j < n; j++)
    {
        iter = 0;

        do
        {
            for (m = j; m < n - 1; m++)
            {
                dd = fabs(d[m]) + fabs(d[m + 1]);

                if (fabs(e[m]) + dd == dd)
                {
                    break;    /* ??? */
                }
            }

            if (m != j)
            {
                if (iter++ == niter)
                {
                    status = -1;
                    break;
                }

                g = (d[j + 1] - d[j]) / (2 * e[j]);

                /* r=pythag(g,1); */

                if ((aa = fabs(g)) > 1)
                {
                    r = aa * sqrt(1 + 1 / (g * g));
                }
                else
                {
                    r = sqrt(1 + g * g);
                }

                g = d[m] - d[j] + e[j] / (g + copysign(r, g));
                s = c = 1;
                p = 0;

                for (i = m - 1; i >= j; i--)
                {
                    f = s * e[i];
                    b = c * e[i];

                    /*  r=pythag(f,g);  */
                    aa = fabs(f);
                    bb = fabs(g);

                    if (aa > bb)
                    {
                        bb /= aa;
                        r = aa * sqrt(1 + bb * bb);
                    }
                    else if (bb == 0)
                    {
                        r = 0;
                    }
                    else
                    {
                        aa /= bb;
                        r = bb * sqrt(1 + aa * aa);
                    }

                    e[i + 1] = r;

                    if (r == 0)
                    {
                        d[i + 1] -= p;
                        e[m] = 0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2 * c * b;
                    d[i + 1] = g + (p = s * r);
                    g = c * r - b;

                    for (k = 0; k < n; k++)
                    {
                        f = z[k * n + i + 1];
                        z[k *n + i + 1] = s * z[k * n + i] + c * f;
                        z[k *n + i] = c * z[k * n + i] - s * f;
                    }
                }

                if (r == 0 && i >= j)
                {
                    continue;
                }

                d[j] -= p;
                e[j] = g;
                e[m] = 0;
            }
        }
        while(m != j);
    }

    return status;
}

void TransitionMatrix::eigenRealSymm(double* aU, int aDim, double* aR, double* aWork)
{
    /* This finds the eigensolution of a real symmetrical matrix aU[aDim*aDim]. In return,
       aU has the right vectors and aR has the eigenvalues.
       aWork[n] is the working space.
       The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
       and then using the QL algorithm with implicit shifts.

       Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
       Ziheng Yang, 23 May 2001
    */
    HouseholderRealSym(aU, aDim, aR, aWork);
    int sts = EigenTridagQLImplicit(aR, aWork, aDim, aU);
	if(sts < 0) throw std::range_error("Error in EigenTridagQLImplicit");
    EigenSort(aR, aU, aDim);

    // Reorder eigenvalues so they are stored in reverse order
    const int mid = aDim/2;
    for(int i=0; i < mid; ++i) { double t = aR[i]; aR[i] = aR[aDim-1-i]; aR[aDim-1-i] = t; }
}

#else

#include "blas.h"
#include "lapack.h"

//
//  Using LAPACK DSYEVR driver routine to compute the eigenvalues and eigenvector of the symmetric input matrix aU[n*n]
//  Reorders the output values so they are ordered as the ones computed by the original eigenRealSym() routine.
//
void TransitionMatrix::eigenRealSymm(double* aU, int aDim, double* aR, double* /* aIgnored */)
{
    int m;
    int info;
    int isuppz[2*N];
    double ALIGN64 tmp_u[N*N64];

#ifndef OPTIMAL_WORKAREAS
	// Allocate fixed workareas
    static const int lwork = 33*N;
    double work[lwork];
    static const int liwork = 10*N;
    int iwork[liwork];
#else
	// Allocate fixed workareas
    static const int lfwork = 33*N;
    double fwork[lfwork];
    static const int lfiwork = 10*N;
    int fiwork[lfiwork];

	// Prepare for getting the optimal sizes
    int lwork=-1, liwork=-1;
    double *work = fwork;
    int *iwork = fiwork;

	// Compute the optimal size of the workareas
	double opt_work;
	int opt_iwork;
    dsyevr_("V", "A", "U", &aDim, aU, &aDim, &D0, &D0, &I0, &I0, &D0, &m, aR, tmp_u, &N64, isuppz, &opt_work, &lwork, &opt_iwork, &liwork, &info);
	if(info != 0) throw FastCodeMLMemoryError("Error sizing workareas");

	// Notice that LAPACK stores an integer value in a double array
    lwork = static_cast<int>(opt_work);
    liwork = opt_iwork;

	if(lwork > lfwork)
	{
		work = new double[lwork];
		//std::cout << "Optimal work:  " << lwork << " (" << lfwork << ")" << std::endl;
		
	}
    if(liwork > lfiwork)
	{
		iwork = new int[liwork];
		//std::cout << "Optimal iwork: " << liwork << " (" << lfiwork << ")" << std::endl;
	}
#endif
    // Compute eigenvalues and eigenvectors for the full symmetric matrix
    dsyevr_("V", "A", "U", &aDim, aU, &aDim, &D0, &D0, &I0, &I0, &D0, &m, aR, tmp_u, &N64, isuppz, work, &lwork, iwork, &liwork, &info);

#ifdef OPTIMAL_WORKAREAS
	// Release workareas, if allocated
	if(lwork > lfwork)   delete [] work;
	if(liwork > lfiwork) delete [] iwork;
#endif

	// Check convergence
	if(info > 0) throw std::range_error("No convergence in dsyevr");
	//if(info < 0) throw std::invalid_argument("Invalid parameter to dsyevr");

	// Reorder eigenvectors (instead the eigenvalues are stored in reverse order)
    for(int c=0; c < aDim; ++c)
    {
        for(int r=0; r < aDim; ++r)
        {
            //aU[r*aDim+c] = tmp_u[(aDim-1-c)*aDim+r];
            aU[r*aDim+c] = tmp_u[(aDim-1-c)*N64+r];
        }
    }
}

#endif

#ifdef USE_LAPACK

void TransitionMatrix::eigenQREV(void)
{
	  /*
       This finds the eigensolution of the rate matrix Q for a time-reversible
       Markov process, using the algorithm for a real symmetric matrix.
       Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
       where S is symmetrical, all elements of pi are positive, and U*V = I.
       space[n] is for storing sqrt(pi).

       [U 0] [Q_0 0] [U^-1 0]    [Root  0]
       [0 I] [0   0] [0    I]  = [0     0]

       Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
    */
    int i, j;

	try {
    if(mNumGoodFreq == N)
    {
		// The S matrix is defined as Q = S*pi
		// Due to the fact that each Q row should sum to zero, the S diagonal values are so adjusted.
		// But also to save multiplications the S diagonal elements are already multiplied by the corresponding codon frequency
		// Also the eigensolver use only half the matrix. So S is filled only for half.
        for(i=0; i < N; ++i)
		{
			//mU[i*N + i] = mS[i*N + i] * mCodonFreq[i];
			mU[i*N + i] = mS[i*N + i];

            for(j=0; j < i; ++j)
			{
                //mU[i*N + j] = mU[j*N + i] = mS[i*N + j] * mSqrtCodonFreq[i] * mSqrtCodonFreq[j];
                mU[i*N + j] = mS[i*N + j] * mSqrtCodonFreq[i] * mSqrtCodonFreq[j];
			}
		}

		// Eigendecomposition of mU into mD (eigenvalues) and mU (eigenvectors), size is N and mV is used as workarea
        eigenRealSymm(mU, N, mD, mV);

		// Construct mV = pi^1/2*mU
		for(j=0; j < N; ++j)
		{
			for(i=0; i < N; ++i)
			{
               mV[j*N + i] = mU[j*N + i] * mSqrtCodonFreq[j];
            }
		}
    }
    else
    {
		int inew, jnew;

        for(i=0, inew=0; i < N; ++i)
        {
            if(mGoodFreq[i])
            {
                for(j=0, jnew=0; j < i; ++j)
				{
                    if(mGoodFreq[j])
                    {
                        mU[inew*mNumGoodFreq + jnew] = mS[i*N + j] * mSqrtCodonFreq[i] * mSqrtCodonFreq[j];
                        ++jnew;
                    }
				}

                mU[inew*mNumGoodFreq + inew] = mS[i*N + i];
                ++inew;
            }
        }

		// Eigendecomposition of mU into mD (eigenvalues) and mU (eigenvectors), size is mNumGoodFreq and mV is used as workarea
		eigenRealSymm(mU, mNumGoodFreq, mD, mV);

		// Construct D (D is stored in reverse order)
        for(i=N-1, inew=mNumGoodFreq-1; i >= 0; --i)
        {
            mD[i] = mGoodFreq[N-1-i] ? mD[inew--] : 0.;
        }

		// Construct R
        for(i=N-1, inew=mNumGoodFreq-1; i >= 0; --i)
        {
            if(mGoodFreq[i])
            {
                for(j=N-1, jnew=mNumGoodFreq-1; j >= 0; --j)
                    if(mGoodFreq[j])
                    {
                        mV[j*N + i] = mU[jnew*mNumGoodFreq + inew] * mSqrtCodonFreq[j];
                        --jnew;
                    }
                    else
                    {
                        mV[j*N + i] = (i == j) ? 1. : 0.;
                    }

                --inew;
            }
            else
                for(j=0; j < N; ++j)
                {
                    mV[i*N + j] = (i == j) ? 1. : 0.;
                }
        }
	}
	}
	catch(std::exception& e)
	{
		std::cout << "Exception in eigensolver: " << e.what() << std::endl;
		throw;
	}
}
#else

void TransitionMatrix::eigenQREV(void)
{
	/*
       This finds the eigensolution of the rate matrix Q for a time-reversible
       Markov process, using the algorithm for a real symmetric matrix.
       Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
       where S is symmetrical, all elements of pi are positive, and U*V = I.
       space[n] is for storing sqrt(pi).

       [U 0] [Q_0 0] [U^-1 0]    [Root  0]
       [0 I] [0   0] [0    I]  = [0     0]

       Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
    */
    int i, j;

	try {
    if(mNumGoodFreq == N)
    {
		// Store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D)
        for(i=0; i < N; ++i)
		{
			mU[i*N + i] = mQ[i*N + i];

            for(j=0; j < i; ++j)
			{
                mU[i*N + j] = mU[j*N + i] = mQ[i*N + j] * mSqrtCodonFreq[i] / mSqrtCodonFreq[j];
			}
		}

		eigenRealSymm(mU, N, mD, mV);

        for(i=0; i < N; ++i)
		{
            for(j=0; j < N; ++j)
            {
                mV[i*N + j] = mU[j*N + i] * mSqrtCodonFreq[j];
            }
		}

        for(i=0; i < N; ++i)
		{
            for(j=0; j < N; ++j)
            {
                mU[i*N + j] /= mSqrtCodonFreq[i];
            }
		}
    }
    else
	{
		int inew, jnew;

        for(i=0, inew=0; i < N; ++i)
        {
            if(mGoodFreq[i])
            {
                for(j=0, jnew=0; j < i; ++j)
				{
                    if(mGoodFreq[j])
                    {
                        mU[inew*mNumGoodFreq + jnew] = mU[jnew*mNumGoodFreq + inew]
                                               = mQ[i*N + j] * mSqrtCodonFreq[i] / mSqrtCodonFreq[j];
                        ++jnew;
                    }
				}

                mU[inew*mNumGoodFreq + inew] = mQ[i*N + i];
                ++inew;
            }
        }

		eigenRealSymm(mU, mNumGoodFreq, mD, mV);

		// Construct D (D is stored in reverse order)
        for(i=N-1, inew=mNumGoodFreq-1; i >= 0; --i)
        {
            mD[i] = mGoodFreq[N-1-i] ? mD[inew--] : 0.;
        }

		// Construct V
        for(i=N-1, inew=mNumGoodFreq-1; i >= 0; --i)
        {
            if(mGoodFreq[i])
            {
                for(j=N-1, jnew=mNumGoodFreq-1; j >= 0; --j)
                    if(mGoodFreq[j])
                    {
                        mV[i*N + j] = mU[jnew*mNumGoodFreq + inew] * mSqrtCodonFreq[j];
                        --jnew;
                    }
                    else
                    {
                        mV[i*N + j] = (i == j) ? 1. : 0.;
                    }

                --inew;
            }
            else
                for(j=0; j < N; ++j)
                {
                    mV[i*N + j] = (i == j) ? 1. : 0.;
                }
        }

		// Construct U
        for(i=N-1, inew=mNumGoodFreq-1; i >= 0; --i)
        {
            if(mGoodFreq[i])
            {
                for(j=N-1, jnew=mNumGoodFreq-1; j >= 0; --j)
                    if(mGoodFreq[j])
                    {
                        mU[i*N + j] = mU[inew*mNumGoodFreq + jnew] / mSqrtCodonFreq[i];
                        --jnew;
                    }
                    else
                    {
                        mU[i*N + j] = (i == j) ? 1. : 0.;
                    }

                --inew;
            }
            else
                for(j=0; j < N; ++j)
                {
                    mU[i*N + j] = (i == j) ? 1. : 0.;
                }
        }
	}
	}
	catch(std::exception& e)
	{
		std::cout << "Exception in eigensolver: " << e.what() << std::endl;
		throw;
	}
}
#endif


void CheckpointableTransitionMatrix::saveCheckpoint(double aScale)
{
	mSavedScale = aScale;
	memcpy(mSavedV, mV, N*N*sizeof(double));
	memcpy(mSavedU, mU, N*N*sizeof(double));
#ifndef USE_LAPACK
	memcpy(mSavedQ, mQ, N*N*sizeof(double));
#else
	memcpy(mSavedS, mS, N*N*sizeof(double));
#endif
	memcpy(mSavedD, mD,   N*sizeof(double));
}


double CheckpointableTransitionMatrix::restoreCheckpoint(void)
{
	memcpy(mV, mSavedV, N*N*sizeof(double));
	memcpy(mU, mSavedU, N*N*sizeof(double));
#ifndef USE_LAPACK
	memcpy(mQ, mSavedQ, N*N*sizeof(double));
#else
	memcpy(mS, mSavedS, N*N*sizeof(double));
#endif
	memcpy(mD, mSavedD,   N*sizeof(double));

	return mSavedScale;
}

