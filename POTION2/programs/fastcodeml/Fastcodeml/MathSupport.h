
#ifndef MATHSUPPORT_H
#define MATHSUPPORT_H

#ifdef _CRAYC
#undef __SSE2__
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif
#include "MatrixSize.h"
#include "CompilerHints.h"

#ifdef USE_LAPACK
#include "blas.h"
#else
#include <cmath>
#endif

//#ifdef USE_MKL_VML
//#include <mkl_vml_functions.h>
//#endif

/// Dot product specialized for 61 elements vectors.
///
/// @param[in] aV1 First vector
/// @param[in] aV2 Second vector
///
/// @return The dot product
///
inline double dot(const double* RESTRICT aV1, const double* RESTRICT aV2)
{
#if 0
	double result;
   __m128d num1, num2, num3, num4;

    num4 = _mm_setzero_pd();

    for(int i=0; i < N-1; i += 2)
    {
        num1 = _mm_load_pd(aV1+i);
        num2 = _mm_load_pd(aV2+i);
        num3 = _mm_mul_pd(num1, num2);
        num4 = _mm_add_pd(num4, num3);
    }
    num4 = _mm_hadd_pd(num4, num4);
    _mm_store_sd(&result, num4);
    result += aV1[60]*aV2[60];
	return result;
#endif
#ifdef USE_LAPACK
	return ddot_(&N, aV1, &I1, aV2, &I1);
#else
	double tot = 0.;
	for(int i=0; i < N; ++i) tot += aV1[i]*aV2[i];
	return tot;
#endif
}

/// Element-wise vector-vector multiplication (specialized to 61 elements vectors)
///
/// @param[in,out] aVres Vector that should be multiplied by the aV one
/// @param[in] aV Multiplicand (that is: for(i=0; i < N; ++i) aVres[i] *= aV[i])
///
inline void elementWiseMult(double* RESTRICT aVres, const double* RESTRICT aV)
{
//#ifdef USE_MKL_VML
//	vdMul(N, aVres, aV, aVres);
//#elif defined(__SSE2__)
//	__m128d num1, num2, num3;
//
//    for(int i=0; i < N-1; )
//    {
//        num1 = _mm_load_pd(aVres+i);
//        num2 = _mm_load_pd(aV+i);
//        num3 = _mm_mul_pd(num1, num2);
//        _mm_store_pd(aVres+i, num3);
//		i += 2;
//
//        num1 = _mm_load_pd(aVres+i);
//        num2 = _mm_load_pd(aV+i);
//        num3 = _mm_mul_pd(num1, num2);
//        _mm_store_pd(aVres+i, num3);
//		i += 2;
//    }
//	aVres[N-1] *= aV[N-1];
//#else
	// Manual unrolling gives the best results here
	for(int i=0; i < 61; ++i) aVres[i] *= aV[i];
#if 0
	for(int i=0; i < 60; )
	{
		aVres[i] *= aV[i]; ++i;
		aVres[i] *= aV[i]; ++i;
		aVres[i] *= aV[i]; ++i;
		aVres[i] *= aV[i]; ++i;
		aVres[i] *= aV[i]; ++i;
		aVres[i] *= aV[i]; ++i;
	}
	aVres[60] *= aV[60];
#endif
//#endif
}


/// Check if two values are sufficiently different
///
/// @param[in] aFirst First number to compare
/// @param[in] aSecond Second term to compare
///
/// @return True if the two parameters differs more than (hardcoded) TOL
///
inline bool isDifferent(double aFirst, double aSecond)
{
	static const double TOL = 1e-8;
	const double diff = aFirst - aSecond;
	return (diff > TOL || diff < -TOL);
}

#ifdef USE_CPV_SCALING

/// Normalize a vector (length 61).
///
/// @param[in,out] aVector The vector to be scaled
///
/// @return The length of the vector
///
inline double normalizeVector(double* RESTRICT aVector)
{
#ifdef USE_LAPACK
	double norm = dnrm2_(&N, aVector, &I1);
	double inv_norm = 1./norm;

	dscal_(&N, &inv_norm, aVector, &I1);
#else
	double norm = 0.;
	for(int i=0; i < N; ++i) norm += aVector[i]*aVector[i];
	norm = sqrt(norm);
	for(int i=0; i < N; ++i) aVector[i] /= norm;
#endif
	return norm;
}

#endif


#endif

