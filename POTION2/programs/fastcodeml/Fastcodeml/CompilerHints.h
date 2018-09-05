
#ifndef COMPILERHINTS_H
#define COMPILERHINTS_H


#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
//  Intel
#define ALIGN64 __declspec(align(64))
#define RESTRICT restrict
//#define PURE    __declspec(const)

#elif defined(__GNUC__)
//  GNU C++
#define ALIGN64 __attribute__ ((aligned (64)))
#define RESTRICT __restrict__
//#define PURE    __attribute__ ((pure))

#elif defined(_MSC_VER)
// Microsoft Visual C++
#define ALIGN64 __declspec(align(64))
#define RESTRICT __restrict
//#define PURE

#elif defined(__PGI)
//  PGI C++
#define ALIGN64 //__attribute__ ((aligned (64)))
#define RESTRICT restrict
//#define PURE    __attribute__ ((pure))

#elif defined(_CRAYC)
//  Cray C++
#define ALIGN64 __attribute__ ((aligned (64)))
#define RESTRICT restrict
//#define PURE    __attribute__ ((pure))

#elif defined(__MTA__)
// Cray XMT
#define ALIGN64
#define RESTRICT
//#define PURE

#else
#warning "Unknown compiler detected"
#define ALIGN64
#define RESTRICT
//#define PURE

#endif


#endif

