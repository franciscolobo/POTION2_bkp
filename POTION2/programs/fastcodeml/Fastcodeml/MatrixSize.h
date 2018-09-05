
#ifndef MATRIXSIZE_H
#define MATRIXSIZE_H

/// The matrix side (64 possible three-letters codons less 3 STOP codons).
///
static const int N = 61;

/// Number of possible codons (64).
///
static const int N64 = 64;

/// Max number of possible tree traversals for each cycle (4).
/// This is also the number of codon classes.
///
static const int Nt = 4;

/// Slot size for a matrix. It should be equal or larger than N*N. Filler based on measurement.
///
static const int MATRIX_SLOT = N*N+3;

/// Slot size for a vector. It should be equal or larger than N+1. The slot size value has been identified by measurement.
/// In the vectors used as CPV the location after the end will contain the CVP norm.
///
#ifdef NEW_LIKELIHOOD
#ifdef USE_CPV_SCALING
#error "Cannot use CPV scaling with new likelihood!"
#endif
static const int VECTOR_SLOT = N;
#else
static const int VECTOR_SLOT = N+1+4;
#endif

/// Alignment to avoid cache line false sharing.
///
static const int CACHE_LINE_ALIGN = 64;

#endif

