
#ifndef TYPES_H
#define TYPES_H

#include <vector>
#include "AlignedAllocator.h"
#include "MatrixSize.h"

/// Array of doubles aligned on a cache line
///
typedef std::vector<double, AlignedAllocator<double, CACHE_LINE_ALIGN> > CacheAlignedDoubleVector;

/// Array of doubles to be used by SSE instructions
///
typedef std::vector<double, AlignedAllocator<double, 16> > SSEAlignedDoubleVector;

#endif

