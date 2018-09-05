
#ifndef ALIGNEDMALLOC_H
#define ALIGNEDMALLOC_H

/// Allocate and return a block of memory aligned as requested.
/// Alignment must be power of 2 (1,2,4,8,16...)
///
/// @param[in] aSize The size (in bytes) of the block to be allocated
/// @param[in] aAlignment The alignment requested (must be a power of 2)
///
/// @return Pointer to the allocated area
///
extern void* alignedMalloc(size_t aSize, size_t aAlignment);

/// Free the memory allocated with alignedMalloc.
///
/// @param[in] aPtr The pointer to the memory (allocated with alignedMalloc) to be freed
///
extern void alignedFree(void* aPtr);

#endif
