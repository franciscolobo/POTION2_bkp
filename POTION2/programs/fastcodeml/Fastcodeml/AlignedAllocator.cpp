// The following headers are required for all allocators.
#include <cstddef>  // Required for size_t and ptrdiff_t and NULL
#include <stdexcept> // Required for std::length_error

// The following headers contain stuff that AlignedAllocator uses.
#include <cstdlib>  // For malloc() and free()

// For XMT the following function should be supplied till the bug is fixed by Cray
#ifdef __MTA__
//extern "C" int posix_memalign(void **memptr, size_t alignment, size_t size);
static int posix_memalign(void **memptr, size_t alignment, size_t size)
{
	*memptr = malloc(size);
	return 0;
}
#endif

#include "AlignedMalloc.h"

#ifndef _MSC_VER
//#include <stdint.h>  // for uintptr_t
#include <malloc.h>
#endif

// Alignment must be power of 2 (1,2,4,8,16...)
void* alignedMalloc(size_t aSize, size_t aAlignment)
{
#ifdef _MSC_VER
#if 0
	--aAlignment;
    uintptr_t r = reinterpret_cast<uintptr_t>(malloc(aSize + aAlignment + sizeof(uintptr_t)));
    if(!r) return NULL;
    uintptr_t t = r + sizeof(uintptr_t);
    uintptr_t o = (t + aAlignment) & ~static_cast<uintptr_t>(aAlignment);
    reinterpret_cast<uintptr_t*>(o)[-1] = r;
    return reinterpret_cast<void*>(o);
#endif
	return _aligned_malloc(aSize, aAlignment);
#else
	void* ptr = NULL;
	if(posix_memalign(&ptr, aAlignment, aSize)) return NULL;
	return ptr;
#endif
}

void alignedFree(void* aPtr)
{
    if(!aPtr) return;
#ifdef _MSC_VER
#if 0
    free(reinterpret_cast<void*>(reinterpret_cast<uintptr_t*>(aPtr)[-1]));
#endif
	_aligned_free(aPtr);
#else
    free(aPtr);
#endif
}


#if 0

// The following headers contain stuff that main() uses.
#include <iostream>  // For std::cout
#include <ostream>   // For std::endl
#include <vector>    // For std::vector
#include "AlignedAllocator.h"

int main()
{
    using namespace std;

    cout << "Constructing l:" << endl;

    vector<double, AlignedAllocator<double, 8> > l;
	l.reserve(10);
    cout << endl << "l.push_back(1729):" << endl;

    l.push_back(1729.);

    cout << endl << "l.push_back(2161):" << endl;

    l.push_back(2161.);

    cout << endl;
	double* p = &l[0];
	int x = reinterpret_cast<int>(p);
	cout << "Aligned on 16: " << x%16 << endl;
	cout << "Aligned on 8:  " << x%8 << endl;
	cout << "Aligned on 4:  " << x%4 << endl;
    cout << endl;

    for (vector<double, AlignedAllocator<double, 8> >::const_iterator i = l.begin(); i != l.end(); ++i) {
        cout << "Element: " << *i << endl;
    }

    cout << endl << "Destroying l:" << endl;
}
#endif

