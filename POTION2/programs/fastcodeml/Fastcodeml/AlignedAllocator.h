
#ifndef ALIGNED_ALLOCATOR_H
#define ALIGNED_ALLOCATOR_H

// The following headers are required for all allocators.
#include <cstddef>  // Required for size_t and ptrdiff_t and NULL
//#include <new>       // Required for placement new and std::bad_alloc
#include <stdexcept> // Required for std::length_error

#include "AlignedMalloc.h"

/// Aligned allocator definition.
/// It will be used to obtain a vector aligned to a given power of 2.
/// Example allocation aligned to 64: std::vector<double, AlignedAllocator<double, 64> > aligned_vector;
///
///  @author Mario Valle - Swiss National Supercomputing Centre (CSCS)
///  @date 2010-12-22 (initial version)
///  @version 1.0
///
template <typename T, size_t A> class AlignedAllocator
{
public:

    // The following will be the same for virtually all allocators.
    typedef T * pointer;
    typedef const T * const_pointer;
    typedef T& reference;
    typedef const T& const_reference;
    typedef T value_type;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    T * address(T& r) const
	{
        return &r;
    }

    const T * address(const T& s) const
	{
        return &s;
    }

    size_t max_size() const
	{
        // The following has been carefully written to be independent of
        // the definition of size_t and to avoid signed/unsigned warnings.
        return (static_cast<size_t>(0) - static_cast<size_t>(1)) / sizeof(T);
    } 

	/// Internal definition for AlignedAllocator.
    /// The following must be the same for all allocators.
	///
    template <typename U> struct rebind
	{
        typedef AlignedAllocator<U, A> other;
    };

    bool operator!=(const AlignedAllocator& other) const
	{
        return !(*this == other);
    }

    void construct(T * const p, const T& t) const
	{
        void * const pv = static_cast<void *>(p);
        new (pv) T(t);
    }

    void destroy(T * const p) const; // Defined below.

    /// Returns true if and only if storage allocated from *this
    /// can be deallocated from other, and vice versa.
    /// Always returns true for stateless allocators.
	///
	/// @param[in] other The other allocator to be compared
	///
	/// @return Always true, this is a stateless allocator.
	///
    bool operator==(const AlignedAllocator& other) const
	{
        return true;
    }

    /// Default constructor. It should be empty for stateless allocators.
	///
    AlignedAllocator() { }

    /// Default copy constructor. It should be empty for stateless allocators.
	///
    AlignedAllocator(const AlignedAllocator&) { }

    /// Default rebinding constructor. It should be empty for stateless allocators.
	///
    template <typename U> AlignedAllocator(const AlignedAllocator<U, A>&) { }

    /// Default destructor. It should be empty for stateless allocators.
	///
    ~AlignedAllocator() { }

	/// The main allocator routine.
	/// The following will be different for each allocator.
	///
	/// @param[in] n Number of objects of type T to allocate
	///
	/// @return The allocated memory
	///
	/// @exception std::length_error Integer overflow
	/// @exception std::bad_alloc Memory allocation failure
	///
    T * allocate(const size_t n) const
	{
        // AlignedAllocator prints a diagnostic message to demonstrate
        // what it's doing. Real allocators won't do this.
        //std::cout << "Allocating " << n << (n == 1 ? " object" : " objects")
        //    << " of size " << sizeof(T) << " aligned on " << A << std::endl;

        // The return value of allocate(0) is unspecified.
        // AlignedAllocator returns NULL in order to avoid depending
        // on malloc(0)'s implementation-defined behavior
        // (the implementation can define malloc(0) to return NULL,
        // in which case the bad_alloc check below would fire).
        // All allocators can return NULL in this case.
        if (n == 0) return NULL;

        // All allocators should contain an integer overflow check.
        // The Standardization Committee recommends that std::length_error
        // be thrown in the case of integer overflow.
        if (n > max_size()) throw std::length_error("AlignedAllocator<T>::allocate() - Integer overflow.");

        // AlignedAllocator wraps aligned malloc().
        void * const pv = alignedMalloc(n * sizeof(T), A);

        // Allocators should throw std::bad_alloc in the case of memory allocation failure.
        if(pv == NULL) throw std::bad_alloc();

        return static_cast<T *>(pv);
    }

    void deallocate(T * const p, const size_t /*n*/) const
	{
        // AlignedAllocator prints a diagnostic message to demonstrate
        // what it's doing. Real allocators won't do this.
        //std::cout << "Deallocating " << n << (n == 1 ? " object" : " objects")
        //    << " of size " << sizeof(T) << "." << std::endl;

        // AlignedAllocator wraps aligned free().
        alignedFree(p);
    }

    /// The following will be the same for all allocators that ignore hints.
	///
    template <typename U> T * allocate(const size_t n, const U * /* const hint */) const
	{
        return allocate(n);
    }


private:
    /// Allocators are not required to be assignable, so
    /// all allocators should have a private unimplemented
    /// assignment operator. Note that this will trigger the
    /// off-by-default (enabled under /Wall) warning C4626
    /// "assignment operator could not be generated because a
    /// base class assignment operator is inaccessible" within
    /// the STL headers, but that warning is useless.
    //AlignedAllocator& operator=(const AlignedAllocator& a) {return this;}
    //AlignedAllocator& operator=(const AlignedAllocator& a) {return const_cast<AlignedAllocator&>(a);}
    AlignedAllocator& operator=(const AlignedAllocator&);
};


// A compiler bug causes it to believe that p->~T() doesn't reference p.

#ifdef _MSC_VER
    #pragma warning(push)
    #pragma warning(disable: 4100) // unreferenced formal parameter
#endif

/// The definition of destroy() must be the same for all allocators.
template <typename T, size_t A> inline void AlignedAllocator<T, A>::destroy(T * const p) const
{
    p->~T();
}

#ifdef _MSC_VER
    #pragma warning(pop)
#endif

#endif

