/*
 * $Revision: 3395 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-15 17:58:20 +0530 (Mon, 15 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration and implementation of Array class and
 * Array algorithms
 *
 * \author Carsten Gutwenger
 *
 * \par License:
 * This file is part of the Open Graph Drawing Framework (OGDF).
 *
 * \par
 * Copyright (C)<br>
 * See README.txt in the root directory of the OGDF installation for details.
 *
 * \par
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * Version 2 or 3 as published by the Free Software Foundation;
 * see the file LICENSE.txt included in the packaging of this file
 * for details.
 *
 * \par
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * \par
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/


#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_ARRAY_H
#define OGDF_ARRAY_H


#include <ogdf/basic/basic.h>


namespace ogdf {

template<class E, class INDEX> class ArrayBuffer;

//! Iteration over all indices \a i of an array \a A.
/**
 * Note that the index variable \a i has to be defined prior to this macro
 * (just as for \c #forall_edges, etc.).
 * <h3>Example</h3>
 *
 *   \code
 *   Array<double> A;
 *   ...
 *   int i;
 *   forall_arrayindices(i, A) {
 *     cout << A[i] << endl;
 *   }
 *   \endcode
 *
 *   Note that this code is equivalent to the following tedious long version
 *
 *   \code
 *   Array<double> A;
 *   ...
 *   int i;
 *   for(i = A.low(); i <= A.high(); ++i) {
 *     cout << A[i] << endl;
 *   }
 *   \endcode
 */
#define forall_arrayindices(i, A) \
	for(i = (A).low(); i<=(A).high(); ++i)

//! Iteration over all indices \a i of an array \a A, in reverse order.
/**
 * Note that the index variable \a i has to be defined prior to this macro
 * (just as for \c #forall_edges, etc.).
 * See \c #forall_arrayindices for an example
 */
#define forall_rev_arrayindices(i, A) \
	for(i = (A).high(); i>=(A).low(); --i)



//! The parameterized class \a Array<E,INDEX> implements dynamic arrays of type \a E.
/**
 * @tparam E     denotes the element type.
 * @tparam INDEX denotes the index type. The index type must be chosen such that it can
 *               express the whole index range of the array instance, as well as its size.
 *               The default index type is \c int, other possible types are \c short and
 *               <code>long long</code> (on 64-bit systems).
 */
template<class E, class INDEX = int> class Array {
public:
	//! Threshold used by \a quicksort() such that insertion sort is
	//! called for instances smaller than \a maxSizeInsertionSort.
	enum { maxSizeInsertionSort = 40 };


	//! Creates an array with empty index set.
	Array() { construct(0,-1); }

	//! Creates an array with index set [0..\a s-1].
	explicit Array(INDEX s) {
		construct(0,s-1); initialize();
	}

	//! Creates an array with index set [\a a..\a b].
	Array(INDEX a, INDEX b) {
		construct(a,b); initialize();
	}

	//! Creates an array with index set [\a a..\a b] and initializes each element with \a x.
	Array(INDEX a, INDEX b, const E &x) {
		construct(a,b); initialize(x);
	}

	//! Creates an array that is a copy of \a A.
	Array(const Array<E,INDEX> &A) {
		copy(A);
	}

	//! Creates an array that is a copy of \a A. The array-size is set to be the number of elements (not the capacity) of the buffer.
	Array(const ArrayBuffer<E,INDEX> &A);

	// destruction
	~Array() {
		deconstruct();
	}

	//! Returns the minimal array index.
	INDEX low() const { return m_low; }

	//! Returns the maximal array index.
	INDEX high() const { return m_high; }

	//! Returns the size (number of elements) of the array.
	INDEX size() const { return m_high - m_low + 1; }

	//! Returns a pointer to the first element.
	E *begin() { return m_pStart; }

	//! Returns a pointer to the first element.
	const E *begin() const { return m_pStart; }

	//! Returns a pointer to one past the last element.
	E *end() { return m_pStop; }

	//! Returns a pointer to one past the last element.
	const E *end() const { return m_pStop; }

	//! Returns a pointer to the last element.
	E *rbegin() { return m_pStop-1; }

	//! Returns a pointer to the last element.
	const E *rbegin() const { return m_pStop-1; }

	//! Returns a pointer to one before the first element.
	E *rend() { return m_pStart-1; }

	//! Returns a pointer to one before the first element.
	const E *rend() const { return m_pStart-1; }

	//! Returns a reference to the element at position \a i.
	const E &operator[](INDEX i) const {
		OGDF_ASSERT(m_low <= i && i <= m_high)
		return m_vpStart[i];
	}

	//! Returns a reference to the element at position \a i.
	E &operator[](INDEX i) {
		OGDF_ASSERT(m_low <= i && i <= m_high)
		return m_vpStart[i];
	}

	//! Swaps the elements at position \a i and \a j.
	void swap(INDEX i, INDEX j) {
		OGDF_ASSERT(m_low <= i && i <= m_high)
		OGDF_ASSERT(m_low <= j && j <= m_high)

		std::swap(m_vpStart[i], m_vpStart[j]);
	}

	//! Reinitializes the array to an array with empty index set.
	void init() {
		//init(0,-1);
		deconstruct();
		construct(0,-1);
	}

	//! Reinitializes the array to an array with index set [0..\a s-1].
	/**
	 * Notice that the elements contained in the array get discarded!
	 */
	void init(INDEX s) { init(0,s-1); }

	//! Reinitializes the array to an array with index set [\a a..\a b].
	/**
	 * Notice that the elements contained in the array get discarded!
	 */
	void init(INDEX a, INDEX b) {
		deconstruct();
		construct(a,b);
		initialize();
	}

	//! Reinitializes the array to an array with index set [\a a..\a b] and sets all entries to \a x.
	void init(INDEX a, INDEX b, const E &x) {
		deconstruct();
		construct(a,b);
		initialize(x);
	}

	//! Assignment operator.
	Array<E,INDEX> &operator=(const Array<E,INDEX> &array2) {
		deconstruct();
		copy(array2);
		return *this;
	}

	//! Sets all elements to \a x.
	void fill(const E &x) {
		E *pDest = m_pStop;
		while(pDest > m_pStart)
			*--pDest = x;
	}

	//! Sets elements in the intervall [\a i..\a j] to \a x.
	void fill(INDEX i, INDEX j, const E &x) {
		OGDF_ASSERT(m_low <= i && i <= m_high)
		OGDF_ASSERT(m_low <= j && j <= m_high)

		E *pI = m_vpStart + i, *pJ = m_vpStart + j+1;
		while(pJ > pI)
			*--pJ = x;
	}

	//! Enlarges the array by \a add elements and sets new elements to \a x.
	/**
	 *  Note: address of array entries in memory may change!
	 * @param add is the number of additional elements; \a add can be negative in order to shrink the array.
	 * @param x is the inital value of all new elements.
	 */
	void grow(INDEX add, const E &x);

	//! Enlarges the array by \a add elements.
	/**
	 *  Note: address of array entries in memory may change!
	 * @param add is the number of additional elements; \a add can be negative in order to shrink the array.
	 */
	void grow(INDEX add);

	//! Resizes (enlarges or shrinks) the array to hold \a newSize elements and sets new elements to \a x.
	/**
	 *  Note: address of array entries in memory may change!
	 * @param newSize is new size of the array
	 * @param x is the inital value of all new elements.
	 */
	void resize(INDEX newSize, const E &x) { grow(newSize - size(), x); }

	//! Resizes (enlarges or shrinks) the array to hold \a newSize elements.
	/**
	 *  Note: address of array entries in memory may change!
	 * @param newSize is new size of the array
	 */
	void resize(INDEX newSize) { grow(newSize - size()); }

	//! Randomly permutes the subarray with index set [\a l..\a r].
	void permute(INDEX l, INDEX r);

	//! Randomly permutes the array.
	void permute() {
		permute(low(), high());
	}

#ifdef OGDF_HAVE_CPP11
	//! Randomly permutes the subarray with index set [\a l..\a r].
	template<class RNG>
	void permute(INDEX l, INDEX r, RNG &rng);

	//! Randomly permutes the array.
	template<class RNG>
	void permute(RNG &rng) {
		permute(low(), high(), rng);
	}
#endif

	//! Performs a binary search for element \a x.
	/**
	 * \pre The array must be sorted!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	inline int binarySearch (const E& e) const {
		return binarySearch(low(), high(), e, StdComparer<E>());
	}

	//! Performs a binary search for element \a x within the array section [l..r] .
	/**
	 * \pre The array must be sorted!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	inline int binarySearch (INDEX l, INDEX r, const E& e) const {
		return binarySearch(l, r, e, StdComparer<E>());
	}

	//! Performs a binary search for element \a x with comparer \a comp.
	/**
	 * \pre The array must be sorted according to \a comp!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	template<class COMPARER>
	inline int binarySearch(const E& e, const COMPARER &comp) const {
		return binarySearch(low(), high(), e, comp);
	}

	//! Performs a binary search for element \a x within the array section [l..r] with comparer \a comp.
	/**
	 * \pre The array must be sorted according to \a comp!
	 * \return the index of the found element, and low()-1 if not found.
	 */
	template<class COMPARER>
	INDEX binarySearch(INDEX l, INDEX r, const E& e, const COMPARER &comp) const {
		if(r<l) return low()-1;
		while(r>l) {
			INDEX m = (r + l)/2;
			if(comp.greater(e, m_vpStart[m]))
				l = m+1;
			else
				r = m;
		}
		return comp.equal(e, m_vpStart[l]) ? l : low()-1;
	}
	//! Performs a linear search for element \a x.
	/**
	 * Warning: This method has linear running time!
	 * Note that the linear search runs from back to front.
	 * \return the index of the found element, and low()-1 if not found.
	 */
	inline INDEX linearSearch (const E& e) const {
		int i;
		for(i = size(); i-->0;)
			if(e == m_pStart[i]) break;
		return i+low();	}

	//! Performs a linear search for element \a x with comparer \a comp.
	/**
	 * Warning: This method has linear running time!
	 * Note that the linear search runs from back to front.
	 * \return the index of the found element, and low()-1 if not found.
	 */
	template<class COMPARER>
	INDEX linearSearch(const E& e, const COMPARER &comp) const {
		int i;
		for(i = size(); i-->0;)
			if(comp.equal(e, m_pStart[i])) break;
		return i+low();
	}

	//! Sorts array using Quicksort.
	inline void quicksort() {
		quicksort(StdComparer<E>());
	}

	//! Sorts subarray with index set [\a l..\a r] using Quicksort.
	inline void quicksort(INDEX l, INDEX r) {
		quicksort(l, r, StdComparer<E>());
	}

	//! Sorts array using Quicksort and a user-defined comparer \a comp.
	/**
	 * @param comp is a user-defined comparer; \a C must be a class providing a \a less(x,y) method.
	 */
	template<class COMPARER>
	inline void quicksort(const COMPARER &comp) {
		if(low() < high())
			quicksortInt(m_pStart,m_pStop-1,comp);
	}

	//! Sorts the subarray with index set [\a l..\a r] using Quicksort and a user-defined comparer \a comp.
	/**
	 * @param l is the left-most position in the range to be sorted.
	 * @param r is the right-most position in the range to be sorted.
	 * @param comp is a user-defined comparer; \a C must be a class providing a \a less(x,y) method.
	 */
	template<class COMPARER>
	void quicksort(INDEX l, INDEX r, const COMPARER &comp) {
		OGDF_ASSERT(low() <= l && l <= high())
		OGDF_ASSERT(low() <= r && r <= high())
		if(l < r)
			quicksortInt(m_vpStart+l,m_vpStart+r,comp);
	}

	//! Removes the components listed in \a ind by shifting the remaining components to the left.
	/**
	 * The "free" positions in the array at the end remain as they are.
	 *
	 * This function is mainly used by Abacus. Other uses are supposed to be rare.
	 *
	 * Memory management of the removed components must be
	 * carefully implemented by the user of this function to avoid
	 * memory leaks.
	 *
	 * @param ind The compenents being removed from the array.
	 */
	void leftShift(ArrayBuffer<INDEX, INDEX> &ind);

	//! Removes the components listed in \a ind by shifting the remaining components to the left.
	/**
	 * The "free" positions in the array at the end are filled with \a val.
	 *
	 * Memory management of the removed components must be
	 * carefully implemented by the user of this function to avoid
	 * memory leaks.
	 *
	 * @param ind specifies the components that are removed from the array.
	 * @param val is the value used to fill the positions that become free.
	 */
	void leftShift(ArrayBuffer<INDEX, INDEX> &ind, const E& val) {
		leftShift(ind);
		fill(high()-ind.size(),high(),val);
	}

	template<class F, class I> friend class ArrayBuffer; // for efficient ArrayBuffer::compact-method

private:
	E *m_vpStart; //!< The virtual start of the array (address of A[0]).
	E *m_pStart;  //!< The real start of the array (address of A[m_low]).
	E *m_pStop;   //!< Successor of last element (address of A[m_high+1]).
	INDEX m_low;    //!< The lowest index.
	INDEX m_high;   //!< The highest index.

	//! Allocates new array with index set [\a a..\a b].
	void construct(INDEX a, INDEX b);

	//! Initializes elements with default constructor.
	void initialize();

	//! Initializes elements with \a x.
	void initialize(const E &x);

	//! Deallocates array.
	void deconstruct();

	//! Constructs a new array which is a copy of \a A.
	void copy(const Array<E,INDEX> &A);

	//! Internal Quicksort implementation with comparer template.
	template<class COMPARER>
	static void quicksortInt(E *pL, E *pR, const COMPARER &comp) {
		size_t s = pR-pL;

		// use insertion sort for small instances
		if (s < maxSizeInsertionSort) {
			for (E *pI = pL+1; pI <= pR; pI++) {
				E v = *pI;
				E *pJ = pI;
				while (--pJ >= pL && comp.less(v,*pJ)) {
					*(pJ+1) = *pJ;
				}
				*(pJ+1) = v;
			}
			return;
		}

		E *pI = pL, *pJ = pR;
		E x = *(pL+(s>>1));

		do {
			while (comp.less(*pI,x)) pI++;
			while (comp.less(x,*pJ)) pJ--;
			if (pI <= pJ) std::swap(*pI++,*pJ--);
		} while (pI <= pJ);

		if (pL < pJ) quicksortInt(pL,pJ,comp);
		if (pI < pR) quicksortInt(pI,pR,comp);
	}

	OGDF_NEW_DELETE
}; // class Array



// enlarges array by add elements and sets new elements to x
template<class E, class INDEX>
void Array<E,INDEX>::grow(INDEX add, const E &x)
{
	if(add==0) return;

	INDEX sOld = size(), sNew = sOld + add;

	// expand allocated memory block
	if(m_pStart != 0) {
		E *p = (E *)realloc(m_pStart, sNew*sizeof(E));
		if(p == 0) OGDF_THROW(InsufficientMemoryException);
		m_pStart = p;
	} else {
		m_pStart = (E *)malloc(sNew*sizeof(E));
		if (m_pStart == 0) OGDF_THROW(InsufficientMemoryException);
	}

	m_vpStart = m_pStart-m_low;
	m_pStop   = m_pStart+sNew;
	m_high   += add;

	// initialize new array entries
	for (E *pDest = m_pStart+sOld; pDest < m_pStop; pDest++)
		new (pDest) E(x);
}

// enlarges array by add elements (initialized with default constructor)
template<class E, class INDEX>
void Array<E,INDEX>::grow(INDEX add)
{
	if(add==0) return;

	INDEX sOld = size(), sNew = sOld + add;

	// expand allocated memory block
	if(m_pStart != 0) {
		E *p = (E *)realloc(m_pStart, sNew*sizeof(E));
		if(p == 0) OGDF_THROW(InsufficientMemoryException);
		m_pStart = p;
	} else {
		m_pStart = (E *)malloc(sNew*sizeof(E));
		if (m_pStart == 0) OGDF_THROW(InsufficientMemoryException);
	}

	m_vpStart = m_pStart-m_low;
	m_pStop   = m_pStart+sNew;
	m_high   += add;

	// initialize new array entries
	for (E *pDest = m_pStart+sOld; pDest < m_pStop; pDest++)
		new (pDest) E;
}

template<class E, class INDEX>
void Array<E,INDEX>::construct(INDEX a, INDEX b)
{
	m_low = a; m_high = b;
	INDEX s = b-a+1;

	if (s < 1) {
		m_pStart = m_vpStart = m_pStop = 0;

	} else {
		m_pStart = (E *)malloc(s*sizeof(E));
		if (m_pStart == 0) OGDF_THROW(InsufficientMemoryException);

		m_vpStart = m_pStart - a;
		m_pStop = m_pStart + s;
	}
}


template<class E, class INDEX>
void Array<E,INDEX>::initialize()
{
	E *pDest = m_pStart;
	try {
		for (; pDest < m_pStop; pDest++)
			new(pDest) E;
	} catch (...) {
		while(--pDest >= m_pStart)
			pDest->~E();
		free(m_pStart);
		throw;
	}
}


template<class E, class INDEX>
void Array<E,INDEX>::initialize(const E &x)
{
	E *pDest = m_pStart;
	try {
		for (; pDest < m_pStop; pDest++)
			new(pDest) E(x);
	} catch (...) {
		while(--pDest >= m_pStart)
			pDest->~E();
		free(m_pStart);
		throw;
	}
}


template<class E, class INDEX>
void Array<E,INDEX>::deconstruct()
{
	if (doDestruction((E*)0)) {
		for (E *pDest = m_pStart; pDest < m_pStop; pDest++)
			pDest->~E();
	}
	free(m_pStart);
}


template<class E, class INDEX>
void Array<E,INDEX>::copy(const Array<E,INDEX> &array2)
{
	construct(array2.m_low, array2.m_high);

	if (m_pStart != 0) {
		E *pSrc = array2.m_pStop;
		E *pDest = m_pStop;
		while(pDest > m_pStart)
			//*--pDest = *--pSrc;
			new (--pDest) E(*--pSrc);
	}
}


// permutes array a from a[l] to a[r] randomly
template<class E, class INDEX>
void Array<E,INDEX>::permute (INDEX l, INDEX r)
{
	OGDF_ASSERT(low() <= l && l <= high())
	OGDF_ASSERT(low() <= r && r <= high())

	E *pI = m_vpStart+l, *pStart = m_vpStart+l, *pStop = m_vpStart+r;
	while(pI <= pStop)
		std::swap(*pI++,*(pStart+randomNumber(0,r-l)));
}


#ifdef OGDF_HAVE_CPP11

// permutes array a from a[l] to a[r] randomly
template<class E, class INDEX>
template<class RNG>
void Array<E,INDEX>::permute (INDEX l, INDEX r, RNG &rng)
{
	OGDF_ASSERT(low() <= l && l <= high())
	OGDF_ASSERT(low() <= r && r <= high())

	std::uniform_int_distribution<int> dist(0,r-l);

	E *pI = m_vpStart+l, *pStart = m_vpStart+l, *pStop = m_vpStart+r;
	while(pI <= pStop)
		std::swap( *pI++, *(pStart + dist(rng)) );
}

#endif


// prints array a to output stream os using delimiter delim
template<class E, class INDEX>
void print(ostream &os, const Array<E,INDEX> &a, char delim = ' ')
{
	for (int i = a.low(); i <= a.high(); i++) {
		if (i > a.low()) os << delim;
		os << a[i];
	}
}


// output operator
template<class E, class INDEX>
ostream &operator<<(ostream &os, const ogdf::Array<E,INDEX> &a)
{
	print(os,a);
	return os;
}

}

#include <ogdf/basic/ArrayBuffer.h>

namespace ogdf {

template<class E, class INDEX>
void Array<E,INDEX>::leftShift(ArrayBuffer<INDEX, INDEX> &ind) {
	const INDEX nInd = ind.size();
	if (nInd == 0) return;

	//! shift all items up to the last element of \a ind to the left
#ifdef OGDF_DEBUG
	if(ind[0] < low() || ind[0] > high())
		OGDF_THROW_PARAM(AlgorithmFailureException, afcIndexOutOfBounds);
#endif

	INDEX j, current = ind[0];
	for (INDEX i = 0; i < nInd - 1; i++) {
#ifdef OGDF_DEBUG
		if(ind[i+1] < low() || ind[i+1] > high())
			OGDF_THROW_PARAM(AlgorithmFailureException, afcIndexOutOfBounds);
#endif

		const INDEX last = ind[i+1];
		for(j = ind[i]+1; j < last; j++)
			m_vpStart[current++] = m_vpStart[j];
	}

	//! copy the rest of the buffer
	for (j = ind[nInd - 1]+1; j < size(); j++)
		m_vpStart[current++] = m_vpStart[j];
}

template<class E, class INDEX>
Array<E,INDEX>::Array(const ArrayBuffer<E, INDEX> &A) {
	construct(0,-1);
	A.compactCopy(*this);
}

} // end namespace ogdf


#endif
