/*
 * $Revision: 2963 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-11-05 18:47:50 +0530 (Mon, 05 Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Declaration of memory manager for allocating small
 *        pieces of memory
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

#ifndef OGDF_MEMORY_H
#define OGDF_MEMORY_H


#include <stdlib.h>
#include <new>


#include <ogdf/internal/basic/PoolMemoryAllocator.h>
#include <ogdf/internal/basic/MallocMemoryAllocator.h>


namespace ogdf {

#define OGDF_MM(Alloc) \
public: \
static void *operator new(size_t nBytes) { \
	if(OGDF_LIKELY(Alloc::checkSize(nBytes))) \
		return Alloc::allocate(nBytes); \
	else \
	return ogdf::MallocMemoryAllocator::allocate(nBytes); \
} \
\
static void operator delete(void *p, size_t nBytes) { \
	if(OGDF_LIKELY(p != 0)) { \
		if(OGDF_LIKELY(Alloc::checkSize(nBytes))) \
			Alloc::deallocate(nBytes, p); \
		else \
			ogdf::MallocMemoryAllocator::deallocate(nBytes, p); \
	} \
} \
static void *operator new(size_t, void *p) { return p; } \
static void operator delete(void *, void *) { }


#define OGDF_NEW new

#ifdef OGDF_MEMORY_MALLOC_TS
#define OGDF_ALLOCATOR ogdf::MallocMemoryAllocator
#else
#define OGDF_ALLOCATOR ogdf::PoolMemoryAllocator
#endif

//! Creates new and delete operators in a class using ogdf's memory allocator.
#define OGDF_NEW_DELETE OGDF_MM(OGDF_ALLOCATOR)

//! Creates new and delete operators in a class using the malloc memory allocator.
#define OGDF_MALLOC_NEW_DELETE OGDF_MM(ogdf::MallocMemoryAllocator)

} // namespace ogdf


#endif
