/*
 * $Revision: 2571 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-07-10 20:55:20 +0530 (Tue, 10 Jul 2012) $
 ***************************************************************/

/** \file
 * \brief Declares class FixedUpwardEmbeddingInserter which inserts
 * edges into an upward planar graph
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


#ifndef OGDF_FIXED_UPWARD_EMBEDDING_INSERTER_H
#define OGDF_FIXED_UPWARD_EMBEDDING_INSERTER_H


#include <ogdf/basic/GraphCopy.h>


namespace ogdf {

//! inserts edges into an upward planar graph
class OGDF_EXPORT FixedUpwardEmbeddingInserter
{
public:
	//! constructor
	FixedUpwardEmbeddingInserter() { }
	//! destructor
	~FixedUpwardEmbeddingInserter() { }


	void call(GraphCopy &GC, const List<edge> &origEdges);


private:
	class Network;

}; // class FixedUpwardEmbeddingInserter



} // end namespace ogdf


#endif
