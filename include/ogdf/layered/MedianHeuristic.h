/*
 * $Revision: 3210 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-01-15 16:28:53 +0530 (Tue, 15 Jan 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class MedianHeuristic
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

#ifndef OGDF_MEDIAN_HEURISTIC_H
#define OGDF_MEDIAN_HEURISTIC_H



#include <ogdf/module/TwoLayerCrossMin.h>


namespace ogdf {


//! The median heuristic for 2-layer crossing minimization.
class OGDF_EXPORT MedianHeuristic : public TwoLayerCrossMin
{
public:
	//! Creates a new instance of the median heuristic.
	MedianHeuristic() { }

	//! Creates a new instance of the median heuristic.
	MedianHeuristic(const MedianHeuristic &crossMin) { }

	//! Returns a new instance of the median heuristic with the same option settings.
	TwoLayerCrossMin *clone() const { return new MedianHeuristic; }

	//! Initializes crossing minimization for hierarchy \a H.
	void init (const HierarchyLevels &levels) { m_weight.init(levels.hierarchy()); }

	//! Calls the median heuristic for level \a L.
	void call (Level &L);

	//! Does some clean-up after calls.
	void cleanup () { m_weight.init(); }

private:
	NodeArray<int> m_weight; //!< The median weight of a node.
};


} // end namespace ogdf


#endif
