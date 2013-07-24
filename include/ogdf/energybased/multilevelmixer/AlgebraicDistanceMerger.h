/*
 * $Revision: 2523 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-07-03 00:29:27 +0530 (Tue, 03 Jul 2012) $
 ***************************************************************/

/** \file
 * \brief Declares class AlgebraicDistanceMerger, a merger for
 *        multilevel-layout based on the AlgebraicDistance.
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

#ifndef OGDF_ALGEBRAIC_DISTANCE_MERGER_H
#define OGDF_ALGEBRAIC_DISTANCE_MERGER_H


#include <ogdf/energybased/multilevelmixer/MultilevelBuilder.h>


namespace ogdf {


class OGDF_EXPORT AlgebraicDistanceMerger : public MultilevelBuilder
{
	double m_levelSizeFactor;

public:
	AlgebraicDistanceMerger() : m_levelSizeFactor(2.0) { }

	void setLevelSizeFactor(double f) {
		m_levelSizeFactor = f;
	}

	double levelSizeFactor() const { return m_levelSizeFactor; }

protected:
	bool buildOneLevel(MultilevelGraph &MLG);

private:
	void computeAlgDist(const Graph &G, EdgeArray<double> &dist);

	template <class T = double>
	class EdgeWeightComparerAsc {
		const EdgeArray<T> *m_pWeight;

	public:
		EdgeWeightComparerAsc(const EdgeArray<T> *pWeight) : m_pWeight(pWeight) { }

		bool less(edge e1, edge e2) const { return (*m_pWeight)[e1] < (*m_pWeight)[e2]; }
		bool operator()(edge e1, edge e2) const { return (*m_pWeight)[e1] < (*m_pWeight)[e2]; }
	};

	template <class T = double>
	class EdgeWeightComparerDesc {
		const EdgeArray<T> *m_pWeight;

	public:
		EdgeWeightComparerDesc(const EdgeArray<T> *pWeight) : m_pWeight(pWeight) { }

		bool less(edge e1, edge e2) const { return (*m_pWeight)[e1] > (*m_pWeight)[e2]; }
		bool operator()(edge e1, edge e2) const { return (*m_pWeight)[e1] > (*m_pWeight)[e2]; }
	};

};


} // end namespace ogdf


#endif
