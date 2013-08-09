/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Definition of a Triple used in contraction-based approximation
 *        algorithm for the minimum steiner tree problem
 *
 * \author Matthias Woste
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

#ifndef TRIPLELC_OGDF_H_
#define TRIPLELC_OGDF_H_

#include <ogdf/basic/Graph.h>

namespace ogdf {

/*!
 * \brief This class represents a triple used by various contraction-based minimum steiner tree approximations.
 */
template<typename T>
class Triple {
public:
	Triple() {
	}
	;
	~Triple() {
	}
	;

	node s0() const {
		return m_s0;
	}
	node s1() const {
		return m_s1;
	}
	node s2() const {
		return m_s2;
	}
	node z() const {
		return m_z;
	}
	T cost() const {
		return m_cost;
	}

	void s0(node u) {
		m_s0 = u;
	}
	void s1(node u) {
		m_s1 = u;
	}
	void s2(node u) {
		m_s2 = u;
	}
	void z(node u) {
		m_z = u;
	}
	void cost(T c) {
		m_cost = c;
	}

private:
	node m_s0, m_s1, m_s2; //!< terminal nodes
	node m_z; //!< center node of the triple
	T m_cost; //!< edge costs of the triple in the original graph
};

} // end namespace ogdf

#endif /* TRIPLELC_OGDF_H_ */
