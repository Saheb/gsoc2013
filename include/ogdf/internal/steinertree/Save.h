/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Interface for various LCA methods
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

#ifndef SAVE_OGDF_H_
#define SAVE_OGDF_H_

#include <ogdf/internal/steinertree/Triple.h>

namespace ogdf {

template<typename T>

/*!
 * \brief This class serves as an interface for different approaches concerning
 * the calculation of save edges
 */
class Save {
public:
	Save() {
	}
	virtual ~Save() {
	}

	/*!
	 * \brief Returns the gain (sum of the save edges) of a node triple
	 * @param u First triple node
	 * @param v Second triple node
	 * @param w Third triple node
	 * @return Sum of the save edges between the three nodes
	 */
	virtual T gain(node u, node v, node w) const = 0;

	/*!
	 * \brief Returns the weight of the save edge between two nodes
	 * @param u First node
	 * @param v Second node
	 * @return The weight of the save edge
	 */
	virtual T saveWeight(node u, node v) const = 0;

	/*!
	 * \brief Returns the save edge between two nodes
	 * @param u First node
	 * @param v Second node
	 * @return The save edge
	 */
	virtual edge saveEdge(node u, node v) const = 0;

	/*!
	 * \brief Updates the weighted tree data structure given a contracted triple
	 * @param t The contracted triple
	 */
	virtual void update(const Triple<T> &t) = 0;

	/*!
	 * \brief Checks, if components of a given triple have already been contracted
	 * @param t Triple to be checked. This is only necessary for the dynamic tree approach.
	 * @return True, if parts of the triple have been contracted before, false otherwise
	 */
	virtual bool alreadyContracted(const Triple<T> &t) const {
		return false;
	}
	;

};
} // end namespace ogdf

#endif /* SAVE_OGDF_H_ */
