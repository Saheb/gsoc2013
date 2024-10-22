/*
 * $Revision: 2963 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-11-05 18:47:50 +0530 (Mon, 05 Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Implementation of class PlanarSubgraphModule.
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


#include <ogdf/module/PlanarSubgraphModule.h>


namespace ogdf {


Module::ReturnType PlanarSubgraphModule::callAndDelete(
	GraphCopy &PG,
	const List<edge> &preferedEdges,
	List<edge> &delOrigEdges,
	bool preferedImplyPlanar)
{
	List<edge> delEdges;

	ReturnType retValue = call(PG, preferedEdges, delEdges, preferedImplyPlanar);

	if(isSolution(retValue))
	{
		ListConstIterator<edge> it;
		for(it = delEdges.begin(); it.valid(); ++it) {
			edge eCopy = *it;

			delOrigEdges.pushBack(PG.original(eCopy));
			PG.delEdge(eCopy);
		}
	}

	return retValue;
}


} // end namespace ogdf
