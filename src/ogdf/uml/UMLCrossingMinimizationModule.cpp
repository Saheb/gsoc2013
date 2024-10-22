/*
 * $Revision: 3135 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2012-12-11 15:41:27 +0530 (Tue, 11 Dec 2012) $
 ***************************************************************/

/** \file
 * \brief implementation of class EdgeInsertionModuleOld
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


#include <ogdf/module/UMLCrossingMinimizationModule.h>


namespace ogdf {


bool UMLCrossingMinimizationModule::checkCrossingGens(const PlanRepUML &prUML)
{
	edge e;
	forall_edges(e,prUML) {
		Graph::EdgeType et = prUML.typeOf(e);
		if (et != Graph::generalization && et != Graph::association)
			return false;
	}

	node v;
	forall_nodes(v,prUML)
	{
		if (prUML.typeOf(v) == PlanRepUML::dummy && v->degree() == 4) {
			adjEntry adj = v->firstAdj();

			edge e1 = adj->theEdge();
			edge e2 = adj->succ()->theEdge();

			if (prUML.typeOf(e1) == Graph::generalization &&
				prUML.typeOf(e2) == Graph::generalization)
				return false;
		}
	}

	return true;
}


} // end namespace ogdf

