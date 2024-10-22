/*
 * $Revision: 3008 $
 *
 * last checkin:
 *   $Author: chimani $
 *   $Date: 2012-11-12 21:04:59 +0530 (Mon, 12 Nov 2012) $
 ***************************************************************/

/** \file
 * \brief Declaration of the variable class for the Branch&Cut algorithm
 * for the Maximum C-Planar SubGraph problem
 *
 * \author Mathias Jansen
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

#ifndef OGDF_CPLANAR_EDGE_H
#define OGDF_CPLANAR_EDGE_H

#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/internal/cluster/Cluster_EdgeVar.h>

#include <ogdf/abacus/variable.h>

namespace ogdf {


class CPlanarEdgeVar : public EdgeVar {
	friend class MaxCPlanarSub;
public:
	CPlanarEdgeVar(abacus::Master *master, double obj, node source, node target) :
		EdgeVar(master, obj, source, target)
	{

	}
	CPlanarEdgeVar(abacus::Master *master, double obj, double lbound, node source, node target) :
			EdgeVar(master, obj, lbound, source, target)
	{

	}

	virtual ~CPlanarEdgeVar() {}

	void printMe(ostream& out) {
		out << "[Var: " << sourceNode() << "->" << targetNode() << " (" << "connect" << ") ZF=" << obj() << "]";
	}

};

}

#endif

