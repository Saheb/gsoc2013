/*
 * $Revision: 3526 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-31 20:00:06 +0530 (Fri, 31 May 2013) $
 ***************************************************************/

/** \file
 * \brief KandinskyLayout for orthogonal planar drawings of
 *          graphs with arbitary degree using special bends
 *          at nodes.
 *
 * \author Moritz Schallab&ouml;ck
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

#ifndef OGDF_KANDINSKY_LAYOUT_H
#define OGDF_KANDINSKY_LAYOUT_H


#include <ogdf/basic/Graph_d.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/module/LayoutPlanRepUMLModule.h>
#include <ogdf/external/coin.h>


#ifdef USE_COIN

namespace ogdf {

//! The class KandinskyLayout creates an orthogonal, planar drawing of an input PlanRepUML. For nodes with a degree > 4, bends at the node
//! are generated. The orthogonal representation is computed using an Integer Program. The OrthoRep is compacted using FlowCompaction.
//! Requires the SYMPHONY IP solver.
class KandinskyLayout : public LayoutPlanRepUMLModule
{
public:
	// Constructor
	KandinskyLayout();

	//! Calls the layout algorithm. Input is a planarized representation \a PG, output is the layout \a drawing. Uses an IP to compute an orthogonal representation of \a PG.
	void call(PlanRepUML &PG, adjEntry adjExternal, Layout &drawing);

	//! Returns the minimal allowed distance between edges and vertices.
	double separation() const;
	//! Sets the minimal allowed distance between edges and vertices to \a sep.
	void separation(double sep);

protected:

	//! Creates the extension graph \a H for the input PlanRep, which has the property of having a degree &le; 4 required for the compaction step.
	void createExtension(CombinatorialEmbedding& emb, const double * sol, const int ncols,
			Graph & H, CombinatorialEmbedding& embH,  AdjEntryArray<int> &angles, AdjEntryArray<BendString> &bends,
			NodeArray<node> & nodeRefs, NodeArray<node> & nodeRefsBack, EdgeArray<edge> & edgeRefs, EdgeArray<edge> & edgeRefsBack, AdjEntryArray<double>& adjEntryOffsets, AdjEntryArray<int>& adjEntryNumNeighbours);

	//! Assigns coordinates to the nodes and edges of the extension graph using FlowCompaction. The resulting \a drawing needs to be adjusted to arrive at a drawing for the input graph of the main call.
	void draw(PlanRep &PG, OrthoRep &OR, CombinatorialEmbedding &E, adjEntry adjExternal, Layout &drawing);

	//! Generate an IP yielding the data for an orthogonal representation for \a emb.
	OsiSolverInterface * getIP(CombinatorialEmbedding &emb);

private:
	inline BendString reverseBendString(const BendString & bs);
	inline void adjustPackedVector(CoinPackedVector &vector, const int index, const double adjust);
	inline int getAngleFor(adjEntry ae, const double * sol, const int E);
	inline int getBendFor(adjEntry ae, const double * sol, const int E);
	inline BendString getBendString(adjEntry ae, const double * sol, const int E);

	bool opt1;
	bool alwaysCenterStraight;
	double sep;

};


} // end namespace ogdf
#endif // OSI_SYM

#endif	/* OGDF_KANDINSKY_LAYOUT_H */
