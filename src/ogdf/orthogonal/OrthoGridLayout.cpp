/*
 * $Revision: 3418 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-18 14:36:44 +0530 (Thu, 18 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implements planar orthogonal drawing algorithm for
 * mixed-upward embedded graphs
 *
 * \author Carsten Gutwenger, Sebastian Leipert
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


#include <ogdf/orthogonal/OrthoGridLayout.h>
#include <ogdf/orthogonal/LongestPathCompaction.h>
#include <ogdf/orthogonal/GridFlowCompaction.h>
#include <ogdf/orthogonal/EdgeRouter.h>
#include <ogdf/internal/orthogonal/RoutingChannel.h>
#include <ogdf/orthogonal/MinimumEdgeDistances.h>
#include <ogdf/orthogonal/OrthoShaper.h>
#include <ogdf/planarity/SimpleEmbedder.h>

namespace ogdf {


OrthoGridLayout::OrthoGridLayout()
{
	//drawing object distances
	separation(1.0);
	m_margin     = 0;
	m_optionProfile = 0;
	//scale layout while improving it
	//m_useScalingCompaction = false;
	//m_scalingSteps = 0;
	m_bendBound = 2; //bounds number of bends per edge in ortho shaper
	m_embedder.set(new SimpleEmbedder);
	m_conCompactor.set(new GridFlowCompaction);
	m_impCompactor.set(new GridFlowCompaction);

	m_orthoStyle = 0;//0; //traditional 0, progressive 1
}

void OrthoGridLayout::doCall(
		PlanRep &PG,
		adjEntry adjExternal,
		GridLayout &drawing,
		IPoint &boundingBox,
		bool fixEmbedding)
{

	// if we have only one vertex in PG ...
	if(PG.numberOfNodes() == 1) {
		node v1 = PG.firstNode();
		//node vOrig = PG.original(v1);

		drawing.x(v1) = m_margin;
		drawing.y(v1) = m_margin;
		m_gridBoundingBox = IPoint(2*m_margin, 2*m_margin);
		return;
	}



	//compaction with scaling: help node cages to pass by each other
	int l_orsep = (int)separation();
	/*if (m_useScalingCompaction)
	{
		m_scalingSteps = 6;
		int scaleFactor = int(1 << m_scalingSteps);
		separation(scaleFactor*separation()); //reduce this step by step in compaction
	}*///if scaling


	//***********************************
	// PHASE 1: determine orthogonal shape

	// expand high-degree vertices and generalization mergers
	//PG.expand();

	//check preconditions, currently not necessary
	//assureDrawability(PG);


	// embed augmented graph (if required)
	if(fixEmbedding) {
		OGDF_ASSERT(PG.representsCombEmbedding());
		if (adjExternal == 0)
			adjExternal = PG.firstEdge()->adjSource();
	}
	else
		m_embedder.get().call(PG,adjExternal);
	// get combinatorial embedding
	CombinatorialEmbedding E(PG);
	E.setExternalFace(E.rightFace(adjExternal));

	// determine orthogonal shape
	OrthoRep OR;

	//OrthoFormerUML OF;
	OrthoShaper OFG;

	//set some options
	OFG.traditional(m_orthoStyle > 0 ? false : true); //prefer 90/270 degree angles over 180/180

	// New Call
	OFG.setBendBound(m_bendBound);
	OFG.call(PG,E,OR,true); //only 4-planar!

	//******************************************************************
	// PHASE 2: construction of a feasible drawing of the expanded graph

	// expand low degree vertices
	// We only treat grid graphs here, i.e. we don't need to expand to
	// consider vertex sizes
	//PG.expandLowDegreeVertices(OR);

	OGDF_ASSERT(PG.representsCombEmbedding());

	// restore embedding
	E.computeFaces();
	E.setExternalFace(E.rightFace(adjExternal));

	// apply constructive compaction heuristics

	OR.normalize();
	OR.gridDissect(&PG); //OR.dissect();

	OR.orientate(PG,odNorth);

	// compute cage information and routing channels
	OR.computeCageInfoUML(PG);

	//double dummyov = 0.0; //overhang
	bool dummyal = false; //align
	//double dummycost = 1.0; //edge cost

	//temporary grid layout
	//GridLayoutMapped gridDrawing(PG,OR,m_separation,dummyov,2);

	//RoutingChannel<int> rcGrid(PG,gridDrawing.toGrid(m_separation),dummyov);
	//rcGrid.computeRoutingChannels(OR, dummyal);


	node v;
	const OrthoRep::VertexInfoUML *pInfoExp;
	forall_nodes(v,PG) {
		pInfoExp = OR.cageInfo(v);

		if (pInfoExp) break;
	}

	m_conCompactor.get().callConstructive(PG,OR,drawing);

	OR.undissect(dummyal);

	// call flow compaction on grid
	GridFlowCompaction fc;

	//fc.scalingSteps(m_scalingSteps);

	m_impCompactor.get().callImprovement(PG,OR,drawing);

	// PHASE 3: routing of edges
	//

	//EdgeRouter router;
	//MinimumEdgeDistances<int> minDistGrid(PG, gridDrawing.toGrid(m_separation));
	//router.setOrSep(int(gridDrawing.toGrid(l_orsep))); //scaling test
	//router.call(PG,OR,gridDrawing,E,rcGrid,minDistGrid, gridDrawing.width(),
	//	gridDrawing.height(), dummyal);


	string msg;
	OGDF_ASSERT(OR.check(msg) == true);

	OR.orientate(PG.firstEdge()->adjSource(),odNorth);

	//*************************************************
	// PHASE 4: apply improvement compaction heuristics

	// call flow compaction on grid
	//fc.improvementHeuristics(PG,OR,drawing,l_orsep);
							 //int(gridDrawing.toGrid(l_orsep)));


	// re-map result
	//gridDrawing.remapGrid(drawing);

	// collapse all expanded vertices by introducing a new node in the center
	// of each cage representing the original vertex
	//PG.collapseVertices(OR,drawing);

	// finally set the bounding box
	computeBoundingBox(PG,drawing);

	separation(l_orsep);
}//call



//-----------------------------------------------------------------------------
//Helpers
//-----------------------------------------------------------------------------
void OrthoGridLayout::classifyEdges(PlanRepUML &PG, adjEntry &adjExternal)
{
	//classify brother-to-brother hierarchy edges to allow alignment
	//when shifting this to planrep, guarantee edge type correction in planarization
	//save external face entry

	//PG.classifyEdges
	//potential direct connection are all non-gen. edges that are alignUpward
	edge e, eSucc;
	for(e = PG.firstEdge(); e; e = eSucc)
	{
		eSucc = e->succ();
		if (PG.typeOf(e) != Graph::generalization)
		{
			adjEntry as = e->adjSource();
			node v = e->source();
			if ( (PG.alignUpward(as))
				&& (PG.typeOf(e->target()) != Graph::dummy)//TODO: crossings ?
				&& (PG.typeOf(v) != Graph::dummy)
			   )
			{
				edge gen1, gen2;
				int stop = 0;
				adjEntry runAE = as->cyclicSucc();
				edge run = runAE->theEdge();
				while ( (stop < v->degree()) &&  //only once
						((PG.typeOf(run) != Graph::generalization) || //search generalization
						(run->source() != v) //outgoing gen
						)
						)
				{
				  stop++;
				  runAE = runAE->cyclicSucc();
				  run = runAE->theEdge();
				}//while
				OGDF_ASSERT(stop <= v->degree());

				//now we have the outgoing generalization (to the merger) at v
				gen1 = run;

				node w = e->target(); //crossings ?
				adjEntry asTwin = as->twin();

				stop = 0;
				runAE = asTwin->cyclicSucc();
				run = runAE->theEdge();
				while ( (stop < w->degree()) &&
						((PG.typeOf(run) != Graph::generalization) ||
						(run->source() != w)
						)
						)
				{
				  stop++;
				  runAE = runAE->cyclicSucc();
				  run = runAE->theEdge();
				}//while
				OGDF_ASSERT(stop <= w->degree());

				//now we have the outgoing generalization (to the merger) at w
				gen2 = run;

				//two possible orientations
				//left to right
				bool ltr = ( gen1->adjSource()->faceCycleSucc() == gen2->adjTarget() );
				//right to left
				bool rtl = ( gen2->adjSource()->faceCycleSucc() == gen1->adjTarget() );
				if (ltr || rtl) //should be disjoint cases because of merger node
				{
				  PG.setBrother(e);

				  //now check if the embedding does include unnecessary nodes in the slope
				  if (ltr)
				  {
					  //there are edges between e and gen2 at target
					  if (!(e->adjTarget()->faceCyclePred() == gen2->adjTarget()))
					  {
						OGDF_ASSERT(v != e->target());
						PG.moveAdj(e->adjTarget(), before, gen2->adjTarget()->twin());
					  }
					  //there are edges between e and gen1 at source
					  if (!(e->adjTarget()->faceCycleSucc() == gen1->adjSource()))
					  {
						//test if we discard the outer face entry
						if (adjExternal == e->adjSource())
						{
							adjExternal = e->adjSource()->faceCyclePred();
						}
						PG.moveAdj(e->adjSource(), after, gen1->adjSource());
					  }
				  }//if gen 1 left of gen 2
				  if (rtl)
				  {
					  //there are edges between e and gen2 at target
					  if (!(e->adjSource()->faceCycleSucc() == gen2->adjSource()))
					  {
						//test if we discard the outer face entry
						if (adjExternal == e->adjTarget())
						{
							adjExternal = e->adjTarget()->faceCycleSucc();
						}
						PG.moveAdj(e->adjTarget(), after, gen2->adjSource());
					  }
					  //there are edges between e and gen1 at source
					  if (!(e->adjSource()->faceCyclePred() == gen1->adjTarget()))
					  {
						PG.moveAdj(e->adjSource(), before, gen1->adjSource());
					  }

				  }//if gen 2 left of gen 1
				}//if
				else PG.setHalfBrother(e);

			}//if upward edge
		}//if not generalization
	}//for


}//classifyedges



// compute bounding box and move final drawing such that it is 0 aligned
// respecting margins
void OrthoGridLayout::computeBoundingBox(
	const PlanRep &PG,
	GridLayout &drawing)
{
	int minX, maxX, minY, maxY;

	minX = maxX = drawing.x(PG.firstNode());
	minY = maxY = drawing.y(PG.firstNode());

	node v;
	forall_nodes(v,PG)
	{
		int x = drawing.x(v);
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;

		int y = drawing.y(v);
		if (y < minY) minY = y;
		if (y > maxY) maxY = y;
	}

	int deltaX = m_margin - minX;
	int deltaY = m_margin - minY;

	forall_nodes(v,PG)
	{
		drawing.x(v) += deltaX;
		drawing.y(v) += deltaY;
	}

	m_gridBoundingBox = IPoint(maxX+deltaX+m_margin, maxY+deltaY+m_margin);
}


} // end namespace ogdf

