/*
 * $Revision: 3368 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-04 23:37:31 +0530 (Thu, 04 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Handling of clique replacement in planarization layout.
 *
 * \author Carsten Gutwenger, Karsten Klein
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


#include <ogdf/internal/planarity/CliqueReplacer.h>
#include <ogdf/basic/Math.h>
#include <ogdf/misclayout/CircularLayout.h>


namespace ogdf {

	CliqueReplacer::CliqueReplacer(GraphAttributes &ga, Graph &G) : m_G(G), m_ga(ga)
	{
	}


	// replace each node set in cliques by a star connecting
	// a new center node with all nodes in set, deletes all
	// edges between nodes in set, lists need to be disjoint
	// TODO: think about directly using the cliquenum array
	// output of findcliques here
	void CliqueReplacer::replaceByStar(List< List<node> > &cliques)
	{
		m_cliqueCircleSize.init(m_G);
		m_cliqueCirclePos.init(m_G);
		m_replacementEdge.init(m_G, false);

		if (cliques.empty()) return;
		//we save membership of nodes in each list
		NodeArray<int> cliqueNum(m_G, -1);
		ListIterator< List<node> > it = cliques.begin();

		int num = 0;
		while (it.valid())
		{
			ListIterator<node> itNode = (*it).begin();
			while (itNode.valid())
			{
				cliqueNum[(*itNode)] = num;
				itNode++;
			}//while

			num++;
			it++;
		}//while

		//now replace each list
		it = cliques.begin();
		while (it.valid())
		{
			node newCenter = replaceByStar((*it), cliqueNum);
			OGDF_ASSERT(newCenter)
				m_centerNodes.pushBack(newCenter);
			//now we compute a circular drawing of the replacement
			//and save its size and the node positions
			m_cliqueCircleSize[newCenter] = circularBound(newCenter);
			it++;
		}//while

	}//replacebystar


	node CliqueReplacer::replaceByStar(List<node> &clique, NodeArray<int> &cliqueNum)
	{
		if (clique.empty()) return 0;
		//insert an additional center node

		node center = m_G.newNode();
		m_ga.width(center) = m_cliqueCenterSize;
		m_ga.height(center) = m_cliqueCenterSize;
#ifdef OGDF_DEBUG
		//should ask for attributes
		if(m_ga.attributes() & GraphAttributes::nodeStyle)
			m_ga.fillColor(center) = Color(0x55,0x55,0x55);
#endif
		//we delete all edges inzident to two clique nodes
		//store all of them first in delEdges
		List<edge> delEdges;
		ListIterator<node> it = clique.begin();
		while (it.valid())
		{
			node v = (*it);

			adjEntry ad;
			int numIt = cliqueNum[v];

			forall_adj(ad, v)
			{
				if (cliqueNum[ad->twinNode()] == numIt)
				{
					if (ad->theEdge()->source() == v)
					{
						//m_cliqueEdges[v].pushBack(new CliqueInfo(ad->theEdge()->target(), ad->theEdge()->index()));
						delEdges.pushBack(ad->theEdge());
					}
				}//if
			}//foralladj

			//connect center node to clique node
			edge inserted = m_G.newEdge(center, v);
			m_replacementEdge[inserted] = true;

			it++;
		}//while

		//now delete all edges
		ListIterator<edge>	itEdge = delEdges.begin();
		while (itEdge.valid())
		{
			//m_pG->delEdge((*itEdge));
			m_G.hideEdge((*itEdge));
			itEdge++;
		}//while

		return center;
	}//replaceByStar


	//compute a drawing of the clique around node center and save its size
	//the call to circular will later be replaced by an dedicated computation
	DRect CliqueReplacer::circularBound(node center)
	{
		//TODO: hier computecliqueposition(0,...) benutzen, rest weglassen
		DRect bb;
		CircularLayout cl;
		Graph G;
		GraphAttributes AG(G);
		NodeArray<node> umlOriginal(G);

		//TODO: we need to assure that the circular drawing
		//parameters fit the drawing parameters of the whole graph
		//umlgraph clique parameter members?

		OGDF_ASSERT(center->degree() > 0)
			node lastNode = 0;
		node firstNode = 0;
		node v;

		adjEntry ae = center->firstAdj();
		do {
			node w = ae->twinNode();
			v = G.newNode();
			umlOriginal[v] = w;

			if (!firstNode) firstNode = v;
			AG.width(v) = m_ga.width(w);
			AG.height(v) = m_ga.height(w);
			ae = ae->cyclicSucc();
			if (lastNode != 0) G.newEdge(lastNode, v);
			lastNode = v;
		} while (ae != center->firstAdj());
		G.newEdge(lastNode, firstNode);

		cl.call(AG);

		forall_nodes(v, G)
		{
			m_cliqueCirclePos[umlOriginal[v]] = DPoint(AG.x(v), AG.y(v));
		}//forallnodes
		bb = AG.boundingBox();

		return bb;
	}//circularBound


	void CliqueReplacer::undoStars()
	{
		SListIterator<node> it = m_centerNodes.begin();
		while (it.valid())
		{
			undoStar(*it, false);
			it++;
		}//while

		m_G.restoreAllEdges();
		m_centerNodes.clear();
		m_replacementEdge.init();

	}//undostars


	//remove the center node and reinsert the deleted edges
	void CliqueReplacer::undoStar(node center, bool restoreAllEdges)
	{
		OGDF_ASSERT(center)

			//TODO: we should only restore the hidden clique edges, maybe there were
			//already hidden edges and we call this for all cliques, but it is global
			if (restoreAllEdges) m_G.restoreAllEdges();

		//remove center node
		m_G.delNode(center);

	}//undostar


	//computes relative positions of all nodes around center on a circle
	//keeping the topological ordering of the nodes, letting the opposite
	//node (to center) of the edge at firstAdj having the position at
	//three o'clock (TODO: unter umstaenden hier aus gegebenen coords in
	//this den Punkt auswaehlen. Aber: Was, wenn keine coords=>optional?
	void CliqueReplacer::computeCliquePosition(node center, double rectMin)
	{
		List<node> adjNodes;
		adjEntry ae = center->firstAdj();
		do
		{
			adjNodes.pushBack(ae->twinNode());
			ae = ae->cyclicPred();
		} while (ae != center->firstAdj());
		computeCliquePosition(adjNodes, center, rectMin);
	}//computeCliquePosition


	// computes relative positions of all nodes in List cList on a minimum size
	// circle (needed to compute positions with different ordering than given in *this).
	// Precondition: nodes in adjNodes are adjacent to center
	// first node in adjNodes is positioned to the right
	void CliqueReplacer::computeCliquePosition(List<node> &adjNodes, node center, double rectMin)//, const adjEntry &startAdj)
	{
		DRect boundingBox;
		OGDF_ASSERT(center->degree() > 0)
			OGDF_ASSERT(center->degree() == adjNodes.size())

			node v;
		double radius = 0.0;
		//TODO: member, parameter
		//const
		double minDist = 1.0;
		//TODO: necessary?
		double minCCDist = 20.0;

		ListIterator<node> itNode = adjNodes.begin();

		//--------------------------------------------------------------------------
		//for the temporary solution (scale clique to fixed rect if possible instead
		//of guaranteeing the rect size in compaction) we check in advance if the sum
		//of diameters plus dists fits into the given rect by heuristic estimate (biggest
		//node size + radius)
		if (rectMin > 0.0)
		{
			double rectDist = m_cliqueCenterSize; //dist to rect border todo: parameter
			double rectBound = rectMin - 2.0*rectDist;
			double maxSize = 0.0;
			double pureSumDiameters = 0.0;
			while (itNode.valid())
			{
				node q  =(*itNode);
				double d = sqrt(
					m_ga.width(q)*m_ga.width(q) + m_ga.height(q)*m_ga.height(q));
				pureSumDiameters += d;

				if (d > maxSize) maxSize = d;

				itNode++;
			}//while
			double totalSum = pureSumDiameters+(center->degree()-1)*minDist;
			//TODO: scling, not just counting
			while (totalSum/Math::pi < rectBound*0.75)
			{
				minDist = minDist + 1.0;
				totalSum += (center->degree()-1.0);
			}//while
			if (minDist > 1.1) minDist -= 1.0;
			//do not use larger value than cliquecentersize (used with separation)
			//if (minDist > m_cliqueCenterSize) minDist = m_cliqueCenterSize;
			itNode = adjNodes.begin();
		}
		//temporary part ends-------------------------------------------------------
		//------------------------------------------
		//first, we compute the radius of the circle

		const int n = center->degree();
		//sum of all diameters around the nodes and the max diameter radius
		double sumDiameters = 0.0, maxR = 0;
		//list of angles for all nodes
		List<double> angles; //node at startAdj gets 0.0
		double lastDiameter = 0.0; //temporary storage of space needed for previous node
		bool first = true;

		while (itNode.valid())
		{
			v  =(*itNode);
			double d = sqrt(
				m_ga.width(v)*m_ga.width(v) + m_ga.height(v)*m_ga.height(v));

			sumDiameters += d;

			if (d/2.0 > maxR) maxR = d/2.0;

			//save current position relative to startadj
			//later on, compute angle out of these values
			if (first)
			{
				angles.pushBack(0.0);
				first = false;
			}
			else
			{
				angles.pushBack(lastDiameter+d/2.0+minDist+angles.back());
			}
			lastDiameter = d/2.0; //its only half diameter...

			itNode++;
		}//while

		OGDF_ASSERT(adjNodes.size() == angles.size())

			if(n == 1) {
				radius      = 0;

			} else if (n == 2) {
				radius      = 0.5*minDist + sumDiameters / 4;

			} else {
				double perimeter = (n*minDist + sumDiameters);
				radius      = perimeter / (2*Math::pi);

				ListIterator<double> it = angles.begin();
				itNode = adjNodes.begin();
				while (it.valid())
				{
					(*it) = (*it)*360.0/perimeter;
					node w = *itNode;
					double angle = Math::pi*(*it)/180.0;
					m_cliqueCirclePos[w].m_x = radius*cos(angle);
					m_cliqueCirclePos[w].m_y = radius*sin(angle);
					itNode++;
					it++;
				}//while
			}//if n>2

			//now we normalize the values (start with 0.0) and
			//derive the bounding box
			v = adjNodes.front();
			double minX = m_cliqueCirclePos[v].m_x,
				maxX = m_cliqueCirclePos[v].m_x;
			double minY = m_cliqueCirclePos[v].m_y,
				maxY = m_cliqueCirclePos[v].m_y;
			itNode = adjNodes.begin();
			while (itNode.valid())
			{
				node w = *itNode;
				double wx = m_cliqueCirclePos[w].m_x;
				double wy = m_cliqueCirclePos[w].m_y;
				if(wx-m_ga.width (w)/2.0 < minX) minX = wx-m_ga.width(w)/2.0;
				if(wx+m_ga.width (w)/2.0 > maxX) maxX = wx+m_ga.width(w)/2.0;
				if(wy-m_ga.height(w)/2.0 < minY) minY = wy-m_ga.height(w)/2.0;
				if(wy+m_ga.height(w)/2.0 > maxY) maxY = wy+m_ga.height(w)/2.0;
				itNode++;
			}
			//allow distance
			minX -= minCCDist;
			minY -= minCCDist;
			//normalize
			//cout<<"\n";

			itNode = adjNodes.begin();
			while (itNode.valid())
			{
				node w = *itNode;
				//cout<<"x1:"<<m_cliqueCirclePos[w].m_x<<":y:"<<m_cliqueCirclePos[w].m_y<<"\n";
				m_cliqueCirclePos[w].m_x -= minX;
				m_cliqueCirclePos[w].m_y -= minY;
				//cout<<"x:"<<m_cliqueCirclePos[w].m_x<<":y:"<<m_cliqueCirclePos[w].m_y<<"\n";
				itNode++;
			}

			//reassign the size, this time it is the final value
			m_cliqueCircleSize[center] = DRect(0.0, 0.0, maxX-minX, maxY-minY);
	}//computecliqueposition

}
