/*
 * $Revision: 3526 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-05-31 20:00:06 +0530 (Fri, 31 May 2013) $
 ***************************************************************/

/** \file
 * \brief Implements class AlgebraicDistanceMerger, a merger for
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

#include <ogdf/energybased/multilevelmixer/AlgebraicDistanceMerger.h>
#include <ogdf/basic/SList.h>


namespace ogdf {


int loss(edge e, NodeArray<bool> tmpMark)
{
	node v = e->source();
	node w = e->target();
	int counter = 1;

	adjEntry adj;
	forall_adj(adj,v)
		tmpMark[adj->twinNode()] = true;

	forall_adj(adj,w) {
		node x = adj->twinNode();
		if(tmpMark[x])
			counter++;
	}

	forall_adj(adj,v)
		tmpMark[adj->twinNode()] = false;

	return counter;
}

bool AlgebraicDistanceMerger::buildOneLevel(MultilevelGraph &MLG)
{
	Graph &G = MLG.getGraph();
	int level = MLG.getLevel() + 1;

	// special case: only three nodes
	int n = G.numberOfNodes();
	if(n <= 3)
		return false;

	// compute algebraic distances
	EdgeArray<double> weight(G);
	computeAlgDist(G,weight);

	double minDist, maxDist;
	minDist = maxDist = weight[G.firstEdge()];
	edge e;
	forall_edges(e,G) {
		double w = weight[e];
		if(w < minDist) minDist = w;
		if(w > maxDist) maxDist = w;
	}

	int minDeg, maxDeg;
	minDeg = maxDeg = G.firstNode()->degree();
	node v;
	forall_nodes(v,G) {
		int d = v->degree();
		if(d < minDeg) minDeg = d;
		if(d > maxDeg) maxDeg = d;
	}

	NodeArray<bool> tmpMark(G,false);

	forall_edges(e,G) {
		node v = e->source();
		node w = e->target();

		double nwe  = (weight[e] - minDist) / (maxDist-minDist);
		//double ndv = double(v->degree() - minDeg) / (maxDeg-minDeg);
		//double ndw = double(w->degree() - minDeg) / (maxDeg-minDeg);

		//weight[e] = 1.0 / ( nwe + 0.5*ndv + 0.5*ndw );
		//weight[e] = 1.0 / ((ndv+ndw) * nwe);
		weight[e] = 1.0 / (nwe * MLG.radius(v) * MLG.radius(w) );
	}

	// create array of all edges sorted by increasing algebraic distance
	Array<edge> edges(G.numberOfEdges());
	int m = 0;
	forall_edges(e,G)
		edges[m++] = e;

	EdgeWeightComparerDesc<> cmp(&weight);
	std::sort(&edges[0], &edges[0] + edges.size(), cmp);

	// compute matching
	NodeArray<bool> marked(G,false);
	SListPure<edge> matching, rest;

	for(int i = 0; i < m; ++i) {
		e = edges[i];

		node src = e->source(), tgt = e->target();
		if(!marked[src] && !marked[tgt]) {
			matching.pushBack(e);
			marked[src] = marked[tgt] = true;
		} else
			rest.pushBack(e);
	}

	// compute edge cover
	SListPure<edge> edgeCover;
	while(!rest.empty()) {
		e = rest.popFrontRet();
		node src = e->source(), tgt = e->target();
		if(!marked[src] || !marked[tgt]) {
			edgeCover.pushBack(e);
			marked[src] = marked[tgt] = true;
		}
	}

	// merge nodes
	bool retVal = false;
	NodeArray<node> parent(G,0);

	while((!matching.empty() || !edgeCover.empty()) && G.numberOfNodes() > n / m_levelSizeFactor)
	{
		e = (!matching.empty()) ? matching.popFrontRet() : edgeCover.popFrontRet();

		node mergeNode  = e->source();
		node parentNode = e->target();
		if(mergeNode->degree() > parentNode->degree())
			swap(mergeNode,parentNode);

		while(parent[parentNode] != 0)
			parentNode = parent[parentNode];

		while(parent[mergeNode] != 0)
			mergeNode = parent[mergeNode];

		if(MLG.getNode(parentNode->index()) == parentNode && MLG.getNode(mergeNode->index()) == mergeNode
			&& parentNode != mergeNode)
		{
			NodeMerge *nm = new NodeMerge(level);
			bool ret = MLG.changeNode(nm, parentNode, MLG.radius(parentNode), mergeNode);
			OGDF_ASSERT(ret);
			MLG.moveEdgesToParent(nm, mergeNode, parentNode, true, m_adjustEdgeLengths);

			if(MLG.postMerge(nm, mergeNode)) {
				parent[mergeNode] = parentNode;
				retVal = true;
			} else {
				delete nm;
				retVal = false;
			}
		}
	}

	return retVal;
}


void AlgebraicDistanceMerger::computeAlgDist(const Graph &G, EdgeArray<double> &dist)
{
	edge e;
	forall_edges(e,G)
		dist[e] = 0.00001;

	NodeArray<double> prev(G), next(G);

	for(int R = 0; R < 3; ++R)
	{
		node v;
		forall_nodes(v,G)
			prev[v] = randomDouble(-0.5,0.5);

		const double w = 0.5;

		for(int k = 0; k < 7; ++k) {
			forall_nodes(v,G) {
				next[v] = 0;

				adjEntry adj;
				forall_adj(adj,v)
					next[v] += prev[adj->twinNode()];

				if(v->degree() > 1)
					next[v] /= v->degree();
			}

			forall_nodes(v,G)
				prev[v] = (1-w)*prev[v] + w*next[v];
		}

		edge e;
		forall_edges(e,G)
			dist[e] += fabs(prev[e->source()] - prev[e->target()]) / 7.0;
	}
}


} // end namespace ogdf

