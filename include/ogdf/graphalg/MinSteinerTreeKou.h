/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration and implementation of the class computing a
 * 		  2(1-1/l) minimum steiner tree approximation according
 * 		  to the algorithm of Kou et al.
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

#ifndef OGDF_MIN_STEINER_TREE_KOU_H_
#define OGDF_MIN_STEINER_TREE_KOU_H_

#include <ogdf/basic/List.h>
#include <ogdf/module/MinSteinerTreeModule.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>

namespace ogdf {

/*!
 * \brief This class implements the Minimum Steiner Tree 2-approximation algorithm by Kou et al.
 *
 * This implementation is based on:
 *
 * (L. Kou, G. Markowsky, L. Berman, A fast algorithm for steiner trees, Acta Informatica,
 * volumne 15, number 2, pages 141-145, 1981)
 */
template<typename T>
class MinSteinerTreeKou: public MinSteinerTreeModule<T> {
public:
	MinSteinerTreeKou() {
	}
	;
	virtual ~MinSteinerTreeKou() {
	}
	;
	/*!
	 * \brief Builds a minimum steiner tree given a weighted graph and a list of terminals
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal,
			EdgeWeightedGraphCopy<T> *&finalSteinerTree);

protected:
	/*!
	 * \brief Builds a complete terminal graph
	 * @param wG the original graph
	 * @param terminals list of terminals
	 * @param predecessor stores for each edge in the complete terminal graph the according path in the original graph
	 * @param completeTerminalGraph the resulting complete terminal graph
	 */
	void calculateCompleteGraph(const EdgeWeightedGraph<T> &wG, const List<node> &terminals, EdgeArray<List<edge> > &predecessor,
			EdgeWeightedGraphCopy<T> &completeTerminalGraph);

	/*!
	 * \brief Swaps an edge in the complete terminal graph with the corresponding shortest path in the original graph
	 * @param completeTerminalGraph the complete terminal graph
	 * @param ssspPred contains for each edge in the complete terminal graph the corresponding path in the original graph
	 * @param mstPred predecessor data structure of a minimum terminal spanning tree
	 * @param finalSteinerTree the resulting steiner tree
	 * @param wG the original graph
	 */
	void reinsertShortestPaths(const EdgeWeightedGraphCopy<T> &completeTerminalGraph, const EdgeArray<List<edge> > &ssspPred,
			const NodeArray<edge> &mstPred, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG);

	/*!
	 * \brief Inserts a shortest path corresponding to an edge in the complete terminal graph
	 * @param ssspPred contains for each edge in the complete terminal graph the corresponding path in the original graph
	 * @param finalSteinerTree the resulting steiner tree
	 * @param wG the original graph
	 */
	void insertPath(const List<edge> &ssspPred, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG);
};

} // end namespace ogdf

// ============= Implementation =================

#include <ogdf/graphalg/Dijkstra.h>

namespace ogdf {

template<typename T>
T MinSteinerTreeKou<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	EdgeWeightedGraphCopy<T> completeTerminalGraph;
	completeTerminalGraph.createEmpty(G);

	edge e;

	forall_listiterators(node, it, terminals) {
		completeTerminalGraph.newNode(*it);
	}

	List<edge> steinerTreeEdges;
	EdgeArray<List<edge> > ssspPred(completeTerminalGraph);

	calculateCompleteGraph(G, terminals, ssspPred, completeTerminalGraph);

	NodeArray<edge> mstPred(completeTerminalGraph);
	EdgeArray<bool> isInTree(completeTerminalGraph);
	T mstWeight = computeMinST(completeTerminalGraph, completeTerminalGraph.edgeWeights(), mstPred, isInTree);

	finalSteinerTree = new EdgeWeightedGraphCopy<T>();
	finalSteinerTree->createEmpty(G);

	reinsertShortestPaths(completeTerminalGraph, ssspPred, mstPred, *finalSteinerTree, G);

	EdgeArray<bool> isInStTree(*finalSteinerTree);
	mstWeight = computeMinST(*finalSteinerTree, finalSteinerTree->edgeWeights(), isInStTree);

	edge nextEdge;
	for (e = finalSteinerTree->firstEdge(); e; e = nextEdge) {
		nextEdge = e->succ();
		if (!isInStTree[e]) {
			finalSteinerTree->delEdge(e);
		}
	}

	mstWeight -= MinSteinerTreeModule<T>::eliminate1DegreeSteinerNodes(*finalSteinerTree, isTerminal);
	OGDF_ASSERT(MinSteinerTreeModule<T>::isSteinerTree(G, terminals, *finalSteinerTree));

	return mstWeight;
}

template<typename T>
void MinSteinerTreeKou<T>::calculateCompleteGraph(const EdgeWeightedGraph<T> &wG, const List<node> &terminals, EdgeArray<List<edge> > &predecessor, EdgeWeightedGraphCopy<T> &completeTerminalGraph)
{
	Dijkstra<T> sssp;
	for (node u = completeTerminalGraph.firstNode(); u->succ(); u = u->succ()) {
		NodeArray<T> d(wG);
		NodeArray<edge> pi(wG);
		sssp.call(wG, wG.edgeWeights(), completeTerminalGraph.original(u), pi, d);
		for (node v = u->succ(); v; v = v->succ()) {
			edge e = completeTerminalGraph.newEdge(u, v, d[completeTerminalGraph.original(v)]);
			predecessor[e].clear();
			for (node t = completeTerminalGraph.original(v); pi[t]; t = pi[t]->opposite(t)) {
				predecessor[e].pushBack(pi[t]);
			}
		}
	}
}

template<typename T>
void MinSteinerTreeKou<T>::reinsertShortestPaths(const EdgeWeightedGraphCopy<T> &completeTerminalGraph, const EdgeArray<List<edge> > &ssspPred, const NodeArray<edge> &mstPred, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG)
{
	node u;
	forall_nodes(u, completeTerminalGraph) {
		if (mstPred[u]) {
			insertPath(ssspPred[mstPred[u]], finalSteinerTree, wG);
		}
	}
}

template<typename T>
void MinSteinerTreeKou<T>::insertPath(const List<edge> &ssspPred, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG)
{
	node edgeSource, edgeTarget, stSource, stTarget;
	edge newE;

	forall_listiterators(edge, it, ssspPred) {
		if ((*it) != 0 && finalSteinerTree.chain(*it).size() == 0) {
			edgeSource = (*it)->source();
			edgeTarget = (*it)->target();

			if ((stSource = finalSteinerTree.copy(edgeSource)) == 0) {
				stSource = finalSteinerTree.newNode(edgeSource);
			}
			if ((stTarget = finalSteinerTree.copy(edgeTarget)) == 0) {
				stTarget = finalSteinerTree.newNode(edgeTarget);
			}
			if ((*it)->source() == finalSteinerTree.original(stSource)) {
				newE = finalSteinerTree.newEdge(stSource, stTarget, wG.weight(*it));
				finalSteinerTree.setEdge((*it), newE);
			}
		}
	}
}

} // end namespace ogdf

#endif /* OGDF_MINIMUM_STEINER_TREE_KOU_H_ */
