/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Mehlhorn's minimum steiner tree
 * 		  2(1-1/l)-approximation algorithm
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

#ifndef OGDF_MIN_STEINER_TREE_MEHLHORN_H_
#define OGDF_MIN_STEINER_TREE_MEHLHORN_H_

#include <ogdf/basic/List.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/internal/steinertree/Voronoi.h>
#include <ogdf/module/MinSteinerTreeModule.h>

namespace ogdf {

/*!
 * \brief Represents a triple as specified in the algorithms description (see paper)
 */
template<typename T>
struct MehlhornTriple {
	node u;
	node v;
	T value;
	edge bridge;
};

/*!
 * \brief Helper class to sort MehlhornTriples lexicographically
 */
template<typename T>
class MehlhornTripleBucketMaxFunc: public BucketFunc<MehlhornTriple<T> > {
public:
	MehlhornTripleBucketMaxFunc() {
	}
	;

	int getBucket(const MehlhornTriple<T>& MT) {
		int source_index = MT.u->index();
		int target_index = MT.v->index();
		OGDF_ASSERT(source_index != target_index);
		if (source_index < target_index)
			return target_index;
		else
			return source_index;
	}
	;
};

/*!
 * \brief Helper class to sort MehlhornTriples lexicographically
 */
template<typename T>
class MehlhornTripleBucketMinFunc: public BucketFunc<MehlhornTriple<T> > {
public:
	MehlhornTripleBucketMinFunc() {
	}
	;

	int getBucket(const MehlhornTriple<T>& MT) {
		int source_index = MT.u->index();
		int target_index = MT.v->index();
		OGDF_ASSERT(source_index != target_index);
		if (source_index < target_index)
			return source_index;
		else
			return target_index;
	}
	;
};

/*!
 * \brief This class implements the Minimum Steiner Tree 2-approximation algorithm by Mehlhorn.
 *
 * This implementation is based on:
 *
 * (K. Mehlhorn, A faster approximation algorithm for the steiner problem in graphs,
 * Information Processing Letters, volume 27, number 3, pages 125-128, 1998)
 */
template<typename T>
class MinSteinerTreeMehlhorn: public MinSteinerTreeModule<T> {
public:
	MinSteinerTreeMehlhorn() {
	}
	;
	virtual ~MinSteinerTreeMehlhorn() {
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
	 * @param voronoi Voronoi regions (providing a mapping from each edge in the complete terminal graph to the shortest path in the original graph)
	 * @param bridges list of edges connecting terminal nodes by voronoi regions
	 * @param completeTerminalGraph the resulting complete terminal graph
	 */
	void calculateCompleteGraph(const EdgeWeightedGraph<T> &wG, const List<node> &terminals, const Voronoi<T> &voronoi,
			EdgeArray<edge> &bridges, EdgeWeightedGraphCopy<T> &completeTerminalGraph);

	/*!
	 * \brief Swaps an edge in the complete terminal graph with the corresponding shortest path in the original graph
	 * @param completeTerminalGraph the complete terminal graph
	 * @param voronoi Voronoi regions (providing a mapping from each edge in the complete terminal graph to the shortest path in the original graph)
	 * @param mstPred predecessor data structure of a minimum terminal spanning tree
	 * @param bridges list of edges connecting terminal nodes by voronoi regions
	 * @param finalSteinerTree the resulting steiner tree
	 * @param wG the original graph
	 */
	void reinsertShortestPaths(EdgeWeightedGraphCopy<T> &completeTerminalGraph, const Voronoi<T> &voronoi,
			const NodeArray<edge> &mstPred, const EdgeArray<edge> &bridges, EdgeWeightedGraphCopy<T> &finalSteinerTree,
			const EdgeWeightedGraph<T> &wG);

	/*!
	 * \brief Inserts a shortest path corresponding to an edge in the complete terminal graph
	 * @param u terminal node needed to access the according predecessor edge in the minimum terminal spanning tree
	 * @param voronoi Voronoi regions (contains shortest paths)
	 * @param finalSteinerTree the resulting steiner tree
	 * @param wG the original graph
	 */
	void insertPath(node u, const Voronoi<T> &voronoi,
			EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG);
};

} // end namespace ogdf

// ============= Implementation =================

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/List.h>
#include <ogdf/graphalg/Dijkstra.h>
#include <sstream>

namespace ogdf {

template<typename T>
T MinSteinerTreeMehlhorn<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	EdgeWeightedGraphCopy<T> completeTerminalGraph;
	completeTerminalGraph.createEmpty(G);

	forall_listiterators(node, it, terminals) {
		completeTerminalGraph.newNode(*it);
	}

	EdgeArray<edge> bridges(completeTerminalGraph);
	Voronoi<T> voronoi(G, terminals);

	calculateCompleteGraph(G, terminals, voronoi, bridges, completeTerminalGraph);

	NodeArray<edge> mstPred(completeTerminalGraph);
	EdgeArray<bool> isInTree(completeTerminalGraph);
	T mstWeight = computeMinST(completeTerminalGraph, completeTerminalGraph.edgeWeights(), mstPred, isInTree);

	finalSteinerTree = new EdgeWeightedGraphCopy<T>;
	finalSteinerTree->createEmpty(G);

	reinsertShortestPaths(completeTerminalGraph, voronoi, mstPred, bridges, *finalSteinerTree, G);

	EdgeArray<bool> isInStTree(*finalSteinerTree);
	mstWeight = computeMinST(*finalSteinerTree, finalSteinerTree->edgeWeights(), isInStTree);
	mstWeight -= MinSteinerTreeModule<T>::eliminate1DegreeSteinerNodes(*finalSteinerTree, isTerminal);
	OGDF_ASSERT(MinSteinerTreeModule<T>::isSteinerTree(G, terminals, *finalSteinerTree));

	return mstWeight;
}

template<typename T>
void MinSteinerTreeMehlhorn<T>::calculateCompleteGraph(const EdgeWeightedGraph<T> &wG, const List<node> &terminals, const Voronoi<T> &voronoi, EdgeArray<edge> &bridges, EdgeWeightedGraphCopy<T> &completeTerminalGraph)
{
	// extract complete graph edges
	List<MehlhornTriple<T> > triples;
	edge e;
	forall_edges(e, wG) {
		MehlhornTriple<T> triple;
		triple.u = voronoi.seed(e->source());
		triple.v = voronoi.seed(e->target());
		if (triple.u != triple.v) {
			triple.value = voronoi.distance(e->source()) + voronoi.distance(e->target()) + wG.weight(e);
			triple.bridge = e;
			triples.pushBack(triple);
		}
	}

	MehlhornTripleBucketMaxFunc<T> mtbMax;
	MehlhornTripleBucketMinFunc<T> mtbMin;

	triples.bucketSort(0, wG.numberOfNodes() - 1, mtbMax);
	triples.bucketSort(0, wG.numberOfNodes() - 1, mtbMin);

	int currentSource = triples.front().u->index();
	int currentTarget = triples.front().v->index();
	ListConstIterator<MehlhornTriple<T> > minTriple = triples.begin();

	for (ListConstIterator<MehlhornTriple<T> > mtIt = triples.begin().succ(); mtIt.valid(); ++mtIt) {
		if (((*mtIt).u->index() == currentSource && (*mtIt).v->index() == currentTarget)
		 || ((*mtIt).u->index() == currentTarget && (*mtIt).v->index() == currentSource)) {
			if ((*mtIt).value < (*minTriple).value) {
				minTriple = mtIt;
			}
		} else {
			// add new direct edge
			edge tmp = completeTerminalGraph.newEdge(completeTerminalGraph.copy((*minTriple).u),
					completeTerminalGraph.copy((*minTriple).v), (*minTriple).value);
			bridges[tmp] = (*minTriple).bridge;

			currentSource = (*mtIt).u->index();
			currentTarget = (*mtIt).v->index();
			minTriple = mtIt;
		}
	}
	// insert last triple
	if (minTriple.valid()) {
		edge tmp = completeTerminalGraph.newEdge(completeTerminalGraph.copy((*minTriple).u),
				completeTerminalGraph.copy((*minTriple).v), (*minTriple).value);
		bridges[tmp] = (*minTriple).bridge;
	}
}

template<typename T>
void MinSteinerTreeMehlhorn<T>::reinsertShortestPaths(EdgeWeightedGraphCopy<T> &completeTerminalGraph, const Voronoi<T> &voronoi, const NodeArray<edge> &mstPred, const EdgeArray<edge> &bridges, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG)
{
	node u;
	forall_nodes(u, completeTerminalGraph) {
		if (mstPred[u] != 0) {
			edge bridge = bridges[mstPred[u]];
			node v = bridge->source();
			node w = bridge->target();
			insertPath(v, voronoi, finalSteinerTree, wG);
			insertPath(w, voronoi, finalSteinerTree, wG);
			edge e = finalSteinerTree.newEdge(finalSteinerTree.copy(v), finalSteinerTree.copy(w), wG.weight(bridge));
			finalSteinerTree.setEdge(bridge, e);
		}
	}
}

template<typename T>
void MinSteinerTreeMehlhorn<T>::insertPath(node u, const Voronoi<T> &voronoi, EdgeWeightedGraphCopy<T> &finalSteinerTree, const EdgeWeightedGraph<T> &wG)
{
	node currentSource;
	node currentTarget = finalSteinerTree.copy(u);
	if (!currentTarget) {
		currentTarget = finalSteinerTree.newNode(u);
	}
	edge e = voronoi.predecessorEdge(u);
	edge newE;
	while (e && finalSteinerTree.chain(e).empty()) { // e is not in ST yet
		if ((currentSource = finalSteinerTree.copy(e->opposite(finalSteinerTree.original(currentTarget)))) == 0) {
			currentSource = finalSteinerTree.newNode(e->opposite(finalSteinerTree.original(currentTarget)));
		}
		if (finalSteinerTree.original(currentSource) == e->source()) {
			newE = finalSteinerTree.newEdge(currentSource, currentTarget, wG.weight(e));
		} else {
			newE = finalSteinerTree.newEdge(currentTarget, currentSource, wG.weight(e));
		}
		finalSteinerTree.setEdge(e, newE);
		currentTarget = currentSource;
		e = voronoi.predecessorEdge(finalSteinerTree.original(currentTarget));
	}
}

} // end namespace ogdf

#endif /* OGDF_MIN_STEINER_TREE_MEHLHORN_H_ */
