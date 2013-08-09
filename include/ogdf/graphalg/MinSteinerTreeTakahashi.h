/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the 2(1-1/l)-approximation algorithm for
 * 		  the minimum steiner tree problem by Matsuyama and Takahashi
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

#ifndef OGDF_MIN_STEINER_TREE_TAKAHASHI_H_
#define OGDF_MIN_STEINER_TREE_TAKAHASHI_H_

#include <ogdf/basic/List.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/module/MinSteinerTreeModule.h>


namespace ogdf {

/*!
 * \brief This class implements the minimum steiner tree 2-approximation algorithm by Takahashi and Matsuyama with improvements proposed by Poggi de Aragao et al..
 *
 * This implementation is based on:
 *
 * (H. Takahashi and A. Matsuyama, An approximate solution for the steiner problem in graphs, Math. Japonica,
 * volume 24, number 6, pages 573-577, 1980)
 *
 * (M. Poggi de Aragao, C. Riberiro, E. Uchoa, R. Werneck, Hybrid Local Search for the Steiner Problem in Graphs,
 * MIC 2001, pages 429-433, 2001)
 */
template<typename T>
class MinSteinerTreeTakahashi: public MinSteinerTreeModule<T> {
public:
	MinSteinerTreeTakahashi() {
	}
	;
	virtual ~MinSteinerTreeTakahashi() {
	}
	;

	/*!
	 * \brief An extended call method with specific start node and additional steiner nodes
	 * that have to be considered terminals \see MinSteinerTreeModule::call
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param newTerminals List of steiner nodes that shall become terminals
	 * @param isNewTerminal Incidence vector for the steiner nodes that shall become terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @param startNode The source node for the special Dijkstra call
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, const List<node> &newTerminals, const NodeArray<bool> &isNewTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree, const node startNode);

	/*!
	 * \brief An extended call method with specific start node \see MinSteinerTreeModule::call
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @param startNode The source node for the special Dijkstra call
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree, const node startNode);

	/*!
	 * \brief Builds a minimum steiner tree given a weighted graph and a list of terminals \see MinSteinerTreeModule::call
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree);

	/*!
	 * \brief Returns the overall running time
	 * @return the overall running time
	 */
	double runningTime();

protected:
	/*!
	 * Modified Dijkstra algorithm to solve the minimum steiner tree problem
	 * @param wG the original graph
	 * @param intermediateTerminalSpanningTree intermediate terminal spanning tree
	 * @param s source node to start from
	 * @param numberOfTerminals number of terminal nodes
	 * @param predecessor predecessor edge in shortest path tree
	 * @param distance //!< distances to nodes
	 * @param isTerminal terminal incivende vector
	 * @return the weight of the intermediateTerminalSpanningTree
	 */
	T terminalDijkstra(const EdgeWeightedGraph<T> &wG, EdgeWeightedGraphCopy<T> &intermediateTerminalSpanningTree, const node s, int numberOfTerminals, NodeArray<edge> &predecessor, NodeArray<T> &distance, const NodeArray<bool> &isTerminal);

private:
	double m_totalTime; //!< Stores the total running time
};

} // end namespace ogdf

// ============= Implementation =================

#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/List.h>
#include <ogdf/graphalg/Dijkstra.h>

namespace ogdf {

template<typename T>
T MinSteinerTreeTakahashi<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	return call(G, terminals, isTerminal, terminals, isTerminal, finalSteinerTree, terminals.front());
}

template<typename T>
T MinSteinerTreeTakahashi<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree, const node startNode)
{
	return call(G, terminals, isTerminal, terminals, isTerminal, finalSteinerTree, startNode);
}

template<typename T>
T MinSteinerTreeTakahashi<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, const List<node> &newTerminals, const NodeArray<bool> &isNewTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree, const node startNode)
{
	double t_total;
	usedTime(t_total);

	EdgeWeightedGraphCopy<T> terminalSpanningTree;
	terminalSpanningTree.createEmpty(G);
	NodeArray<edge> predecessor(G);
	NodeArray<T> distance(G);
	T mstWeight = terminalDijkstra(G, terminalSpanningTree, startNode, newTerminals.size(), predecessor, distance,
			isNewTerminal);

	finalSteinerTree = new EdgeWeightedGraphCopy<T>(G);
	node u;
	List<node> toDelete;
	forall_nodes(u, G) {
		if (terminalSpanningTree.copy(u) == 0) {
			finalSteinerTree->delNode(finalSteinerTree->copy(u));
		}
	}

	EdgeArray<bool> isInTree(*finalSteinerTree);
	computeMinST(*finalSteinerTree, finalSteinerTree->edgeWeights(), isInTree);

	mstWeight = 0;
	edge nextEdge;
	for (edge e = finalSteinerTree->firstEdge(); e; e = nextEdge) {
		nextEdge = e->succ();
		if (!isInTree[e]) {
			finalSteinerTree->delEdge(e);
		} else {
			mstWeight += finalSteinerTree->weight(e);
		}
	}

	mstWeight -= MinSteinerTreeModule<T>::eliminate1DegreeSteinerNodes(*finalSteinerTree, isTerminal);

	m_totalTime = usedTime(t_total);

	OGDF_ASSERT(MinSteinerTreeModule<T>::isSteinerTree(G, terminals, *finalSteinerTree));

	return mstWeight;
}

template<typename T>
T MinSteinerTreeTakahashi<T>::terminalDijkstra(const EdgeWeightedGraph<T> &wG,
		EdgeWeightedGraphCopy<T> &intermediateTerminalSpanningTree, const node s, int numberOfTerminals,
		NodeArray<edge> &predecessor, NodeArray<T> &distance, const NodeArray<bool> &isTerminal)
{
	T mstWeight = 0;
	node v;
	BinaryHeap2<T, node> queue(wG.numberOfNodes()); //priority queue
	int* qpos = new int[wG.numberOfNodes()];
	NodeArray<int> vIndex(wG);
	NodeArray<T> bestDistance(wG);
	NodeArray<bool> isInQueue(wG);
	int i = 0;

	forall_nodes(v, wG) {
		vIndex[v] = i;
		distance[v] = numeric_limits<T>::max();
		bestDistance[v] = numeric_limits<T>::max();
		predecessor[v] = 0;
		queue.insert(v, distance[v], &(qpos[i++]));
		isInQueue[v] = true;
	}

	distance[s] = 0;
	bestDistance[s] = 0;
	queue.decreaseKey(qpos[vIndex[s]], (T) 0);
	node tmpS, tmpT;
	edge e, tmpE;
	int terminalsFound = 1;
	while (!queue.empty() && terminalsFound < numberOfTerminals) {
		v = queue.extractMin();
		isInQueue[v] = false;
		bestDistance[v] = distance[v];
		if (isTerminal[v] && distance[v] > 0) { // new terminal found
			terminalsFound++;
			// insert path from new node to old tree
			tmpT = intermediateTerminalSpanningTree.newNode(v);
			while (distance[v] > 0) {
				distance[v] = 0;
				//				bestDistance[v] = 0;
				queue.insert(v, distance[v], &(qpos[vIndex[v]]));
				isInQueue[v] = true;
				if ((tmpS = intermediateTerminalSpanningTree.copy(predecessor[v]->opposite(v))) == 0) {
					tmpS = intermediateTerminalSpanningTree.newNode(predecessor[v]->opposite(v));
				}
				if (predecessor[v]->target() == v) {
					tmpE = intermediateTerminalSpanningTree.newEdge(tmpS, tmpT, wG.weight(predecessor[v]));
				} else {
					tmpE = intermediateTerminalSpanningTree.newEdge(tmpT, tmpS, wG.weight(predecessor[v]));
				}
				mstWeight += wG.weight(predecessor[v]);
				intermediateTerminalSpanningTree.setEdge(predecessor[v], tmpE);
				tmpT = tmpS;
				v = predecessor[v]->opposite(v);
			}
		} else {
			forall_adj_edges(e, v) {
				node w = e->opposite(v);
				if (distance[w] > distance[v] + wG.weight(e) && bestDistance[w] >= distance[w]) {
					distance[w] = distance[v] + wG.weight(e);
					if (!isInQueue[w]) {
						queue.insert(w, distance[w], &(qpos[vIndex[w]]));
						isInQueue[w] = true;
					} else {
						queue.decreaseKey(qpos[vIndex[w]], distance[w]);
					}
					predecessor[w] = e;
				}
			}
		}
	}
	delete[] qpos;
	return mstWeight;
}

template<typename T>
double MinSteinerTreeTakahashi<T>::runningTime() {
	return m_totalTime;
}

} // end namespace ogdf

#endif /* OGDF_MIN_STEINER_TREE_TAKAHASHI_H_ */
