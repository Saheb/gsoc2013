/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the 1.55-approximation algorithm for the minimum
 * 		  steiner tree problem by Robins and Zelikovsky
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

#ifndef MIN_STEINER_TREE_RZ_LOSS_OGDF_H_
#define MIN_STEINER_TREE_RZ_LOSS_OGDF_H_

#include <ogdf/basic/List.h>
#include <ogdf/basic/HashArray.h>
#include <ogdf/basic/Queue.h>
#include <ogdf/internal/steinertree/FullComponent.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/module/MinSteinerTreeModule.h>


namespace ogdf {

/*!
 * \brief This class implements the loss-contracting approximation algorithm
 * for the steiner tree problem by Robins and Zelikovsky.
 *
 * This implementation is based on:
 *
 * (G. Robins, A. Zelikovsky, Improved Steiner Tree Approximation in Graphs,
 * SODA 2000, pages 770-779, SIAM, 2000)
 */
template<typename T>
class MinSteinerTreeRZLoss: public MinSteinerTreeModule<T> {
public:
	MinSteinerTreeRZLoss() {
		setK(3);
	}
	;
	MinSteinerTreeRZLoss(int v) {
		setK(v);
	}
	;
	virtual ~MinSteinerTreeRZLoss() {
	}
	;
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
	 * \brief Sets the maximal number of terminals in a full component
	 * @param v the maximal number of terminals in a full component
	 */
	void setK(int v) {
		k = v;
	}

	//! Returns the number of generated components
	long numberOfGeneratedComponents() {
		return m_componentsGenerated;
	}

	//! Returns the number of contracted components
	long numberOfContractedComponents() {
		return m_componentsContracted;
	}

	//! Returns the number of components lookups during execution time
	long numberOfComponentLookUps() {
		return m_componentsLookUps;
	}

	//! Returns the total running time
	double runningTime() {
		return m_totalTime;
	}

private:

	/*!
	 * \brief Finds all full components up to currentK terminals per full component.
	 * @param currentK Find full component with this number of terminals
	 * @param terminals List of all terminals
	 * @param isTerminal Terminal node incidence vector
	 * @param distance Distance matrix for the original graph
	 * @param path Shortest paths between nodes in the original graph
	 * @param fullComponents List of all found full components
	 * @param tree Current minimal terminal spanning tree
	 * @param save Data structure for calculating save edges
	 */
	void findFullComponents(int currentK, const List<node> &terminals, const NodeArray<bool> &isTerminal, NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path, List<FullComponent<T>*> &fullComponents, EdgeWeightedGraphCopy<T> &tree, HashArray<int, edge> &save);

	/*!
	 * \brief Builds a minimum terminal spanning tree via a MST call on the complete terminal graph
	 * @param steinerTree The resulting minimum terminal spanning tree
	 * @return The sum of all edge costs of the returned minimum terminal spanning tree
	 */
	T generateMinimumSpanningTree(EdgeWeightedGraphCopy<T> &steinerTree);

	/*!
	 * \brief Builds the complete terminal graph as well as the minimum terminal spanning tree
	 * @param steinerTree The build minimum terminal spanning tree
	 * @param distance Distance matrix for the original graph
	 * @param path Shortest paths between nodes in the original graph
	 */
	void createSteinerTreeAndCompleteTerminalGraph(EdgeWeightedGraphCopy<T> &steinerTree, NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path);

	/*!
	 * \brief Traverses over all full components and finds the one with the highest win-objective.
	 * @param fullComponents List of full Components
	 * @param steinerTree Current steiner tree
	 * @param maxComponent The returned full component with highest win-objective
	 * @param save Data structure for calculating save values
	 * @return the win-objective of the returned full component
	 */
	double findMaxComponent(List<FullComponent<T>*> &fullComponents, EdgeWeightedGraphCopy<T> &steinerTree, FullComponent<T>* &maxComponent, HashArray<int, edge> &save);

	/*!
	 * \brief Calculates the gain for a full component
	 * @param fullComponent List of full Components
	 * @param steinerTree Current steiner tree
	 * @param save Data structure for calculating save values
	 * @return the gain (cost of save edges) of the returned full component
	 */
	double gain(FullComponent<T> &fullComponent, const EdgeWeightedGraphCopy<T> &steinerTree, HashArray<int, edge> &save);

	/*!
	 * \brief Integrates the given full component into the complete terminal graph
	 * @param fullComponent Full component to be integrated
	 * @param newTerminal The list of steiner nodes to be become terminals
	 * @param isNewTerminal Incidence vector for terminal becoming steiner nodes
	 */
	void augmentCompleteTerminalGraph(FullComponent<T> &fullComponent, List<node> &newTerminal, NodeArray<bool> &isNewTerminal);

	/*!
	 * \brief Contracts a full component and integrates it into the given steiner tree
	 * @param steinerTree The steiner tree into which the full component is integrated
	 * @param fullComponent The full component to be contracted and inserted
	 */
	void contractAndAugmentComponent(EdgeWeightedGraphCopy<T> &steinerTree, FullComponent<T> &fullComponent);

	/*!
	 * \brief This function enumerates all possible k-sized subsets of terminals in lexicographical order.
	 *
	 * A function call returns the next one in the order.
	 * @param tuple This array contains the next subset of terminal IDs in the lexicographic order
	 * @param n This is the number of possible elements to chose from (e.g. number of terminals)
	 * @return true if there are further subsets and false if this is the last subset
	 */
	bool nextTuple(Array<int> &tuple, int n);

	/*!
	 * \brief Does an APSP on the original graph.
	 * @param distance The distance matrix
	 * @param path All shortest path between all nodes
	 */
	void allPairsShortestPaths(NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path);

	/*!
	 * \brief Creates a graph representing a full component
	 * @param steinerTuple Array with IDs of steiner nodes to integrate into the component
	 * @param steinerNode Array of all steiner nodes
	 * @param fullSteinerComponentGraph The resulting full component graph
	 * @param distance The complete distance matrix
	 * @param currentK The number of terminals in the full component
	 */
	void createFullSteinerComponentGraph(
		const Array<int> &steinerTuple,
		const Array<node> &steinerNode,
		EdgeWeightedGraphCopy<T> &fullSteinerComponentGraph,
		const NodeArray<NodeArray<T> > &distance,
		const int currentK);

	/*!
	 * \brief Inserts steiner points and their shortest path into a full component
	 * @param s First steiner node
	 * @param t Second steiner node
	 * @param fullSteinerComponentTree Full component to be inserted into
	 * @param distance Complete distance matrix
	 * @param path All shortest paths in the original graph
	 */
	void appendSteinerPointsToFullSteinerComponentTree(node s, node t, EdgeWeightedGraphCopy<T> &fullSteinerComponentTree,
			NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path);

	/*!
	 * \brief Based on a full component graph this method creates a full component by calculation a MST
	 * @param fullSteinerGraph The basic full component graph
	 * @param distance Complete distance matrix
	 * @param path All shortest path
	 * @param fullSteinerComponentTree The resulting full steiner component tree
	 */
	void createFullSteinerComponentTree(EdgeWeightedGraphCopy<T> &fullSteinerGraph, NodeArray<NodeArray<T> > &distance,
			NodeArray<NodeArray<List<edge> > > &path, EdgeWeightedGraphCopy<T> &fullSteinerComponentTree);

	/*!
	 * \brief Updates the current steiner tree my calculating a MST
	 * @param steinerTree The steiner tree to be updated
	 */
	void updateTree(EdgeWeightedGraphCopy<T> &steinerTree);

	/*!
	 * \brief Builds the table for looking up save edges recursively
	 * @param mst Minimum terminal spanning tree
	 * @param u "Starting node" for the calculation
	 * @param processedNodes List of nodes processed in this function call
	 * @param saveArray The resulting map for save edges
	 */
	void buildSaveArray(EdgeWeightedGraphCopy<T> &mst, node u, List<node> &processedNodes, HashArray<int, edge> &saveArray);

	const EdgeWeightedGraph<T> *m_originalGraph; //!< The original edge-weighted graph
	EdgeWeightedGraphCopy<T> *m_completeTerminalGraph; //!< The complete terminal graph
	const NodeArray<bool> *m_isTerminal; //!< Incidence vector for terminal nodes
	const List<node> *m_terminals; //!< List of terminal nodes
	int k; //!< Parameter for the number of terminals in a full component

	double m_startTime; //!< Stores the start time of the algorithm
	double m_totalTime; //!< Stores the overall running time
	long m_componentsGenerated; //!< Number of generated components
	long m_componentsContracted; //!< Number of contracted components
	long m_componentsLookUps; //!< Number of components lookups
};

}

// Implementation

#include <ogdf/basic/Math.h>
#include <ogdf/graphalg/Dijkstra.h>
#include <ogdf/graphalg/MinSteinerTreeTakahashi.h>

namespace ogdf {

template<typename T>
T MinSteinerTreeRZLoss<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	usedTime(m_startTime);

	m_originalGraph = &G;
	m_terminals = &terminals;
	m_isTerminal = &isTerminal;
	m_completeTerminalGraph = new EdgeWeightedGraphCopy<T>();

	if (k > m_terminals->size()) {
		cerr << "ERROR: k is greater than the number of terminals.\n";
		OGDF_THROW_PARAM(AlgorithmFailureException, afcIllegalParameter);
	}

	List<node> newTerminals;
	NodeArray<bool> isNewTerminal(G);

	NodeArray<NodeArray<T> > distance(G);
	NodeArray<NodeArray<List<edge> > > path(G);
	node u;
	forall_nodes(u, G) {
		distance[u] = NodeArray<T>(G);
		path[u] = NodeArray<List<edge> >(G);
		if (isTerminal[u]) {
			isNewTerminal[u] = true;
			newTerminals.pushBack(u);
		} else {
			isNewTerminal[u] = false;
		}
	}
	allPairsShortestPaths(distance, path);

	EdgeWeightedGraphCopy<T> steinerTree;
	createSteinerTreeAndCompleteTerminalGraph(steinerTree, distance, path);

	HashArray<int, edge> save;
	List<node> processedNodes;
	buildSaveArray(steinerTree, steinerTree.firstNode(), processedNodes, save);

	m_componentsGenerated = 0;
	m_componentsLookUps = 0;
	m_componentsContracted = 0;

	List<FullComponent<T>*> fullComponents;
	for (int currentK = 3; currentK <= k; currentK++) {
		findFullComponents(currentK, terminals, isTerminal, distance, path, fullComponents, steinerTree, save);
	}

	m_componentsGenerated = fullComponents.size();

	double r;
	while (fullComponents.size() > 0) {
		FullComponent<T> *maxComponent;
		r = findMaxComponent(fullComponents, steinerTree, maxComponent, save);
		if (r <= 0.0)
			break;
		m_componentsContracted++;
		augmentCompleteTerminalGraph(*maxComponent, newTerminals, isNewTerminal);
		contractAndAugmentComponent(steinerTree, *maxComponent);
		updateTree(steinerTree);
	}

	T tmpMstWeight = 0, bestMstWeight = numeric_limits<T>::max();

	MinSteinerTreeTakahashi<T> mstt;

	OGDF_ASSERT(!finalSteinerTree);

	forall_listiterators(node, it, *m_terminals) {
		EdgeWeightedGraphCopy<T> *tmpSteinerTree = new EdgeWeightedGraphCopy<T>(*m_originalGraph);
		tmpMstWeight = mstt.call(*m_originalGraph, *m_terminals, *m_isTerminal, newTerminals, isNewTerminal, tmpSteinerTree, *it);
		if (tmpMstWeight < bestMstWeight) {
			bestMstWeight = tmpMstWeight;
			if (finalSteinerTree) {
				delete finalSteinerTree;
			}
			finalSteinerTree = tmpSteinerTree;
		} else {
			delete tmpSteinerTree;
		}
	}

	forall_listiterators(FullComponent<T>*, it, fullComponents) {
		delete *it;
	}

	m_totalTime = usedTime(m_startTime);

	return bestMstWeight;
}

template<typename T>
void MinSteinerTreeRZLoss<T>::updateTree(EdgeWeightedGraphCopy<T> &steinerTree)
{
	EdgeArray<bool> isTreeEdge(steinerTree);
	computeMinST(steinerTree, steinerTree.edgeWeights(), isTreeEdge);

	edge nextEdge;
	for (edge e = steinerTree.firstEdge(); e; e = nextEdge) {
		nextEdge = e->succ();
		if (!isTreeEdge[e]) {
			steinerTree.delEdge(e);
		}
	}
}

template<typename T>
void MinSteinerTreeRZLoss<T>::createSteinerTreeAndCompleteTerminalGraph(EdgeWeightedGraphCopy<T> &steinerTree,
		NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path)
{
	T tmp;

	m_completeTerminalGraph->createEmpty(*m_originalGraph);
	steinerTree.createEmpty(*m_originalGraph);

	forall_listiterators(node, it, *m_terminals) {
		m_completeTerminalGraph->newNode(*it);
	}

	node u, v;
	u = m_completeTerminalGraph->firstNode();

	while (u != 0) {
		v = u->succ();
		while (v != 0) {
			tmp = distance[m_completeTerminalGraph->original(u)][m_completeTerminalGraph->original(v)];
			m_completeTerminalGraph->newEdge(u, v, tmp);
			v = v->succ();
		}
		u = u->succ();
	}

	generateMinimumSpanningTree(steinerTree);
}

template<typename T>
T MinSteinerTreeRZLoss<T>::generateMinimumSpanningTree(EdgeWeightedGraphCopy<T> &steinerTree)
{
	node u;
	edge e;
	T tmp;

	forall_nodes(u, *m_completeTerminalGraph) {
		steinerTree.newNode(m_completeTerminalGraph->original(u));
	}

	NodeArray<edge> pred(*m_completeTerminalGraph);
	EdgeArray<bool> isTreeEdge(*m_completeTerminalGraph);
	T mstF = computeMinST(*m_completeTerminalGraph, m_completeTerminalGraph->edgeWeights(), pred, isTreeEdge);

	forall_nodes(u, *m_completeTerminalGraph) {
		e = pred[u];
		if (e != 0) {
			tmp = m_completeTerminalGraph->weight(e);
			steinerTree.newEdge(steinerTree.copy(m_completeTerminalGraph->original(e->source())), steinerTree.copy(m_completeTerminalGraph->original(e->target())), tmp);
		}
	}

	return mstF;
}

template<typename T>
void MinSteinerTreeRZLoss<T>::findFullComponents(int currentK, const List<node> &terminals,
		const NodeArray<bool> &isTerminal, NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path,
		List<FullComponent<T>*> &fullComponents, EdgeWeightedGraphCopy<T> &tree, HashArray<int, edge> &save)
{
	FullComponent<T> *c;
	T minDist;
	node minDistSteinerNode;
	List<node> currentTerminals;
	List<node> steinerNodeList;
	double saveTime = m_startTime;
	node u;
	forall_nodes(u, *m_originalGraph) {
		if (!isTerminal[u]) {
			steinerNodeList.pushBack(u);
		}
	}

	Array<node> steinerNode(steinerNodeList.size());
	Array<node> terminalNode(terminals.size());
	Array<int> steinerTuple(currentK - 2);
	Array<int> terminalTuple(currentK);
	Array<node> terminalToSteinerNodeTuple(currentK);
	int i = 0;
	forall_listiterators(node, it, steinerNodeList) {
		steinerNode[i++] = *it;
	}
	i = 0;
	forall_listiterators(node, it, terminals) {
		terminalNode[i++] = *it;
	}
	for (i = 0; i < currentK - 2; i++) {
		steinerTuple[i] = i;
	}

	node s, t;
	EdgeWeightedGraphCopy<T> *fullSteinerComponentGraph;
	EdgeWeightedGraphCopy<T> *fullSteinerComponentTree;

	do { // iterate over steinerTuples
		saveTime = m_startTime;
		for (i = 0; i < currentK; i++) {
			terminalTuple[i] = i;
			terminalToSteinerNodeTuple[i] = 0;
		}

		do { // iterate over terminalTuples
			bool fail = false;
			for (int p = 0; p < currentK && !fail; p++) {
				t = terminalNode[terminalTuple[p]];
				minDist = numeric_limits<T>::max();
				for (int q = 0; q < currentK - 2 && !fail; q++) {
					s = steinerNode[steinerTuple[q]];
					if (distance[s][t] < numeric_limits<T>::max()) {
						if (distance[s][t] < minDist) {
							minDist = distance[s][t];
							minDistSteinerNode = s;
						}
					} else {
						fail = true;
					}
				}
				if (!fail) {
					terminalToSteinerNodeTuple[p] = minDistSteinerNode;
				}
			}

			if (!fail) {
				fullSteinerComponentGraph = new EdgeWeightedGraphCopy<T>();
				fullSteinerComponentTree = new EdgeWeightedGraphCopy<T>();
				fullSteinerComponentGraph->createEmpty(*m_originalGraph);
				fullSteinerComponentTree->createEmpty(*m_originalGraph);

				createFullSteinerComponentGraph(steinerTuple, steinerNode, *fullSteinerComponentGraph, distance, currentK);
				createFullSteinerComponentTree(*fullSteinerComponentGraph, distance, path, *fullSteinerComponentTree);
				c = new FullComponent<T>(fullSteinerComponentTree);
				for (int p = 0; p < currentK && !fail; p++) {
					(*c).appendTerminalToSteinerTree(terminalToSteinerNodeTuple[p], terminalNode[terminalTuple[p]],
							distance, path, fullSteinerComponentTree);
					(*c).addTerminalToList(terminalNode[terminalTuple[p]]);
				}
				if (gain(*c, tree, save) - (*c).cost() > 0) {
					(*c).calculateLoss(fullSteinerComponentTree);
					fullComponents.pushBack(c);
				} else {
					delete c;
				}
				delete fullSteinerComponentGraph;
				delete fullSteinerComponentTree;
			}
		} while (nextTuple(terminalTuple, terminals.size()));
	} while (nextTuple(steinerTuple, steinerNodeList.size()) && usedTime(saveTime) < 3600);
}

template<typename T>
void MinSteinerTreeRZLoss<T>::createFullSteinerComponentTree(EdgeWeightedGraphCopy<T> &fullSteinerGraph,
		NodeArray<NodeArray<T> > &distance, NodeArray<NodeArray<List<edge> > > &path,
		EdgeWeightedGraphCopy<T> &fullSteinerComponentTree)
{
	node u;

	forall_nodes(u, fullSteinerGraph) {
		fullSteinerComponentTree.newNode(fullSteinerGraph.original(u));
	}

	NodeArray<edge> pred(fullSteinerGraph);
	EdgeArray<bool> isTreeEdge(fullSteinerGraph);
	computeMinST(fullSteinerGraph, fullSteinerGraph.edgeWeights(), pred, isTreeEdge);

	forall_nodes(u, fullSteinerGraph) {
		if (pred[u] != 0) {
			appendSteinerPointsToFullSteinerComponentTree(fullSteinerGraph.original(u),
					fullSteinerGraph.original(pred[u]->opposite(u)), fullSteinerComponentTree, distance, path);
		}
	}
}

template<typename T>
void MinSteinerTreeRZLoss<T>::appendSteinerPointsToFullSteinerComponentTree(node s, node t,
		EdgeWeightedGraphCopy<T> &fullSteinerComponentTree, NodeArray<NodeArray<T> > &distance,
		NodeArray<NodeArray<List<edge> > > &path)
{
	node u, v;
	edge e, f;
	forall_listiterators(edge, it, path[s][t]) {
		e = *it;
		if (fullSteinerComponentTree.chain(e).size() == 0) {

			if ((u = fullSteinerComponentTree.copy(e->source())) == 0) {
				u = fullSteinerComponentTree.newNode(e->source());
			}
			if ((v = fullSteinerComponentTree.copy(e->target())) == 0) {
				v = fullSteinerComponentTree.newNode(e->target());
			}
			f = fullSteinerComponentTree.newEdge(u, v,
					distance[fullSteinerComponentTree.original(u)][fullSteinerComponentTree.original(v)]);
			fullSteinerComponentTree.setEdge(e, f);
		}
	}
}

template<typename T>
void MinSteinerTreeRZLoss<T>::createFullSteinerComponentGraph(
		const Array<int> &steinerTuple,
		const Array<node> &steinerNode,
		EdgeWeightedGraphCopy<T> &fullSteinerComponentGraph,
		const NodeArray<NodeArray<T> > &distance,
		const int currentK)
{
	node u, v, s, t;
	for (int i = 0; i < currentK - 2; i++) {
		u = steinerNode[steinerTuple[i]];
		s = fullSteinerComponentGraph.newNode(u);
		for (int j = i - 1; j >= 0; j--) {
			v = steinerNode[steinerTuple[j]];
			t = fullSteinerComponentGraph.copy(v);
			fullSteinerComponentGraph.newEdge(s, t, distance[u][v]);
		}
	}
}

template<typename T>
bool MinSteinerTreeRZLoss<T>::nextTuple(Array<int> &tuple, int n) {
	const int elements = tuple.size();
	int i = elements - 1;

	OGDF_ASSERT(tuple[i] < n);
	while (tuple[i] == i + n - elements) { // while last possible item...
		if (i == 0) { // if last possible item at first position:
			return false; // the last subset has been found
		}
		--i; // ...check previous item
	}

	// next subset
	tuple[i]++;
	while (i < elements - 1) { // reset next items
		tuple[i+1] = tuple[i] + 1;
		++i;
	}
	return true;
}

template<typename T>
void MinSteinerTreeRZLoss<T>::allPairsShortestPaths(NodeArray<NodeArray<T> > &distance,
		NodeArray<NodeArray<List<edge> > > &path)
{
	node u, v, w;
	edge e;
	double saveTime = m_startTime;

// init
	forall_nodes(u, *m_originalGraph) {
		forall_nodes(v, *m_originalGraph) {
			distance[u][v] = numeric_limits<T>::max();
			distance[v][u] = numeric_limits<T>::max();
		}
	}

	forall_edges(e, *m_originalGraph) {
		distance[e->source()][e->target()] = m_originalGraph->weight(e);
		distance[e->target()][e->source()] = m_originalGraph->weight(e);
		path[e->source()][e->target()].pushBack(e);
		path[e->target()][e->source()].pushBack(e);
	}

	forall_nodes(w, *m_originalGraph) {
		if (usedTime(saveTime) < 3600) {
			saveTime = m_startTime;
			if (!(*m_isTerminal)[w]) {
				forall_nodes(u, *m_originalGraph) {
					forall_nodes(v, *m_originalGraph) {
						if (distance[u][w] < numeric_limits<T>::max()
								&& distance[w][v] < numeric_limits<T>::max()
								&& distance[u][w] + distance[w][v] < distance[u][v]) {
							distance[u][v] = distance[u][w] + distance[w][v];
							path[u][v].clear();
							forall_listiterators(edge, it, path[u][w]) {
								path[u][v].pushBack(*it);
							}
							forall_listiterators(edge, it, path[w][v]) {
								path[u][v].pushBack(*it);
							}
						}
					}
				}
			}
		}
	}
}

template<typename T>
double MinSteinerTreeRZLoss<T>::findMaxComponent(List<FullComponent<T>*> &fullComponents, EdgeWeightedGraphCopy<T> &steinerTree, FullComponent<T> *&maxComponent, HashArray<int, edge> &save)
{
	// build save data structure
	List<node> processedNodes;
	save.clear();
	buildSaveArray(steinerTree, steinerTree.firstNode(), processedNodes, save);

	double cGain;
	double cLoss;
	double r;
	double max = 0;
	ListIterator<FullComponent<T>*> itMaxComponent = 0;
	List<ListIterator<FullComponent<T>*> > toDelete;

	// find max component
	for (ListIterator<FullComponent<T>*> it = fullComponents.begin(); it.valid(); it++) {
		m_componentsLookUps++;
		cGain = gain(*(*it), steinerTree, save) - (*(*it)).cost();
		if (cGain > 0) {
			cLoss = double((*(*it)).loss());
			r = cGain / cLoss;
			if (r > max) {
				max = r;
				maxComponent = (*it);
				itMaxComponent = it;
			}
		} else {
			toDelete.pushBack(it);
		}
	}

	forall_listiterators(ListIterator<FullComponent<T>*>, it, toDelete) {
		delete *(*it);
		fullComponents.del(*it);
	}
	toDelete.clear();
	if (itMaxComponent.valid()) {
		fullComponents.del(itMaxComponent);
	}
	return max;
}

template<typename T>
double MinSteinerTreeRZLoss<T>::gain(FullComponent<T> &fullComponent, const EdgeWeightedGraphCopy<T> &steinerTree, HashArray<int, edge> &save)
{
	List<edge> saveEdges;
	T gain = 0;
	List<node> *terminals = fullComponent.terminals();

// extract edges
	edge e;
	node u, v;
	int edgeIndex;
	for (ListConstIterator<node> it1 = (*terminals).begin(); it1.valid(); it1++) {
		for (ListConstIterator<node> it2 = (*terminals).begin(); it2.valid(); it2++) {
			if (it1 != it2) {
				u = (*it1);
				v = (*it2);
				edgeIndex = min(u->index(), v->index()) * m_originalGraph->numberOfNodes()
						+ max(u->index(), v->index());
				e = save[edgeIndex];
				if (saveEdges.search(e) == -1) {
					saveEdges.pushBack(e);
				}
			}
		}
	}

	forall_listiterators(edge, it, saveEdges) {
		gain += steinerTree.weight(*it);
	}

	return gain;
}

template<typename T>
void MinSteinerTreeRZLoss<T>::augmentCompleteTerminalGraph(FullComponent<T> &fullComponent, List<node> &newTerminal, NodeArray<bool> &isNewTerminal)
{
	List<node> *nodes = fullComponent.nodes();
	List<edge> *edges = fullComponent.edges();

	node u, v;
	forall_listiterators(node, it, *nodes) {
		if (m_completeTerminalGraph->copy(*it) == 0) {
			m_completeTerminalGraph->newNode(*it);
			newTerminal.pushBack(*it);
			isNewTerminal[*it] = true;
		}
	}

	edge e, f;
	forall_listiterators(edge, it, *edges) {
		e = *it;
		u = e->source();
		v = e->target();
		if ((f = m_completeTerminalGraph->searchEdge(m_completeTerminalGraph->copy(u), m_completeTerminalGraph->copy(v)))
				== 0) {
			m_completeTerminalGraph->newEdge(m_completeTerminalGraph->copy(u), m_completeTerminalGraph->copy(v),
					m_originalGraph->weight(e));
		} else {
			m_completeTerminalGraph->setWeight(f, min(m_originalGraph->weight(e), m_completeTerminalGraph->weight(f)));
		}
	}
}

template<typename T>
void MinSteinerTreeRZLoss<T>::contractAndAugmentComponent(EdgeWeightedGraphCopy<T> &steinerTree, FullComponent<T> &fullComponent)
{
	node u, v;
	edge e, f;
	List<edge> *lossEdges = fullComponent.lossEdges();
	HashArray<node, node> *pairedTerminal = fullComponent.pairedTerminal();
	List<edge> *edges = fullComponent.edges();

// set new edge weights for contracted component
	forall_listiterators(edge, it, *edges) {
		e = *it;
		if (lossEdges->search(e) == -1) {
			u = (*pairedTerminal)[e->source()];
			v = (*pairedTerminal)[e->target()];
			f = steinerTree.searchEdge(steinerTree.copy(u), steinerTree.copy(v));
			if (f == 0) {
				steinerTree.newEdge(steinerTree.copy(u), steinerTree.copy(v), m_originalGraph->weight(e));
			} else {
				steinerTree.setWeight(f, min(m_originalGraph->weight(e), steinerTree.weight(f)));
			}
		}
	}
}

template<typename T>
void MinSteinerTreeRZLoss<T>::buildSaveArray(EdgeWeightedGraphCopy<T> &mst, node u, List<node> &processedNodes, HashArray<int, edge> &saveArray)
{
	// traverse through tree and find heaviest edge
	Queue<node> q;
	q.append(u);
	node v, w;
	NodeArray<bool> processed(mst);
	edge e, maxE = 0;
	T maxV = -1;
	T tmp;

	forall_nodes(v, mst) {
		processed[v] = false;
	}
	processed[u] = true;

	while (!q.empty()) {
		v = q.pop();
		processedNodes.pushBack(v);
		forall_adj_edges(e, v) {
			w = e->opposite(v);
			if (!processed[w]) {
				q.append(w);
				processed[w] = true;
				tmp = mst.weight(e);
				if (tmp > maxV) {
					maxE = e;
					maxV = mst.weight(e);
				}
			}
		}
	}

	if (maxE != 0 && maxV > -1) {
		mst.hideEdge(maxE);

		List<node> processedNodes1;
		List<node> processedNodes2;

		buildSaveArray(mst, maxE->source(), processedNodes1, saveArray);
		buildSaveArray(mst, maxE->target(), processedNodes2, saveArray);

		mst.restoreEdge(maxE);

		int maxIndex, minIndex, index;
		node f, g;
		forall_listiterators(node, it1, processedNodes1) {
			f = (*it1);
			forall_listiterators(node, it2, processedNodes2) {
				g = (*it2);
				minIndex = min(mst.original(f)->index(), mst.original(g)->index());
				maxIndex = max(mst.original(f)->index(), mst.original(g)->index());
				index = minIndex * m_originalGraph->numberOfNodes() + maxIndex;
				saveArray[index] = maxE;
			}
		}
	}
}

}

#endif /* MIN_STEINER_TREE_RZ_LOSS_OGDF_H_ */
