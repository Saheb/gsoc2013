/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the staticTree option for calculating
 *        save edges in Zelikovsky's 11/6-approximation
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

#ifndef STATICTREE_OGDF_H_
#define STATICTREE_OGDF_H_

#include <ogdf/basic/HashArray.h>
#include <ogdf/internal/steinertree/Triple.h>

namespace ogdf {

/*!
 * \brief This class computes save edges recursively and stores for every node pair their save edge in a HashArray.
 */
template<typename T>
class StaticTree: public Save<T> {

public:
	/*!
	 * \brief Initializes the data structures and calculates a MST of the given complete terminal graph
	 * @param completeTerminalGraph The complete terminal graph
	 * @param G The original graph
	 */
	StaticTree(const EdgeWeightedGraphCopy<T> &completeTerminalGraph, const EdgeWeightedGraph<T> &G) :
			Save<T>() {
		m_completeTerminalGraph = &completeTerminalGraph;
		m_numberOfNodes = G.numberOfNodes();
		m_steinerTree = new EdgeWeightedGraphCopy<T>(G);
		List<node> m_processedNodes;
		m_save = new HashArray<int, edge>();
		generateSteinerTree();
		buildSaveArray(m_steinerTree->firstNode(), m_processedNodes);
	}
	virtual ~StaticTree() {
		delete m_steinerTree;
		delete m_save;
	}

	/*!
	 * \brief Determines the weight of the save edge between two nodes by a table lookup
	 * @param u First node
	 * @param v First node
	 * @return Weight of the save egde between two given nodes
	 */
	virtual T saveWeight(node u, node v) const {
		int lca = min(u->index(), v->index()) * m_numberOfNodes + max(u->index(), v->index());
		return m_steinerTree->weight((*m_save)[lca]);
	}

	/*!
	 * \brief Determines the save edge between two nodes by a table lookup
	 * @param u First node
	 * @param v First node
	 * @return The save egde between two given nodes
	 */
	virtual edge saveEdge(node u, node v) const {
		int lca = min(u->index(), v->index()) * m_numberOfNodes + max(u->index(), v->index());
		return (*m_save)[lca];
	}

	/*!
	 * \brief Returns the gain (sum of the save edges) of a node triple by an table lookup
	 * @param u First triple node
	 * @param v Second triple node
	 * @param w Third triple node
	 * @return Sum of the save edges between the three nodes
	 */
	virtual T gain(node u, node v, node w) const {
		int uIndex = u->index();
		int vIndex = v->index();
		int wIndex = w->index();
		int save0 = min(uIndex, vIndex) * m_numberOfNodes + max(uIndex, vIndex);
		int save1 = min(vIndex, wIndex) * m_numberOfNodes + max(vIndex, wIndex);
		int save2 = min(uIndex, wIndex) * m_numberOfNodes + max(uIndex, wIndex);

		T saveMax = max(m_steinerTree->weight((*m_save)[save2]),
				max(m_steinerTree->weight((*m_save)[save0]), m_steinerTree->weight((*m_save)[save1])));
		T saveMin = min(m_steinerTree->weight((*m_save)[save2]),
				min(m_steinerTree->weight((*m_save)[save0]), m_steinerTree->weight((*m_save)[save1])));
		return saveMax + saveMin;
	}

	/*!
	 * \brief Updates the weighted tree data structure given a contracted triple
	 *
	 * The update first inserts two 0-cost edges into the complete terminal graph and removes
	 * the two save edges. Afterward the lookup table is rebuild.
	 * @param t The contracted triple
	 */
	virtual void update(const Triple<T> &t) {
		updateSteinerTree(t);

		m_save->clear();
		List<node> m_processedNodes;
		buildSaveArray(m_steinerTree->firstNode(), m_processedNodes);
	}

	/*!
	 * \brief Updates the steiner tree by deleting save edges, removing all direct connections
	 * between the terminals of the contracted triple and connecting them through 0-cost edges
	 * @param t The contracted triple
	 */
	void updateSteinerTree(const Triple<T> &t) {
		int uIndex = t.s0()->index();
		int vIndex = t.s1()->index();
		int wIndex = t.s2()->index();
		int save0 = min(uIndex, vIndex) * m_numberOfNodes + max(uIndex, vIndex);
		int save1 = min(vIndex, wIndex) * m_numberOfNodes + max(vIndex, wIndex);
		int save2 = min(uIndex, wIndex) * m_numberOfNodes + max(uIndex, wIndex);

		edge e0 = (*m_save)[save0];
		edge e1 = (*m_save)[save1];
		edge e2 = (*m_save)[save2];

		if (e0 == e1) {
			m_steinerTree->delEdge(e1);
			m_steinerTree->delEdge(e2);
		} else {
			m_steinerTree->delEdge(e0);
			m_steinerTree->delEdge(e1);
		}

		edge formerConnection0 = m_steinerTree->searchEdge(m_steinerTree->copy(t.s0()), m_steinerTree->copy(t.s1()));
		edge formerConnection1 = m_steinerTree->searchEdge(m_steinerTree->copy(t.s0()), m_steinerTree->copy(t.s2()));
		edge formerConnection2 = m_steinerTree->searchEdge(m_steinerTree->copy(t.s1()), m_steinerTree->copy(t.s2()));

		if (formerConnection0 != 0) {
			m_steinerTree->delEdge(formerConnection0);
		}
		if (formerConnection1 != 0) {
			m_steinerTree->delEdge(formerConnection1);
		}
		if (formerConnection2 != 0) {
			m_steinerTree->delEdge(formerConnection2);
		}

		m_steinerTree->newEdge(m_steinerTree->copy(t.s0()), m_steinerTree->copy(t.s1()), 0);
		m_steinerTree->newEdge(m_steinerTree->copy(t.s0()), m_steinerTree->copy(t.s2()), 0);

	}

	/*!
	 * \brief Builds the lookup table for teh save edges
	 *
	 * This is done by first finding the maximum weighted edge in the graph. This edge
	 * partitions the graph and is the save edge for each pair of nodes such that not both
	 * of the nodes come from the same partition. This procedure is then applied to both
	 * partitions recursively
	 * @param u "Starting node" for traversing a partition in order to find a maximum weighted edge
	 * @param processedNodes List of seen nodes during the traversing
	 */
	void buildSaveArray(node u, List<node> &processedNodes);

	/*!
	 * \brief Builds a steiner tree based on the complete terminal graph by calculating a MST
	 */
	void generateSteinerTree();

protected:

private:
	const EdgeWeightedGraphCopy<T> *m_completeTerminalGraph; //!< The complete terminal graph
	EdgeWeightedGraphCopy<T> *m_steinerTree; //!< The current terminal spanning tree
	HashArray<int, edge> *m_save; //!< Data structure for the lookup table
	int m_numberOfNodes; //!< Number of nodes in the original graph

};
} // end namespace ogdf

// ============= Implementation =================

namespace ogdf {

template<typename T>
void StaticTree<T>::generateSteinerTree() {
	NodeArray<edge> stPred(*m_completeTerminalGraph);
	EdgeArray<bool> isInStTree(*m_completeTerminalGraph);
	node u;
	edge e;

	m_steinerTree->clear();

	forall_nodes(u, *m_completeTerminalGraph)
	{
		m_steinerTree->newNode(m_completeTerminalGraph->original(u));
	}

	computeMinST(*m_completeTerminalGraph, m_completeTerminalGraph->edgeWeights(), stPred, isInStTree);

	forall_nodes(u, *m_completeTerminalGraph)
	{
		if ((e = stPred[u]) != 0) {
			m_steinerTree->newEdge(m_steinerTree->copy(m_completeTerminalGraph->original(e->source())),
					m_steinerTree->copy(m_completeTerminalGraph->original(e->target())),
					m_completeTerminalGraph->weight(e));
		}
	}
}

template<typename T>
void StaticTree<T>::buildSaveArray(node u, List<node> &processedNodes) {
// traverse through tree and find heaviest edge
	Queue<node> q;
	q.append(u);
	node v, w;
	NodeArray<bool> processed(*m_steinerTree);
	edge e, maxE = 0;
	T maxV = -1;

	forall_nodes(v, *m_steinerTree)
	{
		processed[v] = false;
	}
	processed[u] = true;

	while (!q.empty()) {
		v = q.pop();
		processedNodes.pushBack(v);
		forall_adj_edges(e, v)
		{
			w = e->opposite(v);
			if (!processed[w]) {
				q.append(w);
				processed[w] = true;
				if (m_steinerTree->weight(e) > maxV) {
					maxE = e;
					maxV = m_steinerTree->weight(e);
				}
			}
		}
	}

	if (maxE != 0 && maxV > -1) {
		m_steinerTree->hideEdge(maxE);

		List<node> processedNodes1;
		List<node> processedNodes2;

		buildSaveArray(maxE->source(), processedNodes1);
		buildSaveArray(maxE->target(), processedNodes2);

		m_steinerTree->restoreEdge(maxE);

		int maxIndex, minIndex, index;
		node f, g;
		forall_listiterators(node, it1, processedNodes1)
		{
			f = (*it1);
			forall_listiterators(node, it2, processedNodes2)
			{
				g = (*it2);
				minIndex = min(m_steinerTree->original(f)->index(), m_steinerTree->original(g)->index());
				maxIndex = max(m_steinerTree->original(f)->index(), m_steinerTree->original(g)->index());
				index = minIndex * m_numberOfNodes + maxIndex;
				(*m_save)[index] = maxE;
			}
		}
	}
}

} // end namespace ogdf

#endif /* STATICTREE_OGDF_H_ */
