/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the staticLCATree option for calculating
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

#ifndef STATICLCATREE_OGDF_H_
#define STATICLCATREE_OGDF_H_

#include <ogdf/tree/LCA.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/HeaviestEdgeInComponentTree.h>

namespace ogdf {

/*!
 * \brief This class behaves basically the same as CTree instead that the update of the weighted graph is not done dynamically here
 */
template<typename T>
class StaticLCATree: public Save<T> {
public:
	StaticLCATree(const EdgeWeightedGraphCopy<T> &completeTerminalGraph, const EdgeWeightedGraph<T> &G)
	  : Save<T>()
	  , m_tree()
	  , m_treeEdge(m_tree, NULL)
	{
		m_completeTerminalGraph = &completeTerminalGraph;
		m_numberOfNodes = G.numberOfNodes();
		m_steinerTree = new EdgeWeightedGraphCopy<T>(G);
		generateSteinerTree();
		m_cTerminals = new NodeArray<node>(*m_steinerTree, NULL);
		m_root = buildHeaviestEdgeInComponentTree(*m_steinerTree, *m_cTerminals, m_treeEdge, m_tree);
		m_lca = new LCA(m_tree, m_root);
	}
	virtual ~StaticLCATree()
	{
		delete m_steinerTree;
		delete m_lca;
	}

	/*!
	 * \brief Determines the weight of the save edge between two nodes by a LCA query
	 * @param u First node
	 * @param v First node
	 * @return Weight of the save egde between two given nodes
	 */
	virtual T saveWeight(node u, node v) const
	{
		return m_steinerTree->weight(saveEdge(u, v));
	}

	/*!
	 * \brief Determines the save edge between two nodes by a LCA query
	 * @param u First node
	 * @param v First node
	 * @return The save egde between two given nodes
	 */
	virtual edge saveEdge(node u, node v) const
	{
		return m_treeEdge[lca(m_steinerTree->copy(u), m_steinerTree->copy(v))];
	}

	/*!
	 * \brief Returns the gain (sum of the save edges) of a node triple by an LCA query
	 * @param u First triple node
	 * @param v Second triple node
	 * @param w Third triple node
	 * @return Sum of the save edges between the three nodes
	 */
	virtual T gain(node u, node v, node w) const
	{
		node save1 = lca(m_steinerTree->copy(u), m_steinerTree->copy(v));
		node save2 = lca(m_steinerTree->copy(u), m_steinerTree->copy(w));
		if (save1 == save2) {
			save2 = lca(m_steinerTree->copy(v), m_steinerTree->copy(w));
		}
		return (save1 ? m_steinerTree->weight(m_treeEdge[save1]) : 0)
		     + (save2 ? m_steinerTree->weight(m_treeEdge[save2]) : 0);
	}

	/*!
	 * \brief Updates the weighted tree data structure given a contracted triple
	 *
	 * The update first identifies the save edges and removes them. After removing links between the triple
	 * terminals and replacing them with 0-cost edges a new weighted tree is created and the LCA data
	 * structure is rebuild as well.
	 * @param t The contracted triple
	 */
	virtual void update(const Triple<T> &t)
	{
		edge save0 = saveEdge(t.s0(), t.s1());
		edge save1 = saveEdge(t.s0(), t.s2());
		edge save2 = saveEdge(t.s1(), t.s2());

		if (save0 == save1) {
			m_steinerTree->delEdge(save1);
			m_steinerTree->delEdge(save2);
		} else {
			m_steinerTree->delEdge(save0);
			m_steinerTree->delEdge(save1);
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

		m_tree.clear();
		m_cTerminals->fill(NULL);
		m_root = buildHeaviestEdgeInComponentTree(*m_steinerTree, *m_cTerminals, m_treeEdge, m_tree);
		delete m_lca;
		m_lca = new LCA(m_tree, m_root);
	}

	/*!
	 * \brief Returns the node corresponding to the LCA between two given nodes.
	 * @param u First node
	 * @param v Second node
	 * @return The LCA of u and v
	 */
	node lca(node u, node v) const
	{
		return m_lca->call((*m_cTerminals)[u], (*m_cTerminals)[v]);
	}

	/*!
	 * \brief Builds the binary weighted tree recursively top-down based on the terminal spanning tree.
	 * @param u "starting node" in the terminal spanning graph
	 * @param currentNode node for the heaviest edge in the terminal spanning graph containing the "starting node"
	 */
	void buildTree(node u, node currentNode);

	/*!
	 * \brief Builds a minimum terminal spanning tree based on a complete terminal graph.
	 */
	void generateSteinerTree();

protected:

private:
	Graph m_tree; //!< The weighted binary tree to represent the edge weight hierarchy
	NodeArray<edge> m_treeEdge; //!< Maps each inner node of m_tree to an edge in m_steinerTree
	const EdgeWeightedGraphCopy<T> *m_completeTerminalGraph;
	EdgeWeightedGraphCopy<T> *m_steinerTree;
	node m_root;
	int m_numberOfNodes;
	LCA *m_lca;
	NodeArray<node> *m_cTerminals;

};
} // end namespace ogdf

// ============= Implementation =================

namespace ogdf {

template<typename T>
void StaticLCATree<T>::generateSteinerTree()
{
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

} // end namespace ogdf

#endif /* STATICLCATREE_OGDF_H_ */
