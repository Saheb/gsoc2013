/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief A weighted tree as auxiliary data structure for
 * contraction based algorithms
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

#ifndef CTREE_OGDF_H_
#define CTREE_OGDF_H_

#include <ogdf/tree/LCA.h>
#include <ogdf/internal/steinertree/Save.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/internal/steinertree/HeaviestEdgeInComponentTree.h>

namespace ogdf {

/*!
 * \brief Dynamically updatable weighted Tree for determining save edges via LCA computation.
 */
template<typename T>
class CTree: public Save<T> {
public:
	/*!
	 * \brief Builds a weighted binary tree based on the given terminal spanning tree.
	 *
	 * Additionally the LCA data structure is initialized.
	 * @param terminalTree the given terminal spanning tree
	 */
	CTree(EdgeWeightedGraphCopy<T> *terminalTree)
	  : Save<T>()
	  , m_tree()
	  , m_treeEdge(m_tree, NULL)
	  , m_terminalTree(terminalTree)
	  , m_cTerminals(*m_terminalTree, NULL)
	{
		m_root = buildHeaviestEdgeInComponentTree(*m_terminalTree, m_cTerminals, m_treeEdge, m_tree);
		m_lca = new LCA(m_tree, m_root);
	}

	virtual ~CTree()
	{
		delete m_lca;
	}

	/*!
	 * \brief Returns the gain (sum of save edge costs) of the given triple, calculated by an LCA query
	 *
	 * @param u First terminal
	 * @param v Second terminal
	 * @param w Third terminal
	 * @return The gain (sum of save edge costs) of the given triple
	 */
	virtual T gain(node u, node v, node w) const
	{
		node save1 = lca(m_terminalTree->copy(u), m_terminalTree->copy(v));
		node save2 = lca(m_terminalTree->copy(u), m_terminalTree->copy(w));
		if (save1 == save2) {
			save2 = lca(m_terminalTree->copy(v), m_terminalTree->copy(w));
		}
		return weight(save1) + weight(save2);
	}

	/*!
	 * \brief Determines the weight of the save edge between two nodes by a LCA query
	 *
	 * @param u First node
	 * @param v First node
	 * @return Weight of the save edge between two given nodes
	 */
	virtual T saveWeight(node u, node v) const
	{
		return weight(saveEdge(u, v));
	}

	/*!
	 * \brief Determines the save edge between two nodes by a LCA query
	 *
	 * @param u First node
	 * @param v First node
	 * @return The save edge between two given nodes
	 */
	virtual edge saveEdge(node u, node v) const
	{
		OGDF_ASSERT(m_terminalTree->copy(u) && m_terminalTree->copy(v));
		node anc = lca(m_terminalTree->copy(u), m_terminalTree->copy(v));
		return (m_treeEdge[anc] ? m_treeEdge[anc] : NULL);
	}

	/*!
	 * \brief Updates the weighted tree data structure given a contracted triple
	 *
	 * The update of the weighted tree is performed dynamically. To achieve this,
	 * the weighted tree is traversed bottom-up, starting at the three leaves corresponding
	 * to the terminals in the triple. It takes time O(height of m_tree).
	 * @param t The contracted triple
	 */
	virtual void update(const Triple<T> &t)
	{
		// initialize the three terminal nodes
		node s0 = m_terminalTree->copy(t.s0());
		node s1 = m_terminalTree->copy(t.s1());
		node s2 = m_terminalTree->copy(t.s2());

		// determine the save edges
		node save1 = lca(s0, s1);
		node save2 = lca(s0, s2);
		if (save1 == save2) {
			save2 = lca(s1, s2);
		}

		node v0 = m_cTerminals[s0];
		node v1 = m_cTerminals[s1];
		node v2 = m_cTerminals[s2];

		int v0level = m_lca->level(v0);
		int v1level = m_lca->level(v1);
		int v2level = m_lca->level(v2);

		delete m_lca;

		node v = m_tree.newNode();
		node currentNode = m_tree.newNode();
		OGDF_ASSERT(!m_treeEdge[currentNode]);
		m_tree.newEdge(currentNode, v);
		m_cTerminals[s0] = v;
		m_cTerminals[s1] = v;
		m_cTerminals[s2] = currentNode;

		while (v0) {
			// bubble pointers such that level(v0) >= level(v1) >= level(v2)
			if (v1level < v2level) {
				swap(v1, v2);
				swap(v1level, v2level);
			}
			if (v0level < v1level) {
				swap(v0, v1);
				swap(v0level, v1level);
			}
			if (v1level < v2level) {
				swap(v1, v2);
				swap(v1level, v2level);
			}
			// bubble pointers such that weight(v0) <= weight(v1), weight(v2)
			if (weight(v1) > weight(v2)) {
				swap(v1, v2);
				swap(v1level, v2level);
			}
			if (weight(v0) > weight(v1)) {
				swap(v0, v1);
				swap(v0level, v1level);
			}
			// now v0 is the node with the least weight... if equal, with the highest level.

			if (v0 != save1
			 && v0 != save2) {
				edge e;
				forall_adj_edges(e, currentNode) {
					if (e->target() == currentNode) {
						m_tree.delEdge(e);
						break;
					}
				}
				m_tree.newEdge(v0, currentNode);
				currentNode = v0;
			} // else: nothing to do since m_tree is binary and save1/save2 are LCAs

			// set v0 to its parent or to NULL if there is no parent
			v = v0;
			v0 = NULL;
			edge e;
			forall_adj_edges(e, v) {
				if (e->target() == v) {
					v0 = e->source();
					--v0level;
					break;
				}
			}
			if (v1 == e->target()) {
				v1 = e->source();
				--v1level;
			}
			if (v2 == e->target()) {
				v2 = e->source();
				--v2level;
			}
		}
		m_root = currentNode;
		m_tree.delNode(save1);
		m_tree.delNode(save2);

		m_lca = new LCA(m_tree, m_root);
	}

	/*!
	 * \brief Checks, if components of a given triple have already been contracted
	 * @param t Triple to be checked. This is only necessary for the dynamic tree approach.
	 * @return True, if parts of the triple have been contracted before, false otherwise
	 */
	virtual bool alreadyContracted(const Triple<T> &t) const
	{
		return (saveWeight(t.s0(), t.s1()) == 0
		     || saveWeight(t.s0(), t.s2()) == 0
		     || saveWeight(t.s1(), t.s2()) == 0);
	}

protected:

	/*!
	 * \brief Returns the node in m_tree that is the LCA of two nodes
	 *
	 * @param u first node
	 * @param v second node
	 * @return the LCA of u and v
	 */
	node lca(node u, node v) const
	{
		OGDF_ASSERT(u && v);
		return m_lca->call(m_cTerminals[u], m_cTerminals[v]);
	}

	//! Returns the weight of an edge in the terminal tree or 0
	T weight(edge e) const
	{
		return (e ? m_terminalTree->weight(e) : 0);
	}

	//! Returns the associated weight of a node v in m_tree, or 0 if it is not associated.
	T weight(node v) const
	{
		return weight(m_treeEdge[v]);
	}

private:
	Graph m_tree; //!< The weighted binary tree to represent the edge weight hierarchy
	NodeArray<edge> m_treeEdge; //!< Maps each inner node of m_tree to an edge in m_terminalTree
	node m_root; //!< The root node of the weighted binary tree
	const EdgeWeightedGraphCopy<T> *m_terminalTree; //!< The underlying terminal spanning tree this weighted tree instance represents
	NodeArray<node> m_cTerminals; //!< Connects terminal nodes in the terminal spanning tree to their leafs in the weighted tree
	LCA *m_lca; //!< Data structure for calculating the LCAs
};

} // end namespace ogdf

// ============ Implementation ===============



namespace ogdf {


} // end namespace ogdf

#endif /* CTREE_OGDF_H_ */
