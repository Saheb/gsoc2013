/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of a full component used in Robins and Zelikovsky's
 * 		  loss contracting algorithm
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

#ifndef FULLCOMPONENT_OGDF_H_
#define FULLCOMPONENT_OGDF_H_

#include <ogdf/basic/HashArray.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Queue.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/basic/extended_graph_alg.h>

namespace ogdf {

/*!
 * \brief This class represents a full component for the Loss Contraction Algorithm
 */
template<typename T>
class FullComponent {
public:

	/*!
	 * \brief Initializes a full component representing by the given graph
	 * @param wGC Graph representing the full component
	 */
	FullComponent(EdgeWeightedGraphCopy<T> *wGC) {
		m_pairedTerminal = new HashArray<node, node>();
		m_terminals = new List<node>();
		m_lossEdges = new List<edge>();
		m_nodes = new List<node>();
		m_edges = new List<edge>();
		m_loss = 0;
		m_cost = 0;

		edge e;
		forall_edges(e, *wGC)
		{
			m_cost += wGC->weight(e);
			m_edges->pushBack(wGC->original(e));
		}

		node u;
		forall_nodes(u, *wGC)
		{
			m_nodes->pushBack(wGC->original(u));
		}
	}
	;

	virtual ~FullComponent() {
		delete m_pairedTerminal;
		delete m_terminals;
		delete m_lossEdges;
		delete m_nodes;
		delete m_edges;
	}
	;

	/*!
	 * \brief Returns the list of terminals in this full component
	 * @return List of terminals
	 */
	List<node>* terminals() {
		return m_terminals;
	}
	;

	/*!
	 * \brief Returns the list of edges that form the loss forest
	 * @return List of loss edges
	 */
	List<edge>* lossEdges() {
		return m_lossEdges;
	}
	;

	/*!
	 * \brief Returns the relation between terminals and steiner nodes according to the loss calculation
	 *
	 * A terminal and a steiner node are linked if the terminal is the first one on the shortest loss path
	 * starting from the steiner node.
	 * @return HashArray representing the loss relationship between terminals and steiner nodes
	 */
	HashArray<node, node>* pairedTerminal() {
		return m_pairedTerminal;
	}
	;

	/*!
	 * \brief Returns the loss value of this full component
	 * @return The loss value
	 */
	T loss() const {
		return m_loss;
	}
	;

	/*!
	 * \brief Returns the sum of edge costs of this full component
	 * @return The cost value
	 */
	T cost() const {
		return m_cost;
	}
	;

	/*!
	 * \brief Returns a list of all nodes in this full component
	 * @return List of nodes
	 */
	List<node>* nodes() {
		return m_nodes;
	}
	;

	/*!
	 * \brief Returns a list of all edges in this full component
	 * @return List of edges
	 */
	List<edge>* edges() {
		return m_edges;
	}
	;

	/*!
	 * \brief Calculates the loss, both edge set and value, of this full component
	 * @param wgc Full component
	 */
	void calculateLoss(EdgeWeightedGraphCopy<T> *wgc) {
		node s, u, v;
		edge e;

		NodeArray<bool> isTerminal(*wgc);
		forall_nodes(u, *wgc)
		{
			isTerminal[u] = false;
			(*m_pairedTerminal)[wgc->original(u)] = wgc->original(u);
		}

		List<edge> zeroEdges;
		s = m_terminals->front();
		forall_listiterators(node, it, *m_terminals)
		{
			isTerminal[wgc->copy(*it)] = true;
			if (s != (*it)) {
				u = wgc->copy(s);
				v = wgc->copy(*it);
				e = wgc->newEdge(u, v, 0);
				zeroEdges.pushBack(e);
			}
		}

		NodeArray<edge> pred(*wgc);
		EdgeArray<bool> isTreeEdge(*wgc);
		computeMinST(wgc->copy(s), *wgc, wgc->edgeWeights(), pred, isTreeEdge);

		m_loss = 0;
		forall_nodes(u, *wgc)
		{
			e = pred[u];
			if (e != 0 && wgc->weight(e) > 0) {
				m_loss += wgc->weight(e);
				m_lossEdges->pushBack(wgc->original(e));

				findTerminal(u, isTerminal, pred, wgc);
			}
		}

		forall_listiterators(edge, it, zeroEdges)
		{
			wgc->delEdge(*it);
		}
	}
	;

	/*!
	 * \brief Starting from an steiner node find the nearest terminal along a shortest path
	 * @param u Steiner node to start from
	 * @param isTerminal Terminal node incidence vector
	 * @param pred The shortest path predecessor data structure
	 * @param wgc The graph structure representing this full component
	 * @return first terminal on a shortest path starting from a steiner node
	 */
	node findTerminal(node u, NodeArray<bool> &isTerminal, NodeArray<edge> &pred, EdgeWeightedGraphCopy<T> *wgc) {
		if (!isTerminal[u] && (*m_pairedTerminal)[wgc->original(u)] == wgc->original(u) && pred[u] != 0) {
			(*m_pairedTerminal)[wgc->original(u)] = findTerminal(pred[u]->opposite(u), isTerminal, pred, wgc);
		}
		return (*m_pairedTerminal)[wgc->original(u)];
	}

	/*!
	 * \brief Adds a new terminal node to this component.
	 * @param s
	 * @param t
	 * @param distance
	 * @param path
	 * @param wgc
	 */
	void appendTerminalToSteinerTree(node s, node t, NodeArray<NodeArray<T> > &distance,
			NodeArray<NodeArray<List<edge> > > &path, EdgeWeightedGraphCopy<T> *wgc) {
		node u, v;
		edge e, f;
		T edgeCost;
		forall_listiterators(edge, it, path[s][t])
		{
			e = *it;
			if (wgc->chain(e).size() == 0) {

				if ((u = wgc->copy(e->source())) == 0) {
					u = wgc->newNode(e->source());
					m_nodes->pushBack(e->source());
				}
				if ((v = wgc->copy(e->target())) == 0) {
					v = wgc->newNode(e->target());
					m_nodes->pushBack(e->target());
				}
				edgeCost = distance[wgc->original(u)][wgc->original(v)];
				f = wgc->newEdge(u, v, edgeCost);
				wgc->setEdge(e, f);
				m_cost += edgeCost;
				m_edges->pushBack(e);

			}
		}
	}

	/*!
	 * \brief Adds a terminal to the list of terminals
	 * @param t terminal to add
	 */
	void addTerminalToList(node t) {
		m_terminals->pushBack(t);
	}

private:
	List<node> *m_terminals; //!< List of terminal nodes in this full component
	T m_loss; //!< Value of loss for this component
	T m_cost; //!< Value of cost for this component
	List<edge> *m_lossEdges; //!< List of edges that build the loss connection to the terminals
	HashArray<node, node> *m_pairedTerminal; //!< Indicates which steiner node is connected to which terminal through the loss edges
	List<node> *m_nodes; //!< List of nodes of this full component
	List<edge> *m_edges; //!< List of edges of this full component
};

}

#endif /* FULLCOMPONENT_OGDF_H_ */
