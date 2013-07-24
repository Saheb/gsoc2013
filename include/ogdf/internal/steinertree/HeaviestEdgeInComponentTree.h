/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief A function, given an edge-weighted tree, that builds an
 *  auxiliary aborescence where each arc of the input tree
 *  is a node in the aborescence. The weight of each node is at
 *  least the weight of its children.
 *  The construction algorithm takes time O(n log n).
 *
 * \author Stephan Beyer
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

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_HEAVIESTEDGEINCOMPONENTTREE_H
#define OGDF_HEAVIESTEDGEINCOMPONENTTREE_H



namespace ogdf {

/*! Builds the aborescence.
 * @return root node of the aborescence
 */
template<typename T>
node buildHeaviestEdgeInComponentTree(
		const EdgeWeightedGraphCopy<T> &inputTree, //!< the input tree
		NodeArray<node> &externalNodes, //!< the resulting mapping from input nodes to aborescence leaves (must be NULL'ed in input!)
		NodeArray<edge> &treeEdge, //!< the resulting mapping from each (inner) node of the aborescence to an edge in the input tree
		Graph &outputTree) //!< the output aborescence
{
	// sort edges by weight
	Array< Prioritized<edge, T> > sortEdges(inputTree.numberOfEdges());
	int i = 0;
	edge e;
	forall_edges(e, inputTree) {
		sortEdges[i] = Prioritized<edge,T>(e, inputTree.weight(e));
		++i;
	}
	sortEdges.quicksort();

	// insert edges into forest (which in the end makes up a tree)
	NodeArray<node *> root(outputTree);
	List<node *> garbage;
	node edgeNode;
	for (i = 0; i < inputTree.numberOfEdges(); ++i) {
		edgeNode = outputTree.newNode();
		e = sortEdges[i].item();
		treeEdge[edgeNode] = e;

		node u = e->source();
		node v = e->target();
		if (externalNodes[u]) {
			node *uRoot = root[externalNodes[u]];
			OGDF_ASSERT(uRoot);
			while (root[*uRoot] != uRoot) {
				*uRoot = *root[*uRoot];
				uRoot = root[*uRoot];
				OGDF_ASSERT(uRoot);
			}
			outputTree.newEdge(edgeNode, *uRoot);
			root[edgeNode] = uRoot;
			if (externalNodes[v]) {
				node *vRoot = root[externalNodes[v]];
				OGDF_ASSERT(vRoot);
				while (root[*vRoot] != vRoot) {
					*vRoot = *root[*vRoot];
					vRoot = root[*vRoot];
					OGDF_ASSERT(vRoot);
				}
				outputTree.newEdge(edgeNode, *vRoot);
				*vRoot = edgeNode;
			} else {
				externalNodes[v] = edgeNode;
			}
		} else {
			externalNodes[u] = edgeNode;
			if (externalNodes[v]) {
				node *vRoot = root[externalNodes[v]];
				OGDF_ASSERT(vRoot);
				while (root[*vRoot] != vRoot) {
					*vRoot = *root[*vRoot];
					vRoot = root[*vRoot];
					OGDF_ASSERT(vRoot);
				}
				outputTree.newEdge(edgeNode, *vRoot);
				root[edgeNode] = vRoot;
			} else {
				root[edgeNode] = new node;
				garbage.pushBack(root[edgeNode]);
				externalNodes[v] = edgeNode;
			}
		}
		*root[edgeNode] = edgeNode;
	}

	// garbage collect
	forall_listiterators(node *, it, garbage) {
		delete *it;
	}

	return edgeNode;
}


} // end namespace ogdf

#endif
