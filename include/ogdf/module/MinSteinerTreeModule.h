/** \file
 * \brief Interface of Minimum Steiner Tree Algorithms
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
 * Boston, MA 02110-1301, USA.treeEdges
 *
 * \see  http://www.gnu.org/copyleft/gpl.html
 ***************************************************************/

#ifdef _MSC_VER
#pragma once
#endif

#ifndef OGDF_MIN_STEINER_TREE_MODULE_H
#define OGDF_MIN_STEINER_TREE_MODULE_H

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/internal/steinertree/FullComponent.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/fileformats/GraphIO.h>
#include <sstream>

namespace ogdf {

/*!
 * \brief This class serves as an interface for various methods to calculate minimum
 * Steiner trees.
 *
 * Furthermore it supplies some helping methods.
 */
template<typename T>
class MinSteinerTreeModule {
public:
	MinSteinerTreeModule() {
	}
	virtual ~MinSteinerTreeModule() {
	}

	/*!
	 * \brief Builds a minimum steiner tree given a weighted graph and a list of terminals
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	virtual T call(const EdgeWeightedGraph<T> &G,
			const List<node> &terminals,
			const NodeArray<bool> &isTerminal,
			EdgeWeightedGraphCopy<T> *&finalSteinerTree
			) = 0;

	/*!
	 * \brief Eliminates steiner nodes with degree 1 and returns the achieved improvement
	 * @param steinerTree The given steiner tree
	 * @param isTerminal Incidence vector indicating terminal nodes
	 * @return
	 */
	T eliminate1DegreeSteinerNodes(EdgeWeightedGraphCopy<T> &steinerTree, const NodeArray<bool> &isTerminal);

	/*!
	 * \brief Writes an SVG file of a minimum steiner tree in the original graph
	 * @param G The original weighted graph
	 * @param terminals The list of terminal nodes
	 * @param steinerTree The steiner tree of the given graph
	 * @param filename The name of the output file
	 */
	void drawSVG(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename);

	/*!
	 * \brief Writes a SVG that shows only the given steiner tree
	 * @param terminals The list of terminal nodes
	 * @param steinerTree The steiner tree to be drawn
	 * @param filename The name of the output file
	 */
	void drawSteinerTreeSVG(const List<node> &terminals, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename);

	/*!
	 * \brief Checks if a given tree is acually a steiner tree
	 *
	 * - |E| = |V| - 1
	 * - all terminal nodes in the graph
	 * - all terminal nodes with degree >= 1
	 * - all steiner nodes in the graph with degree > 1
	 * @param G The original graph
	 * @param terminals The list of terminal nodes
	 * @param steinerTree The steiner tree to be checked
	 * @return true, if the given steiner tree is actually one, false otherwise
	 */
	bool isSteinerTree(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const EdgeWeightedGraphCopy<T> &steinerTree);

};

template<typename T>
bool MinSteinerTreeModule<T>::isSteinerTree(
	const EdgeWeightedGraph<T> &G,
	const List<node> &terminals,
	const EdgeWeightedGraphCopy<T> &steinerTree)
{
	node u;

	if (steinerTree.numberOfEdges() != steinerTree.numberOfNodes() - 1) {
		return false;
	}

	forall_listiterators(node, it, terminals)
	{
		u = steinerTree.copy(*it);
		if (u == 0 || u->degree() < 1) {
			return false;
		}
	}

	forall_nodes(u, steinerTree)
	{
		if (terminals.search(steinerTree.original(u)) == -1) {
			if (u->degree() <= 1) {
				return false;
			}
		}
	}

	return true;
}

template<typename T>
T MinSteinerTreeModule<T>::eliminate1DegreeSteinerNodes(EdgeWeightedGraphCopy<T> &steinerTree,
		const NodeArray<bool> &isTerminal)
{
	node u, v;
	edge e;
	T delWeights(0);
	List<node> toProcess;
	forall_nodes(u, steinerTree)
	{
		v = steinerTree.original(u);
		if (u->degree() == 1 && !isTerminal[v]) {
			toProcess.pushBack(u);
		}
	}
	forall_listiterators(node, it, toProcess)
	{
		u = *it;
		while (u->degree() == 1 && !isTerminal[steinerTree.original(u)]) {
			e = u->firstAdj()->theEdge();
			delWeights += steinerTree.weight(e);
			v = e->opposite(u);
			steinerTree.delNode(u);
			u = v;
		}
	}
	return delWeights;
}

template<typename T>
void MinSteinerTreeModule<T>::drawSteinerTreeSVG(const List<node> &terminals, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename) {
	GraphAttributes GAorig(steinerTree,
	  GraphAttributes::nodeGraphics |
	  GraphAttributes::nodeStyle |
	  GraphAttributes::nodeLabel |
	  GraphAttributes::edgeGraphics |
	  GraphAttributes::edgeStyle |
	  GraphAttributes::edgeLabel);

	node v;
	string s;

	forall_nodes(v, steinerTree)
	{
		std::stringstream out;
		GAorig.x(v) = 10;
		GAorig.y(v) = 10;
		GAorig.width(v) = 5;
		GAorig.height(v) = 5;
		if (terminals.search(steinerTree.original(v)) < 0) {
			out << "S" << steinerTree.original(v);
			GAorig.fillColor(v) = "#00ff00";
			GAorig.label(v) = out.str();
		} else {
			out << "T" << steinerTree.original(v);
			GAorig.fillColor(v) = "#ff0000";
			GAorig.label(v) = out.str();
		}
	}

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(15.0);
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	fmmm.call(GAorig);
	GraphIO::drawSVG(GAorig, filename);
}

template<typename T>
void MinSteinerTreeModule<T>::drawSVG(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const EdgeWeightedGraphCopy<T> &steinerTree, const char *filename) {
	GraphAttributes GA(G,
	  GraphAttributes::nodeGraphics |
	  GraphAttributes::nodeStyle |
	  GraphAttributes::nodeLabel |
	  GraphAttributes::edgeGraphics |
	  GraphAttributes::edgeStyle |
	  GraphAttributes::edgeLabel);

	node v;
	edge e;
	string s;

	EdgeArray<bool> isSteinerTreeEdge(G);
	forall_edges(e, G)
	{
		isSteinerTreeEdge[e] = false;
	}
	forall_edges(e, steinerTree)
	{
		isSteinerTreeEdge[steinerTree.original(e)] = true;
	}

	forall_nodes(v, G)
	{
		std::stringstream out;
		GA.x(v) = 10;
		GA.y(v) = 10;
		GA.width(v) = 5;
		GA.height(v) = 5;
		if (terminals.search(v) < 0) {
			out << "S" << v;
			GA.fillColor(v) = "#ffffff";
			GA.label(v) = out.str();
		} else {
			out << "T" << v;
			GA.fillColor(v) = "#ff0000";
			GA.label(v) = out.str();
		}
	}

	forall_edges(e, G)
	{
		if (isSteinerTreeEdge[e]) {
			GA.strokeColor(e) = "#ff0000";
		}
	}

	FMMMLayout fmmm;

	fmmm.useHighLevelOptions(true);
	fmmm.unitEdgeLength(15.0);
	fmmm.newInitialPlacement(true);
	fmmm.qualityVersusSpeed(FMMMLayout::qvsGorgeousAndEfficient);

	fmmm.call(GA);

	GraphIO::drawSVG(GA, filename);

	return;
}

} // end namespace ogdf

#endif
