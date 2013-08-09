/*
 * $Revision: 3379 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-06 15:57:54 +0530 (Sat, 06 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Zelikovsky's 11/6-approximation algorithm
 * 	      for the minimum steiner tree problem.
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

#ifndef MIN_STEINER_TREE_ZELIKOVSKY_OGDF_H_
#define MIN_STEINER_TREE_ZELIKOVSKY_OGDF_H_

#include <ogdf/basic/List.h>
#include <ogdf/internal/steinertree/CTree.h>
#include <ogdf/internal/steinertree/Triple.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraphCopy.h>
#include <ogdf/internal/steinertree/Voronoi.h>
#include <ogdf/module/MinSteinerTreeModule.h>

namespace ogdf {

/*!
 * \brief Helper class for pass=one option (sorts Triples descending after gain)
 */
template<typename T>
class TripleComparer {
public:
	static T compare(const Triple<T>& x1, const Triple<T>& x2) {
		return (x2.cost() - x1.cost());
	}
	OGDF_AUGMENT_STATICCOMPARER(Triple<T> )
};

/*!
 * \brief This class implements the 11/6-approximation algorithm by Zelikovsky
 * for the minimum steiner tree problem along with variants and practical improvements.
 *
 * This implementation is based on:
 *
 * (A. Zelikovsky, An 11/6-Approximation Algorithm for the Network Steiner Problem,
 * Algorithmica, volume 9, number 5, pages 463-470, Springer, 1993)
 *
 * (A. Zelikovsky, A faster approximation algorithm for the steiner problem in graphs,
 * Information Processing Letters, volume 46, number 2, pages 79-83, 1993)
 *
 * (A. Zelikovsky, Better approximation bound for the network and euclidean steiner
 * tree problems, Technical Report, 2006)
 */
template<typename T>
class MinSteinerTreeZelikovsky: public MinSteinerTreeModule<T> {
public:

	//! Choice of objective function
	enum WinCalculation {
		absolute, //!< win=gain-cost

		relative
	//!< win=gain/cost
	};

	//! Choice of triple generation
	enum TripleGeneration {
		exhaustive, //!< try all possibilities

		voronoi, //!< use voronoi regions

		none
	//!< generate triples "on the fly"
	};

	//! Switches immediate triple dropping
	enum TripleReducing {
		on, //!< removes triples as soon as their gain is known to be non positive

		off
	//!< keeps triples all the time
	};

	//! Different methods for obtaining save edges
	enum SaveCalculation {
		staticTree, //!< stores explicitly the save edge for every pair of terminals
					//!< needs O(n^2) space but has fast query times

		staticLCATree,	//!< builds a "weight tree" (save edges are inner nodes, terminals are leafes
						//!< and searches save edges via LCA calculation of two nodes

		dynamicLCATree, //!< same as staticLCATree but each time a triple has been contracted
						//!< the "weight tree" is updated dynamically rather than completely new from scratch.
						//!< Has the fastest update time

		hybrid
	//!< uses staticTree for the triple generation phase (many queries)
	//!< and dynamicLCATree during the contraction phase (few updates)
	};

	//! Enables a heuristic version
	enum Pass {
		one, //!< heuristic: evaluate all triples, sort them descending by gain, traverse sorted triples
			 //!< once, contract when possible

		multi
	//!< normal, greedy version
	};

	MinSteinerTreeZelikovsky(WinCalculation wc = absolute, TripleGeneration tg = voronoi, TripleReducing tr = on, SaveCalculation sc = hybrid, Pass pass = multi)
	 : m_winCalculation(wc),
	   m_tripleGeneration(tg),
	   m_tripleReducing(tr),
	   m_saveCalculation(sc),
	   m_pass(pass)
	{
	}

	virtual ~MinSteinerTreeZelikovsky() {
	}

	/*!
	 * \brief Builds a minimum steiner tree given a weighted graph and a list of terminals \see MinSteinerTreeModule::call
	 * @param G The weighted input graph
	 * @param terminals The list of terminal nodes
	 * @param isTerminal A bool array of terminals
	 * @param finalSteinerTree The finals steiner tree
	 * @return The objective value (sum of edge costs) of the final steiner tree
	 */
	T call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree);

	//! Sets type of gain calculation \see MinSteinerTreeZelikovsky::WinCalculation
	void winCalculation(WinCalculation wc)
	{
		m_winCalculation = wc;
	}

	//! Returns type of gain calculation currently in use \see MinSteinerTreeZelikovsky::WinCalculation
	WinCalculation winCalculation() const
	{
		return m_winCalculation;
	}

	//! Sets type of triple generation \see MinSteinerTreeZelikovsky::TripleGeneration
	void tripleGeneration(TripleGeneration tg)
	{
		m_tripleGeneration = tg;
	}

	//! Returns type of triple generation currently in use \see MinSteinerTreeZelikovsky::TripleGeneration
	TripleGeneration tripleGeneration() const
	{
		return m_tripleGeneration;
	}

	//! Sets type of triple reducing \see MinSteinerTreeZelikovsky::TripleReducing
	void tripleReducing(TripleReducing tr)
	{
		m_tripleReducing = tr;
	}

	//! Returns type of triple reducing currently in use \see MinSteinerTreeZelikovsky::TripleReducing
	TripleReducing tripleReducing() const
	{
		return m_tripleReducing;
	}

	//! Sets type of save calculation \see MinSteinerTreeZelikovsky::SaveCalculation
	void saveCalculation(SaveCalculation sv)
	{
		m_saveCalculation = sv;
	}

	//! Returns type of save calculation currently in use \see MinSteinerTreeZelikovsky::SaveCalculation
	SaveCalculation saveCalculation() const
	{
		return m_saveCalculation;
	}

	//! Sets type of pass \see MinSteinerTreeZelikovsky::Pass
	void pass(Pass p) {
		m_pass = p;
	}

	//! Returns type of pass currently in use \see MinSteinerTreeZelikovsky::Pass
	Pass pass() const
	{
		return m_pass;
	}

	//! Returns the number of generated triples
	long numberOfGeneratedTriples() const
	{
		return m_triplesGenerated;
	}

	//! Returns the number of contracted triples
	long numberOfContractedTriples() const
	{
		return m_triplesContracted;
	}

	//! Returns the number of triple lookups during execution time
	long numberOfTripleLookUps() const
	{
		return m_tripleLookUps;
	}

	//! Returns the running time for the Zelikovsky algorithm, i.e. without the eventual heuristic
	double zelikovskyRunningTime() const
	{
		return m_totalZelikovskyTime;
	}

	//! Returns the total running time
	double runningTime() const
	{
		return m_totalTime;
	}

protected:
	/*!
	 * \brief Generates a minimum terminal spanning tree based on a complete terminal graph
	 * @param steinerTree is assigned the resulting minimum terminal spanning tree
	 */
	void generateMinimumSpanningTree(EdgeWeightedGraphCopy<T> *steinerTree);

	/*!
	 * \brief Initializes data structures for the complete terminal graph
	 * @param newTerminals data structure for saving chosen center nodes as terminals
	 * @param isNewTerminal incidence vector for newTerminals
	 * @param d complete distance matrix for the original graph
	 */
	void initCompleteGraph(List<node> &newTerminals, NodeArray<bool> &isNewTerminal, Array<NodeArray<T> > &d);

	/*!
	 * \brief Builds a complete terminal graph
	 * @param d complete distance matrix for original graph
	 */
	void calculateCompleteGraph(Array<NodeArray<T> > &d);

	/*!
	 * \brief Generates triples according to the chosen option \see TripleGeneration
	 * @param triples list of generated triples
	 * @param d full distance matrix of the original graph
	 * @param save data structure for calculation save edges
	 * @return triple with highest gain
	 */
	Triple<T> generateTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d, const Save<T> &save);

	/*!
	 * \brief Generates triples according to voronoi regions
	 * @param triples list of generated triples
	 * @param d full distance matrix of the original graph
	 * @param voronoi Voronoi regions around terminals
	 * @param save data structure for calculation save edges
	 * @return triple with highest gain
	 */
	Triple<T> generateVoronoiTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d,
			const Voronoi<T> &voronoi, const Save<T> &save);

	/*!
	 * \brief Generates all possible triples
	 * @param triples list of generated triples
	 * @param d full distance matrix of the original graph
	 * @param isTerminal incidence array indicating terminal nodes
	 * @param save data structure for determining save edges
	 * @return triple with highest gain
	 */
	Triple<T> generateExhaustiveTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d, const NodeArray<bool> *isTerminal, const Save<T> &save);

	/*!
	 * \brief Contracts a triple and updates the save data structure
	 * @param triple triple to be contracted
	 * @param save save data structure
	 * @param newTerminals list storing center nodes chosen to be new terminals
	 * @param isNewTerminal incidence vector for newTerminals
	 * @param contractedTriples list of all contracted triples
	 */
	void contractTriple(const Triple<T> &triple, Save<T> &save, List<node> &newTerminals, NodeArray<bool> &isNewTerminal,
			List<Triple<T> > &contractedTriples);

	/*!
	 * \brief Contraction phase for algorithm generating triples on demand \see MinSteinerTreeZelikovsky::none
	 * @param save save data structure
	 * @param d complete distance matrix for original graph
	 * @param newTerminals list storing center nodes chosen to be new terminals
	 * @param isNewTerminal incidence vector for newTerminals
	 * @param contractedTriples list of all contracted triples
	 */
	void tripleOnDemand(Save<T> &save, const Array<NodeArray<T> > &d, List<node> &newTerminals,
			NodeArray<bool> &isNewTerminal, List<Triple<T> > &contractedTriples);

	/*!
	 * \brief Contraction phase for the original version of the algorithm \see MinSteinerTreeZelikovsky::multi
	 * @param triples list of all generated triples
	 * @param save save data structure
	 * @param newTerminals list storing center nodes chosen to be new terminals
	 * @param isNewTerminal incidence vector for newTerminals
	 * @param contractedTriples list of all contracted triples
	 */
	void multiPass(List<Triple<T> > &triples, Save<T> &save, List<node> &newTerminals, NodeArray<bool> &isNewTerminal,
			List<Triple<T> > &contractedTriples);

	/*!
	 * \brief Contraction phase for the one pass heuristic \see MinSteinerTreeZelikovsky::one
	 * @param triples list of all generated triples
	 * @param save save data structure
	 * @param newTerminals list storing center nodes chosen to be new terminals
	 * @param isNewTerminal incidence vector for newTerminals
	 * @param contractedTriples list of all contracted triples
	 */
	void onePass(List<Triple<T> > &triples, Save<T> &save, List<node> &newTerminals, NodeArray<bool> &isNewTerminal,
			List<Triple<T> > &contractedTriples);

	/*!
	 * \brief Calculate the win
	 */
	double calcWin(double gain, const Triple<T> &triple) const
	{
		switch (winCalculation()) {
		case relative:
			return gain / triple.cost();
		case absolute:
		default:
			return gain - triple.cost();
		}
	}

private:
	WinCalculation m_winCalculation; //!< Chosen option for gain calculation \see WinCalculation
	TripleGeneration m_tripleGeneration; //!< Chosen option for triple generation \see TripleGeneration
	TripleReducing m_tripleReducing; //!< Chosen option for triple reducing \see TripleReducing
	SaveCalculation m_saveCalculation; //!< Chosen option for save calculation \see SaveCalculation
	Pass m_pass; //!< Chosen option for pass \see Pass

	const EdgeWeightedGraph<T> *m_originalGraph; //!< The original edge-weighted graph
	EdgeWeightedGraphCopy<T> *m_completeTerminalGraph; //!< The complete terminal graph
	const NodeArray<bool> *m_isTerminal; //!< Incidence vector for terminal nodes
	const List<node> *m_terminals; //!< List of terminal nodes

	double m_startTime; //!< Stores the start time of the algorithm
	double m_totalZelikovskyTime; //!< Stores the time needed by the Zelikovsky part, i.e. without the eventual heuristic
	double m_totalTime; //!< Stores the overall running time
	long m_triplesGenerated; //!< Number of generated triples
	long m_triplesContracted; //!< Number of contracted triples
	long m_tripleLookUps; //!< Number of triple lookups

};

} // end namespace ogdf

// ============= Implementation =================

#include <ogdf/graphalg/Dijkstra.h>
#include <ogdf/graphalg/MinSteinerTreeTakahashi.h>
#include <ogdf/internal/steinertree/Save.h>
#include <ogdf/internal/steinertree/StaticLCATree.h>
#include <ogdf/internal/steinertree/StaticTree.h>

namespace ogdf {

template<typename T>
T MinSteinerTreeZelikovsky<T>::call(const EdgeWeightedGraph<T> &G, const List<node> &terminals, const NodeArray<bool> &isTerminal, EdgeWeightedGraphCopy<T> *&finalSteinerTree)
{
	usedTime(m_startTime);

	m_originalGraph = &G;
	m_terminals = &terminals;
	m_isTerminal = &isTerminal;

	NodeArray<bool> isNewTerminal(*m_originalGraph, false);
	List<node> newTerminals;
	List<Triple<T> > triples;
	List<Triple<T> > contractedTriples;
	Array<NodeArray<T> > d(0, m_terminals->size() - 1);

	initCompleteGraph(newTerminals, isNewTerminal, d);

	Save<T> *save;
	finalSteinerTree = new EdgeWeightedGraphCopy<T>(); // compute an initial Steiner tree first

	switch (saveCalculation()) {
	case staticTree:
		save = new StaticTree<T>(*m_completeTerminalGraph, *m_originalGraph);
		break;
	case dynamicLCATree:
		generateMinimumSpanningTree(finalSteinerTree);
		save = new CTree<T>(finalSteinerTree);
		break;
	case staticLCATree:
		save = new StaticLCATree<T>(*m_completeTerminalGraph, *m_originalGraph);
		break;
	case hybrid:
		generateMinimumSpanningTree(finalSteinerTree);
		save = new CTree<T>(finalSteinerTree);
		break;
	}

	m_triplesGenerated = 0;
	m_tripleLookUps = 0;
	m_triplesContracted = 0;

	switch (tripleGeneration()) {
	case none:
		tripleOnDemand(*save, d, newTerminals, isNewTerminal, contractedTriples);
		break;
	case exhaustive:
	case voronoi:
		Triple<T> genMaxTriple;
		if (saveCalculation() == hybrid) {
			Save<T> *saveTriple = new StaticTree<T>(*m_completeTerminalGraph, *m_originalGraph);
			genMaxTriple = generateTriples(triples, d, *saveTriple);
			delete saveTriple;
		} else {
			genMaxTriple = generateTriples(triples, d, *save);
		}
		if (genMaxTriple.cost() > 0) {
			contractTriple(genMaxTriple, *save, newTerminals, isNewTerminal, contractedTriples);
			if (pass() == multi) {
				multiPass(triples, *save, newTerminals, isNewTerminal, contractedTriples);
			} else {
				onePass(triples, *save, newTerminals, isNewTerminal, contractedTriples);
			}
		}
		break;
	}

	delete m_completeTerminalGraph;
	delete save;

	m_totalZelikovskyTime = usedTime(m_startTime);

	MinSteinerTreeTakahashi<T> mstt;
	T tmpMstWeight = 0, bestMstWeight = numeric_limits<T>::max();

	forall_listiterators(node, it, *m_terminals) {
		EdgeWeightedGraphCopy<T> *tmpSteinerTree = new EdgeWeightedGraphCopy<T>(*m_originalGraph);
		tmpMstWeight = mstt.call(*m_originalGraph, *m_terminals, *m_isTerminal, newTerminals, isNewTerminal, tmpSteinerTree, *it);
		if (tmpMstWeight < bestMstWeight) {
			bestMstWeight = tmpMstWeight;
			delete finalSteinerTree;
			finalSteinerTree = tmpSteinerTree;
		} else {
			delete tmpSteinerTree;
		}
	}

	m_totalTime = usedTime(m_startTime) + m_totalZelikovskyTime;

	return bestMstWeight;
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::initCompleteGraph(List<node> &newTerminals, NodeArray<bool> &isNewTerminal, Array<NodeArray<T> > &d)
{
	m_completeTerminalGraph = new EdgeWeightedGraphCopy<T>();
	m_completeTerminalGraph->createEmpty(*m_originalGraph);

	node u;

	forall_nodes(u, *m_originalGraph) {
		isNewTerminal[u] = false;
	}

	forall_listiterators(node, it, *m_terminals) {
		m_completeTerminalGraph->newNode(*it);
		newTerminals.pushBack(*it);
		isNewTerminal[*it] = true;
	}

	calculateCompleteGraph(d);
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::calculateCompleteGraph(Array<NodeArray<T> > &d)
{
	Dijkstra<T> sssp;
	NodeArray<edge> pi(*m_originalGraph);

	for (node u = m_completeTerminalGraph->firstNode(); u; u = u->succ()) {
		d[u->index()].init(*m_originalGraph);
		pi.init(*m_originalGraph);
		sssp.call(*m_originalGraph, m_originalGraph->edgeWeights(), m_completeTerminalGraph->original(u), pi,
				d[u->index()]);
		for (node v = u->succ(); v; v = v->succ()) {
			m_completeTerminalGraph->newEdge(u, v, d[u->index()][m_completeTerminalGraph->original(v)]);
		}
	}
}

template<typename T>
Triple<T> MinSteinerTreeZelikovsky<T>::generateTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d,
		const Save<T> &save)
{
	Triple<T> genMaxTriple;
	if (tripleGeneration() == voronoi) {
		Voronoi<T> voronoi(*m_originalGraph, *m_terminals);
		genMaxTriple = generateVoronoiTriples(triples, d, voronoi, save);
	} else
	if (tripleGeneration() == exhaustive) {
		genMaxTriple = generateExhaustiveTriples(triples, d, m_isTerminal, save);
	}
	return genMaxTriple;
}

template<typename T>
Triple<T> MinSteinerTreeZelikovsky<T>::generateVoronoiTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d,
		const Voronoi<T> &voronoi, const Save<T> &save)
{
	double maxWin = 0;
	Triple<T> maxTriple;
	maxTriple.cost(0);
	double saveTime = m_startTime;

	for (node u = m_completeTerminalGraph->firstNode(); u && usedTime(saveTime) < 3600; u = u->succ()) {
		saveTime = m_startTime;
		const NodeArray<T> *uDistance = &(d[u->index()]);
		for (node v = u->succ(); v; v = v->succ()) {
			const NodeArray<T> *vDistance = &(d[v->index()]);
			for (node w = v->succ(); w; w = w->succ()) {
				const NodeArray<T> *wDistance = &(d[w->index()]);
				const node uO = m_completeTerminalGraph->original(u);
				const node vO = m_completeTerminalGraph->original(v);
				const node wO = m_completeTerminalGraph->original(w);
				const double gain = save.gain(uO, vO, wO);

				if (tripleReducing() == off || gain > 0) {
					node center;
					T min = numeric_limits<T>::max();
					forall_listiterators(node, it, voronoi.nodesInRegion(uO)) {
						T tmp = (*uDistance)[*it] + (*vDistance)[*it] + (*wDistance)[*it];
						if (min > tmp) {
							center = *it;
							min = tmp;
						}
					}
					forall_listiterators(node, it, voronoi.nodesInRegion(vO)) {
						T tmp = (*uDistance)[*it] + (*vDistance)[*it] + (*wDistance)[*it];
						if (min > tmp) {
							center = *it;
							min = tmp;
						}
					}
					forall_listiterators(node, it, voronoi.nodesInRegion(wO)) {
						T tmp = (*uDistance)[*it] + (*vDistance)[*it] + (*wDistance)[*it];
						if (min > tmp) {
							center = *it;
							min = tmp;
						}
					}

					Triple<T> triple;
					triple.cost(min);
					double win = calcWin(gain, triple);
					if (tripleReducing() == off
					 || win > 0) {
						++m_triplesGenerated;
						triple.s0(uO);
						triple.s1(vO);
						triple.s2(wO);
						triple.z(center);
						triples.pushBack(triple);
						if (win > maxWin) {
							maxTriple = triple;
							maxWin = win;
						}
					}
				}
			}
		}
	}
	// XXX: hm, what happens if maxTriple is not set?
	return maxTriple;
}

template<typename T>
Triple<T> MinSteinerTreeZelikovsky<T>::generateExhaustiveTriples(List<Triple<T> > &triples, const Array<NodeArray<T> > &d,
		const NodeArray<bool> *isTerminal, const Save<T> &save)
{
	double maxWin = 0;
	Triple<T> maxTriple;
	maxTriple.cost(0);
	double saveTime = m_startTime;

	for (node u = m_completeTerminalGraph->firstNode(); u != 0 && usedTime(saveTime) < 3600; u = u->succ()) {
		saveTime = m_startTime;
		for (node v = u->succ(); v; v = v->succ()) {
			for (node w = v->succ(); w; w = w->succ()) {
				const node uO = m_completeTerminalGraph->original(u);
				const node vO = m_completeTerminalGraph->original(v);
				const node wO = m_completeTerminalGraph->original(w);
				Triple<T> triple;
				triple.cost(numeric_limits<T>::max());

				node x;
				forall_nodes(x, *m_originalGraph) {
					T tmp = d[u->index()][x] + d[v->index()][x] + d[w->index()][x];
					if (triple.cost() > tmp) {
						triple.z(x);
						triple.cost(tmp);
					}
				}

				const double win = calcWin(double(save.gain(uO, vO, wO)), triple);
				if (tripleReducing() == off
				 || win > 0) {
					++m_triplesGenerated;
					triple.s0(uO);
					triple.s1(vO);
					triple.s2(wO);
					triples.pushBack(triple);
					if (win > maxWin) {
						if (maxWin != 0) {
							triples.pushBack(maxTriple);
						}
						maxTriple = triple;
						maxWin = win;
					}
				}
			}
		}
	}
	if (tripleReducing() == off && maxTriple.cost() > 0) {
		triples.pushBack(maxTriple);
	}
	return maxTriple;
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::tripleOnDemand(Save<T> &save, const Array<NodeArray<T> > &d, List<node> &newTerminals,
		NodeArray<bool> &isNewTerminal, List<Triple<T> > &contractedTriples)
{
	T currentWin, win = 0;
	node u, v, s0, s1, s2;
	edge e, save1;
	T min, max, tmp, save1V, save2V;
	Triple<T> maxTriple;
	double elapsedTime;

	do {
		win = 0;
		elapsedTime = m_startTime;

		forall_nodes(u, *m_originalGraph) {
			if (!(*m_isTerminal)[u]) {
				min = numeric_limits<T>::max();
				s0 = 0;
				forall_nodes(v, *m_completeTerminalGraph) {
					if (min > (tmp = d[v->index()][u])) {
						min = tmp;
						s0 = v;
					}
				}

				s1 = 0;
				save1V = 0;
				bool maxSet = false;
				forall_nodes(v, *m_completeTerminalGraph) {
					if (v == s0) {
						continue;
					}
					tmp = save.saveWeight(m_completeTerminalGraph->original(v), m_completeTerminalGraph->original(s0))
					    - d[v->index()][u];
					if (!maxSet || max < tmp) {
						max = tmp;
						s1 = v;
						save1 = save.saveEdge(m_completeTerminalGraph->original(s0),
								m_completeTerminalGraph->original(s1));
						save1V = save.saveWeight(m_completeTerminalGraph->original(s0),
								m_completeTerminalGraph->original(s1));
						maxSet = true;
					}
				}

				currentWin = 0;
				s2 = 0;
				save2V = 0;
				forall_nodes(v, *m_completeTerminalGraph) {
					if (v != s0 && v != s1) {
						if ((save1 =
						  (e = save.saveEdge(m_completeTerminalGraph->original(s0),
						    m_completeTerminalGraph->original(v))))) {
							e = save.saveEdge(m_completeTerminalGraph->original(s1),
									m_completeTerminalGraph->original(v));
							save2V = save.saveWeight(m_completeTerminalGraph->original(s1),
									m_completeTerminalGraph->original(v));
						} else {
							save2V = save.saveWeight(m_completeTerminalGraph->original(s0),
									m_completeTerminalGraph->original(v));
						}
						tmp = save1V + save2V - d[s0->index()][u] - d[s1->index()][u] - d[v->index()][u];
						if (currentWin < tmp) {
							currentWin = tmp;
							s2 = v;
						}
					}
				}

				if (currentWin > 0 && currentWin > win) {
					win = currentWin;
					maxTriple.s0(m_completeTerminalGraph->original(s0));
					maxTriple.s1(m_completeTerminalGraph->original(s1));
					maxTriple.s2(m_completeTerminalGraph->original(s2));
					maxTriple.z(u);
					maxTriple.cost(save1V + save2V - currentWin);
				}
			}
		}

		if (win > 0) {
			contractTriple(maxTriple, save, newTerminals, isNewTerminal, contractedTriples);
		}
	} while (win > 0 && usedTime(elapsedTime) < 3600);
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::onePass(List<Triple<T> > &triples, Save<T> &save, List<node> &newTerminals,
		NodeArray<bool> &isNewTerminal, List<Triple<T> > &contractedTriples)
{
	ListIterator< Triple<T> > maxTriple;
	TripleComparer<T> tc;

	triples.quicksort(tc);
	for (ListIterator<Triple<T> > it = triples.begin(); it.valid(); ++it) {
		++m_tripleLookUps;
		double tmp = calcWin(double(save.gain((*it).s0(), (*it).s1(), (*it).s2())), *it);
		if (tmp > 0) {
			contractTriple(*it, save, newTerminals, isNewTerminal, contractedTriples);
		}
	}
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::multiPass(List<Triple<T> > &triples, Save<T> &save, List<node> &newTerminals,
		NodeArray<bool> &isNewTerminal, List<Triple<T> > &contractedTriples)
{
	double win = 0;
	ListIterator< Triple<T> > maxTriple;
	double elapsedTime;

	do {
		elapsedTime = m_startTime;
		win = 0;
		ListIterator< Triple<T> > nextIt;
		for (ListIterator< Triple<T> > it = triples.begin(); it.valid(); it = nextIt) {
			nextIt = it.succ();
			++m_tripleLookUps;
			if (tripleReducing() == on && save.alreadyContracted(*it)) {
				triples.del(it);
			} else {
				double tmp = calcWin(double(save.gain((*it).s0(), (*it).s1(), (*it).s2())), *it);
				if (win < tmp) {
					win = tmp;
					maxTriple = it;
				} else {
					if (tripleReducing() == on && tmp <= 0) {
						triples.del(it);
					}
				}
			}
		}

		if (win > 0) {
			contractTriple(*maxTriple, save, newTerminals, isNewTerminal, contractedTriples);
			if (tripleReducing() == on) {
				triples.del(maxTriple);
			}
		}
	} while (win > 0 && usedTime(elapsedTime) < 3600);
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::contractTriple(const Triple<T> &triple, Save<T> &save, List<node> &newTerminals,
		NodeArray<bool> &isNewTerminal, List<Triple<T> > &contractedTriples)
{
	m_triplesContracted++;
	contractedTriples.pushBack(triple);
	save.update(triple);
	if (!(*m_isTerminal)[triple.z()]) {
		newTerminals.pushBack(triple.z());
		isNewTerminal[triple.z()] = true;
	}
}

template<typename T>
void MinSteinerTreeZelikovsky<T>::generateMinimumSpanningTree(EdgeWeightedGraphCopy<T> *steinerTree)
{
	NodeArray<edge> stPred(*m_completeTerminalGraph);
	EdgeArray<bool> isInStTree(*m_completeTerminalGraph);
	edge e;

	steinerTree->createEmpty(*m_originalGraph);

	forall_listiterators(node, it, *m_terminals) {
		steinerTree->newNode(*it);
	}

	computeMinST(*m_completeTerminalGraph, m_completeTerminalGraph->edgeWeights(), stPred, isInStTree);

	node u;
	forall_nodes(u, *m_completeTerminalGraph) {
		if ((e = stPred[u]) != 0) {
			steinerTree->newEdge(steinerTree->copy(m_completeTerminalGraph->original(e->source())),
					steinerTree->copy(m_completeTerminalGraph->original(e->target())),
					m_completeTerminalGraph->weight(e));
		}
	}
}

}
// end namespace ogdf

#endif /* MIN_STEINER_TREE_ZELIKOVSKY_OGDF_H_ */
