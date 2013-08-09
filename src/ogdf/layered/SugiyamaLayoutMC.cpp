/*
 * $Revision: 3400 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-16 12:40:39 +0530 (Tue, 16 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of Sugiyama algorithm for multicore architecturres
 *
 * \author Carsten Gutwenger
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


#include <ogdf/layered/SugiyamaLayoutMC.h>
#include <ogdf/layered/LongestPathRanking.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/packing/TileToRowsCCPacker.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>


namespace ogdf {


//---------------------------------------------------------
// SugiyamaLayoutMC
// Sugiyama drawing algorithm for hierarchical graphs
//---------------------------------------------------------

SugiyamaLayoutMC::SugiyamaLayoutMC()
{
	m_ranking.set(new LongestPathRanking);
	m_layout .set(new FastHierarchyLayout);
	m_packer .set(new TileToRowsCCPacker);

	m_fails = 4;
	m_runs = 15;
	m_transpose = true;
	m_permuteFirst = false;

	m_arrangeCCs = true;
	m_minDistCC = 20;
	m_pageRatio = 1.0;

	m_evenOdd = false;
	m_sortAfterLayers = 10;

	m_maxLevelSize = -1;
	m_numLevels = -1;
	m_timeReduceCrossings = 0.0;
}


void SugiyamaLayoutMC::call(GraphAttributes &AG)
{
	doCall(AG,false);
}


void SugiyamaLayoutMC::call(GraphAttributes &AG, NodeArray<int> &rank)
{
	doCall(AG,false,rank);
}


void SugiyamaLayoutMC::doCall(GraphAttributes &AG, bool umlCall)
{
	NodeArray<int> rank;
	doCall(AG, umlCall, rank);
}


void SugiyamaLayoutMC::doCall(GraphAttributes &AG, bool umlCall, NodeArray<int> &rank)
{
	const Graph &G = AG.constGraph();
	if (G.numberOfNodes() == 0)
		return;

	// compute connected component of G
	NodeArray<int> component(G);
	m_numCC = connectedComponents(G,component);

	const bool optimizeHorizEdges = (umlCall || rank.valid());
	if(!rank.valid())
		m_ranking.get().call(AG.constGraph(),rank);

	if(m_arrangeCCs) {
		// intialize the array of lists of nodes contained in a CC
		Array<List<node> > nodesInCC(m_numCC);

		node v;
		forall_nodes(v,G)
			nodesInCC[component[v]].pushBack(v);

		Hierarchy H;
		H.createEmpty(G);
		const GraphCopy &GC = H;

		EdgeArray<edge> auxCopy(G);
		Array<DPoint> boundingBox(m_numCC);
		Array<DPoint> offset1(m_numCC);
		NodeArray<bool> mark(GC);

		m_numLevels = m_maxLevelSize = 0;

		int totalCrossings = 0;
		for(int i = 0; i < m_numCC; ++i)
		{
			// adjust ranks in cc to start with 0
			int minRank = numeric_limits<int>::max();
			ListConstIterator<node> it;
			for(it = nodesInCC[i].begin(); it.valid(); ++it)
				if(rank[*it] < minRank)
					minRank = rank[*it];

			if(minRank != 0) {
				for(it = nodesInCC[i].begin(); it.valid(); ++it)
					rank[*it] -= minRank;
			}

			H.initByNodes(nodesInCC[i],auxCopy,rank);
			HierarchyLevels levels(H);
			//cout << "GC: " << GC.numberOfNodes() << ", " << GC.numberOfEdges() << endl;

			if(m_evenOdd) reduceCrossingsEvenOdd(levels); else reduceCrossings(levels);
			totalCrossings += m_nCrossings;

			m_layout.get().call(levels,AG);

			double
				minX =  numeric_limits<double>::max(),
				maxX = -numeric_limits<double>::max(),
				minY =  numeric_limits<double>::max(),
				maxY = -numeric_limits<double>::max();

			node vCopy;
			forall_nodes(vCopy,GC)
			{
				mark[vCopy] = false;
				node v = GC.original(vCopy);
				if(v == 0) continue;

				if(AG.x(v)-AG.width (v)/2 < minX) minX = AG.x(v)-AG.width(v) /2;
				if(AG.x(v)+AG.width (v)/2 > maxX) maxX = AG.x(v)+AG.width(v) /2;
				if(AG.y(v)-AG.height(v)/2 < minY) minY = AG.y(v)-AG.height(v)/2;
				if(AG.y(v)+AG.height(v)/2 > maxY) maxY = AG.y(v)+AG.height(v)/2;
			}

			if(optimizeHorizEdges)
			{
				for(int i = 0; i < levels.size(); ++i) {
					const Level &L = levels[i];
					for(int j = 0; j < L.size(); ++j) {
						node v = L[j];
						if(!GC.isDummy(v)) continue;
						edge e = GC.original(v->firstAdj()->theEdge());
						if(e == 0) continue;
						node src = GC.copy(e->source());
						node tgt = GC.copy(e->target());

						if(H.rank(src) == H.rank(tgt)) {
							int minPos = levels.pos(src), maxPos = levels.pos(tgt);
							if(minPos > maxPos) std::swap(minPos,maxPos);

							bool straight = true;
							const Level &L_e = levels[H.rank(src)];
							for(int p = minPos+1; p < maxPos; ++p) {
								if(!H.isLongEdgeDummy(L_e[p]) && mark[L_e[p]] == false) {
									straight = false; break;
								}
							}
							if(straight) {
								AG.bends(e).clear();
								mark[v] = true;
							}
						}
					}
				}
			}

			edge eCopy;
			forall_edges(eCopy,GC)
			{
				edge e = GC.original(eCopy);
				if(e == 0 || eCopy != GC.chain(e).front()) continue;

				const DPolyline &dpl = AG.bends(e);
				ListConstIterator<DPoint> it;
				for(it = dpl.begin(); it.valid(); ++it)
				{
					if((*it).m_x < minX) minX = (*it).m_x;
					if((*it).m_x > maxX) maxX = (*it).m_x;
					if((*it).m_y < minY) minY = (*it).m_y;
					if((*it).m_y > maxY) maxY = (*it).m_y;
				}
			}

			minX -= m_minDistCC;
			minY -= m_minDistCC;

			boundingBox[i] = DPoint(maxX - minX, maxY - minY);
			offset1    [i] = DPoint(minX,minY);

			m_numLevels = max(m_numLevels, levels.size());
			for(int i = 0; i <= levels.high(); i++) {
				Level &l = levels[i];
				m_maxLevelSize = max(m_maxLevelSize, l.size());
			}
		}

		m_nCrossings = totalCrossings;

		// call packer
		Array<DPoint> offset(m_numCC);
		m_packer.get().call(boundingBox,offset,m_pageRatio);

		// The arrangement is given by offset to the origin of the coordinate
		// system. We still have to shift each node and edge by the offset
		// of its connected component.

		for(int i = 0; i < m_numCC; ++i)
		{
			const List<node> &nodes = nodesInCC[i];

			const double dx = offset[i].m_x - offset1[i].m_x;
			const double dy = offset[i].m_y - offset1[i].m_y;

			// iterate over all nodes in ith CC
			ListConstIterator<node> it;
			for(it = nodes.begin(); it.valid(); ++it)
			{
				node v = *it;

				AG.x(v) += dx;
				AG.y(v) += dy;

				edge e;
				forall_adj_edges(e,v)
				{
					if(e->isSelfLoop() || e->source() != v) continue;

					DPolyline &dpl = AG.bends(e);
					ListIterator<DPoint> it;
					for(it = dpl.begin(); it.valid(); ++it)
					{
						(*it).m_x += dx;
						(*it).m_y += dy;
					}
				}
			}
		}

	} else {
		int minRank = numeric_limits<int>::max();
		node v;
		forall_nodes(v,G)
			if(rank[v] < minRank)
				minRank = rank[v];

		if(minRank != 0) {
			forall_nodes(v,G)
				rank[v] -= minRank;
		}

		Hierarchy H(G,rank);
		HierarchyLevels levels(H);
		const GraphCopy &GC = H;

		m_compGC.init(GC);
		forall_nodes(v,GC) {
			node vOrig = GC.original(v);
			if(vOrig == 0)
				vOrig = GC.original(v->firstAdj()->theEdge())->source();

			m_compGC[v] = component[vOrig];
		}

		if(m_evenOdd) reduceCrossingsEvenOdd(levels); else reduceCrossings(levels);
		m_compGC.init();

		m_layout.get().call(levels,AG);

		if(optimizeHorizEdges)
		{
			NodeArray<bool> mark(GC,false);
			for(int i = 0; i < levels.size(); ++i) {
				const Level &L = levels[i];
				for(int j = 0; j < L.size(); ++j) {
					node v = L[j];
					if(!GC.isDummy(v)) continue;
					edge e = GC.original(v->firstAdj()->theEdge());
					if(e == 0) continue;
					node src = GC.copy(e->source());
					node tgt = GC.copy(e->target());

					if(H.rank(src) == H.rank(tgt)) {
						int minPos = levels.pos(src), maxPos = levels.pos(tgt);
						if(minPos > maxPos) std::swap(minPos,maxPos);

						bool straight = true;
						const Level &L_e = levels[H.rank(src)];
						for(int p = minPos+1; p < maxPos; ++p) {
							if(!H.isLongEdgeDummy(L_e[p]) && mark[L_e[p]] == false) {
								straight = false; break;
							}
						}
						if(straight) {
							AG.bends(e).clear();
							mark[v] = true;
						}
					}
				}
			}
		}

		m_numLevels = levels.size();
		m_maxLevelSize = 0;
		for(int i = 0; i <= levels.high(); i++) {
			Level &l = levels[i];
			if (l.size() > m_maxLevelSize)
				m_maxLevelSize = l.size();
		}

	}

}


bool SugiyamaLayoutMC::transposeLevel(int i, HierarchyLevels &levels)
{
	bool improved = false;

	if (m_levelChanged[i] || m_levelChanged[i-1] || m_levelChanged[i+1]) {
		Level &L = levels[i];

		for (int j = 0; j < L.high(); j++) {
			if (levels.transpose(L[j])) improved = true;
		}
	}

	if (improved) levels.buildAdjNodes(i);
	return (m_levelChanged[i] = improved);
}


void SugiyamaLayoutMC::doTranspose(HierarchyLevels &levels)
{
	m_levelChanged.fill(true);

	bool improved;
	do {
		improved = false;

		for (int i = 0; i <= levels.high(); ++i)
			improved |= transposeLevel(i,levels);
	} while (improved);
}


void SugiyamaLayoutMC::doTransposeRev(HierarchyLevels &levels)
{
	m_levelChanged.fill(true);

	bool improved;
	do {
		improved = false;

		for (int i = levels.high(); i >= 0 ; --i)
			improved |= transposeLevel(i,levels);
	} while (improved);
}


//const int sugiSortAfterLayers = 5;

void SugiyamaLayoutMC::barycenter(Level &L, bool doSorting)
{
	const double sugiEps = 0.1;

	// special case: just one node on layer
	if(L.high() == 0) {
		m_weight[L[0]] = 0;
		return;
	}

	double M = -numeric_limits<double>::infinity();
	for (int j = 0; j <= L.high(); ++j) {
		node v = L[j];
		const Array<node> &adjNodes = L.adjNodes(v);

		double sumWeights = 0;
		for(int k = 0; k <= adjNodes.high(); ++k)
			sumWeights += m_weight[adjNodes[k]];

		m_weight[v] = (adjNodes.high() < 0) ? M :
			sumWeights / double(adjNodes.size());
	}

	// handle "isolated" nodes
	for(int j = 0; j <= L.high(); ++j)
	{
		if(m_weight[L[j]] <= M) {
			int b = j+1;
			while(b <= L.high() && m_weight[L[b]] <= M)
				++b;

			if(j == 0) {
				double wb = (b <= L.high()) ? m_weight[L[b]]-1.0 : double(b-j-1);
				for(int i = j-1; i >= 0; --i, wb -= 1.0)
					m_weight[L[i]] = wb;

			} else if(b > L.high()) {
				double wa = m_weight[L[j-1]] + 1.0;
				for(int i = j; i < b; ++i, wa += 1.0)
					m_weight[L[i]] = wa;

			} else {
				double wa = m_weight[L[j-1]];
				double wb = m_weight[L[b]];
				if(wa > wb) swap(wa,wb);

				double delta = (wb-wa) / (b-j+1);
				double offset = delta;
				for(int i = j; i < b; ++i, offset += delta)
					m_weight[L[i]] = wa + offset;
			}
			j = b;
		}
	}

	if(doSorting) {
		L.sortByWeightOnly(m_weight);
		for(int j = 0; j <= L.high(); ++j)
			m_weight[L[j]] = j;

	} else {
		for(int j = 0; j <= L.high(); ++j) {
			m_weight[L[j]] += sugiEps*j;
		}
		//L.sortByWeightOnly(m_weight);
	}
}


int SugiyamaLayoutMC::traverseTopDown(HierarchyLevels &levels)
{
	levels.direction(HierarchyLevels::downward);

	int t = 0;
	//if(levels[0].high() == 0) t++;

	Level &Lt = levels[t];
	for(int j = 0; j <= Lt.high(); ++j)
		m_weight[Lt[j]] = j;

	int counter = 0;
	for (int i = t+1; i <= levels.high(); ++i)
	{
		//Level &L = levels[i];
		barycenter(levels[i], (++counter % m_sortAfterLayers) == 0);
	}

	counter = 0;
	for (int i = t+1; i <= levels.high(); ++i) {
		if((++counter % m_sortAfterLayers) != 0)
			levels[i].sortByWeightOnly(m_weight);
	//	cout << "Level " << i << ": ";
	//	Level &L = levels[i];
	//	for (int j = 0; j <= L.high(); ++j) {
	//		cout << L[j] << " [" << m_weight[L[j]] << "]  ";
	//	}
	//	cout << endl;
	}

	if (m_transpose)
		doTranspose(levels);

	if(m_arrangeCCs == false)
		levels.separateCCs(m_numCC, m_compGC);

	return levels.calculateCrossings();
}


int SugiyamaLayoutMC::traverseBottomUp(HierarchyLevels &levels)
{
	levels.direction(HierarchyLevels::upward);

	Level &Lb = levels[levels.high()];
	for(int j = 0; j <= Lb.high(); ++j)
		m_weight[Lb[j]] = j;

	int counter = 0;
	for (int i = levels.high()-1; i >= 0; --i) {
		Level &L = levels[i];
		barycenter(L, (++counter % m_sortAfterLayers) == 0);

	}

	counter = 0;
	for (int i = levels.high()-1; i >= 0; --i) {
		if((++counter % m_sortAfterLayers) != 0)
			levels[i].sortByWeightOnly(m_weight);
	//	cout << "Level " << i << ": ";
	//	Level &L = levels[i];
	//	for (int j = 0; j <= L.high(); ++j) {
	//		cout << L[j] << " [" << m_weight[L[j]] << "]  ";
	//	}
	//	cout << endl;
	}

	if (m_transpose)
		doTransposeRev(levels);

	if(m_arrangeCCs == false)
		levels.separateCCs(m_numCC, m_compGC);

	return levels.calculateCrossings();
}


void SugiyamaLayoutMC::reduceCrossings(HierarchyLevels &levels)
{
	const GraphCopy &GC = levels.hierarchy();
	GraphAttributes GA(GC,  GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::nodeLabel);

	__int64 t;
	System::usedRealTime(t);

	//node x;
	//forall_nodes(x,GC) {
	//	GA.label(x) = to_string(x->index());
	//}
	//GraphIO::writeGML(GA, "gc.gml");

	if(m_permuteFirst)
		levels.permute();

	int nCrossingsOld, nCrossingsNew;
	NodeArray<int> bestPos;

	m_nCrossings = nCrossingsOld = levels.calculateCrossings();
	levels.storePos(bestPos);

	if (m_nCrossings == 0) {
		t = System::usedRealTime(t);
		m_timeReduceCrossings = double(t) / 1000;
		return;
	}

	m_weight.init(GC);

	if (m_transpose) {
		m_levelChanged.init(-1,levels.size());
		m_levelChanged[-1] = m_levelChanged[levels.size()] = false;
	}


	for (int i = 1; ;i++ ) {
		int nFails = m_fails+1;

		do {
			nCrossingsNew = traverseTopDown(levels);
			//cout << "cr[t->d] = " << nCrossingsNew << endl;

			if (nCrossingsNew < nCrossingsOld) {
				if (nCrossingsNew < m_nCrossings) {
					levels.storePos(bestPos);

					if ((m_nCrossings = nCrossingsNew) == 0)
						break;
				}
				nCrossingsOld = nCrossingsNew;
				nFails = m_fails+1;
			} else
				--nFails;

			nCrossingsNew = traverseBottomUp(levels);
			//cout << "cr[b->u] = " << nCrossingsNew << endl;

			if (nCrossingsNew < nCrossingsOld) {
				if (nCrossingsNew < m_nCrossings) {
					levels.storePos(bestPos);

					if ((m_nCrossings = nCrossingsNew) == 0)
						break;
				}
				nCrossingsOld = nCrossingsNew;
				nFails = m_fails+1;
			} else
				--nFails;

		} while (nFails > 0);

		if (m_nCrossings == 0 || i >= m_runs)
			break;

		levels.permute();

		nCrossingsOld = levels.calculateCrossings();
		if (nCrossingsOld < m_nCrossings) {
			levels.storePos(bestPos);

			m_nCrossings = nCrossingsOld;
		}
	}

	levels.restorePos(bestPos);

	//for(int i = 0; i < levels.size(); i++)
	//	cout << i << ": " << levels[i].size() << endl;

	m_weight.init();
	m_levelChanged.init();

	t = System::usedRealTime(t);
	m_timeReduceCrossings = double(t) / 1000;
}


int SugiyamaLayoutMC::traverseEvenOdd(HierarchyLevels &levels)
{
	levels.direction(HierarchyLevels::downward);

	// process odd layers
	//cout << "process odd" << endl;
	for(int i = 1; i <= levels.high(); i += 2)
	{
		//cout << "Level " << i << endl;
		Level &L = levels[i];

		for(int j = 0; j <= L.high(); ++j) {
			node v = L[j];
			double sum_in  = 0;
			double sum_out = 0;

			adjEntry adj;
			forall_adj(adj,v) {
				node x = adj->twinNode();
				if(adj->theEdge()->source() == v)
					sum_out += levels.pos(x);
				else
					sum_in  += levels.pos(x);
			}

			m_weight[v] = sum_in / v->indeg() + sum_out/v->outdeg();
		}

		L.sort(m_weight);
	}

	// process even layers
	//cout << "process even" << endl;
	for(int i = 0; i <= levels.high(); i += 2)
	{
		//cout << "Level " << i << endl;
		Level &L = levels[i];

		for(int j = 0; j <= L.high(); ++j) {
			node v = L[j];
			double sum_in  = 0;
			double sum_out = 0;

			adjEntry adj;
			forall_adj(adj,v) {
				node x = adj->twinNode();
				if(adj->theEdge()->source() == v)
					sum_out += levels.pos(x);
				else
					sum_in  += levels.pos(x);
			}

			m_weight[v] = sum_in / v->indeg() + sum_out/v->outdeg();
		}

		L.sort(m_weight);
	}

	if(m_arrangeCCs == false)
		levels.separateCCs(m_numCC, m_compGC);

	return levels.calculateCrossings();
}


void SugiyamaLayoutMC::reduceCrossingsEvenOdd(HierarchyLevels &levels)
{
	int nCrossingsOld, nCrossingsNew;
	NodeArray<int> bestPos;

	m_nCrossings = nCrossingsOld = levels.calculateCrossings();
	levels.storePos(bestPos);

	if (m_nCrossings == 0) return;
	cout << nCrossingsOld << endl;

	m_weight.init(levels.hierarchy());

	if (m_transpose) {
		m_levelChanged.init(-1,levels.size());
		m_levelChanged[-1] = m_levelChanged[levels.size()] = false;
	}

	for (int i = 1; ; i++ ) {
		int nFails = m_fails+1;

		do {
			nCrossingsNew = traverseEvenOdd(levels);
			cout << nCrossingsNew << endl;

			if (nCrossingsNew < nCrossingsOld) {
				if (nCrossingsNew < m_nCrossings) {
					levels.storePos(bestPos);

					if ((m_nCrossings = nCrossingsNew) == 0)
						break;
				}
				nCrossingsOld = nCrossingsNew;
				nFails = m_fails+1;
			} else
				--nFails;

		} while (nFails > 0);

		if (m_nCrossings == 0 || i >= m_runs)
			break;

		cout << "============" << endl;
		levels.permute();

		nCrossingsOld = levels.calculateCrossings();
		if (nCrossingsOld < m_nCrossings) {
			levels.storePos(bestPos);

			m_nCrossings = nCrossingsOld;
		}
	}

	levels.restorePos(bestPos);

	m_weight.init();
	m_levelChanged.init();
}


} // end namespace ogdf
