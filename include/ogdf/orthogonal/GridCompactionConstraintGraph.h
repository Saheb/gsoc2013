/*
 * $Revision: 3425 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-22 13:49:37 +0530 (Mon, 22 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declares GridCompactionConstraintGraph.
 *
 * I.e. a representation of constraint graphs (dependency graphs)
 * used in compaction algorithms.
 *
 * \author Karsten Klein
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


#ifndef OGDF_GRID_COMP_CONSTR_GRAPH_H
#define OGDF_GRID_COMP_CONSTR_GRAPH_H


#include <ogdf/orthogonal/OrthoRep.h>
#include <ogdf/internal/orthogonal/RoutingChannel.h>
#include <ogdf/orthogonal/MinimumEdgeDistances.h>


namespace ogdf {


	// types of edges in the constraint graph
	enum ConstraintEdgeType {
		cetBasicArc,
		cetVertexSizeArc,
		cetVisibilityArc,
		cetFixToZeroArc, //can be compacted to zero length, can be fixed
		cetReducibleArc, //can be compacted to zero length
		cetMedianArc //inserted to replace some reducible in fixzerolength
	};

//---------------------------------------------------------
// GridCompactionConstraintGraph
// base class implementing common behaviour of all parameterized
// GridCompactionConstraintGraph<ATYPE>
//---------------------------------------------------------
class GridCompactionConstraintGraphBase : protected Graph
{
public:

	// output for debugging only
	void writeGML(const char *fileName) const ;
	void writeGML(ostream &os) const;
	//output edges on external face
	void writeGML(const char *fileName, NodeArray<bool> one) const ;
	void writeGML(ostream &os, NodeArray<bool> one) const;

	//return constraint arc representing input edge e in constraint graph
	edge basicArc(edge e) const {
		return m_edgeToBasicArc[e];
	}

	//***************************
	//return some edge properties
	//***************************
	//edge lies on cage border
	bool onBorder(edge e) const {return m_border[e]>0;}
	//these are subject to length fixation if length < sep
	bool fixOnBorder(edge e) const {return (m_border[e] == 2);}

	const PlanRep& getPlanRep() const {return *m_pPR;}

	edge pathToOriginal(node v) {return m_pathToEdge[v];}

protected:
	// construction
	GridCompactionConstraintGraphBase(const OrthoRep &OR,
		const PlanRep &PG,
		OrthoDir arcDir);

	// computes topological numbering on the segments of the constraint graph.
	void computeTopologicalSegmentNum(NodeArray<int> &topNum);

	// remove "arcs" from visibArcs which we already have in the constraint graph
	// (as basic arcs)
	void removeRedundantVisibArcs(SListPure<Tuple2<node,node> > &visibArcs);

	const OrthoRep *m_pOR;
	const PlanRep *m_pPR;
	OrthoDir m_arcDir;
	OrthoDir m_oppArcDir;

	NodeArray<SListPure<node> > m_path; // list of nodes contained in a segment
	NodeArray<node> m_pathNode;         // segment containing a node in PG
	EdgeArray<edge> m_edgeToBasicArc;   // basic arc representing an edge in PG

	EdgeArray<int>    m_cost;    // cost of an edge
	EdgeArray<ConstraintEdgeType> m_type;

	EdgeArray<int> m_border; //only used for cage precompaction in flowcompaction computecoords

	NodeArray<edge> m_pathToEdge; //save the (single!) edge (segment) for a pathNode
	NodeArray<edge> m_originalEdge; //save edge for the basic arcs

	// embeds constraint graph such that all sources and sinks lie in a common
	// face
	void embed();

	virtual void writeLength(ostream &os, edge e) const = 0;

private:

	void insertPathVertices(const PlanRep &PG);
	void dfsInsertPathVertex(
		node v,
		node pathVertex,
		NodeArray<bool> &visited,
		const NodeArray<node> &genOpposite);

	void insertBasicArcs(const PlanRep &PG);

	SList<node> m_sources;
	SList<node> m_sinks;
	int m_edgeCost;

};


//---------------------------------------------------------
// GridCompactionConstraintGraph
// represents a constraint graph used for compaction
//  vertices: maximally connected horiz. (resp. vert.) paths
//  basic arcs: paths connected by edges of opposite direction
//  vertex size arcs: care for minimum size of cages
//  visibility arcs: paths seeing each other
// Each edge has a (minimum) length and cost.
//---------------------------------------------------------
template<class ATYPE>
class GridCompactionConstraintGraph : public GridCompactionConstraintGraphBase
{
public:
	// construction
	GridCompactionConstraintGraph(const OrthoRep &OR,
		const PlanRep &PG,
		OrthoDir arcDir,
		ATYPE sep) :
			GridCompactionConstraintGraphBase(OR, PG, arcDir)
	{
		OGDF_ASSERT(&(const Graph &)PG == &(const Graph &)OR);
		m_sep       = sep;
		m_length.init((Graph&)*this, sep);
		m_extraNode.init((Graph&)*this, false);
		m_extraOfs.init((Graph&)*this, 0);
		m_extraRep.init((Graph&)*this, 0);

		initializeCosts();

	}//constructor

	// output for debugging only
	void writeGML(const char *fileName) const ;
	void writeGML(ostream &os) const;


	// returns underlying graph
	const Graph &getGraph() const { return (const Graph&)*this; }
	Graph &getGraph() { return (Graph&)*this; }


	const OrthoRep &getOrthoRep() const {
		return *m_pOR;
	}


	// returns list of nodes contained in segment v
	// Precodn.: v is in the constraint graph
	const SListPure<node> &nodesIn(node v) const {
		return m_path[v];
	}

	// returns the segment (path node in constraint graph) containing v
	// Precond.: v is a node in the associated planarized representation
	node pathNodeOf(node v) const {
		return m_pathNode[v];
	}

	// returns length of edge e
	// Precond.: e is an edge in the constraint graph
	ATYPE length(edge e) const {
		return m_length[e];
	}

	// returns cost of edge e
	// Precond.: e is an edge in the constraint graph
	int cost(edge e) const {
		return m_cost[e];
	}

	// returns type of edge e
	// Precond.: e is an edge in the constraint graph
	ConstraintEdgeType typeOf(edge e) const {
		return m_type[e];
	}

	//returns node status
	bool extraNode(node v) const {
		return m_extraNode[v];
	}
	//returns extraNode position, change to save mem, only need some entries
	ATYPE extraOfs(node v) const {
		return m_extraOfs[v];
	}
	//returns extraNode existing anchor representant
	node extraRep(node v) const {
		return m_extraRep[v];
	}

	// computes the total costs for coordintes given by pos, i.e.,
	// the sum of the weighted lengths of edges in the constraint graph.
	ATYPE computeTotalCosts(const NodeArray<ATYPE> &pos) const;

	// inserts arcs connecting segments which can see each other in a drawing
	// of the associated planarized representation PG which is given by
	// posDir and posOppDir.
	void insertVisibilityArcs(
		const PlanRep &PG,
		const NodeArray<ATYPE> &posDir,
		const NodeArray<ATYPE> &posOrthDir);


	//set min sep for multi edge original
	void setMinimumSeparation(const PlanRep &PG,
		const NodeArray<int> coord,
		const MinimumEdgeDistances<ATYPE> &minDist);


	// embeds constraint graph such that all sources and sinks lie in a common
	// face
	void embed() {
		GridCompactionConstraintGraphBase::embed();
	}


	// performs feasibility test for position assignment pos, i.e., checks if
	// the segment positions given by pos fulfill the constraints in the
	// compaction constraint graph
	// (for debuging only)
	bool isFeasible(const NodeArray<ATYPE> &pos);

	//returns the separation value
	ATYPE separation() const {return m_sep;}

	//return PG result for flowcompaction
	bool areMulti(edge e1, edge e2) const;


private:
	//---------------------------------------------------------
	// GridCompactionConstraintGraph::Interval
	// represents an interval on the sweep line
	//---------------------------------------------------------
	struct Interval
	{
		Interval(node v, ATYPE low, ATYPE high) {
			m_low = low;
			m_high = high;
			m_pathNode = v;
		}

		ATYPE m_low, m_high; // lower and upper bound
		node m_pathNode;     // corresponding segment

		// output operator
		friend ostream &operator<<(ostream &os,
			const Interval &interval)
		{
			os << "[" << interval.m_low << "," << interval.m_high <<
				";" << interval.m_pathNode << "]";
			return os;
		}

	};

	//---------------------------------------------------------
	// GridCompactionConstraintGraph::SegmentComparer
	// comparer class used for sorting segments by increasing position
	// (given by segPos) such that two overlapping segments come in the
	// order imposed by the embedding (given by secSort: segment which
	// comes first has secSort = 0, the other has 1)
	//---------------------------------------------------------
	class SegmentComparer
	{
	public:
		SegmentComparer(const NodeArray<ATYPE> &segPos,
			const NodeArray<int> &secSort) {
			m_pPos = &segPos;
			m_pSec = &secSort;
		}

		int compare(const node &x, const node &y) const {
			ATYPE d = (*m_pPos)[x] - (*m_pPos)[y];
			if (d < 0)
				return -1;
			else if (d > 0)
				return 1;
			else
				return (*m_pSec)[x] - (*m_pSec)[y];
		}

		OGDF_AUGMENT_COMPARER(node)
	private:
		const NodeArray<ATYPE> *m_pPos;
		const NodeArray<int>    *m_pSec;
	};

	virtual void writeLength(ostream &os, edge e) const {
		os << m_length[e];
	}

	void setBasicArcsZeroLength(const PlanRep &PG);
	void resetGenMergerLengths(const PlanRep &PG, adjEntry adjFirst);
	void setBoundaryCosts(adjEntry cornerDir,adjEntry cornerOppDir);

	bool checkSweepLine(const List<Interval> &sweepLine);

	ATYPE m_sep;

	EdgeArray<ATYPE> m_length;  // length of an edge

	NodeArray<bool> m_extraNode; //node does not represent drawing node
	//as we dont have positions, we save a drawing representant and an offset
	NodeArray<ATYPE> m_extraOfs; //offset of extra node to its rep, should change this
	NodeArray<node> m_extraRep; //existing representant of extranodes position anchor

	//**********************
	//COST SETTINGS SECTION

	// we make vertex size arcs more expensive than basic arcs in order
	// to get small cages
	// should be replaced by option/value dependent on e.g. degree
	int m_vertexArcCost;  //get small cages
	int m_doubleBendCost; //try to minimize double bends

	//factor of costs relative to generalization
	static const int c_vertexArcFactor;
	static const int c_doubleBendFactor; // = 20; //double bends cost factor*vertexArcCost

protected:
	//node v has no representation in drawing, only internal representation
	void setExtra(node v, node rep, ATYPE ofs)
		{m_extraNode[v] = true; m_extraRep[v] = rep; m_extraOfs[v] = ofs;}

	void initializeCosts()
	{
		// we make vertex size arcs more expensive than basic arcs in order
		// to get small cages; not necessary if cage size fixed in improvement
		// cost should be dependend on degree
		// Z.B. DURCH OPTION ODER WERT; DER VON DER ZAHL ADJAZENTER KANTEN ABHAENGIG IST
		// should be derived by number of edges times something
		int costGen = 1;

		m_vertexArcCost = costGen;

		m_doubleBendCost = c_doubleBendFactor*m_vertexArcCost;
	}//initializeCosts
};

//********************************
//initialization of static members
template<class ATYPE>
const int GridCompactionConstraintGraph<ATYPE>::c_vertexArcFactor = 20;
template<class ATYPE>
const int GridCompactionConstraintGraph<ATYPE>::c_doubleBendFactor = 20; //double bends cost mxxx*vertexArcCost


//************************************
//
// implementation of member functions
//
//************************************
}

#include <ogdf/planarity/PlanRep.h>

namespace ogdf {


// set cost of edges on the cage boundary to 0
template<class ATYPE>
void GridCompactionConstraintGraph<ATYPE>::setBoundaryCosts(
	adjEntry cornerDir,
	adjEntry cornerOppDir)
{
	/*
	adjEntry adj;
	for (adj = cornerDir; m_pOR->direction(adj) == m_arcDir; adj = adj->faceCycleSucc())
		m_cost[m_edgeToBasicArc[adj]] = 0;
	for (adj = cornerOppDir; m_pOR->direction(adj) == m_oppArcDir; adj = adj->faceCycleSucc())
		m_cost[m_edgeToBasicArc[adj]] = 0;
	*/
	//test fuer multi separation
	adjEntry adj;
	for (adj = cornerDir; m_pOR->direction(adj) == m_arcDir; adj = adj->faceCycleSucc())
	{
		m_cost[m_edgeToBasicArc[adj]] = 0;

		if (m_pathNode[adj->twin()->cyclicSucc()->theNode()] &&
			(m_pOR->direction(adj->faceCycleSucc()) == m_arcDir)
			)
			m_originalEdge[m_pathNode[adj->twin()->cyclicSucc()->theNode()]] =
				m_pPR->original(adj->twin()->cyclicSucc()->theEdge());

	}
	for (adj = cornerOppDir; m_pOR->direction(adj) == m_oppArcDir; adj = adj->faceCycleSucc())
	{
		m_cost[m_edgeToBasicArc[adj]] = 0;

		if (m_pathNode[adj->twin()->cyclicSucc()->theNode()])
			m_originalEdge[m_pathNode[adj->twin()->cyclicSucc()->theNode()]] =
				m_pPR->original(adj->twin()->cyclicSucc()->theEdge());
	}
}




template<class ATYPE>
void GridCompactionConstraintGraph<ATYPE>::setBasicArcsZeroLength(
	const PlanRep &PG)
{
	edge e;
	forall_edges(e,PG)
	{
		edge arc = m_edgeToBasicArc[e];
		if (arc == 0) continue;

		node v = e->source();
		node w = e->target();
		if ( ((PG.typeOf(v) == Graph::dummy) && (PG.typeOf(w) == Graph::dummy) &&
			(v->degree() == 2) && w->degree() == 2) &&
			(m_pOR->angle(e->adjSource()) == m_pOR->angle(e->adjTarget()) ) && //no uturns
			(PG.typeOf(e) != Graph::generalization)
		   )
		{
			m_length[arc] = 0;
			m_type[arc] = cetFixToZeroArc;
			//we make fixtozero arcs as expensive as possible
			m_cost[arc] = m_doubleBendCost;
		}
	}
}




// computes the total costs for coordinates given by pos, i.e.,
// the sum of the weighted lengths of edges in the constraint graph.
template<class ATYPE>
ATYPE GridCompactionConstraintGraph<ATYPE>::computeTotalCosts(
	const NodeArray<ATYPE> &pos) const
{
	ATYPE c = 0;

	edge e;
	forall_edges(e,*this)
	{
		c += cost(e) * (pos[e->target()] - pos[e->source()]);
	}

	return c;
}


//
// insertion of visibility arcs

// checks if intervals on the sweep line are in correct order
template<class ATYPE>
bool GridCompactionConstraintGraph<ATYPE>::checkSweepLine(const List<Interval> &sweepLine)
{
	if (sweepLine.empty())
		return true;

	ListConstIterator<Interval> it = sweepLine.begin();

	if((*it).m_high < (*it).m_low)
		return false;

	ATYPE x = (*it).m_low;

	for(++it; it.valid(); ++it) {
		if((*it).m_high < (*it).m_low)
			return false;
		if ((*it).m_high > x)
			return false;
		x = (*it).m_low;
	}

	return true;
}


// inserts arcs connecting segments which can see each other in a drawing
// of the associated planarized representation PG which is given by
// posDir and posOppDir.

template<class ATYPE>
void GridCompactionConstraintGraph<ATYPE>::insertVisibilityArcs(
	const PlanRep &PG,
	const NodeArray<ATYPE> &posDir,
	const NodeArray<ATYPE> &posOrthDir)
{
	OrthoDir segDir    = OrthoRep::prevDir(m_arcDir);
	OrthoDir segOppDir = OrthoRep::nextDir(m_arcDir);

	// list of visibility arcs which have to be inserted
	SListPure<Tuple2<node,node> > visibArcs;

	// lower and upper bound of segments
	NodeArray<ATYPE> low(getGraph()), lowReal(getGraph()), high(getGraph());
	NodeArray<ATYPE> segPos(getGraph(), 0); // position of segments
	NodeArray<int>   topNum(getGraph()/*,0*/); // secondary sorting criteria for segments

	// compute position and lower/upper bound of segments
	// We have to take care that segments cannot be shifted one upon the other,
	// e.g., if we have two segments (l1,h1) and (l2,h2) with l2 > h2 and
	// the distance l2-h2 is smaller than separation, the segments can see
	// each other. We do this by enlarging the lower bound of all segments
	// by separation if this lower bound is realized by a bend.
	//
	// Note: Be careful at segments attached at a vertex which are closer
	// than separation to each other. Possible solution: Remove visibility
	// arcs of segments which are connected by orthogonal segments to the
	// same vertex and bend in opposite directions.
	node v;
	forall_nodes(v,*this) {

		//special case nodes
		if (m_path[v].empty()) continue;

		SListConstIterator<node> it = m_path[v].begin();

		segPos[v] = posDir[*it];
		low[v] = high[v] = posOrthDir[*it];
		node nodeLow = *it;
		for(++it; it.valid(); ++it) {
			ATYPE x = posOrthDir[*it];
			if (x < low [v]) {
				low [v] = x;
				nodeLow = *it;
			}
			if (x > high[v]) high[v] = x;
		}
		lowReal[v] = low[v];
		Graph::NodeType typeLow = PG.typeOf(nodeLow);
		if(typeLow == Graph::dummy || typeLow == Graph::generalizationExpander) {
			/*bool subtractSep = true;
			if (nodeLow->degree() == 2) {
				adjEntry adj;
				forall_adj(adj,nodeLow) {
					if(m_pOR->direction(adj) == m_arcDir || m_pOR->direction(adj) == m_oppArcDir)
						break;
				}
				if (adj) {
					for(adjEntry adj2 = adj->faceCycleSucc();
						m_pOR->direction(adj2) == m_pOR->direction(adj);
						adj2 = adj2->twin()->faceCycleSucc()) ;
					if(posDir[adj->theNode()] == posDir[adj2->twinNode()])
					subtractSep = false;
				}
			}
			//if (subtractSep)*/
				low[v] -= m_sep;
		}
	}

	// correct "-= m_sep" ...
	OrthoDir dirMin = OrthoRep::prevDir(m_arcDir);
	OrthoDir dirMax = OrthoRep::nextDir(m_arcDir);
	bool isCaseA = (m_arcDir == odEast || m_arcDir == odSouth);
	const int angleAtMin = (m_arcDir == odEast || m_arcDir == odSouth) ? 3 : 1;
	const int angleAtMax = (m_arcDir == odEast || m_arcDir == odSouth) ? 1 : 3;
	forall_nodes(v,PG)
	{
		if(PG.expandAdj(v) == 0) continue;
		const OrthoRep::VertexInfoUML &vi = *m_pOR->cageInfo(v);

		int i = 0;
		adjEntry adj;

		for (adj = (isCaseA) ? vi.m_corner[dirMin]->faceCycleSucc()->faceCycleSucc() : vi.m_corner[dirMin]->faceCycleSucc();
			m_pOR->direction((isCaseA) ? adj : adj->faceCycleSucc()) == dirMin; //m_pOR->direction(adj) == dirMin;
			adj = adj->faceCycleSucc())
		{
			adjEntry adjCross = adj->cyclicPred();
			adjEntry adjTwin = adjCross->twin();

			adjEntry adjPred = adj->faceCyclePred();
			ATYPE delta = (isCaseA) ?
				min(abs(posOrthDir[adjPred->theNode()] - posOrthDir[adjPred->twinNode()]), m_sep) :
				min(abs(posOrthDir[adj->theNode()] - posOrthDir[adj->twinNode()]), m_sep);
			ATYPE boundary = (isCaseA) ?
				min(posOrthDir[adjPred->theNode()], posOrthDir[adjPred->twinNode()]) :
				min(posOrthDir[adj->theNode()],     posOrthDir[adj->twinNode()]);

			if (PG.typeOf(adjCross->theEdge()) == Graph::generalization)
			{
				if (isCaseA) {
					if(PG.typeOf(adjTwin->theNode()) == Graph::generalizationExpander &&
						m_pOR->angle(adjTwin) == 2)
					{
						node s1 = m_pathNode[adjTwin->theNode()];
						node s2 = m_pathNode[adjTwin->cyclicSucc()->twinNode()];
						low[s1] = lowReal[s1] - delta; // minDist.delta(v,dirMin,i);
						low[s2] = lowReal[s2] - delta; //minDist.delta(v,dirMin,i);
					}
					++i;
				} else {
					++i;
					if(PG.typeOf(adjTwin->theNode()) == Graph::generalizationExpander &&
						m_pOR->angle(adjTwin->cyclicPred()) == 2)
					{
						node s1 = m_pathNode[adjTwin->theNode()];
						node s2 = m_pathNode[adjTwin->cyclicPred()->twinNode()];
						low[s1] = lowReal[s1] - delta; //minDist.delta(v,dirMin,i);
						low[s2] = lowReal[s2] - delta; //minDist.delta(v,dirMin,i);
					}
				}
				continue;
			}

			//we save the current direction and stop if we run in opposite
			OrthoDir runDir = m_pOR->direction(adjCross);
			// if -> while
			while (PG.typeOf(adjTwin->theNode()) == Graph::dummy &&
				adjTwin->theNode()->degree() == 2 &&
				m_pOR->angle(adjTwin) == angleAtMin)
			{
				// We handle the case if an edge segment adjacent to a vertex
				// is separated by less than separation from edge segments above.
				node s = m_edgeToBasicArc[adjCross]->source();
				if(lowReal[s] != low[s])
				{
					if(low[s] >= boundary) // nothing to do?
						break;
					low[s] = boundary;
					//low[s] += m_sep - delta; //minDist.delta(v,dirMin,i);

					// If the compaction has eliminated bends, we can have the situation
					// that segment s has length 0 and the next segment s' (following the
					// edge) is at the same position (the edge arc has length 0).
					// In this case, the low-value of s' must be lowered (low[s'] := lowReal[s']
					// is approproate). The same situation can appear several times in a
					// row.
					//collect chains of segments compacted to zero length
					for( ; ; ) { //while(true/*lowReal[s] == high[s]*/) {
						do {
							adjCross = adjCross->faceCycleSucc();
						} while(m_pOR->direction(adjCross) == segDir ||
							m_pOR->direction(adjCross) == segOppDir);

						if(adjCross->theNode()->degree() != 2) // no longer a bend point?
							break;

						node sNext = m_edgeToBasicArc[adjCross]->opposite(s);

						if(segPos[sNext] != segPos[s])
							break;

						low[sNext] = lowReal[sNext];  //?
						s = sNext;
					}//while
				}//if

				adjTwin = adjCross->twin(); // update of twin for while
				//check if we have to stop
				if (runDir != m_pOR->direction(adjCross)) break;
			}//while dummy bend
		}

		i = 0;
		for (adj = (isCaseA) ? vi.m_corner[dirMax]->faceCycleSucc() : vi.m_corner[dirMax]->faceCycleSucc()->faceCycleSucc();
			m_pOR->direction((isCaseA) ? adj->faceCycleSucc() : adj) == dirMax; // m_pOR->direction(adj) == dirMax;
			adj = adj->faceCycleSucc())
		{
			adjEntry adjCross = adj->cyclicPred();
			adjEntry adjTwin = adjCross->twin();

			//ATYPE delta = -posOrthDir[adj->twinNode()] + posOrthDir[adj->theNode()];
			adjEntry adjPred = adj->faceCyclePred();
			ATYPE delta = (isCaseA) ?
				min(abs(posOrthDir[adj->twinNode()] - posOrthDir[adj->theNode()]), m_sep) :
				min(abs(posOrthDir[adjPred->theNode()] - posOrthDir[adjPred->twinNode()]), m_sep);
			ATYPE boundary = (isCaseA) ?
				min(posOrthDir[adj->twinNode()], posOrthDir[adj->theNode()]) :
				min(posOrthDir[adjPred->theNode()],     posOrthDir[adjPred->twinNode()]);

			if (PG.typeOf(adjCross->theEdge()) == Graph::generalization)
			{
				if (isCaseA) {
					++i;
					if(PG.typeOf(adjTwin->theNode()) == Graph::generalizationExpander &&
						m_pOR->angle(adjTwin->cyclicPred()) == 2)
					{
						node s1 = m_pathNode[adjTwin->theNode()];
						node s2 = m_pathNode[adjTwin->cyclicPred()->twinNode()];
						low[s1] = lowReal[s1] - delta; //minDist.delta(v,dirMax,i);
						low[s2] = lowReal[s2] - delta; //minDist.delta(v,dirMax,i);
					}
				} else {
					if(PG.typeOf(adjTwin->theNode()) == Graph::generalizationExpander &&
						m_pOR->angle(adjTwin) == 2)
					{
						node s1 = m_pathNode[adjTwin->theNode()];
						node s2 = m_pathNode[adjTwin->cyclicSucc()->twinNode()];
						low[s1] = lowReal[s1] - delta; //minDist.delta(v,dirMax,i);
						low[s2] = lowReal[s2] - delta; //minDist.delta(v,dirMax,i);
					}
					++i;
				}
				continue;
			}


			//we save the current direction and stop if we run in opposite
			OrthoDir runDir = m_pOR->direction(adjCross);
			// if -> while
			while (PG.typeOf(adjTwin->theNode()) == Graph::dummy &&
				adjTwin->theNode()->degree() == 2 &&
				m_pOR->angle(adjTwin) == angleAtMax)
			{
				node s = m_edgeToBasicArc[adjCross]->target();
				if(lowReal[s] != low[s])
				{
					if(low[s] >= boundary) // nothing to do?
						break;
					low[s] = boundary;
					//low[s] += m_sep - delta; //minDist.delta(v,dirMax,i);

					// If the compaction has eliminated bends, we can have the situation
					// that segment s has length 0 and the next segment s' (following the
					// edge) is at the same position (the edge arc has length 0).
					// In this case, the low-value of s' must be lowered (low[s'] := lowReal[s']
					// is approproate). The same situation can appear several times in a
					// row.
					//collect chains of segments compacted to zero length
					for( ; ; ) /*lowReal[s] == high[s]*/
					{
						do
						{
							adjCross = adjCross->faceCycleSucc();
						} while(m_pOR->direction(adjCross) == segDir ||
								m_pOR->direction(adjCross) == segOppDir);

						if(adjCross->theNode()->degree() != 2) // no longer a bend point?
							break;

						node sNext = m_edgeToBasicArc[adjCross]->opposite(s);

						if(segPos[sNext] != segPos[s])
							break;

						low[sNext] = lowReal[sNext];//was: low[s]
						s = sNext;
					}
				}//if lowreal != low

				adjTwin = adjCross->twin(); // update of twin for while

				//check if we have to stop
				if (runDir != m_pOR->direction(adjCross)) break;
			}//while dummy
		}
	}

	// compute topological numbering of segments as second sorting criteria
	// in order to process overlapping segments in the order imposed by the
	// embedding
	computeTopologicalSegmentNum(topNum);


	// sort segments
	SegmentComparer cmpBySegPos(segPos,topNum);
	List<node> sortedPathNodes;
	allNodes(sortedPathNodes);
	sortedPathNodes.quicksort(cmpBySegPos);

	// add segments in the order given by sortedPathNodes to sweep line
	List<Interval> sweepLine;

	ListIterator<node> itV;
	for(itV = sortedPathNodes.begin(); itV.valid(); ++itV)
	{
	 //special case nodes
		if (m_path[*itV].empty()) continue;
		OGDF_ASSERT_IF(dlExtendedChecking,checkSweepLine(sweepLine));

		node v = *itV;
		ListIterator<Interval> it;
		for(it = sweepLine.begin(); it.valid(); ++it) {
			if ((*it).m_low < high[v])
				break;
		}

		if (!it.valid()) {
			sweepLine.pushBack(Interval(v,low[v],high[v]));
			continue;
		}

		if((*it).m_high <= low[v]) {
			sweepLine.insertBefore(Interval(v,low[v],high[v]),it);
			continue;
		}

		ListIterator<Interval> itUp = it, itSucc;
		// we store if itUp will be deleted in order not to
		// access the deleted iterator later
		bool isItUpDel = ( ((*itUp).m_low >= low[v]) && ((*itUp).m_high <= high[v]) );

		for(; it.valid() && (*it).m_low >= low[v]; it = itSucc) {
			itSucc = it.succ();
			if ((*it).m_high <= high[v]) {
				visibArcs.pushBack(Tuple2<node,node>((*it).m_pathNode,v));
				sweepLine.del(it);
			}
		}

		if (it == itUp && (*it).m_high > high[v]) {
			node w = (*it).m_pathNode;
			sweepLine.insertAfter(Interval(w,(*it).m_low,low[v]),it);
			(*it).m_low = high[v];
			sweepLine.insertAfter(Interval(v,low[v],high[v]),it);
			visibArcs.pushBack(Tuple2<node,node>(w,v));

		} else {
			if ( (!isItUpDel) && itUp != it && (*itUp).m_low < high[v]) {
				(*itUp).m_low = high[v];
				visibArcs.pushBack(Tuple2<node,node>((*itUp).m_pathNode,v));
			}
			if (it.valid()) {
				if ((*it).m_high > low[v]) {
					(*it).m_high = low[v];
					visibArcs.pushBack(Tuple2<node,node>((*it).m_pathNode,v));
				}
				sweepLine.insertBefore(Interval(v,low[v],high[v]),it);

			} else {
				sweepLine.pushBack(Interval(v,low[v],high[v]));
			}
		}

	}

	// remove all arcs from visibArcs that are already in the constraint graph
	removeRedundantVisibArcs(visibArcs);

	// compute original adjacency entry corresponding to a segment
	// We use this in order to omit visibility arcs between segments which
	// belong to the same edge if they can see each other from the same side
	// of the edge; if they see each other from different sides the arc is
	// essential!
	NodeArray<adjEntry> correspEdge(getGraph(),0);
	forall_nodes(v,PG) {
		node seg = m_pathNode[v];
		adjEntry adj;
		forall_adj(adj,v) {
			if(m_pOR->direction(adj) != segDir) continue;
			edge eAdj = adj->theEdge();
			edge eOrig = PG.original(eAdj);
			if (eOrig == 0) continue;
			if (adj == eAdj->adjSource())
				correspEdge[seg] = eOrig->adjSource();
			else
				correspEdge[seg] = eOrig->adjTarget();
		}
	}

	// remove visibility arcs between ...
	SListIterator<Tuple2<node,node> > itT, itTSucc, itTPred;
	for(itT = visibArcs.begin(); itT.valid(); itT = itTSucc) {
		itTSucc = itT.succ();
		node v = (*itT).x1(), w = (*itT).x2();

		// remove arcs which connect segments belonging to the same edge
		if (correspEdge[v] && (correspEdge[v] == correspEdge[w]))
		{
			if (itTPred.valid())
				visibArcs.delSucc(itTPred);
			else
				visibArcs.popFront();
		}

		else
			itTPred = itT;
	}



	for(itT = visibArcs.begin(); itT.valid(); ++itT) {
		//***********************************CHECK if
		node v = (*itT).x1(), w = (*itT).x2();
		if (!(m_extraNode[v] || m_extraNode[w]))
			{
			//******************************CHECK if
	  node boundRepresentant1 = m_path[v].front();
	  node boundRepresentant2 = m_path[w].front();
	  node en1 = m_pPR->expandedNode(boundRepresentant1);
	  node en2 = m_pPR->expandedNode(boundRepresentant2);
	  //do not insert visibility in cages
	  if (!( ( en1 && en2 ) && ( en1 == en2) ))
		  {
			edge e = newEdge(v,w);

			//hier vielleicht multiedges abfangen: length auf max(min(sep, dists), minDist.sep)

			m_length[e] = m_sep;
			m_cost  [e] = 0;
			m_type  [e] = cetVisibilityArc;

			//writeGML("visibilityinserted.gml");
		  }//special if 2
			}//special if 1
	}

	OGDF_ASSERT_IF(dlExtendedChecking,checkSweepLine(sweepLine));

}//insertvisibilityarcs



// performs feasibility test for position assignment pos, i.e., checks if
// the segment positions given by pos fulfill the constraints in the
// compaction constraint graph
template<class ATYPE>
bool GridCompactionConstraintGraph<ATYPE>::isFeasible(
	const NodeArray<ATYPE> &pos)
{
	edge e;
	forall_edges(e, getGraph()) {
		node v = m_path[e->source()].front();
		node w = m_path[e->target()].front();;
		if (pos[w] - pos[v] < length(e)) {
			cout << "feasibility check failed for edge " << e << endl;
			cout << "  representatives: " << v << ", " << w << endl;
			cout << "  length: " << length(e) << endl;
			cout << "  actual distance: " << pos[w] - pos[v] << endl;
			cout << "  type of " << e << ": ";
			switch(m_type[e]) {
			case cetBasicArc: cout << "basic arc" << endl;
				break;
			case cetVertexSizeArc: cout << "vertex-size arc" << endl;
				break;
			case cetVisibilityArc: cout << "visibility arc" << endl;
				break;
			case cetMedianArc: cout << "median arc" << endl;
				break;
			case cetReducibleArc: cout << "reducible arc" <<endl;
				break;
			case cetFixToZeroArc: cout << "fixtozero arc" <<endl;

			}
			return false;
		}
	}

	return true;
}

//colouring for extranode
template<class ATYPE>
void GridCompactionConstraintGraph<ATYPE>::writeGML(const char *filename) const
{
	ofstream os(filename);
	writeGML(os);
}

template<class ATYPE>
void GridCompactionConstraintGraph<ATYPE>::writeGML(ostream &os) const
{
	const Graph &G = *this;

	NodeArray<int> id(getGraph());
	int nextId = 0;

	os.setf(ios::showpoint);
	os.precision(10);

	os << "Creator \"ogdf::GridCompactionConstraintGraphBase::writeGML\"\n";
	os << "graph [\n";
	os << "  directed 1\n";

	node v;
	forall_nodes(v,G) {
		os << "  node [\n";
		os << "    id " << (id[v] = nextId++) << "\n";

		if (m_extraNode[v]) {
			os << "    label \"" << "0" << "\"\n";
		} else {
			os << "    label \"" << m_pPR->expandedNode(m_path[v].front()) << "\"\n";
		}

		os << "    graphics [\n";
		os << "      x 0.0\n";
		os << "      y 0.0\n";
		os << "      w 30.0\n";
		os << "      h 30.0\n";
		if (m_extraNode[v])
			os << "      fill \"#00FFFF\"\n";
		else
			os << "      fill \"#FFFF00\"\n";
		os << "    ]\n"; // graphics

		os << "  ]\n"; // node
	}

	edge e;
	forall_edges(e,G) {
		os << "  edge [\n";
		os << "    source " << id[e->source()] << "\n";
		os << "    target " << id[e->target()] << "\n";

#if 0
		// show edge lengths as edge lables (not yet supported)
		os << "    label \"";
		writeLength(os, e);
		os << "\"\n";
#endif

		os << "    graphics [\n";

		os << "      type \"line\"\n";
		os << "      arrow \"last\"\n";
		switch(m_type[e])
		{
		case cetBasicArc: // red
			os << "      fill \"#FF0000\"\n";
			break;
		case cetVertexSizeArc: // blue
			os << "      fill \"#0000FF\"\n";
			break;
		case cetVisibilityArc: // green
			os << "      fill \"#00FF00\"\n";
			break;
		case cetReducibleArc: // red
			os << "      fill \"#FF0000\"\n";
			break;
		case cetFixToZeroArc: // violett
			os << "      fill \"#AF00FF\"\n";
			break;
		case cetMedianArc: // rose
			os << "      fill \"#FF00FF\"\n";
			break;
		}

		os << "    ]\n"; // graphics

#if 0
		os << "    LabelGraphics [\n";
		os << "      type \"text\"\n";
		os << "      fill \"#000000\"\n";
		os << "      anchor \"w\"\n";
		os << "    ]\n";
#endif

		os << "  ]\n"; // edge
	}

	os << "]\n"; // graph
}


} // end namespace ogdf


#endif
