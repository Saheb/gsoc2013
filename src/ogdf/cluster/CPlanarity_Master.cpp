/*
 * $Revision: 3504 $
 *
 * last checkin:
 *   $Author: beyer $
 *   $Date: 2013-05-16 18:19:39 +0530 (Thu, 16 May 2013) $
 ***************************************************************/

/** \file
 * \brief Implementation of the CPlanarityMaster class for the Branch-Cut-Price algorithm
 * for c-planarity testing.
 *
 * This class is managing the optimization.
 * Variables and initial constraints are generated and pools are initialized.
 * Since variables correspond to the edges of a complete graph, node pairs
 * are used mostly instead of edges.
 *
 * \author Markus Chimani, Karsten Klein
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

#include <ogdf/basic/basic.h>

#ifdef USE_ABACUS

#include <ogdf/internal/cluster/CPlanarity_Master.h>
#include <ogdf/internal/cluster/CPlanarity_Sub.h>
#include <ogdf/internal/cluster/Cluster_ChunkConnection.h>
#include <ogdf/internal/cluster/Cluster_MaxPlanarEdges.h>
#include <ogdf/planarity/BoyerMyrvold.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/cluster/ClusterArray.h>
#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/fileformats/GraphIO.h>

using namespace ogdf;
using namespace abacus;


#ifdef OGDF_DEBUG
void CPlanarityMaster::printGraph(const Graph &G) {
	edge e;
	int i=0;
	Logger::slout() << "The Given Graph" << endl;
	forall_edges(e,G) {
		Logger::slout() << "Edge " << i++ << ": (" << e->source()->index() << "," << e->target()->index() << ") " << endl;
	}
}
#endif


CPlanarityMaster::CPlanarityMaster(
	const ClusterGraph &C,
		int heuristicLevel,
		int heuristicRuns,
		double heuristicOEdgeBound,
		int heuristicNPermLists,
		int kuratowskiIterations,
		int subdivisions,
		int kSupportGraphs,
		double kHigh,
		double kLow,
		bool perturbation,
		double branchingGap,
		const char *time,
		bool dopricing,
		int numAddVariables,
		double strongConstraintViolation,
		double strongVariableViolation) :
		CP_MasterBase(C, heuristicLevel, heuristicRuns, heuristicOEdgeBound, heuristicNPermLists,
		kuratowskiIterations, subdivisions, kSupportGraphs, kHigh, kLow, perturbation, branchingGap,
		time, dopricing, numAddVariables, strongConstraintViolation, strongVariableViolation), m_ca(0), m_ssg(0)
	//Master("CPlanarity", true, false, OptSense::Min) //no pricing so far
{

	// Reference to the given ClusterGraph and the underlying Graph.
	//m_C = &C;
	//m_G = &(C.getGraph());
	// Create a copy of the graph as we need to modify it
	//m_solutionGraph = new GraphCopy(*m_G);

	// Define the maximum number of variables needed.
	// The actual number needed may be much smaller, so there
	// is room for improvement...
	//ToDo: Just count how many vars are added

	//Max number of edges and extension edges, to be compared
	//with actual number of variables used.
	int nComplete = (m_G->numberOfNodes()*(m_G->numberOfNodes()-1)) / 2;
	m_nMaxVars = nComplete-m_G->numberOfEdges();

	// Initialize the node array to keep track of created variables
	m_varCreated.init(*m_G);
	node v;
	forall_nodes(v, *m_G)
	{
		m_varCreated[v].init(*m_G, false);
	}

	// Setting parameters
	m_nKuratowskiIterations = kuratowskiIterations;
	m_nSubdivisions = subdivisions;
	m_nKuratowskiSupportGraphs = kSupportGraphs;
	m_heuristicLevel = heuristicLevel;
	m_nHeuristicRuns = heuristicRuns;
	m_usePerturbation = perturbation;
	m_kuratowskiBoundHigh = kHigh;
	m_kuratowskiBoundLow = kLow;
	m_branchingGap = branchingGap;
	m_heuristicFractionalBound = heuristicOEdgeBound;
	m_nHeuristicPermutationLists = heuristicNPermLists;
	m_mpHeuristic = true;

	// Further settings
	m_nCConsAdded = 0;
	m_nKConsAdded = 0;
	m_solvesLP = 0;
	m_varsInit = 0;
	m_varsAdded = 0;
	m_varsPotential = 0;
	m_varsMax = 0;
	m_varsCut = 0;
	m_varsKura = 0;
	m_varsPrice = 0;
	m_varsBranch = 0;
	m_activeRepairs = 0;
	m_repairStat.init(100);
	m_shrink = true;

	m_cNodes.init(C);
	cluster c;
	forall_clusters(c, C)
	{
		c->getClusterNodes(m_cNodes[c]);
	}
}


CPlanarityMaster::~CPlanarityMaster() {
	//delete m_maxCpuTime; done in base class
	//delete m_solutionGraph; done in base class
	if (m_ssg != 0) delete m_ssg;
}


Sub *CPlanarityMaster::firstSub() {
	return new CPlanaritySub(this);
}


// Replaces current m_solutionGraph by new GraphCopy based on \a connection list
void CPlanarityMaster::updateBestSubGraph(List<nodePair> &connection) {

	// Creates a new GraphCopy \a m_solutionGraph and deletes all edges
	// TODO: Extend GraphCopySimple to be usable here: Allow
	// edge deletion and add pure node initialization.
	// Is the solutiongraph used during computation anyhow?
	// Otherwise only store the lists
	delete m_solutionGraph;
	m_solutionGraph = new GraphCopy(*m_G);

	// Delete all edges that have been stored previously in edge lists
	m_connectionOneEdges.clear();

	ListConstIterator<nodePair> cit = connection.begin();
	node cv, cw;
	while (cit.valid()) {

		// Add all new connection edges to \a m_solutionGraph
		cv = m_solutionGraph->copy((*cit).v1);
		cw = m_solutionGraph->copy((*cit).v2);
		m_solutionGraph->newEdge(cv,cw);

		m_connectionOneEdges.pushBack(*cit);
		cit++;
	}

#ifdef OGDF_DEBUG
	GraphIO::writeGML(*m_solutionGraph, "UpdateSolutionGraph.gml");
	//Just for special debugging purposes:
	if (true) {
		ClusterArray<cluster> ca(*m_C);
		Graph GG;
		NodeArray<node> na(*m_G);
		ClusterGraph CG(*m_C,GG, ca, na);

		cit = connection.begin();

		List<edge> le;

		while (cit.valid()) {

			// Add all new connection edges to \a m_solutionGraph
			cv = na[(*cit).v1];
			cw = na[(*cit).v2];
			edge e = GG.newEdge(cv,cw);
			le.pushBack(e);

			cit++;
		}

		ClusterGraphAttributes CGA(CG, GraphAttributes::edgeType | GraphAttributes::nodeType |
				GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::edgeStyle);
		ListConstIterator<edge> it = le.begin();
		while (it.valid())
		{
			cout << (*it)->graphOf() << "\n";
			cout << &GG << "\n";
			CGA.strokeColor(*it) = "#FF0000";
			it++;
		}
		GraphIO::writeGML(CGA, "PlanarExtension.gml");
	}
#endif
}


void CPlanarityMaster::getConnectionOptimalSolutionEdges(List<nodePair> &edges) const {
	edges.clear();
	ListConstIterator<nodePair> it;
	for (it=m_connectionOneEdges.begin(); it.valid(); ++it) {
		edges.pushBack(*it);
	}
}


//todo: is called only once, but could be sped up the same way as the co-conn check
//Returns number of edges to be added to achieve cluster connectivity for \a c
//Could use code similar to the one in ClusterAnalysis for speedup
double CPlanarityMaster::clusterConnection(cluster c, GraphCopy &gc) {
	// For better performance, a node array is used to indicate which nodes are contained
	// in the currently considered cluster.
	NodeArray<bool> vInC(gc,false);
	double connectNum = 0.0; //Minimum number of connection edges
	// First check, if the current cluster \a c is a leaf cluster.
	// If so, compute the number of edges that have at least to be added
	// to make the cluster induced graph connected.
	if (c->cCount()==0) { 	//cluster \a c is a leaf cluster
		GraphCopy *inducedC = new GraphCopy((const Graph&)gc);
		node v,w;
		//List<node> clusterNodes;
		//c->getClusterNodes(clusterNodes); // \a clusterNodes now contains all (original) nodes of cluster \a c.
		ListConstIterator<node> it;
		for (it=m_cNodes[c].begin(); it.valid(); ++it) {
			vInC[gc.copy(*it)] = true;
		}

		// Delete all nodes from \a inducedC that do not belong to the cluster,
		// in order to obtain the cluster induced graph.
		v = inducedC->firstNode();
		while (v!=0)  {
			w = v->succ();
			if (!vInC[inducedC->original(v)]) inducedC->delNode(v);
			v = w;
		}

		// Determine number of connected components of cluster induced graph.
		//Todo: check could be skipped
		if (!isConnected(*inducedC)) {

			NodeArray<int> conC(*inducedC);
			//at least #connected components - 1 edges have to be added.
			connectNum = connectedComponents(*inducedC,conC) - 1;
		}
		delete inducedC;
	// Cluster \a c is an "inner" cluster. Process all child clusters first.
	} else {	//c->cCount is != 0, process all child clusters first
		ListConstIterator<cluster> cit;
		for (cit=c->cBegin(); cit.valid(); ++cit) {
			connectNum += clusterConnection(*cit,gc);
		}

		// Create cluster induced graph.
		GraphCopy *inducedC = new GraphCopy((const Graph&)gc);
		node v,w;
		//List<node> clusterNodes;
		//c->getClusterNodes(clusterNodes); //\a clusterNodes now contains all (original) nodes of cluster \a c.
		ListConstIterator<node> it;
		for (it=m_cNodes[c].begin(); it.valid(); ++it) {
			vInC[gc.copy(*it)] = true;
		}
		v = inducedC->firstNode();
		while (v!=0)  {
			w = v->succ();
			if (!vInC[inducedC->original(v)]) inducedC->delNode(v);
			v = w;
		}

		// Now collapse each child cluster to one node and determine #connected components of \a inducedC.
		List<node> oChildClusterNodes;
		List<node> cChildClusterNodes;
		for (cit=c->cBegin(); cit.valid(); ++cit) {
			getClusterNodes((*cit),oChildClusterNodes);
			ListConstIterator<node> it;
			node copy;
			// Compute corresponding nodes of graph \a inducedC.
			for (it=oChildClusterNodes.begin(); it.valid(); ++it) {
				copy = inducedC->copy(gc.copy(*it));
				cChildClusterNodes.pushBack(copy);
			}
			inducedC->collaps(cChildClusterNodes);
			oChildClusterNodes.clear();
			cChildClusterNodes.clear();
		}
		// Now, check \a inducedC for connectivity.
		if (!isConnected(*inducedC)) {

			NodeArray<int> conC(*inducedC);
			//at least #connected components - 1 edges have to added.
			connectNum += connectedComponents(*inducedC,conC) - 1;
		}
		delete inducedC;
	}
	return connectNum;
}//clusterConnection

double CPlanarityMaster::heuristicInitialLowerBound()
{
	//Heuristics?
	/*
		 * Heuristic can be improved by checking, how many additional C-edges have to be added at least.
		 * A first simple approach is the following:
		 * Since the Graph has to be completely connected in the end, all chunks have to be connected.
		 * Thus the numbers of chunks minus 1 summed up over all clusters is a trivial lower bound.

		* We perform a bottom-up search through the cluster-tree, each time checking the cluster
		 * induced Graph for connectivity. If the Graph is not connected, the number of chunks -1 is added to
		 * a counter. For "inner" clusters we have to collapse all child clusters to one node,
		 * in order to obtain a correct result.
		 */

	GraphCopy gcc(*m_G);
	cluster c = m_C->rootCluster();

	double cconn = clusterConnection(c, gcc);
	// Todo Adding a lower bound for the co connectivity is not that simple, as you cannot
	// just sum up the values for each cluster (imagine a circle of clusters).
	// We could at least add the maximum of complement-bags over all clusters.
#ifdef OGDF_ASSERT
	Logger::slout() << "Initial lower bound" << cconn << "\n";
#endif

	return cconn;
}//heuristicInitialLowerBound

double CPlanarityMaster::heuristicInitialUpperBound() {

	//Todo: Nice heuristic
	//Can we just use the number of edges needed
	//to make both the clusters and their complement connected independently?
	return 3*m_G->numberOfNodes() - 6 - m_G->numberOfEdges();//m_nMaxVars;
}//heuristicInitialUpperBound

void CPlanarityMaster::nodeDistances(node u, NodeArray<NodeArray<int> > &dist) {

	// Computing the graphtheoretical distances of node u
	NodeArray<bool> visited(*m_G);
	List<node> queue;
	visited.fill(false);
	visited[u] = true;
	int nodesVisited = 1;
	adjEntry adj;
	node v;
	forall_adj(adj,u) {
		visited[adj->twinNode()] = true;
		nodesVisited++;
		dist[u][adj->twinNode()] += 1;
		queue.pushBack(adj->twinNode());
	}
	while (!queue.empty() || nodesVisited!=m_G->numberOfNodes()) {
		v = queue.front();
		queue.popFront();
		forall_adj(adj,v) {
			if (!visited[adj->twinNode()]) {
				visited[adj->twinNode()] = true;
				nodesVisited++;
				dist[u][adj->twinNode()] += (dist[u][v]+1);
				queue.pushBack(adj->twinNode());
			}
		}
	}
}

bool CPlanarityMaster::goodVar(node a, node b) {
	if (m_varCreated[a][b] || m_varCreated[b][a]) return false;
	return true; //add all variables even if they are bad
	//we may need to add such variables as we could have made the bad
	//decision before and need to detect this...

/*//or do a simple planarity check in advance
	Logger::slout() << "Good Var? " << a << "->" << b << ": ";
	GraphCopy GC(*m_G);
	edge e = GC.newEdge(GC.copy(a),GC.copy(b));
	BoyerMyrvold bm;
	int ret =  bm.planarDestructive(GC);
	Logger::slout() << ret << "\n";
	return ret;
	*/
}

// Create variables for complete connectivity - any solution allowed
void CPlanarityMaster::createCompConnVars(List<CPlanarEdgeVar*>& initVars)
{
	//initVars.clear(); We don't care if there are already vars added
	//We create a copy of the clustergraph and insert connections to
	//make the clusters connected. Afterwards, we check if the complements
	//need to be made connected and add corresponding edges
	Graph G;
	ClusterArray<cluster> oriCluster(*m_C);
	NodeArray<node> copyNode(*m_G);
	ClusterGraph cg(*m_C, G, oriCluster, copyNode);

	NodeArray<node> oriNode(G);
	node w;
	forall_nodes(w, *m_G)
	{
		oriNode[copyNode[w]] = w;
	}
	//*******************************************************
	// First, we add the edges to make the clusters connected
	List<edge> addedEdges;
	//Todo: Use post order traversal to do this without recursion
	//forall_postOrderClusters(c,C)
	makeCConnected(cg, G, addedEdges, true); //use simple cc method
	ListIterator<edge> it = addedEdges.begin();
	while (it.valid())
	{
		node u = (*it)->source();
		node v = (*it)->target();
		initVars.pushBack( createVariable(oriNode[u],oriNode[v]) );
#ifdef OGDF_DEBUG
		cout <<"Added var "<< oriNode[u]->index()<<":" <<oriNode[v]->index()<<"\n";
#endif
		it++;
	}
	//*******************************************************
	// Now complement connnectivity (TODO)
	// TODO: do this optionally, experimentally compare performance
}

//create the variables at start of optimization
void CPlanarityMaster::createInitialVariables(List<CPlanarEdgeVar*>& initVars) {
	// In any case, add a fixed edge in size 2 clusters (can be deleted anyway
	// and helps to handle intermediate and final solutions as graphs).
	cluster c;
	forall_clusters(c,*m_C) {
		if (c->cCount() == 0 && c->nCount() == 2)
		{
			ListIterator<node> it = c->nBegin();
			node v = (*it);
			it++;
			if (!m_G->searchEdge(*it,v))
			{
				initVars.pushBack( createVariable(v, *it, 1.0));
			}
		}

	}
	// In case of pricing, create an initial variable pool allowing
	// connectivity
	if (pricing())
		createCompConnVars(initVars);
}

//! Create variables for external cluster connections in case we search
//! only in the bag-reduced search space. Uses satchel information.
//! Cluster-based adding of external connections may be less efficient
//! than directly adding the necessary variables, as we may check the
//! same node pair several times, and traverse large parts of the graph
//! several times, but there is no clear characterization
//! of "necessary" without running through each cluster's satchels so far.
void CPlanarityMaster::addExternalConnections(cluster c, List<CPlanarEdgeVar*>& connectVars) {
	// We want to traverse the satchels of c, that is minimal bags in the complement of c
	// that are connected to a vertex (outeractive) in c. Satchels may be connected
	// to several outeractive vertices and vice versa.

	// First we create a NodeArray to mark what we have seen already
	NodeArray<int> mark(*m_G, 0); //value 0 means not touched so far
	//List<node> cnodes;
	//c->getClusterNodes(cnodes);
	ListConstIterator<node> it = getClusterNodes(c).begin();//cnodes.begin();
	while (it.valid())
	{
		mark[*it] = 1; // value 1 means part of the cluster, must be skipped
		it++;
	}
	// We also mark the clusters on the path from c to the root,
	// as these are the clusters that we don't add completely when
	// touching them during a satchel detection run.
	// During the same run, we compute the cluster depth of c (root has depth 0).
	Array<bool> notrpath(m_C->maxClusterIndex()+1);
	int cdepth = 0;
	for (int i = 0; i <= m_C->maxClusterIndex(); i++) notrpath[i] = true;
	cluster rc = c->parent(); //we never touch c anyway
	while (rc) {
		notrpath[rc->index()] = false;
		rc = rc->parent();
		cdepth++;
	}
	// We also store all vertices in c's complement that are qualified (wrt activity level)
	List<node> qualifiedComplement;
	node qc;
	forall_nodes(qc,*m_G)
	{
		if ((mark[qc] == 0) && (m_ca->minIOALevel(qc)<= cdepth))
			qualifiedComplement.pushBack(qc);
	}
#ifdef OGDF_DEBUG
	Logger::slout()<<"Qualified complement size: "<<qualifiedComplement.size()<<"\n";
#endif
	// We also store for each cluster the information if its content has
	// already completely been added to the satchel, this helps to
	// shortcut the addition when touching a higher level cluster (otherwise
	// we would need to process all successor clusters' vertices multiple times).
	// This means however that we have a slight overhead for initialization.
	Array<bool> unprocessed(m_C->maxClusterIndex()+1);
	forall_clusters(rc, *m_C)
	{
		unprocessed[rc->index()] = true;
	}
	// In addition we keep the info if a vertex is qualified part of a satchel
	// in order to separate the pool of connection vertices. We use a single
	// NodeArray for all satchels, as it is faster to set and delete the satchel
	// vertex entries each time.
	NodeArray<bool> inActiveSatchel(*m_G, false);
#ifdef OGDF_DEBUG
				Logger::slout()<<"*Searching satchels for next cluster with depth "<<cdepth<<"*\n";
#endif
	// Now we start at each of the outeractive vertices and check if
	// we find a connection to an unmarked vertex. This then is the
	// beginning of a new satchel to explore.
	// The ClusterAnalysis object we have created has outer activity lists:
	// Note that these do include all outeractive vertices, not just direct children.
	List<node>& oaNodes = m_ca->oaNodes(c);
	it = oaNodes.begin();
	OGDF_ASSERT(it.valid()); //May never be empty here
	while (it.valid())
	{
		edge e;
		//Check for edges that lead to external vertices.
		forall_adj_edges(e,*it)
		{
			node w = e->opposite(*it);
			if (mark[w] == 0)
			{   // A vertex that we haven't seen yet. Traverse its CC in complement(c)
				// We maintain two lists (could reduce this to a single one...), one
				// for the whole satchel, the other for the process queue.
#ifdef OGDF_DEBUG
				Logger::slout()<<"New satchel start\n";
#endif
				List<node> queue;   // Stores vertices scheduled for processing.
				List<node> satchel; // Stores satchel vertices qualified for connections.
				mark[w] = 2;
				queue.pushBack(w);
				// We only add vertices to the satchel list that are ia/oa with level <= c level
				// The first vertex is always qualified as it is inneractive wrt c.
				OGDF_ASSERT(m_ca->minIOALevel(w)<= cdepth);
				if (m_ca->minIOALevel(w)<= cdepth)
				{
					satchel.pushBack(w);
					inActiveSatchel[w] = true;
				}
				while (!queue.empty()) {
					w = queue.popFrontRet();
					forall_adj_edges(e,w)
					{
						node u = e->opposite(w);
						if (mark[u] == 0) {
							// A new member of our current satchel run
							mark[u] = 2;
							queue.pushBack(u);
							if (m_ca->minIOALevel(w)<= cdepth)
							{
								satchel.pushBack(u);
								inActiveSatchel[u] = true;
							}
							// we also need to check if this vertex is
							// part of a cluster outside the path from c to root
							// because then we need to add the rest of the cluster, too.
							// forall cluster vertices mark and add to queue
							rc = m_C->clusterOf(u);
							if (notrpath[rc->index()] && unprocessed[rc->index()])
							{
								// Add all cluster vertices, check if child clusters have
								// been added already.
								List<cluster> cqueue;
								cqueue.pushBack(rc);
								while (!cqueue.empty())
								{
									cluster cc = cqueue.popFrontRet();
									unprocessed[cc->index()] = false;
									// Run through all cluster vertices
									ListConstIterator<node> itn = cc->nBegin();
									while (itn.valid())
									{
										node vc = *itn;
										if (mark[vc] == 0)
										{
											mark[vc] = 2;
											queue.pushBack(vc);
										}
										itn++;
									}
									// Run through all children
									ListConstIterator<cluster> itc = cc->cBegin();
									while (itc.valid())
									{
										OGDF_ASSERT(notrpath[(*itc)->index()]);
										if (unprocessed[(*itc)->index()]) {
											cqueue.pushBack(*itc);
										}
										itc++;
									}
								}//while
							}// if cluster qualified
						}// if new satchel vertex
					}//forall adj
				}//while connected component
#ifdef OGDF_DEBUG
				Logger::slout()<<"Found a satchel CC with size "<< satchel.size()<<"\n";
#endif
				// Now we create the connections between all vertices in qualifiedComplement
				// not in the satchel and the qualified vertices in the satchel.
				// Afterwards we reset the status of the satchel vertices in inActiveSatchel to false;
				// (Could also set inActiveSatchel here using a single run over the list)
				ListConstIterator<node> itqc = qualifiedComplement.begin();
				while (itqc.valid())
				{
					node qcn = *itqc;
					if (!inActiveSatchel[qcn])
					{
						//Add vars if necessary to all qualified satchel vertices
						ListConstIterator<node> its = satchel.begin();
						while (its.valid())
						{
							node sn = *its;
							if (goodVar(qcn,sn)) {
								if (!m_G->searchEdge(qcn,sn)){
									if(pricing())
										m_inactiveVariables.pushBack( nodePair(qcn,sn) );
									else
										connectVars.pushBack( createVariable(qcn,sn) );
								}
								++m_varsMax;
							}
							its++;
						}//while satchel vertices
					}
					itqc++;
				}

				//Reset satchel status
				ListConstIterator<node> itsn = satchel.begin();
				while (itsn.valid()){
					inActiveSatchel[*(itsn)] = false;
					itsn++;
				}
			}// if new satchel
		}
		it++;
	}

}//addexternalconnections

//! Create variables for inner cluster connections in case we search
//! only in the bag-reduced search space.
void CPlanarityMaster::addInnerConnections(cluster c, List<CPlanarEdgeVar*>& connectVars) {
	OGDF_ASSERT(m_ca);
	// In case there is only a single outgoing edge, there is nothing
	// to do (each bag has at least one outgoing connection).
	if (m_ca->outerActive(c)<2) {
		OGDF_ASSERT(m_ca->numberOfBags(c) < 2);
		return;
	}
	// Even if there are more outgoing edges all might emerge from
	// the same bag (we separate the tests just for testing purposes).
	if (m_ca->numberOfBags(c) < 2) return;
	// Now we have at least 2 bags, both with outgoing edges.
	// We add all variables that correspond to edges between outeractive
	// vertices among different bags in the cluster.
	// This information is stored in our ClusterAnalysis object.
	// We run through the stored vertices and check if connections are needed.
	ListIterator<node> it = m_ca->oaNodes(c).begin();
	while (it.valid()) {
		int bagindex = m_ca->bagIndex((*it),c);
		ListConstIterator<node> it2 = it.succ();
		while (it2.valid()){
			//Check if vertices are from different bags
			if (bagindex != m_ca->bagIndex((*it2),c))
			{
				 if(!m_G->searchEdge((*it),(*it2))) {
					if(goodVar((*it),(*it2))) {
						if(pricing())
							m_inactiveVariables.pushBack(nodePair((*it),(*it2)));
						else
							connectVars.pushBack( createVariable((*it),(*it2)) );
					}
					++m_varsMax;
				}
			}
			it2++;
		}
		it++;
	}

}


//! Checks which of the inactive vars are needed to cover all chunk connection constraints.
//! Those then are added to the connectVars.
void CPlanarityMaster::generateVariablesForFeasibility(
		const List<ChunkConnection*>& ccons,
		List<CPlanarEdgeVar*>& connectVars)
{
	List<ChunkConnection*> cpy(ccons);
//	forall_listiterators(ChunkConnection*, ccit, cpy) {
//		(*ccit)->printMe();
//	}

	//First we check which of the constraints are already covered by existing
	//connect vars and delete them.
	ListIterator<CPlanarEdgeVar*> itev = connectVars.begin();
	while (itev.valid())
	{
		nodePair np((*itev)->sourceNode(),(*itev)->targetNode());
		ListIterator<ChunkConnection*> ccit = cpy.begin();
		while(ccit.valid()) {
			if((*ccit)->coeff(np)) {
				ListIterator<ChunkConnection*> delme = ccit;
				++ccit;
				cpy.del(delme);
			} else
				++ccit;
		}
		itev++;
	}//while connect vars

	ArrayBuffer<ListIterator<nodePair> > creationBuffer(ccons.size());
	forall_nonconst_listiterators(nodePair, npit, m_inactiveVariables) {
		bool select = false;

		ListIterator<ChunkConnection*> ccit = cpy.begin();
		while(ccit.valid()) {
			if((*ccit)->coeff(*npit)) {
				ListIterator<ChunkConnection*> delme = ccit;
				++ccit;
				cpy.del(delme);
				select = true;
			} else
				++ccit;
		}
		if(select) {
			creationBuffer.push(npit);
		}
		if(cpy.size()==0) break;
	}

	OGDF_ASSERT(cpy.size()==0);
	//Logger::slout() << "Creating " << creationBuffer.size() << " Connect-Variables for feasibility\n";
	m_varsInit = creationBuffer.size();
	// realize creationList
	for(int i = creationBuffer.size(); i-->0;) {
	  connectVars.pushBack( createVariable( creationBuffer[i] ) );
	}
}//generateVariablesForFeasability

void CPlanarityMaster::initializeOptimization() {
	m_nSep = 0;
	m_solState = ss_udf;
	//we don't try heuristic improvement (edge addition)
	heuristicLevel(0);
	//enumerationStrategy(BreadthFirst);
	// Create an analysis object to check for vertex activity state
	// Todo: As a partition into independent bags has to be done
	// externally, we could save work by using these results instead
	// of recomputing parts here. However, in the enclosing module,
	// we might work on a larger input graph that is partitioned,
	// i.e. we might get a copy of a part here only.
	m_ca = new ClusterAnalysis(*m_C, false); //use outer active lists, but no indy bag info

	if (pricing())
		varElimMode(NoVarElim);
	else
		varElimMode(ReducedCost);
	conElimMode(Basic);
	if(pricing())
		pricingFreq(1);

	//----------------------------Creation of Variables--------------------------------//
	// List for connection edges
	List<CPlanarEdgeVar*> connectVars; // MCh: ArrayBuffer would speed this up

	// First, we create variables necessary for an initial solution (abacus needs some vars)
	// For the search space reduction case, this only creates variables that
	// will be necessary in any solution. For the pricing case, we may need to
	// add some more.
	createInitialVariables(connectVars);

#ifdef OGDF_DEBUG
	cout << "Creating "<<connectVars.size()<<" initial variables\n";
#endif

	int nComplete = (m_G->numberOfNodes()*(m_G->numberOfNodes()-1))/2;
	int nConnectionEdges = nComplete - m_G->numberOfEdges();

	m_varsMax = 0;
	node u,v;

	// First we use ClusterAnalysis to identify the edges that are necessary.
	if (m_shrink)
	{
#ifdef OGDF_DEBUG
Logger::slout() << "Starting shrinking\n";
#endif
		// We check for two restrictions:
		// Only edges connecting bags in a cluster
		// between two outeractive vertices are needed for connectivity.
		// Only edges connecting satchels of a cluster.
		//
		// We assume that there is only a single independent bag, i.e. all
		// bags in a cluster c have outgoing edges.
		// Run through all clusters, connection variables for pairs of outeractive
		// vertices from different bags.
		cluster c;
		forall_clusters(c, *m_C)
		{
			addInnerConnections(c, connectVars);
			if (c!=m_C->rootCluster())
				addExternalConnections(c, connectVars);
		}


	}
	else {
		// The edge check here is slow, can be sped up as we only handle
		// planar graphs (see BiconnectedShellingOrder).
		// All node pairs for which no initial variable is generated are
		// classified either as inactive variable or used for connectVar generation.
		forall_nodes(u,*m_G) {
			v = u->succ();
			while (v!=NULL) {
				if(!m_G->searchEdge(u,v)) {
					if(goodVar(u,v)) {
						if(pricing())
							m_inactiveVariables.pushBack(nodePair(u,v));
						else
							connectVars.pushBack( createVariable(u,v) );
					}
					++m_varsMax;
				}
				v = v->succ();
			}
		}
	}
	m_varsPotential = m_inactiveVariables.size();

	//-------------------Creation of ChunkConnection-Constraints------------------------//
	//TODO: This can be used to add good connectivity edges

	//adds all connections, so we either need to reduce this to the
	//shrinked searchspace (todo, compare), or skip.
	//WE could safely skip the whole part, but computation of some minor
	//values is intermingled with other sections, so in order not to mess it up,
	//I leave it for now as it is, wasting computation time. KK

		int nChunks = 0;

	List<ChunkConnection*> constraintsCC;

	// The Graph, in which each cluster-induced Graph is temporarily stored
	Graph subGraph;

	// Since the function inducedSubGraph(..) creates a new Graph \a subGraph, the original
	// nodes have to be mapped to the copies. This mapping is stored in \a orig2new.
	NodeArray<node> orig2new;

	ClusterArray<int> numCEdges(*m_C, 0);//keeps number of edges in induced graph, used
	//below for constraints
	// Iterate over all clusters of the Graph
	ListConstIterator<node> it;
	cluster c;
	forall_clusters(c,*m_C) {

		//List<node> clusterNodes;
		//c->getClusterNodes(clusterNodes);

		// Compute the cluster-induced Subgraph
		it = getClusterNodes(c).begin();
		inducedSubGraph(*m_G, it, subGraph, orig2new);
		numCEdges[c] = subGraph.numberOfEdges();

		// Compute the number of connected components
		NodeArray<int> components(subGraph);
		int nCC = connectedComponents(subGraph,components);
		nChunks+=nCC;
		// If the cluster consists of more than one connected component,
		// ChunkConnection constraints can be created.
		if (nCC > 1) {

			// Determine each connected component (chunk) of the current cluster-induced Graph
			for (int i=0; i<nCC; ++i) {

				ArrayBuffer<node> cC(subGraph.numberOfNodes());
				ArrayBuffer<node> cCComplement(subGraph.numberOfNodes());
				node v;
				forall_nodes(v,*m_G/*subGraph*/) {
					node n = orig2new[v];
					if(n) {
						if (components[n] == i) cC.push(v);
						else cCComplement.push(v);
					}
				}
				// Creating corresponding constraint
				if (!m_shrink) {
					constraintsCC.pushBack(new ChunkConnection(this, cC, cCComplement));
				}
				// Avoiding duplicates if cluster consists of 2 chunks
				if (nCC == 2) break;

			}//connected component has been processed
		}//end if(nCC > 1)
	}//end forall_clusters
	if(pricing())
		generateVariablesForFeasibility(constraintsCC, connectVars);


	/* This part is useless for pure CPlanarity checks, but we might want to
	   test how fast the implementation is when checking non-planar graphs,
	   therefore the code is not deleted. However, in this case first the
	   alternative part for the shrinked search space model needs to be added,
	   as so far all node pairs are added to the constraint.
	//------------Creation of MaxPlanarEdges-Constraints---------------//

	List<MaxPlanarEdgesConstraint*> constraintsMPE;

	int nMaxPlanarEdges = 3*m_G->numberOfNodes() - 6 - m_G->numberOfEdges();
//temp removal of all mpe
	if (m_G->numberOfNodes() > 2)
		constraintsMPE.pushBack(new MaxPlanarEdgesConstraint(this,nMaxPlanarEdges));

	//List<node> clusterNodes;
	List<nodePair> clusterEdges;
	forall_clusters(c,*m_C) {

		if (c == m_C->rootCluster()) continue;
		//clusterNodes.clear();
		clusterEdges.clear();
		const List<node> &clusterNodes = getClusterNodes(c);
		//c->getClusterNodes(clusterNodes);
		if (clusterNodes.size() >= 4) {
			nodePair np;
			ListConstIterator<node> it;
			ListConstIterator<node> it_succ;
			for (it=clusterNodes.begin(); it.valid(); ++it) {
				it_succ = it.succ();
				while (it_succ.valid()) {
					np.v1 = (*it); np.v2 = (*it_succ);
					clusterEdges.pushBack(np);
					it_succ++;
				}
			}
			int maxPlanarEdges = 3*clusterNodes.size() - 6 - numCEdges[c];//TODO: This ignores original edges in induced graph,
			//subtract their number, could use edge number from subgraph computation above
			constraintsMPE.pushBack(new MaxPlanarEdgesConstraint(this,maxPlanarEdges,clusterEdges));
		}
	}*/

	//------------------------Adding Constraints to the Pool---------------------//

	// Adding constraints to the standardpool
	ArrayBuffer<Constraint *> initConstraints(constraintsCC.size()/*+constraintsMCC.size()+constraintsMPE.size()*/,false);

	ListConstIterator<ChunkConnection*> ccIt;
	updateAddedCCons(constraintsCC.size());
	for (ccIt = constraintsCC.begin(); ccIt.valid(); ++ccIt) {
		initConstraints.push(*ccIt);
	}
	//KK: Could be slow however, so comparison would be nice
	//ListConstIterator<MinimalClusterConnection*> mccIt;
	//for (mccIt = constraintsMCC.begin(); mccIt.valid(); ++mccIt) {
	//	initConstraints.push(*mccIt);
	//}
	//ListConstIterator<MaxPlanarEdgesConstraint*> mpeIt;
	//for (mpeIt = constraintsMPE.begin(); mpeIt.valid(); ++mpeIt) {
	//	initConstraints.push(*mpeIt);
	//}

	//---------------------Create Search Space Graph with initial variables-//
	// Currently this graph is only used for mincut in CPlanarity_Sub
	if (m_ssg != 0) delete m_ssg;
	m_ssg = new GraphCopy(*m_G); //note that we don't care about clusters here
#ifdef OGDF_DEBUG
	Logger::slout() << "SSG creation size: "<<m_ssg->numberOfNodes()<<" "<<m_ssg->numberOfEdges()<<"\n";
	cout << "SSG creation size: "<<m_ssg->numberOfNodes()<<" "<<m_ssg->numberOfEdges()<<"\n";
#endif
	//---------------------Adding Variables to the Pool---------------------//

	// Adding variables to the standardpool
	ArrayBuffer<Variable *> edgeVariables(connectVars.size(),false);
	ListConstIterator<CPlanarEdgeVar*> eIt;

//#ifdef OGDF_DEBUG
	Logger::ssout() << "Creating "<<connectVars.size()<<" variables\n";
	Logger::ssout() << "out of a maximum of "<<nConnectionEdges<< " conn vars\n";  // and "<<nMaxPlanarEdges<< " max planar connections\n";
//#endif
	for (eIt = connectVars.begin(); eIt.valid(); ++eIt) {
		edgeVariables.push(*eIt);
		// update the search space graph
		m_ssg->newEdge(m_ssg->copy((*eIt)->sourceNode()),m_ssg->copy((*eIt)->targetNode()));
	}
	//-------------------------------
	// We check if we added enough variables: Is the graph with all potential edges compconn?
	//TODO
#ifdef OGDF_DEBUG
	//Just for special debugging purposes:
	if (true) {
		ClusterArray<cluster> ca(*m_C);
		Graph GG;
		NodeArray<node> na(*m_G);
		ClusterGraph CG(*m_C,GG, ca, na);

		eIt = connectVars.begin();

		List<edge> le;

		while (eIt.valid()) {
			node cv, cw;
			// Add all new connection edges to \a m_solutionGraph
			cv = na[(*eIt)->sourceNode()];
			cw = na[(*eIt)->targetNode()];
			edge e = GG.newEdge(cv,cw);
			le.pushBack(e);

			eIt++;
		}

		ClusterGraphAttributes CGA(CG, GraphAttributes::edgeType | GraphAttributes::nodeType |
				GraphAttributes::nodeGraphics | GraphAttributes::edgeGraphics | GraphAttributes::edgeStyle);
		ListConstIterator<edge> it = le.begin();
		while (it.valid())
		{
			//cout << (*it)->graphOf() << "\n";
			//cout << &GG << "\n";
			CGA.strokeColor(*it) = "#FF0000";
			it++;
		}
		cluster cc;
		forall_clusters(cc,*m_C)
		{
			CGA.height(cc) = 10;
			CGA.width(cc) = 10;
			CGA.strokeWidth(cc) = 1.0;

		}
		GraphIO::writeGML(CGA, "CompleteExtension.gml");
		GraphIO::writeGML(CG, "CompleteExtensionCG.gml");
	}
#endif
	//-------------------------------


	//---------------------Initializing the Pools---------------------------//

	int poolsize = (getGraph()->numberOfNodes() * getGraph()->numberOfNodes());
	if (useDefaultCutPool())
		initializePools(initConstraints, edgeVariables, m_nMaxVars, poolsize, true);
	else
	{
		initializePools(initConstraints, edgeVariables, m_nMaxVars, 0, false);
		//TODO: How many of them?
		m_cutConnPool = new StandardPool<Constraint, Variable>(this, poolsize, true);
		m_cutKuraPool = new StandardPool<Constraint, Variable>(this, poolsize, true);
	}


	//---------------------Initialize Lower Bound---------------------------//

	//if we check only for c-planarity, we cannot set bounds
	//if (!m_checkCPlanar)
	//{
		//This seems to be no good practice with ABACUS
		//dualBound(heuristicInitialLowerBound()); // TODO-TESTING

#ifdef OGDF_DEBUG
		cout << "Dualbound: "<<dualBound()<<"\n";
		cout << "Infinity:"  <<infinity()<<"\n";
#endif
	//---------------------Initialize Upper Bound---------------------------//
	//TODO: Should be working here but abacus wont work if set
		//primalBound(heuristicInitialUpperBound());
	//TODO: Check bound initialization here

	//----------------------Setting Parameters------------------------------//

//	conElimMode(Master::NonBinding);
	maxCpuTime(*m_maxCpuTime);

	Logger::ssout() << "#Nodes: " << m_G->numberOfNodes() << "\n";
	Logger::ssout() << "#Edges: " << m_G->numberOfEdges() << "\n";
	Logger::ssout() << "#Clusters: " << m_C->numberOfClusters() << "\n";
	Logger::ssout() << "#Chunks: " << nChunks << "\n";


}

// returns coefficients of all variables in connect in constraint con
// as list coeffs
void CPlanarityMaster::getCoefficients(Constraint* con,
	const List<CPlanarEdgeVar* > & connect, List<double> & coeffs)
{
	coeffs.clear();
	ListConstIterator<CPlanarEdgeVar*> cIt = connect.begin();
	while (cIt.valid())
	{
		coeffs.pushBack(con->coeff(*cIt));
		cIt++;
	}//while

}


//output statistics
//and change the list of deleted edges in case only c-planarity is tested
//(to guarantee that the list is non-empty if input is not c-planar)
void CPlanarityMaster::terminateOptimization() {
	delete m_ca;

	char prefix[4] = "CP-";
	char fprefix[3] = "F-";
	char* pre;
	if (m_shrink) pre = prefix;
	else pre = fprefix;
	if (isCP()) m_solState = ss_cp;
	else m_solState = ss_ncp;
	Logger::slout() << "=================================================\n";
	Logger::slout() << "Terminate Optimization:\n";
	Logger::slout() << "(primal Bound: " << primalBound() << ")\n";
	Logger::slout() << "(dual Bound: " << dualBound() << ")\n";
	//if(m_checkCPlanar2) {
	Logger::slout() << "*** " << (isCP() ? "" : "NON ") << "C-PLANAR ***\n";
	//} else {
	//	Logger::slout() << "*** " << (feasibleFound() ? "" : "NON ") << "C-PLANAR ***\n";
	//}
	Logger::slout() << "=================================================\n";

	Logger::ssout() << "\n";

	Logger::ssout() << pre<<"C-Planar: " << isCP() << "\n";
	Logger::ssout() << pre<<"Time: "<< getDoubleTime(totalTime()) << "\n";
	Logger::ssout() << pre<<"LP-Time: " << getDoubleTime(lpSolverTime()) << "\n";
	Logger::ssout() << "Search space: " << (m_shrink ? " reduced " : " complete ") << "\n";

	Logger::ssout() << "\n";

	Logger::ssout() << pre<<"#BB-nodes: " << nSub() << "\n";
	Logger::ssout() << pre<<"#LP-relax: " << m_solvesLP << "\n";
	Logger::ssout() << pre<<"#Separations: " << m_nSep <<"\n";

	Logger::ssout() << pre<<"#Cut-Constraints: " << m_nCConsAdded << "\n";
	Logger::ssout() << pre<<"#Kura-Constraints: " << m_nKConsAdded << "\n";
	Logger::ssout() << pre<<"#Vars-init: " << m_varsInit << "\n";
	Logger::ssout() << pre<<"#Vars-used: " << m_varsAdded << "\n";
	Logger::ssout() << pre<<"#Vars-potential: " << m_varsPotential << "\n";
	Logger::ssout() << pre<<"#Vars-max: " << m_varsMax << "\n";
	Logger::ssout() << pre<<"#Vars-cut: " << m_varsCut << "\n";
	Logger::ssout() << pre<<"#Vars-kurarepair: " << m_varsKura << "\n";
	Logger::ssout() << pre<<"#Vars-price: " << m_varsPrice << "\n";
	Logger::ssout() << pre<<"#Vars-branch: " << m_varsBranch << "\n";
	Logger::ssout() << pre<<"#Vars-unused: " << m_inactiveVariables.size() << "\n";
	Logger::ssout() << pre<<"KuraRepair-Stat: <";

	for(int i =0; i<m_repairStat.size(); ++i) {
		Logger::ssout() << m_repairStat[i] << ",";
	}
	Logger::ssout() << ">\n";

	node n,m;
	edge e;
	forall_nodes(n, *m_G) {
		forall_nodes(m, *m_G) {
			if(m->index()<=n->index()) continue;
			forall_adj_edges(e, n) {
				if(e->opposite(n)==m) {
#ifdef OGDF_DEBUG
					Logger::slout() << "ORIG: " << n << "-" << m << "\n";
#endif
					continue;
				}
			}
		}
	}
	forall_nodes(n, *m_G) {
		forall_nodes(m, *m_G) {
			if(m->index()<=n->index()) continue;
			forall_adj_edges(e, n) {
				if(e->opposite(n)==m) {
					goto wup;
				}
			}
			forall_listiterators(nodePair, it, m_inactiveVariables) {
				if( ((*it).v1==n && (*it).v2==m) || ((*it).v2==n && (*it).v1==m)) {
					goto wup;
				}
			}
#ifdef OGDF_DEBUG
			Logger::slout() << "CONN: " << n << "-" << m << "\n";
#endif
			wup:;
		}
	}

	globalPrimalBound = primalBound();
	globalDualBound = dualBound();
}


#endif // USE_ABACUS
