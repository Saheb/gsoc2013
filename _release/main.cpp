#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/fileformats/GraphIO.h>
//#include <emscripten/bind.h>
using namespace ogdf;

int my_main(int j)
{
	Graph G;
	randomSimpleGraph(G, 10, 20);
	DfsAcyclicSubgraph DAS;
	DAS.callAndReverse(G);
	GraphIO::writeGML(G, "test.gml");

	return G.numberOfNodes()*5;
}
//using namespace emscripten;
//EMSCRIPTEN_BINDINGS(my_module) {
  //  function("my_main", &my_main);
//}
