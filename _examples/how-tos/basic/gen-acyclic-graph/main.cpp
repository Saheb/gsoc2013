#include <ogdf/basic/Graph.h>
#include <ogdf/basic/graph_generators.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <emscripten/bind.h>
using namespace emscripten;

int main()
{
	ogdf::Graph G;
	ogdf::randomSimpleGraph(G, 10, 20);
	ogdf::DfsAcyclicSubgraph DAS;
	DAS.callAndReverse(G);
	ogdf::GraphIO::writeGML(G, "test.gml");

	return 0;
}

EMSCRIPTEN_BINDINGS(my_module) {
    function("main", &main);
}
