gsoc2013
========

a) To test and use ogdf.js directly

	1) Clone this repository.
	2) Go to src/ogdf.js
	3) ogdf.js is the file which contains bindings for some OGDF functions.
	4) Read main.html to know how to use ogdf.js and how to call functions.

b) To get it running so that you can add bindings for more functions

	1) Clone this repository.
	2) Get emscripten from this link - https://github.com/kripken/emscripten
	3) Add emscripten to PATH variable. export PATH=/path/to/emscripten:$PATH
	4) Now you can use emcc to compile cpp files.
	5) EMSCRIPTEN Bindings needs to be created for the functions we want to use. To see example go to end of file src/ogdf/basic/Graph.cpp
	6) For more details related to embind functions - https://github.com/kripken/emscripten/wiki/embind
	7) To compile and get your bindings in ogdf.js go to folder src/ogdf/basic 
		 emcc --bind -o ../../ogdf.js/ogdf.js Graph.cpp PoolMemoryAllocator.cpp graph_generators.cpp graphics.cpp GraphAttributes.cpp ../	fileformats/GraphIO.cpp ../layered/acyclic_subgraph.cpp simple_graph_alg.cpp -I/path/to/gsoc2013/include
    8) Your new ogdf.js will be present in src/ogdf.js

List of function which can be tried are:

Graph::Graph()
Graph::numberOfNodes()
Graph::numberOfEdges()
Graph::maxNodeIndex()
Graph::maxEdgeIndex()
Graph::firstNode()
Graph::lastNode()
Graph::firstEdge()
Graph::lastEdge()	
Graph::chooseNode()
Graph::chooseEdge()
Graph::empty()
Graph::newNode()
Graph::newNode(int)
Graph::newEdge(node,node)
Graph::newEdge(node,node,int)

GraphAttributes()
GraphAttributes(Graph&,long)
GraphAttributes::x
GraphAttributes::y
GraphAttributes::width	
GraphAttributes::height
GraphAttributes::bends
GraphAttributes::strokeColor
GraphAttributes::fillColor
setX
setY
setWidth
setHeight
setEdgeColor
setEdgeColor
	
List<edge>()
List<edge>::size
List<edge>::empty

List<node>()
List<node>::size
List<node>::empty

GraphList<EdgeElement>()
GraphList<EdgeElement>::size
GraphList<EdgeElement>::empty

GraphList<NodeElement>()
GraphList<NodeElement>::size
GraphList<NodeElement>::empty

Color::Color()
Color::Color(Color::Name)
Color::toString()

enum_  Color::Name::Red Color::Name::Blue Color::Name::Green)
DPolyline::DPolyline()
		
randomGraph
randomSimpleGraph
completeGraph
completeBipartiteGraph
planarConnectedGraph

DfsAcyclicSubgraph()
DfsAcyclicSubgraph::call
