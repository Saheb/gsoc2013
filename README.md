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
