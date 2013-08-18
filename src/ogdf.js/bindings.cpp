#include <ogdf/basic/Array.h>
#include <ogdf/basic/AdjEntryArray.h>
#include <ogdf/fileformats/GmlParser.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/GraphObserver.h>

#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/fileformats/GmlParser.h>

#include <ogdf/basic/graph_generators.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/CombinatorialEmbedding.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/NodeArray.h>
#include <ogdf/basic/FaceArray.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/basic/Array2D.h>
#include <ogdf/planarity/PlanarizationGridLayout.h>
#include <ogdf/planarlayout/SchnyderLayout.h>

#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/layered/GreedyCycleRemoval.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/Queue.h>

#include <ogdf/basic/Logger.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/fileformats/GmlParser.h>
#include <ogdf/fileformats/OgmlParser.h>
#include <sstream>
#include <map>

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/Queue.h>

#include <ogdf/basic/GridLayout.h>
#include <ogdf/cluster/ClusterGraphAttributes.h>
#include <ogdf/internal/steinertree/EdgeWeightedGraph.h>
#include <sstream>

#include <ogdf/basic/graphics.h>

#include <ogdf/layered/Hierarchy.h>
#include <ogdf/layered/SugiyamaLayout.h>
#include <ogdf/layered/LongestPathRanking.h>
#include <ogdf/layered/BarycenterHeuristic.h>
#include <ogdf/layered/SplitHeuristic.h>
#include <ogdf/layered/FastHierarchyLayout.h>
#include <ogdf/layered/OptimalHierarchyClusterLayout.h>
#include <ogdf/packing/TileToRowsCCPacker.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/basic/Thread.h>
#include <ogdf/basic/CriticalSection.h>

#include <ogdf/layered/OptimalRanking.h>
#include <ogdf/layered/DfsAcyclicSubgraph.h>
#include <ogdf/basic/simple_graph_alg.h>
#include <ogdf/graphalg/MinCostFlowReinelt.h>
#include <ogdf/basic/GraphCopy.h>

#include <emscripten/bind.h>

namespace ogdf{
//--------------------------Functions only for JS--------------------------//

void setX(GraphAttributes& GA,node n, double val)
{
	GA.x(n) = val;
}
void setY(GraphAttributes& GA,node n, double val)
{
	GA.y(n) = val;
}

void setWidth(GraphAttributes& GA,node n, double val)
{
	GA.width(n) = val;
}

void setHeight(GraphAttributes& GA,node n, double val)
{
	GA.height(n) = val;
}

void setEdgeColor(GraphAttributes& GA,edge e, Color c)
{
	GA.strokeColor(e) = c;
}
void setNodeColor(GraphAttributes& GA,node n, Color c)
{
	GA.fillColor(n) = c;
}
void setStrokeWidth(GraphAttributes &GA,edge e, float f)
{
	GA.strokeWidth(e) = f;
}

std::string getSVG(const GraphAttributes &A)
{
	std::stringstream os;
	ogdf::GraphIO::drawSVG(A,os);
	std::string str;
	str = os.str();
	return str;
}
//--------------------------------------------------------------//

}
using namespace emscripten;


EMSCRIPTEN_BINDINGS(graph) {
    class_<ogdf::Graph>("Graph")
		.constructor()		
		.function("numberOfNodes",&ogdf::Graph::numberOfNodes)
		.function("numberOfEdges",&ogdf::Graph::numberOfEdges)
		.function("maxNodeIndex",&ogdf::Graph::maxNodeIndex)
		.function("maxEdgeIndex",&ogdf::Graph::maxEdgeIndex)
		.function("firstNode",&ogdf::Graph::firstNode,allow_raw_pointers())
		.function("lastNode",&ogdf::Graph::lastNode,allow_raw_pointers())
		.function("firstEdge",&ogdf::Graph::firstEdge,allow_raw_pointers())
		.function("lastEdge",&ogdf::Graph::lastEdge,allow_raw_pointers())	
		.function("chooseNode",&ogdf::Graph::chooseNode,allow_raw_pointers())
		.function("chooseEdge",&ogdf::Graph::chooseEdge,allow_raw_pointers())	
		.function("empty",&ogdf::Graph::empty)
		.function("newNode", select_overload<ogdf::node()>(&ogdf::Graph::newNode),allow_raw_pointers())
		.function("newNode",select_overload<ogdf::node(int)>(&ogdf::Graph::newNode),allow_raw_pointers())
		.function("newEdge", select_overload<ogdf::edge(ogdf::node,ogdf::node)>(&ogdf::Graph::newEdge),allow_raw_pointers())
		//.function("newEdge", select_overload<ogdf::edge(ogdf::node,ogdf::node,int)>(&ogdf::Graph::newEdge),allow_raw_pointers())
		//.function("newEdge", select_overload<ogdf::edge(ogdf::node)>(&ogdf::Graph::newNode),allow_raw_pointers())
		//.function("newEdge", select_overload<ogdf::edge()>(&ogdf::Graph::newNode),allow_raw_pointers())
		;
	class_<ogdf::NodeElement>("NodeElement")
		;
	class_<ogdf::EdgeElement>("EdgeElement")
		;
	//register_vector<std::shared_ptr<ogdf::NodeElement>>("NodeElement");
	//register_vector<std::shared_ptr<ogdf::EdgeElement>>("EdgeElement");
	class_<ogdf::GraphAttributes>("GraphAttributes")
        .constructor()
		.constructor<ogdf::Graph&,long>()
		.smart_ptr<std::shared_ptr<ogdf::GraphAttributes>>()
		.function("x", select_overload<double&(ogdf::node)>(&ogdf::GraphAttributes::x),allow_raw_pointers())
		//.function("x", select_overload<double(ogdf::node)>(&ogdf::GraphAttributes::x))
		.function("y", select_overload<double&(ogdf::node)>(&ogdf::GraphAttributes::y),allow_raw_pointers())
		//.function("y", select_overload<double(ogdf::node)>(&ogdf::GraphAttributes::y))
		.function("width", select_overload<double&(ogdf::node)>(&ogdf::GraphAttributes::width),allow_raw_pointers())
		//.function("width", select_overload<double(ogdf::node)>(&ogdf::GraphAttributes::width))
		.function("height", select_overload<double&(ogdf::node)>(&ogdf::GraphAttributes::height),allow_raw_pointers())
		//.function("height", select_overload<double(ogdf::node)>(&ogdf::GraphAttributes::height))
		.function("bends", select_overload<ogdf::DPolyline&(ogdf::edge)>(&ogdf::GraphAttributes::bends),allow_raw_pointers())
		//.function("bends", select_overload<const ogdf::DPolyline&(ogdf::edge)>(&ogdf::GraphAttributes::bends),allow_raw_pointers())
		.function("strokeColor", select_overload<ogdf::Color&(ogdf::edge)>(&ogdf::GraphAttributes::strokeColor),allow_raw_pointers())
		//.function("strokeColor", select_overload<const ogdf::Color&(ogdf::edge)>(&ogdf::GraphAttributes::strokeColor),allow_raw_pointers())
		.function("fillColor", select_overload<ogdf::Color&(ogdf::node)>(&ogdf::GraphAttributes::fillColor),allow_raw_pointers())
		//.function("fillColor", select_overload<const ogdf::Color&(ogdf::node)>(&ogdf::GraphAttributes::fillColor),allow_raw_pointers())
		.function("strokeWidth", select_overload<float&(ogdf::edge)>(&ogdf::GraphAttributes::strokeWidth),allow_raw_pointers())
		//.function("strokeWidth", select_overload<float(ogdf::edge)>(&ogdf::GraphAttributes::strokeWidth),allow_raw_pointers())
		.function("strokeWidth", select_overload<float&(ogdf::node)>(&ogdf::GraphAttributes::strokeWidth),allow_raw_pointers())
		//.function("strokeWidth", select_overload<float(ogdf::edge)>(&ogdf::GraphAttributes::strokeWidth),allow_raw_pointers())
		.function("setStrokeType", select_overload<void(ogdf::edge,ogdf::StrokeType)>(&ogdf::GraphAttributes::setStrokeType),allow_raw_pointers())
		//.function("setStrokeType", select_overload<void(ogdf::node,ogdf::StrokeType)>(&ogdf::GraphAttributes::setStrokeType),allow_raw_pointers())
		//.function("strokeType", select_overload<ogdf::StrokeType (ogdf::edge)>(const &ogdf::GraphAttributes::strokeType ))
		//.function("strokeType", select_overload<ogdf::StrokeType (ogdf::edge)>(const &ogdf::GraphAttributes::strokeType ))
		.function("setX",&ogdf::setX,allow_raw_pointers())
		.function("setY",&ogdf::setY,allow_raw_pointers())
		.function("setWidth",&ogdf::setWidth,allow_raw_pointers())
		.function("setHeight",&ogdf::setHeight,allow_raw_pointers())
		.function("setEdgeColor",&ogdf::setEdgeColor,allow_raw_pointers())
		.function("setNodeColor",&ogdf::setNodeColor,allow_raw_pointers())	
		.function("setStrokeWidth",&ogdf::setStrokeWidth,allow_raw_pointers())
		;
		constant("nodeGraphics",0x00001);
		constant("edgeGraphics",0x00002);
		constant("nodeStyle",0x00800);
		constant("edgeStyle",0x00400);
		constant("edgeType",0x00040);
		constant("nodeType",0x00080);
		
	class_<ogdf::List<ogdf::edge>>("List<edge>")
		.constructor()
		.function("size",&ogdf::List<ogdf::edge>::size)
		.function("empty",&ogdf::List<ogdf::edge>::empty)
		;
	class_<ogdf::List<ogdf::node>>("List<node>")
		.constructor()
		.function("size",&ogdf::List<ogdf::node>::size)
		.function("empty",&ogdf::List<ogdf::node>::empty)
		;
	class_<ogdf::List<ogdf::DPoint>>("List<DPoint>")
		.constructor()
		;
	class_<ogdf::GraphList<ogdf::EdgeElement>>("GraphList<edge>")
		.constructor()
		.function("size",&ogdf::GraphList<ogdf::EdgeElement>::size)
		.function("empty",&ogdf::GraphList<ogdf::EdgeElement>::empty)
		;
	class_<ogdf::GraphList<ogdf::NodeElement>>("GraphList<node>")
		.constructor()
		.function("size",&ogdf::GraphList<ogdf::NodeElement>::size)
		.function("empty",&ogdf::GraphList<ogdf::NodeElement>::empty)
		;
	class_<ogdf::Color>("Color")
		.constructor()
		.constructor<ogdf::Color::Name>()
		.function("toString",&ogdf::Color::toString)
		;
	enum_<ogdf::StrokeType>("StrokeType")
		.value("stNone",ogdf::StrokeType::stNone)
		.value("stSolid",ogdf::StrokeType::stSolid)
		.value("stDash",ogdf::StrokeType::stDash)
		;
	enum_<ogdf::Color::Name>("Name")
		.value("Red",ogdf::Color::Name::Red)
		.value("Blue",ogdf::Color::Name::Blue)
		.value("Green",ogdf::Color::Name::Green)
		;
	class_<ogdf::DPoint>("DPoint")
		.constructor()
		.constructor<double,double>()
		;

	class_<ogdf::DPolyline/*,base<ogdf::List<ogdf::DPoint>*/>("DPolyline")
		.constructor()
		//.function("pushBack",&ogdf::List<DPoint>::pushBack)
		;

	function("randomGraph", &ogdf::randomGraph);
	function("randomSimpleGraph", &ogdf::randomSimpleGraph);
	function("completeGraph", &ogdf::completeGraph);
	function("completeBipartiteGraph", &ogdf::completeBipartiteGraph);
	function("planarConnectedGraph",&ogdf::planarConnectedGraph);

	class_<ogdf::DfsAcyclicSubgraph/*, base<ogdf::AcyclicSubgraphModule>*/>("DfsAcyclicSubgraph")
		.smart_ptr<std::shared_ptr<ogdf::DfsAcyclicSubgraph>>()
		.constructor()
		.function("call",&ogdf::DfsAcyclicSubgraph::call)
		;
		//class_<ogdf::AcyclicSubgraphModule>("AcyclicSubgraphModule")
		//.smart_ptr<std::shared_ptr<ogdf::AcyclicSubgraphModule>>()
		//.function("callAndReverse",select_overload<void(ogdf::Graph&)>(&ogdf::AcyclicSubgraphModule::callAndReverse))
		//.function("callAndReverse",select_overload<void(ogdf::Graph&)>(&ogdf::AcyclicSubgraphModule::callAndReverse))
		;
	class_<ogdf::OptimalRanking,base<ogdf::RankingModule> >("OptimalRanking")
		.constructor()
		;

	class_<ogdf::LongestPathRanking>("LongestPathRanking")
		.constructor()
		;

	class_<ogdf::SugiyamaLayout>("SugiyamaLayout")
		.constructor()
		.function("call",select_overload<void(ogdf::GraphAttributes&)>(&ogdf::SugiyamaLayout::call))
		//.function("setRanking",&ogdf::SugiyamaLayout::setRanking)
		//.function("setCrossMin",&ogdf::SugiyamaLayout::setCrossMin)
		;



	
	class_<ogdf::GraphIO>("GraphIO")
        .constructor()
		.class_function("readGML", select_overload<bool(ogdf::Graph&,const char*)>(&ogdf::GraphIO::readGML),allow_raw_pointers())
		.class_function("readGML", select_overload<bool(ogdf::Graph&,const string&)>(&ogdf::GraphIO::readGML))
		//.function("readGML", select_overload<bool(ogdf::Graph&,std::istream)>(&ogdf::GraphIO::readGML))
		.class_function("writeGML", select_overload<bool(const ogdf::Graph&,const char*)>(&ogdf::GraphIO::writeGML),allow_raw_pointers())
		.class_function("writeGML", select_overload<char*(const ogdf::Graph&,const string&)>(&ogdf::GraphIO::writeGML),allow_raw_pointers())
		.class_function("getSVG",(&ogdf::getSVG))
		;	
}		
