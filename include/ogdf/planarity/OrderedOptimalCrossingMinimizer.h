/*
 * $Revision: 3440 $
 *
 * last checkin:
 *   $Author: gutwenger $
 *   $Date: 2013-04-22 18:34:09 +0530 (Mon, 22 Apr 2013) $
 ***************************************************************/

/** \file
 * \brief Declaration of class OrderedOptimalCrossingMinimizer
 *
 * \author Markus Chimani
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

#ifndef OGDF_ORDERED_OPTIMAL_CROSSING_MINIMIZER_H
#define OGDF_ORDERED_OPTIMAL_CROSSING_MINIMIZER_H

#include <ogdf/external/abacus.h>

#include <ogdf/basic/basic.h>
#include <ogdf/basic/Logger.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/basic/List.h>
#include <ogdf/basic/Array.h>
#include <ogdf/basic/SList.h>
#include <ogdf/basic/ArrayBuffer.h>
#include <ogdf/basic/EdgeArray.h>
#include <ogdf/basic/GraphCopy.h>
#include <ogdf/graphalg/GraphReduction.h>
#include <ogdf/module/CrossingMinimizationModule.h>
#include <ogdf/basic/ModuleOption.h>
#include <ogdf/planarity/KuratowskiSubdivision.h>
#include <ogdf/basic/MinHeap.h>

#define ATTRIBUTE(T, N) private: T _##N; public: const T& N() const { return _##N; } T& N() { return _##N; }
#define SAFE_DELETE(P) if(P) { delete P; P = 0; }


namespace ogdf {

	// forward declaration
	template<class C> class Compare_Equals;
	template<class C> class Compare_CheckedEqual;

class OGDF_EXPORT OrderedOptimalCrossingMinimizer : public CrossingMinimizationModule {

#ifndef USE_ABACUS
protected:
	// Guess what...
	virtual ReturnType doCall(PlanRep &PG,
		int cc,
		const EdgeArray<int>  *_cost,
		const EdgeArray<bool> *_forbid,
		const EdgeArray<__uint32> *_subgraphs,
		int& crossingNumber)
	{ THROW_NO_ABACUS_EXCEPTION; return retError; };
};
#else // USE_ABACUS

public:

	enum INDSTUFF { // just to make it work. *should* be in master (is used in Master&Subproblem)
		I_INT = 0,
		I_RND = 1,
		INDUCEMENTS = 2,  //  # of indices
		I_USEFORSEPARATE = 1 }; // use this solution for separation


	enum GraphHint { GH_None, GH_Complete, GH_CompleteBipartite, GH_Hypercube, GH_ToroidalGrid, GH_Petersen };
	enum HintEffects { HE_KuratowskisMinusOne =        0x0001,  // 1, 2, 4, 8, ...
					   HE_AllSubKuratowskis =          0x0002,
					   HE_EdgeOrder =                  0x0004,
					   HE_NodeOrder =                  0x0008,
					   HE_IterativeLowerBound =        0x0010,
					   HE_HighKuratowskiCutsStatic =   0x0020,
					   HE_ExpensiveKuratowski =        0x0040,
					   HE_HypercubeMinusOne =          0x0080,
					   HE_ToroidalGridMinusOne =       0x0100,
					   HE_Simplicity =			       0x0200
						};
	enum SeparationMode { SM_Simple       = 0x0001,
						  SM_BoyerMyrvold = 0x0002 };

	enum PricingInit { PI_NoPricing = 0,
					   PI_Normal    = 1 };

	enum BranchingMode { BM_Traditional = 0,
						 BM_CompleteOdd = 1 };

	static int completeGuy(int n) {
		return int((n/2) * int((n-1)/2) * int((n-2)/2) * int((n-3)/2)) / 4;
	}
	static int bipartiteZara(int n, int m) {
		return int((n)/2) * int((n-1)/2) * int((m)/2) * int((m-1)/2);
	}
	static int hypercubeEggGuy(int n) {
		return int((5 * pow(4.0,n))/32) - int(int((n*n+1)/2) * pow(2.0,n-2));
	}
	static int toroidalCr(int n, int m) {
		return (n<m) ? ((n-2)*m) : ((m-2)*n);
	}

	class SimpleSeparationParams {
		ATTRIBUTE( int, runs )
		ATTRIBUTE( int, desperateRuns )
		ATTRIBUTE( int, maxCuts )
	public:
		SimpleSeparationParams() :
			_runs(50),
			_desperateRuns(300),
			_maxCuts(20) {}
	};

	class BoyerMyrvoldSeparationParams {
		ATTRIBUTE( int, runs )
		ATTRIBUTE( int, desperateRuns )
		ATTRIBUTE( int, extractions )
		ATTRIBUTE( int, runCuts )
		ATTRIBUTE( int, maxCuts )
		ATTRIBUTE( bool, bundle )
		ATTRIBUTE( bool, noE2 )
		ATTRIBUTE( bool, veryDifferent )
	public:
		BoyerMyrvoldSeparationParams() :
			_runs(20),
			_desperateRuns(100),
			_extractions(100),
			_runCuts(80),
			_maxCuts(1000),
			_bundle(false),
			_noE2(true),
			_veryDifferent(false) {}
	};

protected:

	static const double EPS /*= 0.00001*/;
	static const double SEG_EPS /*= 0.0001*/;

	 //forwards
	class Master;
	class Subproblem;
	class CrossingVariableBase;
	class SimpleCrossingVariable;
	class OrderedCrossingVariable;
	//class CrossingLocationComparer;
	class AbacusConstraint;
	class OmegaExistenceConstraint;
	class ReverseOmegaExistenceConstraint;
	class CyclicOrderConstraint;
	class MirrorOrderConstraint;
	class StackedOrderConstraint;
	class TriangleConstraint;
	class KuratowskiConstraintBase;
	class BasicKuratowskiConstraint;
	class EdgeOrderConstraint;
	class NodeOrderConstraint;

	class CrossingConfiguration {
	public:
//		CrossingConfiguration(const CrossingConfiguration& CC) : crossingEdges(*(CC.crossingEdges.graphOf())) {
//			edge e;
//			crossingNo = 0;
//			forall_edges(e, *(crossingEdges.graphOf())) {
//				forall_listiterators(edge, it, CC.crossingEdges[e])
//					crossingEdges[e].pushBack(*it);
//			}
//		}
		CrossingConfiguration(const Graph& G) : crossingNo(-1), crossingEdges(G) {}
		CrossingConfiguration(const PlanRep& PG, int crNo, bool direct) : crossingNo(crNo), crossingEdges() {
			if(direct) extractDirect(PG, crNo);
			else extractIndirect(PG, crNo);
		}
		void extractDirect(const PlanRep& PG, int crNo);
		void extractIndirect(const PlanRep& PG, int crNo);
		void paste(PlanRep& PG);
		int getCrossingNo() const { return crossingNo; }
		const List<edge>& getCrossingEdges(edge e) const {
			return crossingEdges[e];
		}
		bool probablyExchangeCrossingEdges(edge e, List<edge>& newList) {
			bool exc = false;
			if(crossingEdges[e].size() != newList.size())
				exc = true;
			else {
				ListConstIterator<edge> ii = crossingEdges[e].begin();
				ListConstIterator<edge> ei = newList.begin();
				for(; ii.valid(); ++ii, ++ei) {
					if(*ii != *ei) {
						exc = true;
						break;
					}
				}
			}
			if(exc)
				crossingEdges[e].exchange(newList);
			return exc;
		}
		//////////////////////
	private:
		int crossingNo;
		EdgeArray< List<edge> > crossingEdges;
	public:
		static edge targetCrossingEdge(const PlanRep& PR, edge x);
		static edge targetSuccEdge(const PlanRep& PR, edge x);
	private:
		edge findCrossingPosition(const PlanRep& PR, edge onEdge, edge crEdge);
	};

	struct ctwoEdge {
		edge x[2];
		ctwoEdge(const edge y0, const edge y1) { x[0]=y0; x[1]=y1; };
		ctwoEdge(const ctwoEdge& y) { x[0]=y[0]; x[1]=y[1]; };
		inline edge& operator[](int i) { return x[i]; }
		inline const edge& operator[](int i) const { return x[i]; }
	};

	class CrossingVariableBase : public abacus::Variable {
	public:
		const ctwoEdge e;
		CrossingVariableBase(abacus::Master* m, double cost, const edge t1, const edge t2, bool tLocal = false) :
			abacus::Variable(m,0,0,tLocal, // master / sub=0 / dynamic=false / local=false
				cost, // obj
				0,1,abacus::VarType::Binary),
			e(t1,t2) {}
	};

	class SimpleCrossingVariable : public CrossingVariableBase {
		friend std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::SimpleCrossingVariable& v);
	public:
		SimpleCrossingVariable(abacus::Master* m, const edge t1, const edge t2) :
			CrossingVariableBase(m, ((OrderedOptimalCrossingMinimizer::Master*)m)->getCost(t1, t2), t1, t2) {}
//		virtual ~SimpleCrossingVariable() {}
		inline bool correspondsTo(const edge& ee) const { return ee == e[0] || ee == e[1]; }
		inline bool correspondsTo(const node& n) const { return e[0]->source() ==  n || e[0]->target() ==  n || e[1]->source() ==  n || e[1]->target() ==  n; }
		inline static bool preferedOrder(const edge e1, const edge e2) { return e1->index() < e2->index(); }
		virtual const char* name() const { return "SimpleCrossingVariable"; }
		inline edge commonEdge(const SimpleCrossingVariable* s) const { return (s->correspondsTo(e[0]))?e[0]:( s->correspondsTo(e[1])?e[1]:0 ); }
		inline edge otherEdge(const edge& ee) const { return (ee==e[0])?e[1]:e[0]; }
		OGDF_NEW_DELETE
	};

	class OrderedCrossingVariable : public CrossingVariableBase {
		friend std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::OrderedCrossingVariable& v);
		const edge bef;
	public:
		OrderedCrossingVariable(abacus::Master* m, const edge _e, const edge _f, const edge _g) :
			CrossingVariableBase(m, 0, _e, _f,((OrderedOptimalCrossingMinimizer*)m)->localVariables()), bef(_g) {}
//		virtual ~OrderedCrossingVariable() {}

//		bool relatesTo(const edge& ee) const { return CrossingVariableBase::correspondsTo(ee) || ee == bef; }
//		bool relatesTo(const node& n) const { return CrossingVariableBase::correspondsTo(n) || bef->source() ==  n || bef->target() ==  n; }
		edge base() const { return e[0]; }
		edge crossedBy() const { return e[1]; }
		edge before() const { return bef; }
//		bool correspondsTo(const CrossingVariableBase* cp) const { return CrossingVariableBase::correspondsTo(cp->e[0]) && CrossingVariableBase::correspondsTo(cp->e[1]); }
		bool isTwin(const OrderedCrossingVariable* t) const { return base()==t->base() && crossedBy()==t->before() && before()==t->crossedBy(); }
		virtual const char* name() const { return "OrderedCrossingVariable"; }
		OGDF_NEW_DELETE
	};

	class CrossingVariableComparer {
	public:
		typedef CrossingVariableBase* CrossingVariablePointer;
		inline static int compare(const CrossingVariablePointer& x, const CrossingVariablePointer& y) {
			OGDF_ASSERT(x);
			OGDF_ASSERT(y);
			int ret;
			//cerr << " x.name=" << x->name() << "/y.name=" << y->name() << " " << flush; // fails because of constness
//			cerr << " x.typeid=" << typeid(*x).name() << "/y.typeid=" << typeid(*y).name() << " " << flush;
			const SimpleCrossingVariable *sx = dynamic_cast<const SimpleCrossingVariable*>(x);
			const SimpleCrossingVariable *sy = dynamic_cast<const SimpleCrossingVariable*>(y);
			const OrderedCrossingVariable *ox = sx ? 0 : dynamic_cast<const OrderedCrossingVariable*>(x);
			const OrderedCrossingVariable *oy = sy ? 0 : dynamic_cast<const OrderedCrossingVariable*>(y);

			if(sx && oy)
				return 1;
			else if(ox && sy)
				return -1;
			else if(sx && sy) {
				if((ret = (sx->e[0]->index() - sy->e[0]->index()))) return ret;
				return sx->e[1]->index() - sy->e[1]->index();
			} else if(ox && oy) {
				if((ret = (ox->base()->index() - oy->base()->index()))) return ret;
				if((ret = (ox->crossedBy()->index() - oy->crossedBy()->index()))) return ret;
				return ox->before()->index() - oy->before()->index();
			}
			OGDF_ASSERT( false ); // SHOULD NEVER BE HERE
			return 0; // avoid compiler warning
		}
		OGDF_AUGMENT_COMPARER(CrossingVariablePointer)
	};
	class EdgeComparer {
	public:
		static int compare(const edge& x, const edge& y) { return x->index() - y->index(); }
		OGDF_AUGMENT_STATICCOMPARER(edge)
	};

	class AbacusConstraint : public abacus::Constraint {
	public:
		AbacusConstraint(abacus::Master *master, const abacus::Sub *sub,
				abacus::CSense::SENSE sense, double rhs,
				bool dynamic, bool local, bool liftable) :
			abacus::Constraint(master, sub, sense, rhs, dynamic, local, liftable) {}
		AbacusConstraint(const AbacusConstraint& c) : abacus::Constraint(c) {};
		static const OrderedCrossingVariable* preferredVariable(const OrderedCrossingVariable* t1, const OrderedCrossingVariable* t2) {
			return ( t1->crossedBy()->index() < t1->before()->index() ) ? t1 : t2;
		}
	protected:
		inline OrderedOptimalCrossingMinimizer::Master* master() const {
			return (OrderedOptimalCrossingMinimizer::Master*)master_;
		}
		inline bool localVariables() const { return master()->localVariables(); }
		OGDF_NEW_DELETE
	};

	class DummyConstraint : public AbacusConstraint {
	public:
		DummyConstraint(abacus::Master* m) : AbacusConstraint(m, 0, abacus::CSense::Greater, 1, 0, 0, 1) {}
		double coeff(const abacus::Variable* cvar) const {
			return dynamic_cast<const SimpleCrossingVariable*>(cvar) ? 1 : 0;
		}
		virtual const char* name() const { return "DummyConstraint"; }
	};

	class OmegaExistenceConstraint : public AbacusConstraint {
		const OrderedCrossingVariable* ocv;
	public:
		OmegaExistenceConstraint(abacus::Master* m, const OrderedCrossingVariable* t) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, 0, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			ocv(t) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* scv;
			return (cvar == (const abacus::Variable*)ocv) ? -1 :
				( ((scv = dynamic_cast<const SimpleCrossingVariable*>(cvar)) && scv->correspondsTo(ocv->base()) && scv->correspondsTo(ocv->crossedBy())) ? 1 : 0);
		}
		virtual const char* name() const { return "OmegaExistenceConstraint"; }
	};

	class ReverseOmegaExistenceConstraint : public AbacusConstraint {
		const OrderedCrossingVariable* ocv;
	public:
		ReverseOmegaExistenceConstraint(abacus::Master* m, const OrderedCrossingVariable* t) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, 0, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			ocv(t) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* scv;
			return (cvar == (const abacus::Variable*)ocv) ? -1 :
				( ((scv = dynamic_cast<const SimpleCrossingVariable*>(cvar)) && scv->correspondsTo(ocv->base()) && scv->correspondsTo(ocv->before())) ? 1 : 0);
		}
		virtual const char* name() const { return "ReverseOmegaExistenceConstraint"; }
	};

	class CyclicOrderConstraint : public AbacusConstraint {
//		friend std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::CyclicOrderConstraint& k);
		const edge e, f, g, h;
	public:
		CyclicOrderConstraint(abacus::Master* m, const edge _e, const edge _f, const edge _g, const edge _h) :
			AbacusConstraint(m, 0, abacus::CSense::Less, 2, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			e(_e), f(_f), g(_g), h(_h) {}
		double coeff(const abacus::Variable* cvar) const {
			const OrderedCrossingVariable* o = dynamic_cast<const OrderedCrossingVariable*>(cvar);
			if(o && o->base() == e) {
				if( (o->crossedBy()==f && o->before()==g ) ||
					(o->crossedBy()==g && o->before()==h ) ||
					(o->crossedBy()==h && o->before()==f ) ) return 1;
			}
			return 0;
		}
		virtual const char* name() const { return "CyclicOrderConstraint"; }
		void print(std::ostream& os) const {
			os << "(e=" << e->index() << ";f=" << f->index() << ";g=" << g->index() << ";h=" << h->index() << ")";
		}
	};

	class MirrorOrderConstraint : public AbacusConstraint {
		const OrderedCrossingVariable* ocv;
	public:
		MirrorOrderConstraint(abacus::Master* m, const OrderedCrossingVariable* t) :
			AbacusConstraint(m, 0, abacus::CSense::Less, 1, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			ocv(t) {
				OGDF_ASSERT( ocv->crossedBy()->index() < ocv->before()->index() );
			}
		double coeff(const abacus::Variable* cvar) const {
			if(const OrderedCrossingVariable* o = dynamic_cast<const OrderedCrossingVariable*>(cvar)) {
				if(o == ocv) return 1;
				if(o->base() == ocv->base() && o->crossedBy()==ocv->before() && o->before()==ocv->crossedBy() ) return 1;
			}
			return 0;
		}
		virtual const char* name() const { return "MirrorOrderConstraint"; }
	};

	class StackedOrderConstraint : public AbacusConstraint {
		const OrderedCrossingVariable* ocv;
	public:
		StackedOrderConstraint(abacus::Master* m, const OrderedCrossingVariable* t) :
			AbacusConstraint(m, 0, abacus::CSense::Less, 1, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			ocv(t) {
				OGDF_ASSERT( ocv->crossedBy()->index() < ocv->before()->index() );
			}
		double coeff(const abacus::Variable* cvar) const {
			if(const OrderedCrossingVariable* o = dynamic_cast<const OrderedCrossingVariable*>(cvar)) {
				if(o == ocv) return -1;
				if(o->base() == ocv->base() && o->crossedBy()==ocv->before() && o->before()==ocv->crossedBy() ) return -1;
			} else {
				const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
				if(s->correspondsTo(ocv->base()) && ( s->correspondsTo(ocv->crossedBy()) || s->correspondsTo(ocv->before()) )) return 1;
			}
			return 0;
		}
		virtual const char* name() const { return "StackedOrderConstraint"; }
	};

	class TriangleConstraint : public AbacusConstraint {
		const edge base, crossedBy, before;
		const node tri;
	public:
		TriangleConstraint(abacus::Master* m, const edge _base, const node _tri, const edge _crossedBy, const edge _before) :
			AbacusConstraint(m, 0, abacus::CSense::Less, 2, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			base(_base), crossedBy(_crossedBy), before(_before), tri(_tri) {
				OGDF_ASSERT( !crossedBy->isIncident(tri) && !before->isIncident(tri) && !base->commonNode(crossedBy) && !base->commonNode(before) );
			}
		double coeff(const abacus::Variable* cvar) const {
			if(const OrderedCrossingVariable* o = dynamic_cast<const OrderedCrossingVariable*>(cvar)) {
				return (o->base() == base && o->crossedBy() == crossedBy && o->before() == before ) ? 1 : 0; // TODO: check correctness: crossdeBy & before (wothout o->) where: e1 & e2
			} else {
				const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
				for(int i = 1; i-->0;) { // will be unrolled (hopefully)
					if(s->e[i]->isIncident(tri)) {
						if((s->e[1-i] == crossedBy && s->e[i]->isIncident(base->target())) ||
						   (s->e[1-i] == before    && s->e[i]->isIncident(base->source()))) {
							return 1;
						}
						if((s->e[1-i] == before    && s->e[i]->isIncident(base->target())) ||
						   (s->e[1-i] == crossedBy && s->e[i]->isIncident(base->source()))) {
							return -1;
						}
					}
				}
				return (s->correspondsTo(crossedBy) && s->correspondsTo(before)) ? -1 : 0;
			}
		}
		virtual const char* name() const {
			return "TriangleConstraint";
		}
		virtual unsigned hashKey() const {
			return 97*base->index() + 73*crossedBy->index() + 51*before->index() + tri->index();
		}
		virtual bool equal(const abacus::ConVar *cv) const {
			return equals(*(const TriangleConstraint*)cv);
		}
		bool equals(const TriangleConstraint& tc) const {
			return base == tc.base && crossedBy == tc.crossedBy && before == tc.before && tri == tc.tri;
		}
	};

	class SimplicityConstraint : public AbacusConstraint {
		const edge e;
	public:
		SimplicityConstraint(abacus::Master* m, const edge te) :
			AbacusConstraint(m, 0, abacus::CSense::Less, 1, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			e(te) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			return ((s) && s->correspondsTo(e)) ? 1 : 0;
		}
		virtual const char* name() const {
			return "SimplicityConstraint";
		}
	};

	class EdgeOrderConstraint : public AbacusConstraint {
		const edge e1, e2;
	public:
		EdgeOrderConstraint(abacus::Master* m, const edge te1, const edge te2) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, 0, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			e1(te1), e2(te2) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			return (s) ? s->correspondsTo(e1) - s->correspondsTo(e2) : 0;
		}
		virtual const char* name() const {
			return "EdgeOrderConstraint";
		}
	};

	class NodeOrderConstraint : public AbacusConstraint {
		const node n1, n2;
	public:
		NodeOrderConstraint(abacus::Master* m, const node tn1, const node tn2) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, 0, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			n1(tn1), n2(tn2) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			return (s) ? s->correspondsTo(n1) - s->correspondsTo(n2) : 0;
		}
		virtual const char* name() const {
			return "NodeOrderConstraint";
		}
	};

	class KuratowskiConstraintBase : public AbacusConstraint, protected Logger {
	public:
		typedef int Restrictiveness;
		typedef int KuratowskiType;
		friend class Subproblem;
	protected:
		KuratowskiType type;
	public:
		static const KuratowskiType KT_K33/* = 0*/; // its on purpose, that K33 differes from KTBip(3,3) and K5 from KTComp(5) !
		static const KuratowskiType KT_K5 /*= -1*/;
		static KuratowskiType KTComplete(int n) { return -n; }
		static KuratowskiType KTBipartite(int n, int m) { return (n << 16) + m; }
		static int allTypeToCr(KuratowskiType kt) {
			return (kt == KT_K33 || kt == KT_K5) ? 1 : typeToCr(kt);
		}
		static int typeToCr(KuratowskiType kt) {
			OGDF_ASSERT(kt != KT_K33 && kt != KT_K5);
			if(kt<0) return completeGuy(-kt);
			return bipartiteZara(kt >> 16, kt & 0xff);
		}
		int typeToCr() const { return typeToCr(this->type); }
		int allTypeToCr() const { return allTypeToCr(this->type); }
		static int typeToPaths(KuratowskiType kt) {
			//TODO
			OGDF_ASSERT(false); return 0;
		}
		int typeToPaths() const { return typeToPaths(this->type); }
		static Restrictiveness RGreaterPlus(int i) { return i; }
		static Restrictiveness REqualPlus(int i) { return -i-1; }
		static int RPlus(Restrictiveness r) { return (r<0) ? -r-1 : r; }
		static abacus::CSense::SENSE RSense(Restrictiveness r) { return (r<0) ? abacus::CSense::Equal : abacus::CSense::Greater; }
		bool isCompleteOdd() const { return type == KT_K5 || (type < 0 && (-type)%2==1); }

		KuratowskiType subdivision2type(const KuratowskiSubdivision& K) const {
			 return (K.size() == 9) ? KT_K33 : KT_K5;
		}

		KuratowskiConstraintBase(abacus::Master* m, KuratowskiType t, bool dynamic, Restrictiveness rness = 0, bool local = false) :
			AbacusConstraint(m, 0, RSense(rness), allTypeToCr(t)+RPlus(rness), dynamic, local, 1), type(t) {} //pseudo-liftable

		KuratowskiConstraintBase(const KuratowskiConstraintBase& k) : AbacusConstraint(k), Logger(), type(k.type) {};

		virtual KuratowskiConstraintBase* clone() const = 0;

		void branchMe(const Subproblem* S, KuratowskiConstraintBase*& k1, KuratowskiConstraintBase*& k2) const {
			OGDF_ASSERT( sense()->sense() ==  abacus::CSense::Greater );
			OGDF_ASSERT( k1 == NULL && k2 == NULL );
			k1 = clone();
			k2 = clone();
			k1->sub_ = k2->sub_= S;

			k1->local_ = k2->local_ = true;
			k1->sense_ = abacus::CSense::Equal;
			k2->rhs_ = rhs_+ ((type == KT_K5 || (type < 0 && (-type)%2==1))?2:1);
		}

		static bool crossingApplicable(KuratowskiType kt, int p1, int p2) {
			if(p1==p2) return false; // never applicable on the same path
			if( kt == KT_K33) {
				return p1/3 != p2/3 && p1%3 != p2%3;
			}
			if( kt == KT_K5 ) {
				if(p1 > p2) { int t = p1; p1=p2; p2=t; }
				switch(p1) {
				case 0: return p2>=7;
				case 1: return p2==5 || p2==6 || p2==9;
				case 2: return p2==4 || p2==6 || p2==8;
				case 3: return p2==4 || p2==5 || p2==7;
				case 4: return p2==9;
				case 5: return p2==8;
				case 6: return p2==7;
				default: return false;
				}
			}
			return true; // TO DO: gosh... i have no idea on the creation... need to check that more thoroughly
		}
		bool crossingApplicable(int p1, int p2) const { return crossingApplicable(type, p1, p2); }

		virtual bool equalConst(const KuratowskiConstraintBase &cv) const = 0;
	};

	class BasicKuratowskiConstraint : public KuratowskiConstraintBase {
		friend class OrderedOptimalCrossingMinimizer::Subproblem;
	protected:
		Array< Array<edge> > newCrossings;
		Array< CrossingVariableBase* > oldCrossings; // only very primitive constraints!
	public:
		virtual const char* name() const {
			return "BasicKuratowskiConstraint";
		}
		friend std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::BasicKuratowskiConstraint& k);

		BasicKuratowskiConstraint(abacus::Master* m, const KuratowskiSubdivision& K, bool dynamic, Restrictiveness rness = 0 ) :
			KuratowskiConstraintBase(m, subdivision2type(K), dynamic, rness) {
			newCrossings.init(K.size());
			for(int i = K.size(); i-->0;) {
				newCrossings[i].init(K[i].size());
				int j = 0;
				forall_listiterators(edge, it, K[i]) {
					newCrossings[i][j++] = *it;
				}
				newCrossings[i].quicksort(master()->edgeComparer);
			}
			// blank initialization for oldCrossings. filling has to be done by Subproblem.
		}

		BasicKuratowskiConstraint(const BasicKuratowskiConstraint& b) : KuratowskiConstraintBase(b), newCrossings(b.newCrossings), oldCrossings(b.oldCrossings) {}

		virtual KuratowskiConstraintBase* clone() const { return new BasicKuratowskiConstraint(*this); }

		virtual unsigned hashKey() const {
			unsigned h = 0;
			for(int i = newCrossings.size(); i-->0;) {
				for(int j = newCrossings[i].size(); j-->0;) {
					h = (83*h+newCrossings[i][j]->index())%99991;
				}
			}
			for(int j = oldCrossings.size(); j-->0;) {
				CrossingVariableBase* cp = (CrossingVariableBase*)oldCrossings[j];
				h = (79*h+cp->e[0]->index()+3*cp->e[1]->index())%99991;
			}
			return h;
		}

		virtual bool equal(const abacus::ConVar *cv) const {
			const BasicKuratowskiConstraint* c = dynamic_cast<const BasicKuratowskiConstraint*>(cv);
			OGDF_ASSERT( c );
			return equals(*c);
		}

		virtual bool equalConst(const KuratowskiConstraintBase &cv) const {
			const BasicKuratowskiConstraint* c = dynamic_cast<const BasicKuratowskiConstraint*>(&cv);
			if(!c) return false;
			return equals(*c);
		}

		virtual bool equals(const BasicKuratowskiConstraint& ki) const {
			if(newCrossings.size() != ki.newCrossings.size()) return false;
			for(int i = newCrossings.size(); i-->0;) {
				if(newCrossings[i].size() != ki.newCrossings[i].size()) return false;
				for(int j = newCrossings[i].size(); j-->0;) {
					if( newCrossings[i][j] != ki.newCrossings[i][j] ) return false;
				}
			}
			if(oldCrossings.size() != ki.oldCrossings.size()) return false;
			for(int j = oldCrossings.size(); j-->0;) {
				if( oldCrossings[j] != ki.oldCrossings[j] ) return false;
			}
			return true; // no missmatch found
		}

		virtual double coeff(const abacus::Variable* cvar) const {
//			if(oldCrossings.size() == 1) {
//				CrossingVariableBase* c = (CrossingVariableBase*) cvar;
//				CrossingVariableBase* o = oldCrossings[0];
//				if(master()->crossingVariableComparer.compare(o,c) == 0) return -1;
//				return find(cvar) ? 1 : 0;
//			}
			// we have to cast-away const from cvar as the array oldCrossings stores non-const pointers!
			 CrossingVariableComparer::CrossingVariablePointer cv = (CrossingVariableComparer::CrossingVariablePointer)const_cast<abacus::Variable*>(cvar);
//			return (find(cv) ? 1 : 0) - ((oldCrossings.binarySearch(cv,master()->crossingVariableComparer)>=0) ? 1 : 0);
			OGDF_ASSERT( !((oldCrossings.binarySearch(cv,master()->crossingVariableComparer)>=0) && find(cv)) )
			return (oldCrossings.binarySearch(cv,master()->crossingVariableComparer)>=0) ? -1 : (find(cv) ? 1 : 0);
		}
	private:
		bool findInPath(const abacus::Variable* cvar) const {
			if(const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar)) return findInPath(s);
			else return findInPath(dynamic_cast<const OrderedCrossingVariable*>(cvar));
		}
		bool findInPath(const SimpleCrossingVariable* cvar) const {
			return findInPath(cvar->e[0],cvar->e[1]);
		}
		bool findInPath(const OrderedCrossingVariable* cvar) const {
			return findInPath(cvar->base(), cvar->crossedBy()) || findInPath(cvar->base(), cvar->before());
		}
		bool findInPath(edge e, edge f) const {
			for(int i = newCrossings.size(); i-->0;) {
				if(newCrossings[i].binarySearch(e,master()->edgeComparer) >= 0
						&& newCrossings[i].binarySearch(f,master()->edgeComparer) >= 0)
					return true;
			}
			return false;
		}

		bool find(const abacus::Variable* cvar) const {
			int path = -1;
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			if(!s) return false;
			edge e = s->e[0];
			edge f = s->e[1];
			for(int i = newCrossings.size(); i-->0;) {
				if(path < 0 || crossingApplicable(path, i)) {
					if(e && newCrossings[i].binarySearch(e,master()->edgeComparer) >= 0) {// found
						if(f) {
							e = 0;
							path = i;
							if(newCrossings[i].binarySearch(f,master()->edgeComparer) >= 0) // found too -> wrong
								return false;
						} else return true; // okay
					} else if(f && newCrossings[i].binarySearch(f,master()->edgeComparer) >= 0) {
						if(e) {
							f = 0;
							path = i;
						} else return true; // okay
					}
				}
			}
			return false; // nope
		}
	};

	class KnSubdivConstraint : public KuratowskiConstraintBase {
	public:
		KnSubdivConstraint(OrderedOptimalCrossingMinimizer::Master* m, edge tReplace, node tVia, bool dynamic, Restrictiveness rness = 0, bool local = false) :
				KuratowskiConstraintBase(m, KTComplete(m->givenGraph->numberOfNodes()-1), dynamic, rness, local), replace(tReplace), via(tVia) {
		}
		KnSubdivConstraint(const KnSubdivConstraint& b) : KuratowskiConstraintBase(b), replace(b.replace), via(b.via) {}
		virtual KuratowskiConstraintBase* clone() const { return new KnSubdivConstraint(*this); }
	private:
		edge replace;
		node via;
	public:
		virtual unsigned hashKey() const {
			unsigned h = replace->index()*83+via->index();
			return h;
		}

		virtual bool equal(const abacus::ConVar *cv) const {
			const KnSubdivConstraint* c = dynamic_cast<const KnSubdivConstraint*>(cv);
			OGDF_ASSERT( c );
			return equals(*c);
		}

		virtual bool equalConst(const KuratowskiConstraintBase &cv) const {
			const KnSubdivConstraint* c = dynamic_cast<const KnSubdivConstraint*>(&cv);
			if(!c) return false;
			return equals(*c);
		}

		virtual bool equals(const KnSubdivConstraint& ki) const {
			return replace==ki.replace && via==ki.via;
		}

		virtual double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* cv = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			if(!cv) return 0;
			if(cv->correspondsTo(replace) ||
					(cv->correspondsTo(via) && cv->correspondsTo(replace->source()) && cv->correspondsTo(replace->target())))
				return 0;
			return 1;
		}
	};

	class KnMinusConstraint : public KuratowskiConstraintBase {
	public:
		KnMinusConstraint(OrderedOptimalCrossingMinimizer::Master* m, int numTrue, const NodeArray<bool>& aktnodes, bool dynamic, Restrictiveness rness = 0, bool local = false) :
				KuratowskiConstraintBase(m, KTComplete(numTrue), dynamic, rness, local) {
			int i = m->givenGraph->numberOfNodes()-numTrue;
			notin.init(i);
			node n;
			forall_nodes(n,*(m->givenGraph)) {
				if(! aktnodes[n]) notin[--i]=n;
			}
			OGDF_ASSERT( i==0 )
		}
		KnMinusConstraint(OrderedOptimalCrossingMinimizer::Master* m, node notnode, bool dynamic, Restrictiveness rness = 0, bool local = false) :
				KuratowskiConstraintBase(m, KTComplete(m->givenGraph->numberOfNodes()-1), dynamic, rness, local) {
			notin.init(1);
			notin[0] = notnode;
		}
		KnMinusConstraint(OrderedOptimalCrossingMinimizer::Master* m, node notnode1, node notnode2, bool dynamic, Restrictiveness rness = 0, bool local = false) :
				KuratowskiConstraintBase(m, KTComplete(m->givenGraph->numberOfNodes()-1), dynamic, rness, local) {
			notin.init(2);
			notin[0] = notnode1;
			notin[1] = notnode2;
		}
		KnMinusConstraint(const KnMinusConstraint& b) : KuratowskiConstraintBase(b), notin(b.notin) {}
		virtual KuratowskiConstraintBase* clone() const { return new KnMinusConstraint(*this); }
	private:
		Array<node> notin;
	public:
		virtual unsigned hashKey() const {
			unsigned h = 0;
			for(int j = notin.size(); j-->0;) {
				h = (83*h+notin[j]->index())%99991;
			}
			return h;
		}

		virtual bool equal(const abacus::ConVar *cv) const {
			const KnMinusConstraint* c = dynamic_cast<const KnMinusConstraint*>(cv);
			OGDF_ASSERT( c );
			return equals(*c);
		}

		virtual bool equalConst(const KuratowskiConstraintBase& cv) const {
			const KnMinusConstraint* c = dynamic_cast<const KnMinusConstraint*>(&cv);
			if(!c) return false;
			return equals(*c);
		}

		virtual bool equals(const KnMinusConstraint& ki) const {
			if(notin.size() != ki.notin.size()) return false;
			for(int j = notin.size(); j-->0;) {
				if( notin[j] != ki.notin[j] ) return false;
			}
			return true; // no missmatch found
		}

		virtual double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* cv = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			if(!cv) return 0;
			for(int j = notin.size(); j-->0;) {
				if(cv->correspondsTo(notin[j])) return 0;
			}
			return 1;
		}
	};

	class HypercubeConstraint : public AbacusConstraint {
		const int bit;
		const int elem;
	public:
		HypercubeConstraint(abacus::Master* m, const int tbit, bool isset, const int trhs) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, trhs, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			bit(tbit), elem(isset ? bit : 0) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			return (s
					&& ((s->e[0]->source()->index() & bit) == elem)
					&& ((s->e[1]->source()->index() & bit) == elem)
					&& ((s->e[0]->target()->index() & bit) == elem)
					&& ((s->e[1]->target()->index() & bit) == elem)
					) ? 1: 0;
		}
		virtual const char* name() const {
			return "HypercubeConstraint";
		}
	};

	class ToroidalSubgridConstraint : public AbacusConstraint {
		const int f,n;
		const bool div;
	public:
		ToroidalSubgridConstraint(abacus::Master* m, const int tf, const int tn, const bool tdiv, const int trhs) :
			AbacusConstraint(m, 0, abacus::CSense::Greater, trhs, 0, 0, 1),  // non-dynamic, non-local, pseudo-liftable
			f(tf), n(tn) , div(tdiv) {}
		double coeff(const abacus::Variable* cvar) const {
			const SimpleCrossingVariable* s = dynamic_cast<const SimpleCrossingVariable*>(cvar);
			return (!s || isRemoved(s->e[0]) || isRemoved(s->e[1])) ? 0 : 1;
		}
		bool isRemoved(edge e) const {
			int u = e->source()->index();
			int v = e->target()->index();
			if(div) {
				return ((u/n==f)&&(v/n==f)) ? true : false;
			} else {
				return ((u%n==f)&&(v%n==f)) ? true : false;
			}
		}
		virtual const char* name() const {
			return "ToroidalSubgridConstraint";
		}
	};

//	class SingletonCrossingsKuratowskiConstraint : public KuratowskiConstraint {
//	protected:
//		Array<>
//		Array< CrossingVariablePointer > oldCrossings;
//	public:
//		virtual const char* name() /*const*/ {
//			return "SingletonCrossingsKuratowskiConstraint";
//		}
//		friend std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::SingletonCrossingsKuratowskiConstraint& k);
//
//		SingletonCrossingsKuratowskiConstraint(abacus::Master* m, KuratowskiSubdivision& K, bool dynamic, Restrictiveness rness) :
//				KuratowskiConstraint(m, dynamic, rness) { // K \subset givenGraph
//			...
//		}
//
//		SingletonCrossingsKuratowskiConstraint(const SingletonCrossingsKuratowskiConstraint& b) :
//				SingletonCrossingsKuratowskiConstraint(b), ... {
//			...
//		}
//
//		virtual KuratowskiConstraintBase* clone() const { return new SingletonCrossingsKuratowskiConstraint(*this); }
//
//		virtual unsigned hashKey() /*const*/ {
//			unsigned h = 0;
//			...
//			return h;
//		}
//
//		virtual bool equal(abacus::ConVar *cv) /*const*/ {
//			return equals(*(SingletonCrossingsKuratowskiConstraint*)cv);
//		}
//
//		bool equals(const SingletonCrossingsKuratowskiConstraint& ki) const {
//			...
//		}
//
//		double coeff(const abacus::Variable* cvar) const {
//			...
//			return 0; // nope
//		}
//	};
//
//	class OrderedCrossingsKuratowskiConstraint : public KuratowskiConstraint {
//	private:
//		...
//	public:
//		virtual const char* name() /*const*/ {
//			return "OrderedCrossingsKuratowskiConstraint";
//		}
//		friend std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::OrderedCrossingsKuratowskiConstraint& k);
//
//		OrderedCrossingsKuratowskiConstraint(abacus::Master* m, KuratowskiSubdivision& K, bool dynamic, Restrictiveness rness) :
//				KuratowskiConstraintBase(m, dynamic, rness) { // K \subset givenGraph
//			...
//		}
//
//		OrderedCrossingsKuratowskiConstraint(const SingletonCrossingsKuratowskiConstraint& b) :
//				SingletonCrossingsKuratowskiConstraint(b), ... {
//			...
//		}
//
//		virtual KuratowskiConstraintBase* clone() const { return new OrderedCrossingsKuratowskiConstraint(*this); }
//
//		virtual unsigned hashKey() /*const*/ {
//			unsigned h = 0;
//			...
//			return h;
//		}
//
//		virtual bool equal(abacus::ConVar *cv) /*const*/ {
//			return equals(*(OrderedCrossingsKuratowskiConstraint*)cv);
//		}
//
//		bool equals(const OrderedCrossingsKuratowskiConstraint& ki) const {
//			...
//		}
//
//		double coeff(const abacus::Variable* cvar) const {
//			...
//			return 0; // nope
//		}
//	};

/*	template<class C, class COMP>
	class KillingTop10Heap : public Top10Heap<Prioritized<C*> > {
	public:
		KillingTop10Heap(int size) : Top10Heap<Prioritized<C*> >(size) {}
		void pushAndKill(double val, C* x) {
			Prioritized<C*> vo;
			Prioritized<C*> nv(x, val);
			if(returnedSomething( Top10Heap<Prioritized<C*> >::push(nv, vo) ))
				delete vo.item();
		}
		void pushAndKillNoRedundancy(double val, C* x) {
			for(int i = Top10Heap<Prioritized<C*> >::size(); i-->0;) {
				C* k = Top10Heap<Prioritized<C*> >::operator[](i).item();
				OGDF_ASSERT( x )
				OGDF_ASSERT( k )
				if(COMP::compare(k,x)) {
					delete x;
					return;
				}
			}
			pushAndKill(val, x);
		}
	};*/

	struct Nunchaku {
		Nunchaku() : a(0), b(0) {} // null constructor to allow use in arrays, arraybuffers, etc.
		Nunchaku(const edge x, const edge y) : a(x), b(y) { OGDF_ASSERT( x->commonNode(y) ); }
		Nunchaku(const Nunchaku& n) : a(n.a), b(n.b) {}
		Nunchaku& operator=(const Nunchaku& n) { a=n.a; b=n.b; return *(this); }
		edge a; // should be const, but then the = operator won't work anymore...
		edge b; // same here
	};

	class Subproblem : public abacus::Sub, protected Logger {

		void clearInduced() {
			SAFE_DELETE(inducedPlanarization[I_INT]);
			SAFE_DELETE(inducedPlanarization[I_RND]);
			inducedSeparationSimpleCrossings.clear();
		}

	public:
		Subproblem(abacus::Master* m, abacus::Sub* father, abacus::BranchRule* br) :
				abacus::Sub(m, father, br),
				realIter(0),
				generatedConVar(0),
				storedCurrentCrossings(false),
				bufferedLinearOrders() {
			lout() << "[Branch] Subproblem " << id() << " generated\n"; // (father: " << father->id() << ")
			inducedPlanarization[I_RND] = 0;
			inducedPlanarization[I_INT] = 0;
		}
		Subproblem(abacus::Master* m,
				double conRes,
				double varRes,
				double nnzRes,
				bool relativeRes,
				ArrayBuffer<abacus::PoolSlot<abacus::Constraint, abacus::Variable> *> *constraints = 0,
				ArrayBuffer<abacus::PoolSlot<abacus::Variable, abacus::Constraint> *> *variables = 0) :
					abacus::Sub(m, conRes, varRes, nnzRes, relativeRes, constraints, variables),
					realIter(0),
					generatedConVar(0),
					storedCurrentCrossings(false),
					bufferedLinearOrders() {
			lout() << "Initial Subproblem " << id() << " generated (no father)\n";
			inducedPlanarization[I_RND] = 0;
			inducedPlanarization[I_INT] = 0;
		}
		~Subproblem() {
			clearInduced();
			OGDF_ASSERT( bufferedLinearOrders.empty() );
		};
		abacus::Sub *generateSon(abacus::BranchRule* br) {
			if(master()->nSub()==1) { // i'm the root subproblem, and i'm gonna die -> let's give a status report...
				sout() << "\t" << ceil(this->dualBound()-0.5) // lower bound
					   << "\t" << master()->primalBound() // upper bound
					   << "\t" << master()->totalTime()->seconds();
			}
			clearInduced();
			OGDF_ASSERT( bufferedLinearOrders.empty() );
			return new Subproblem(master(), this, br);
		}
		bool feasible();
		int separate();
		int pricing();
		int makeFeasible();
		int improve(double& d);
		int solveLp();
		int generateBranchRules(ArrayBuffer<abacus::BranchRule*> &rules);

		int separateLinearOrder();
		int separateTriangles();
		int separateBoyerMyrvold();
		int separateSimple();

		KuratowskiConstraintBase* constructKuratowskiConstraint(const PlanRep& G, const KuratowskiSubdivision& K);

		int checkAmbiguitiesAndPrice();
		void createVarsAndCorrespondingLOCs(edge e, edge f, edge g, ArrayBuffer<OrderedCrossingVariable*>& vars);
		double realizeSolutions();
		void primalHeuristics(int inducement);
		void orderCrossings(const edge e, List<edge>& L);
		double orderingWeight(const edge e, const edge f, const ListIterator<edge>& l, const ListIterator<edge>& r);

		int checkKnHighKuratowskiCutsStatic();
		int checkKnmHighKuratowskiCutsStatic();

	protected:

		int realIter;

		PlanRep* inducedPlanarization[INDUCEMENTS];
		List< SimpleCrossingVariable* > inducedSeparationSimpleCrossings;

		int generatedConVar;
		bool storedCurrentCrossings;
		SList<AbacusConstraint*> bufferedLinearOrders;

		OrderedOptimalCrossingMinimizer::Master* master() const {
			return (OrderedOptimalCrossingMinimizer::Master*)master_;
		}

		void findKuratowski(Graph& R, KuratowskiSubdivision& K);
		OrderedOptimalCrossingMinimizer::CrossingConfiguration* callBoundHeuristic(int integer);
//		inline double generateVariable(edge e1, int e1id, edge e2, int e2id, OrderedOptimalCrossingMinimizer::CrossingVariable** var);

	};

	friend std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::SimpleCrossingVariable& v);
	friend std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::OrderedCrossingVariable& v);
	friend std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::BasicKuratowskiConstraint& k);

	class Master : /*public CrossingMinimizationModule,*/ public abacus::Master, protected Logger {

		static int MAX_PAIRS(int e) {
			return (e*e - e) / 2;
		}

		class ActiveVariables {
			struct VarEntry {
				int round;
				int index;
				VarEntry() : round(-1), index(-1) {}
			};
			struct MultiEntry : VarEntry {
				EdgeArray<VarEntry*>* before;
				MultiEntry() : VarEntry() { before = 0; }
				~MultiEntry() {
					edge e;
					if(before) {
						forall_edges(e, *(before->graphOf())) {
							if((*before)[e]) delete (*before)[e];
						}
						delete before;
					}
				}
				void prepareBefore(const Graph& G) {
					OGDF_ASSERT( !before );
					before = new EdgeArray<VarEntry*>(G);
					edge e;
					forall_edges(e, G) {
						(*before)[e] = 0;
					}
				}
				void registerBefore(const edge e) {
					OGDF_ASSERT( before );
					OGDF_ASSERT( !(*before)[e] );
					(*before)[e] = new VarEntry();
				}
			};
			OrderedOptimalCrossingMinimizer::Master* master;
			EdgeArray<int> count;
		public: // TODO: private
			EdgeArray< EdgeArray< MultiEntry > > vars;
		private:
			int currentRound;
		public:
			ActiveVariables(OrderedOptimalCrossingMinimizer::Master* m) :
					master(m), count(*(m->givenGraph),0), vars(*(m->givenGraph)), currentRound(2) {
				edge e;
				const Graph& mG = *(master->givenGraph);
				forall_edges(e, mG) {
					vars[e].init(mG); // no init necessary: default constructed IS 0-init
				}
			}
			~ActiveVariables() {}

			inline int numOrdered(const edge e) const {
				return count[e];
			}

			void registerNewOrdered(const edge e, const edge f) {
				edge g;
				++count[e];
				vars[e][f].prepareBefore(*(master->givenGraph));
				forall_edges(g, *(master->givenGraph)) {
					if(f == g || e == g) continue;
					if(vars[e][g].before) {
						vars[e][g].registerBefore(f);
						vars[e][f].registerBefore(g);
					}
				}
			}

			inline bool exists(const edge e, const edge f, const edge g) const {
				return (vars[e][f].before != 0) && (*(vars[e][f].before))[g];
			}
			inline bool isOrdered(const edge e, const edge f) const {
				return (vars[e][f].before != 0);
			}
			inline int indexSafe(const edge e, const edge f, const edge g) const {
				return exists(e,f,g) ? index(e,f,g) : -2;
			}
			inline int index(const edge e, const edge f, const edge g) const {
				//return ((*(vars[e][f].before))[g]->round == currentRound) ? (*(vars[e][f].before))[g]->index : -2;
				if( (*(vars[e][f].before))[g]->round == currentRound )
					return  (*(vars[e][f].before))[g]->index;
				else {
					cerr << e->index() << "," << f->index() << "," << g->index() << ": sr="<< (*(vars[e][f].before))[g]->round << " cr=" << currentRound << "\n";
					return -2;
				}
			}
			inline int roundlessindex(const edge e, const edge f, const edge g) const {
				return (*(vars[e][f].before))[g]->index;
			}
			inline int index(const edge e, const edge f) const {
				return (vars[e][f].round == currentRound) ? vars[e][f].index : -2;
			}

			int nextVar(const edge e, const edge f, edge& g) const {
				const EdgeArray<VarEntry*>& A= *(vars[e][f].before);
				g = g ? g->succ() : master->givenGraph->firstEdge();
				int r = -2;
				while(g && (!A[g] || (r=index(e,f,g))<0)) { g=g->succ(); }
				return r;
			}

			void loadIndices(OrderedOptimalCrossingMinimizer::Subproblem* sub) { // TODO: speed this up, if indices are stable, i.e., if only new variables have to be indexed
				if(++currentRound >> 14) currentRound = 5;
				//cout << "\t\tCURRENT ROUND " << currentRound << endl;
				OrderedCrossingVariable* o;
				SimpleCrossingVariable* s;
				for(int i = 0; i < sub->nVar(); ++i) {
					if(( o = dynamic_cast<OrderedCrossingVariable*>(sub->variable(i)) )) {
						OGDF_ASSERT( vars[o->base()][o->crossedBy()].before );
						OGDF_ASSERT( (*(vars[o->base()][o->crossedBy()].before))[o->before()] );

						OGDF_ASSERT( (*(vars[o->base()][o->crossedBy()].before))[o->before()]->index == -1 || (*(vars[o->base()][o->crossedBy()].before))[o->before()]->index == i );

						(*(vars[o->base()][o->crossedBy()].before))[o->before()]->round = currentRound;
						(*(vars[o->base()][o->crossedBy()].before))[o->before()]->index = i;
//						cout << "\t\t" << i << ": " << o->base() << "x" << o->crossedBy() << "b" << o->before() << " [" << currentRound << "]\n";
					} else {
						s = dynamic_cast<SimpleCrossingVariable*>(sub->variable(i));
						OGDF_ASSERT( vars[s->e[0]][s->e[1]].index == -1 || vars[s->e[0]][s->e[1]].index == i );
						OGDF_ASSERT( vars[s->e[1]][s->e[0]].index == -1 || vars[s->e[1]][s->e[0]].index == i );
						vars[s->e[0]][s->e[1]].round = vars[s->e[1]][s->e[0]].round = currentRound;
						vars[s->e[0]][s->e[1]].index = vars[s->e[1]][s->e[0]].index = i;
//						cout << "\t\t" << i << ": " << s->e[0] << "x" <<s->e[1] << " [" << currentRound << "]\n";
					}
				}
			}
		};

		friend class ActiveVariables;
		friend class Subproblem;
		friend class SimpleCrossingVariable;
		friend class OrderedCrossingVariable;
		friend class KuratowskiConstraintBase;
		friend class BasicKuratowskiConstraint;
		friend class KnMinusConstraint;
		friend class KnSubdivConstraint;
		friend class OrderedOptimalCrossingMinimizer;

		void setDefaultSettings();

	public:

		Master() :
				abacus::Master("OrderedOptimalCrossingMinimizer::Master",true, true, abacus::OptSense::Min),
				bestSolution(0),
				m_isTrivial(false) {
			activeVars = 0;
			linearOrderPool = 0;
			trianglePool = 0;
			kuratowskiPool = 0;
			hintedPool = 0;
			branchingPool = 0;
			setDefaultSettings();
			for(int i = INDUCEMENTS; i-->0;)
				inducedCrossingConfiguration[i] = 0; //new CrossingConfiguration();
		}
		~Master() {
		}

	private:
		void clearAfterCall() {
			for(int i = INDUCEMENTS; i-->0;)
				SAFE_DELETE(inducedCrossingConfiguration[i]);
			SAFE_DELETE(bestSolution);
			SAFE_DELETE(activeVars);
			SAFE_DELETE(linearOrderPool);
			SAFE_DELETE(trianglePool);
			SAFE_DELETE(kuratowskiPool);
			SAFE_DELETE(hintedPool);
			SAFE_DELETE(branchingPool);
			nunchakus.init();
		}

	public:
		abacus::Sub* firstSub() {
			int maxVars = (givenGraph->numberOfEdges()+1) * givenGraph->numberOfEdges() * givenGraph->numberOfEdges();
			const int myMaxVars = 5000000;
			if(maxVars > myMaxVars) {
				lout() << "Warning: restricting maxVars to " << myMaxVars << " instead of " << maxVars << ".\n";
				maxVars = myMaxVars;
			}
			return new OrderedOptimalCrossingMinimizer::Subproblem(this, 1000, maxVars, upperbound, false);
		}

		void initializeOptimization();

		int enumerationStrategy(const abacus::Sub* s1, const abacus::Sub* s2);

		CrossingVariableComparer crossingVariableComparer;
		EdgeComparer edgeComparer;

		CrossingConfiguration* inducedCrossingConfiguration[INDUCEMENTS];

	public:
		void useThisUpperBound(int i) { m_useThisUpperBound = i; }

		int numStartHeuristics() const { return m_numStartHeuristics; }
		void numStartHeuristics(int i) { m_numStartHeuristics = i; }

		void setStartHeuristic(CrossingMinimizationModule* p) { m_startHeuristic.set(p); }
		void setBoundHeuristic(CrossingMinimizationModule* p) { m_boundHeuristic.set(p); }

		PricingInit pricingInit() const { return m_pricingInit; }
		void pricingInit(PricingInit p) { m_pricingInit = p; }

		int separationMode() const { return m_separationMode; }
		void separationMode(int i) { m_separationMode=i; }

		BranchingMode branchingMode() const { return m_branchingMode; }
		void branchingMode(BranchingMode p) { m_branchingMode = p; }

		const BoyerMyrvoldSeparationParams& boyerMyrvoldSeparationParams() const { return m_boyerMyrvoldSeparationParams; }
		BoyerMyrvoldSeparationParams& boyerMyrvoldSeparationParams() { return m_boyerMyrvoldSeparationParams; }

		const SimpleSeparationParams& simpleSeparationParams() const { return m_simpleSeparationParams; }
		SimpleSeparationParams& simpleSeparationParams() { return m_simpleSeparationParams; }

		int maxTriangleCuts() const { return m_maxTriangleCuts; }
		void maxTriangleCuts(int i) { m_maxTriangleCuts = i; }

		int maxLinearOrderCuts() const { return m_maxLinearOrderCuts; }
		void maxLinearOrderCuts(int i) { m_maxLinearOrderCuts = i; }

		int maxNewVars() const { return m_maxNewVars; }
		void maxNewVars(int i) { m_maxNewVars = i; }

		int numCutHighKuratowskis() const { return m_numCutHighKuratowskis; }
		void numCutHighKuratowskis(int i) { m_numCutHighKuratowskis = i; }

		int numBaseEdgesForCutHighKuratowskis() const { return m_numBaseEdgesForCutHighKuratowskis; }
		void numBaseEdgesForCutHighKuratowskis(int i) { m_numBaseEdgesForCutHighKuratowskis = i; }

		int maxMinutes() const { return m_maxMinutes; }
		void maxMinutes(int m) { m_maxMinutes = m; }

		double roundUp() const { return m_roundUp; }
		void roundUp(double r) { m_roundUp = r; }

		bool isOptimal() { return m_isTrivial || (status() == Optimal); }

		char* writeResult() { return m_writeResult; }
		void writeResult(char* b) { m_writeResult = b; }

		bool writeIntermediateResultsToo() const { return m_writeIntermediateResultsToo; }
		void writeIntermediateResultsToo(bool b) { m_writeIntermediateResultsToo = b; }

		GraphHint graphHint() const { return m_graphHint; }
		void graphHint(GraphHint h) { m_graphHint = h; }

		int hintEffects() const { return m_hintEffects; }
		void hintEffects(int h) { m_hintEffects = h; }

		bool reduceMemory() const { return m_reduceMemory; }
		void reduceMemory(bool b) { m_reduceMemory = b; }

		bool localVariables() const { return m_localVariables; }
		void localVariables(bool b) { m_localVariables = b; }

	protected:
		virtual ReturnType doCall(PlanRep &PG,
			int cc,
			const EdgeArray<int>  *cost,
			const EdgeArray<bool> *forbid,
			int& crossingNumber);

	private:

		OrderedOptimalCrossingMinimizer::CrossingConfiguration* bestSolution;
		ActiveVariables* activeVars;

		bool m_isTrivial;

		int m_minVariables;
		int m_maxVariables;
		int m_startVariables;
		int m_usedVariables;

		void updateBestSolution(OrderedOptimalCrossingMinimizer::CrossingConfiguration* ncc, bool heur) {
			if(bestSolution == ncc) return; // repetition -> init
			upperBoundSource = heur ? OrderedOptimalCrossingMinimizer::Master::SS_ILP_Heuristic : OrderedOptimalCrossingMinimizer::Master::SS_ILP;
			if(bestSolution) delete bestSolution;
			bestSolution = ncc;
			if(writeIntermediateResultsToo()) doWriteBestSolution();
		}

		int getCost(const edge e1, const edge e2) const { // call with e1, e2 of givenGraph
			return useCost() ? (*cost)[e1]*(*cost)[e2] : 1;
		}
		bool variableAllowed(edge e1, edge e2);

		abacus::StandardPool<abacus::Constraint, abacus::Variable>* linearOrderPool;
		abacus::StandardPool<abacus::Constraint, abacus::Variable>* trianglePool; // separation will only add constraints that are not already in there (otherwise they would have been separated already, and would be in the active ILP, leading to a non-addition)
		abacus::StandardPool<abacus::Constraint, abacus::Variable>* hintedPool;
		abacus::NonDuplPool<abacus::Constraint, abacus::Variable>* kuratowskiPool; // is NONDUP neccessary? yes, in a single iteration, we might add the same thing multiple times -> TODO improve on that, by doing the DUP-check by myself, within that iteration, and add the stuff to a STD instead of DUP afterwards...
		abacus::StandardPool<abacus::Constraint, abacus::Variable>* branchingPool;

		void hintsEdgeOrder(List<abacus::Constraint*>& prelist); // Kn, Qn, Pnm, Tnm
		void hintsSimplicity(List<abacus::Constraint*>& prelist); // Kn, Qn, Pnm, Tnm

		void hintsKnKuratowskiMinusOne(List<abacus::Constraint*>& prelist);
		void hintsKnAllSubKuratowskis(List<abacus::Constraint*>& prelist);
		void hintsKnNodeOrder(List<abacus::Constraint*>& prelist);
		void helperHintsKnAllSubKuratowskis(NodeArray<bool>& aktnodes, node posNode, int num, List<abacus::Constraint*>& prelist);
		void hintsKnExpensiveKuratowski(List<abacus::Constraint*>& prelist);

		void hintsKnmKuratowskiMinusOne(List<abacus::Constraint*>& prelist);
		void hintsKnmAllSubKuratowskis(List<abacus::Constraint*>& prelist);
		void hintsKnmNodeOrder(List<abacus::Constraint*>& prelist);
		void hintsKnmEdgeOrder(List<abacus::Constraint*>& prelist);

		void hintsQnHypercubeMinusOne(List<abacus::Constraint*>& prelist);
		void hintsTnmNodeOrder(List<abacus::Constraint*>& prelist);
		void hintsTnmSubgrids(List<abacus::Constraint*>& prelist);

		void hintsPnmNodeOrder(List<abacus::Constraint*>& prelist);


		EdgeArray< Array<Nunchaku> > nunchakus;
		void initNunchakus();

		enum SolutionSource { SS_Trivial, SS_ILP, SS_ILP_Heuristic, SS_Heuristic, SS_Kn, SS_Knm, SS_NoSolution };

		const Graph*           givenGraph;
		PlanRep*               resultingGraph;
		const EdgeArray<int>*  cost;
		const EdgeArray<bool>* forbid;

		int m_useThisUpperBound;
		int m_numStartHeuristics;
		ModuleOption<CrossingMinimizationModule> m_startHeuristic;
		ModuleOption<CrossingMinimizationModule> m_boundHeuristic;
//		SubgraphPlanarizer m_solutionChecker;

		int m_separationMode;
		PricingInit m_pricingInit;
		BranchingMode m_branchingMode;
		BoyerMyrvoldSeparationParams m_boyerMyrvoldSeparationParams;
		SimpleSeparationParams m_simpleSeparationParams;

		int m_maxLinearOrderCuts;
		int m_maxTriangleCuts;
		int m_maxNewVars;
		int m_numCutHighKuratowskis;
		int m_numBaseEdgesForCutHighKuratowskis;

		int m_maxMinutes;
		double m_roundUp;

		GraphHint m_graphHint;
		int m_hintEffects;

		char* m_writeResult;
		bool m_writeIntermediateResultsToo;

//		int numMinNodes;
//		int numMinEdges;
//		int numMinMaxCrossingPairs;
//		int numExpNodes;
//		int numExpEdges;
//		int numExpMaxCrossingPairs;

		int lowerbound;
		int upperbound;
		SolutionSource upperBoundSource;

		bool m_reduceMemory;
		bool m_localVariables;

		bool probablyUpdateLowerBound(int lb) {
			if(lb <= lowerbound) return false;
			lout(LL_MINOR) << "New Lower Bound! " << lb << "\n";
			lowerbound = lb;
			return true;
		}
		bool probablyUpdateUpperBound(int ub, SolutionSource ubs) {
			--ub; // only something better than that is interesting...
			if(ub >= upperbound) return false;
			lout(/*LL_MINOR*/) << "New Upper Bound: " << (ub+1) << " [" << upperBoundSource << "->" << ubs << "]\n";
			upperbound = ub;
			upperBoundSource = ubs;
			return true;
		}

		void calcLowerBounds();
		void calcUpperBounds();

		//! Reminder: OptimalCrossingMinimizerBase assumes that the given graph is connected
		OrderedOptimalCrossingMinimizer::CrossingConfiguration* initBounds();
		OrderedOptimalCrossingMinimizer::CrossingConfiguration* createHeuristicStartSolution();

		void doWriteBestSolution();

		bool m_useCost;
		bool m_useForbid;
		bool useCost() const { return m_useCost; }
		bool useForbid() const { return m_useForbid; }
	}; // end Master

protected:
	Master blaster; //!< delegate
	friend class Master;

	virtual ReturnType doCall(PlanRep &PG,
			int cc,
			const EdgeArray<int>  *cost,
			const EdgeArray<bool> *forbid,
			const EdgeArray<__uint32>  *subgraphs,
			int& crossingNumber) {
		blaster.m_useCost = (cost != 0);
		blaster.m_useForbid = (forbid != 0);
		if( subgraphs != 0 ) throw 4711;
		return blaster.doCall(PG, cc, cost, forbid, crossingNumber);
	}

public:
	friend class DeletingTop10Heap<OrderedOptimalCrossingMinimizer::CyclicOrderConstraint,double,Compare_Equals<OrderedOptimalCrossingMinimizer::CyclicOrderConstraint> >;
	friend class DeletingTop10Heap<OrderedOptimalCrossingMinimizer::TriangleConstraint,double,Compare_Equals<OrderedOptimalCrossingMinimizer::TriangleConstraint> >;
	friend class DeletingTop10Heap<OrderedOptimalCrossingMinimizer::KnSubdivConstraint,double,Compare_Equals<OrderedOptimalCrossingMinimizer::KnSubdivConstraint> >;
	friend class DeletingTop10Heap<OrderedOptimalCrossingMinimizer::KuratowskiConstraintBase,double,Compare_CheckedEqual<OrderedOptimalCrossingMinimizer::KuratowskiConstraintBase> >;

	OrderedOptimalCrossingMinimizer() : blaster() {}

	CrossingMinimizationModule *clone() const;

	void useThisUpperBound(int i) { blaster.useThisUpperBound(i); }

	int numStartHeuristics() const { return blaster.numStartHeuristics(); }
	void numStartHeuristics(int i) { blaster.numStartHeuristics(i); }

	void setStartHeuristic(CrossingMinimizationModule* p) { blaster.setStartHeuristic(p); }
	void setBoundHeuristic(CrossingMinimizationModule* p) { blaster.setBoundHeuristic(p); }

	PricingInit pricingInit() const { return blaster.pricingInit(); }
	void pricingInit(PricingInit p) { blaster.pricingInit(p); }

	int separationMode() const { return blaster.separationMode(); }
	void separationMode(int i) { blaster.separationMode(i); }

	BranchingMode branchingMode() const { return blaster.branchingMode(); }
	void branchingMode(BranchingMode p) { blaster.branchingMode(p); }

//	const BoyerMyrvoldSeparationParams& startBoyerMyrvoldSeparationParams() const { return blaster.startBoyerMyrvoldSeparationParams(); }
//	BoyerMyrvoldSeparationParams& startBoyerMyrvoldSeparationParams() { return blaster.startBoyerMyrvoldSeparationParams(); }
//
//	const SimpleSeparationParams& startSimpleSeparationParams() const { return blaster.startSimpleSeparationParams(); }
//	SimpleSeparationParams& startSimpleSeparationParams() { return blaster.startSimpleSeparationParams(); }

	const BoyerMyrvoldSeparationParams& boyerMyrvoldSeparationParams() const { return blaster.boyerMyrvoldSeparationParams(); }
	BoyerMyrvoldSeparationParams& boyerMyrvoldSeparationParams() { return blaster.boyerMyrvoldSeparationParams(); }

	const SimpleSeparationParams& simpleSeparationParams() const { return blaster.simpleSeparationParams(); }
	SimpleSeparationParams& simpleSeparationParams() { return blaster.simpleSeparationParams(); }

	int maxTriangleCuts() const { return blaster.maxTriangleCuts(); }
	void maxTriangleCuts(int i) { blaster.maxTriangleCuts(i); }

	int maxLinearOrderCuts() const { return blaster.maxLinearOrderCuts(); }
	void maxLinearOrderCuts(int i) { blaster.maxLinearOrderCuts(i); }

	int maxNewVars() const { return blaster.maxNewVars(); }
	void maxNewVars(int i) { blaster.maxNewVars(i); }

	int numCutHighKuratowskis() const { return blaster.numCutHighKuratowskis(); }
	void numCutHighKuratowskis(int i) { blaster.numCutHighKuratowskis(i); }

	int numBaseEdgesForCutHighKuratowskis() const { return blaster.numBaseEdgesForCutHighKuratowskis(); }
	void numBaseEdgesForCutHighKuratowskis(int i) { blaster.numBaseEdgesForCutHighKuratowskis(i); }

	int maxMinutes() const { return blaster.maxMinutes(); }
	void maxMinutes(int m) { blaster.maxMinutes(m); }

	double roundUp() const { return blaster.roundUp(); }
	void roundUp(double r) { blaster.roundUp(r); }

	bool isOptimal() { return blaster.isOptimal(); }

	char* writeResult() { return blaster.writeResult(); }
	void writeResult(char* b) { blaster.writeResult(b); }

	bool writeIntermediateResultsToo() const { return blaster.writeIntermediateResultsToo(); }
	void writeIntermediateResultsToo(bool b) { blaster.writeIntermediateResultsToo(b); }

	GraphHint graphHint() const { return blaster.graphHint(); }
	void graphHint(GraphHint h) { blaster.graphHint(h); }

	int hintEffects() const { return blaster.hintEffects(); }
	void hintEffects(int h) { blaster.hintEffects(h); }

	int tailOffNLp() const { return blaster.tailOffNLp(); }
	void tailOffNLp(int n) { blaster.tailOffNLp(n); }

	double tailOffPercent() const { return blaster.tailOffPercent(); }
	void tailOffPercent(double p) { blaster.tailOffPercent(p); }

	bool reduceMemory() const { return blaster.reduceMemory(); }
	void reduceMemory(bool b) { blaster.reduceMemory(b); }

	bool localVariables() const { return blaster.localVariables(); }
	void localVariables(bool b) { blaster.localVariables(b); }
};

std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::SimpleCrossingVariable& v);
std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::OrderedCrossingVariable& v);
std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::BasicKuratowskiConstraint& k);
//std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::CyclicOrderConstraint& k);

//std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::SimplicityConstraint& k);
//std::ostream &operator<<(std::ostream &os, /*const*/ OrderedOptimalCrossingMinimizer::KuratowskiConstraint& k);

std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::SimpleSeparationParams& p);
std::ostream &operator<<(std::ostream &os, const OrderedOptimalCrossingMinimizer::BoyerMyrvoldSeparationParams& p);

#endif // USE_ABACUS

} // namespace

#endif // OGDF_ORDERED_OPTIMAL_CROSSING_MINIMIZER_H
