//**************************************************************
//  tiny "regression test" for Steiner approximation algorithms
//
//  Tested:
//    - GraphIO::readSTP()
//    - MinSteinerTreeKou
//    - MinSteinerTreeMehlhorn
//    - MinSteinerTreeTakahashi
//    - MinSteinerTreeRZLoss
//    - MinSteinerTreeZelikovsky in all its variants
//
//  Author: Stephan Beyer
//**************************************************************

#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/graphalg/MinSteinerTreeKou.h>
#include <ogdf/graphalg/MinSteinerTreeMehlhorn.h>
#include <ogdf/graphalg/MinSteinerTreeTakahashi.h>
#include <ogdf/graphalg/MinSteinerTreeRZLoss.h>
#include <ogdf/graphalg/MinSteinerTreeZelikovsky.h>

using namespace ogdf;

template <typename T>
static T
computeMinStT(MinSteinerTreeModule<T> *&alg, const string &desc, const EdgeWeightedGraph<T> &wG, const List<node> &terminals, const NodeArray<bool> &isTerminal)
{
	__int64 time, used_time;
	EdgeWeightedGraphCopy<T> *finalSteinerTree = NULL;

	System::usedRealTime(time);
	T obj = alg->call(wG, terminals, isTerminal, finalSteinerTree);
	used_time = System::usedRealTime(time);

	delete finalSteinerTree;

	cout << "  " << desc;
	for (string::size_type w = desc.size(); w < 36; ++w) {
		cout << " ";
	}
	cout << "result: " << obj << "\tin " << 0.001*double(used_time) << " sec.\n";

	return obj;
}

template <typename T>
static void
zelikovskyVariant(MinSteinerTreeModule<T> *&alg, string &desc, int id)
{
	typename MinSteinerTreeZelikovsky<T>::WinCalculation wc; // absolute, relative
	typename MinSteinerTreeZelikovsky<T>::TripleGeneration tg; // exhaustive, voronoi, none
	typename MinSteinerTreeZelikovsky<T>::TripleReducing tr; // on, off
	typename MinSteinerTreeZelikovsky<T>::SaveCalculation sc; // staticTree, staticLCATree, dynamicLCATree, hybrid
	typename MinSteinerTreeZelikovsky<T>::Pass pass; // one, multi

	desc = "Zelikovsky(";

	switch ((id % (2*4*2*3*2)) / (2*4*2*3)) {
	case 0:
		wc = MinSteinerTreeZelikovsky<T>::absolute;
		desc += "abs,";
		break;
	case 1:
		wc = MinSteinerTreeZelikovsky<T>::relative;
		desc += "rel,";
		break;
	}

	switch ((id % (2*4*2*3)) / (2*4*2)) {
	case 0:
		tg = MinSteinerTreeZelikovsky<T>::exhaustive;
		desc += "exh,";
		break;
	case 1:
		tg = MinSteinerTreeZelikovsky<T>::voronoi;
		desc += "vor,";
		break;
	case 2:
		tg = MinSteinerTreeZelikovsky<T>::none;
		desc += "non,";
		break;
	}

	switch ((id % (2*4*2)) / (2*4)) {
	case 0:
		tr = MinSteinerTreeZelikovsky<T>::on;
		desc += "I,";
		break;
	case 1:
		tr = MinSteinerTreeZelikovsky<T>::off;
		desc += "O,";
		break;
	}

	switch ((id % (2*4)) / 2) {
	case 0:
		sc = MinSteinerTreeZelikovsky<T>::staticTree;
		desc += "stT,";
		break;
	case 1:
		sc = MinSteinerTreeZelikovsky<T>::staticLCATree;
		desc += "sLT,";
		break;
	case 2:
		sc = MinSteinerTreeZelikovsky<T>::dynamicLCATree;
		desc += "dLT,";
		break;
	case 3:
		sc = MinSteinerTreeZelikovsky<T>::hybrid;
		desc += "hyb,";
		break;
	}

	switch (id % 2) {
	case 0:
		pass = MinSteinerTreeZelikovsky<T>::one;
		desc += "one)";
		break;
	case 1:
		pass = MinSteinerTreeZelikovsky<T>::multi;
		desc += "mul)";
		break;
	}

	alg = new MinSteinerTreeZelikovsky<T>(wc, tg, tr, sc, pass);
}

template <typename T>
static bool _regSteinerTree(EdgeWeightedGraph<T> &wG)
{
	List<node> terminals;
	NodeArray<bool> isTerminal(wG);

	if (!GraphIO::readSTP(wG, terminals, isTerminal, "test.stp")) {
		cerr << "Could not read `test.stp'." << endl;
		return false;
	}

	int j, i = 0, state = 1;
	while (1) {
		MinSteinerTreeModule<T> *alg;
		string desc;
		switch (i++) {
		case 0:
			alg = new MinSteinerTreeKou<T>();
			desc = "Kou";
			break;
		case 1:
			alg = new MinSteinerTreeMehlhorn<T>();
			desc = "Mehlhorn";
			break;
		case 2:
			alg = new MinSteinerTreeTakahashi<T>();
			desc = "Takahashi";
			break;
		case 3:
			alg = new MinSteinerTreeRZLoss<T>();
			desc = "RZLoss(def)";
			j = i;
			break;
		default:
			switch (state) {
			case 1:
				alg = new MinSteinerTreeRZLoss<T>(i - j + 2),
				desc = "RZLoss(k=";
				desc += char(i - j + 0x32);
				desc += ")";
				if (i - j + 2 == terminals.size()) {
					++state;
					j = i;
				}
				break;
			case 2:
				zelikovskyVariant(alg, desc, i-j-1);
				if (i - j == 2*4*2*3*2) {
					++state;
					j = i;
				}
				break;
			default:
				return true;
			}
		}
		T obj = computeMinStT(alg, desc, wG, terminals, isTerminal);
		delete alg;

		// 18 is the optimal value of the test instance
		if (obj < 18 || obj > 18*2) {
			return false;
		}
	}
	return true;
}

bool regSteinerTree()
{
	cout << "-> double:\n";
	EdgeWeightedGraph<double> Gd;
	if (!_regSteinerTree(Gd)) {
		return false;
	}

	cout << "-> int:\n";
	EdgeWeightedGraph<int> Gi;
	return _regSteinerTree(Gi);
}
