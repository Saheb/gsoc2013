/*
 * minisat_example.cpp
 *
 *  Created on: 12.11.2012
 *      Author: Robert Zeranski
 *
 */

#include <stdio.h>
#include <ogdf/basic/Graph.h>
#include <ogdf/minisat/Minisat.h>
#include <iostream>

using namespace std;
using namespace ogdf;
using namespace MinisatInternal;
using namespace Minisat;

int main (int argc, char **argv)
{

	Formula F;
	Model model;
	clause cl = F.newClause();

	cl->addMultiple(4,-1,-2,-3,4);
	cl->writeToConsole();

	F.newVars(11);

	for (int i = 1; i < 10; i++) {
			clause c = F.newClause();
			if (i % 2 == 0) c->add(i);
				else c->add(-i);
			c->add(i+1);
			F.finalizeClause(c);
	}

	F.finalizeClause(cl);

	bool val = F.solve( model );

	cout << "#vars = " 			<< F.getVariableCount() << endl;
	cout << "#clauses = " 		<< F.getClauseCount() 	<< endl;
	cout << "F satisfiable = " 	<< val					<< endl;
	if (val) model.printModel();
	F.reset();
#if 0
	F.writeFormulaToDimacs("example.cnf");
#endif

	return 0;
}


