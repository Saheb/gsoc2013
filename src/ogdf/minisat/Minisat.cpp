// disable VC++ warnings when using fopen
#define _CRT_SECURE_NO_WARNINGS

#include <ogdf/minisat/Minisat.h>

using MinisatInternal::Var;
using MinisatInternal::Lit;
using MinisatInternal::lbool;
using MinisatInternal::vec;

namespace Minisat
{

void Clause::addMultiple ( int Amount, ... )
{
	va_list params;
	Var paramVar;
	va_start(params, Amount);
	for ( int i = 0; i < Amount; i++ )
	{
		paramVar = va_arg(params, Var);
		Lit l;
		if ( paramVar >= 0 )
		{
			l = MinisatInternal::mkLit(paramVar-1, true);
		}
		else
		{
			l = MinisatInternal::mkLit( (paramVar + 1) * (-1), false);
		}
		m_ps.push(l);

	}
	va_end(params);
}

Clause *Formula::newClause()
{
	m_Clauses.push_back( new Clause );
	return (m_Clauses.back());

}

void Formula::finalizeClause( const clause cl )
{
	for ( int i = 0; i < cl->m_ps.size(); i++ )
	{
		// if an vairable does not exist, it will be generated (and all between the gap)
		if ( ! ( MinisatInternal::var(cl->m_ps[i]) < Solver::nVars()) )
		{
			int max = Solver::nVars();
			for ( int j = 0; j < MinisatInternal::var(cl->m_ps[i]) + 1 - max; j++ )
				Solver::newVar();
		}
	}
	Solver::addClause( cl->m_ps );
}

bool Formula::finalizeNotExtensibleClause ( const clause cl )
{
	//proofs if variables from clause are valid (still known to formula)
	for ( int i = 0; i < cl->m_ps.size(); i++ )
	{
		if ( ! ( MinisatInternal::var(cl->m_ps[i]) < Solver::nVars()) )
		{
			m_messages << "Variable " << i << " is not present.";
			return false;
		}
	}
	Solver::addClause( cl->m_ps );
	return true;
}

Clause *Formula::getClause( const int pos )
{
	if ( pos < (int)m_Clauses.size() )
		return m_Clauses[pos];
	else
		return NULL;
}

bool Formula::solve( Model &ReturnModel )
{
	bool solv = Solver::solve();

	if ( solv )
		ReturnModel.setModel ( *this );

	return solv;
}

void Formula::removeClause( int i )
{
	MinisatInternal::CRef cr = Solver::clauses[i];
	Solver::removeClause( cr );
	int j, k;
	for ( j = k = 0; j < Solver::clauses.size(); j++ )
	{
		if ( !i == j )
			clauses[k++] = clauses[j] ;
	}
	clauses.shrink( j - k );
	delete &m_Clauses[i];
	m_Clauses.erase( m_Clauses.begin() + i );
}

void Formula::reset()
{
	free();
	Solver::assigns.clear();
	Solver::vardata.clear();
	Solver::activity.clear();
	Solver::seen.clear();
	Solver::polarity.clear();
	Solver::decision.clear();
	Solver::trail.clear();
	Solver::dec_vars = 0;
	Solver::model.clear();
}

void Formula::free()
{
	//int a = Solver::clauses.size();
	for ( int i = 0; i < Solver::clauses.size(); i++ )
	{
		Solver::removeClause(Solver::clauses[i]);
		delete m_Clauses[i];

	}
	Solver::clauses.shrink(Solver::clauses.size());
	m_Clauses.clear();
}

bool Formula::readDimacsFile( const char *file )
{
	std::ifstream in;
	std::string currentString;
	int currentInt;
	int clauseCount = 0;
	bool isCNF = false;

	Clause *clauses = NULL;
	in.open(file);

	for (;;) {
		if (in.eof())
			break;
		if (!isCNF) {
			in >> currentString;
			if (currentString == "p") {
				in >> currentString;
				if (currentString == "cnf") {
					isCNF = true;
				}
			}
		}
		if (isCNF) {
			in >> currentInt;
			newVars(currentInt);
			in >> clauseCount;
			clauses = new Clause[currentInt];
			break;
		}
	}

	if (isCNF) {
		int i = 0;
		while (!in.eof()) {
			in >> currentInt;
			if (currentInt == 0
			 && i < clauseCount) {
				clause c = newClause();
				clauses[i].m_ps.copyTo(c->m_ps);
				finalizeClause(c);
				i++;
			}
			else if (i < clauseCount)
				clauses[i].add(currentInt);
		}
	}

	in.close();
	if (clauses)
		delete[] clauses;
	return true;
}

void Formula::writeFormulaToDimacs( const char *filename)
{
	FILE* f = fopen(filename, "w");
	fprintf(f, "p cnf %d %d\n", getVariableCount(), getClauseCount());
	fprintf(f, "\n");
	for ( int i = 0; i < (int)m_Clauses.size(); i++ )
	{
		for (int j = 0; j < m_Clauses[i]->m_ps.size(); j++)
		{
			std::cout 	<< "Sign : "
						<< MinisatInternal::sign(m_Clauses[i]->m_ps[j])
						<< "Var : "
						<< MinisatInternal::var(m_Clauses[i]->m_ps[j]) + 1
						<< std::endl;
			fprintf(f,
					" %c%d ",
					Clause::convertLitSign(m_Clauses[i]->m_ps[j]),
					MinisatInternal::var(m_Clauses[i]->m_ps[j]) + 1 );
		}
		fprintf(f, "0\n");
	}
	fclose(f);
}

} // MinisatAdvanced
