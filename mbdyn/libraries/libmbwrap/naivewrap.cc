/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * The Naive Solver is copyright (C) 2004 by
 * Paolo Mantegazza <mantegazza@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

#include <set>

#include "spmh.h"
#include "spmapmh.h"
#include "naivewrap.h"
#include "mthrdslv.h"
#include "dgeequ.h"

/* NaiveSolver - begin */
NaiveSolver::NaiveSolver(const integer &size, const doublereal& dMP,
		NaiveMatrixHandler *const a)
: LinearSolver(0),
iSize(size),
dMinPiv(dMP < 0 ? 0 : dMP),
piv(size),
A(a)
{
	NO_OP;
}

NaiveSolver::~NaiveSolver(void)
{
	NO_OP;
}

void
NaiveSolver::SetMat(NaiveMatrixHandler *const a)
{
	A = a;
}

void
NaiveSolver::Reset(void)
{
	bHasBeenReset = true;
}

void
NaiveSolver::Solve(void) const
{
	if (bHasBeenReset) {
      		const_cast<NaiveSolver *>(this)->Factor();
      		bHasBeenReset = false;
	}

	integer rc = naivslv(A->ppdRows, iSize, A->piNzc, A->ppiCols,
			LinearSolver::pdRhs, LinearSolver::pdSol, &piv[0]);
	integer err = (rc & NAIVE_MASK);
	if (err) {
		switch (err) {
		case NAIVE_ERANGE:
			silent_cerr("NaiveSolver: ERANGE"
				<< std::endl);
			break;

		default:
			silent_cerr("NaiveSolver: (" << rc << ")"
				<< std::endl);
			break;
		}

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
NaiveSolver::Factor(void)
throw(LinearSolver::ErrFactor)
{
	integer rc = naivfct(A->ppdRows, iSize,
			A->piNzr, A->ppiRows, 
			A->piNzc, A->ppiCols,
			A->ppnonzero, 
			&piv[0], dMinPiv);

	integer err = (rc & NAIVE_MASK);
	if (err) {
		integer idx = (rc & NAIVE_MAX);
		switch (err) {
		case NAIVE_ENULCOL:
			silent_cerr("NaiveSolver: ENULCOL(" << idx << ")"
				<< std::endl);
			throw LinearSolver::ErrNullColumn(idx, MBDYN_EXCEPT_ARGS);
	
		case NAIVE_ENOPIV:
			silent_cerr("NaiveSolver: ENOPIV(" << idx << ")"
				<< std::endl);
			throw LinearSolver::ErrNoPivot(idx + 1, MBDYN_EXCEPT_ARGS);
	
		case NAIVE_ERANGE:
			silent_cerr("NaiveSolver: ERANGE"
				<< std::endl);
			break;

		default:
			silent_cerr("NaiveSolver: (" << rc << ")"
				<< std::endl);
			break;
		}

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* NaiveSolver - end */

/* NaiveSparseSolutionManager - begin */

NaiveSparseSolutionManager::NaiveSparseSolutionManager(const integer Dim,
	const doublereal dMP, SolutionManager::ScaleWhen ms)
: A(0),
VH(Dim),
ms(ms)
{
	SAFENEWWITHCONSTRUCTOR(A, NaiveMatrixHandler, NaiveMatrixHandler(Dim));
	SAFENEWWITHCONSTRUCTOR(pLS, NaiveSolver, NaiveSolver(Dim, dMP, A));

	pLS->pdSetResVec(VH.pdGetVec());
	pLS->pdSetSolVec(VH.pdGetVec());

	pLS->SetSolutionManager(this);
}

NaiveSparseSolutionManager::~NaiveSparseSolutionManager(void) 
{
	if (A != 0) {
		SAFEDELETE(A);
		A = 0;
	}
}

void
NaiveSparseSolutionManager::MatrReset()
{
	pLS->Reset();
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
void
NaiveSparseSolutionManager::Solve(void)
{
	if (ms != SolutionManager::NEVER) {
		if (pLS->bReset()) {
			if (msr.empty() || ms == SolutionManager::ALWAYS) {
				// (re)compute
				doublereal rowcnd = -1., colcnd = -1., amax = -1.;
				dgeequ<NaiveMatrixHandler, NaiveMatrixHandler::const_iterator>(*A, msr, msc, rowcnd, colcnd, amax);
			}
			// in any case scale matrix and right-hand-side
			dgeequ_scale<NaiveMatrixHandler, NaiveMatrixHandler::const_iterator>(*A, msr, msc);
		}
		dgeequ_scale(VH, &msr[0]);
	}

	pLS->Solve();

	if (ms != SolutionManager::NEVER) {
		// scale solution
		dgeequ_scale(VH, &msc[0]);
	}
}

/* Rende disponibile l'handler per la matrice */
MatrixHandler*
NaiveSparseSolutionManager::pMatHdl(void) const
{
	return A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
NaiveSparseSolutionManager::pResHdl(void) const
{
	return &VH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
NaiveSparseSolutionManager::pSolHdl(void) const
{
	return &VH;
}

/* NaiveSparseSolutionManager - end */

/* NaivePermSparseSolutionManager - begin */


template<class T>
NaiveSparsePermSolutionManager<T>::NaiveSparsePermSolutionManager(
	const integer Dim, 
	const doublereal dMP,
	SolutionManager::ScaleWhen ms)
: NaiveSparseSolutionManager(Dim, dMP, ms),
dMinPiv(dMP < 0 ? 0 : dMP),
TmpH(Dim),
ePermState(PERM_NO)
{
	perm.resize(Dim, 0);
	invperm.resize(Dim, 0);

	// replace matrix handler
	SAFEDELETE(A);
	A = 0;
	SAFENEWWITHCONSTRUCTOR(A, NaivePermMatrixHandler,
			NaivePermMatrixHandler(Dim, perm, invperm));

	dynamic_cast<NaiveSolver *>(pLS)->SetMat(A);

	MatrInitialize();
}

template<class T>
NaiveSparsePermSolutionManager<T>::~NaiveSparsePermSolutionManager(void) 
{
	NO_OP;
}

template<class T>
void
NaiveSparsePermSolutionManager<T>::MatrReset(void)
{
	if (ePermState == PERM_INTERMEDIATE) {
		ePermState = PERM_READY;
	}

	NaiveSparseSolutionManager::MatrReset();
}

template<class T>
void
NaiveSparsePermSolutionManager<T>::BackPerm(void)
{
	/* NOTE: use whatever is stored in pLS - someone could
	 * trick us into using its memory */
	doublereal *pd = pLS->pdGetResVec();

	ASSERT(pd != TmpH.pdGetVec());
	
	for (integer i = 0; i < A->iGetNumCols(); i++) {
		pd[invperm[i]] = TmpH(i + 1);
	}
}


/* Risolve il sistema: Fattorizzazione + Backward Substitution */
template<class T>
void
NaiveSparsePermSolutionManager<T>::Solve(void)
{
	doublereal *pd = 0;

	if (ms != SolutionManager::NEVER) {
		if (pLS->bReset()) {
			if (msr.empty() || ms == SolutionManager::ALWAYS) {
				// (re)compute
				doublereal rowcnd, colcnd, amax;
				dgeequ<NaivePermMatrixHandler, NaivePermMatrixHandler::const_iterator>(dynamic_cast<NaivePermMatrixHandler&>(*A), msr, msc, rowcnd, colcnd, amax);
			}
			// in any case scale matrix and right-hand-side
			dgeequ_scale<NaivePermMatrixHandler, NaivePermMatrixHandler::const_iterator>(dynamic_cast<NaivePermMatrixHandler&>(*A), msr, msc);
		}
		dgeequ_scale(VH, &msr[0]);
	}

	if (ePermState == PERM_NO) {
		ComputePermutation();

	} else if (ePermState == PERM_READY) {
		/* We need to use local storage to allow BackPerm();
		 * save and restore original pointer */
		pd = pLS->pdSetSolVec(TmpH.pdGetVec());
	}

	pLS->Solve();

	if (ePermState == PERM_READY) {
		BackPerm();

		ASSERT(pd != 0);
		pLS->pdSetSolVec(pd);
	}

	if (ms != SolutionManager::NEVER) {
		// scale solution
		dgeequ_scale(VH, &msc[0]);
	}
}

/* Inizializzatore "speciale" */
template<class T>
void
NaiveSparsePermSolutionManager<T>::MatrInitialize()
{
	ePermState = PERM_NO;
	for (integer i = 0; i < A->iGetNumRows(); i++) {
		perm[i] = i;
		invperm[i] = i;
	}

	MatrReset();
}

//explicit specializations

extern "C" {
#include "colamd.h"
}

template<>
void
NaiveSparsePermSolutionManager<Colamd_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	A->MakeCCStructure(Ai, invperm);
	doublereal knobs[COLAMD_KNOBS];
	integer stats[COLAMD_STATS];
	integer Alen = mbdyn_colamd_recommended(Ai.size(), A->iGetNumRows(),
			A->iGetNumCols());
	Ai.resize(Alen);
	mbdyn_colamd_set_defaults(knobs);
	if (!mbdyn_colamd(A->iGetNumRows(), A->iGetNumCols(), Alen,
		&Ai[0], &invperm[0], knobs, stats))
	{
		silent_cerr("colamd permutation failed" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	for (integer i = 0; i < A->iGetNumRows(); i++) {
		perm[invperm[i]] = i;
	}
	ePermState = PERM_INTERMEDIATE;
}


#ifdef USE_BOOST
/* NaivePermSparseSolutionManager - end */
#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/bandwidth.hpp"
#include <boost/graph/wavefront.hpp>

#ifdef HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP
#include "boost/graph/cuthill_mckee_ordering.hpp"
#endif /* HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_KING_ORDERING_HPP
#include "boost/graph/king_ordering.hpp"
#endif /* HAVE_BOOST_GRAPH_KING_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP
#include <boost/graph/sloan_ordering.hpp>
#endif /* HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP
#include "boost/graph/minimum_degree_ordering.hpp"
#endif /* HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP
template<>
void
NaiveSparsePermSolutionManager<rcmk_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	std::vector<integer> Ac;
	

	invperm.resize(A->iGetNumCols());
	A->MakeCCStructure(Ai, Ac);


/* boost */

	typedef boost::adjacency_list<
			boost::setS,
			boost::vecS,
			boost::undirectedS,
			boost::property<
				boost::vertex_color_t,
				boost::default_color_type,
				boost::property<
					boost::vertex_degree_t,
					integer
// 					,
// 					boost::property<
// 						boost::vertex_priority_t,
// 						double
// 					>
				>
			>
		> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::vertices_size_type size_type;


	Graph G(A->iGetNumRows());
	for (int col=0; col<A->iGetNumCols(); col++) {
		for (int i = Ac[col]; i < Ac[col + 1]; i++) {
			int row = Ai[i];
			if (row != col) {
				boost::add_edge(row, col, G);
				boost::add_edge(col, row, G);
			}
		}
	}

	
	std::vector<Vertex> inv_perm(num_vertices(G));
	boost::cuthill_mckee_ordering(G, inv_perm.rbegin());
#if 0
	boost::sloan_ordering(G, inv_perm.rbegin(),
		boost::get(boost::vertex_color, G),
		boost::make_degree_map(G),
		boost::get(boost::vertex_priority, G));
#endif

	for (integer i = 0; i < A->iGetNumRows(); i++) {
		invperm[i] = inv_perm[i];
		perm[invperm[i]] = i;
	}
	ePermState = PERM_INTERMEDIATE;
}
#endif /* HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP
template<>
void
NaiveSparsePermSolutionManager<sloan_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	std::vector<integer> Ac;
	

	invperm.resize(A->iGetNumCols());
	A->MakeCCStructure(Ai, Ac);


/* boost */

	typedef boost::adjacency_list<
			boost::setS,
			boost::vecS,
			boost::undirectedS,
			boost::property<
				boost::vertex_color_t,
				boost::default_color_type,
				boost::property<
					boost::vertex_degree_t,
					integer,
					boost::property<
						boost::vertex_priority_t,
						double
					>
				>
			>
		> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::vertices_size_type size_type;


	Graph G(A->iGetNumRows());
	for (int col=0; col<A->iGetNumCols(); col++) {
		for (int i = Ac[col]; i < Ac[col + 1]; i++) {
			int row = Ai[i];
			if (row != col) {
				boost::add_edge(row, col, G);
				boost::add_edge(col, row, G);
			}
		}
	}

	
	std::vector<Vertex> inv_perm(num_vertices(G));
	std::vector<Vertex> _perm(num_vertices(G));
	boost::sloan_ordering(G, _perm.begin(),
		boost::get(boost::vertex_color, G), 
		boost::make_degree_map(G),
		boost::get(boost::vertex_priority, G));

	for (integer i = 0; i < A->iGetNumRows(); i++) {
// 		invperm[i] = inv_perm[i];
// 		perm[invperm[i]] = i;
		perm[i] = _perm[i];
		invperm[perm[i]] = i;
	}
	ePermState = PERM_INTERMEDIATE;
}
#endif /* HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_KING_ORDERING_HPP
template<>
void
NaiveSparsePermSolutionManager<king_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	std::vector<integer> Ac;
	

	invperm.resize(A->iGetNumCols());
	A->MakeCCStructure(Ai, Ac);


/* boost */

	typedef boost::adjacency_list<
			boost::setS,
			boost::vecS,
			boost::undirectedS,
			boost::property<
				boost::vertex_color_t,
				boost::default_color_type,
				boost::property<
					boost::vertex_degree_t,
					integer
				>
			>
		> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::vertices_size_type size_type;


	Graph G(A->iGetNumRows());
	for (int col=0; col<A->iGetNumCols(); col++) {
		for (int i = Ac[col]; i < Ac[col + 1]; i++) {
			int row = Ai[i];
			if (row != col) {
				boost::add_edge(row, col, G);
				boost::add_edge(col, row, G);
			}
		}
	}

	std::vector<Vertex> inv_perm(num_vertices(G));
	boost::king_ordering(G, inv_perm.rbegin());

	for (integer i = 0; i < A->iGetNumRows(); i++) {
		invperm[i] = inv_perm[i];
		perm[invperm[i]] = i;
	}
	ePermState = PERM_INTERMEDIATE;
}
#endif /* HAVE_BOOST_GRAPH_KING_ORDERING_HPP */

#ifdef HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP
template<>
void
NaiveSparsePermSolutionManager<md_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	std::vector<integer> Ac;
	

	invperm.resize(A->iGetNumCols());
	A->MakeCCStructure(Ai, Ac);

	typedef boost::adjacency_list<
			boost::setS,
			boost::vecS,
			boost::directedS
		> Graph;
	typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
	typedef boost::graph_traits<Graph>::vertices_size_type size_type;


	Graph G(A->iGetNumRows());
	for (int col=0; col<A->iGetNumCols(); col++) {
		for (int i = Ac[col]; i < Ac[col + 1]; i++) {
			int row = Ai[i];
			if (row != col) {
				boost::add_edge(row, col, G);
				boost::add_edge(col, row, G);
			}
		}
	}

	boost::property_map<Graph, boost::vertex_index_t>::type 
		id = boost::get(boost::vertex_index, G);
	std::vector<Vertex> inv_perm(num_vertices(G));
	std::vector<Vertex> _perm(num_vertices(G));
	std::vector<integer> degree(A->iGetNumRows(), 0);
	std::vector<integer> supernode_sizes(A->iGetNumRows(), 1); // init has to be 1
	int delta = 0;

	boost::minimum_degree_ordering(
		G,
		make_iterator_property_map(&degree[0], id, degree[0]),
		&invperm[0],
		&perm[0],
		boost::make_iterator_property_map(
			&supernode_sizes[0],
			id,
			supernode_sizes[0]
		), 
		delta,
		id
	);

// 	for (integer i = 0; i < A->iGetNumRows(); i++) {
// 		invperm[i] = inv_perm[i];
// 		perm[invperm[i]] = i;
// 	}
	ePermState = PERM_INTERMEDIATE;
	


}
#endif /* HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP */
#endif /* USE_BOOST */

#ifdef USE_METIS
extern "C" {
#include "metiswrap.h"
}

template<>
void
NaiveSparsePermSolutionManager<metis_ordering>::ComputePermutation(void)
{
	std::vector<integer> Ai;
	std::vector<integer> Ac;
	

	invperm.resize(A->iGetNumCols());
	A->MakeCCStructure(Ai, Ac);


/* metis */

	typedef std::set<integer> row_cont_type;
	std::vector<row_cont_type> col_indices(A->iGetNumCols());
	
	integer NZ = 0;
	{
		integer ncols = A->iGetNumCols();
		for (integer col_id=0; col_id<ncols; col_id++) {
			for (integer i=Ac[col_id]; i<Ac[col_id+1]; i++) {
				int row_id = Ai[i];
				if (row_id != col_id) {
					//A
					row_cont_type& row1 = col_indices[col_id];
					if (row1.find(row_id) == row1.end()) {
						NZ++;
						row1.insert(row_id);
					}
					// + A^T
					row_cont_type& row2 = col_indices[row_id];
					if (row2.find(col_id) == row2.end()) {
						NZ++;
						row2.insert(col_id);
					}
				}
			}
		}
	}
	Ai.resize(NZ);
	{
		integer x_ptr = 0;
		row_cont_type::iterator ri;
		row_cont_type::const_iterator re;
		integer NCols = A->iGetNumCols();
		for (integer col = 0; col < NCols; col++) {
			Ac[col] = x_ptr;
			re = col_indices[col].end();
			for (ri = col_indices[col].begin(); ri != re; ri++) {
				Ai[x_ptr] = *ri;
				x_ptr++;
			}
		}
		Ac[NCols] = x_ptr;
	}
	int numflag = 0;
	int n = A->iGetNumCols();
	int options[8] = {1, 3, 1, 3, 0, 1, 200, 3};
//	METIS_EdgeND(&n, &(Ac[0]), &Ai[0], &numflag, options, &(perm[0]), &(invperm[0]));
 	METIS_NodeND(&n, &(Ac[0]), &Ai[0], &numflag, options, &perm[0], &invperm[0]);

	ePermState = PERM_INTERMEDIATE;

}

#endif //USE_METIS

// #ifdef HAVE_UMFPACK
// extern "C" {
// #include "amd.h"
// }
// template<>
// void
// NaiveSparsePermSolutionManager<amd_ordering>::ComputePermutation(void)
// {
// 	std::vector<integer> Ai;
// 	std::vector<integer> Ac;
// 	
// 
// 	invperm.resize(A->iGetNumCols());
// 	A->MakeCCStructure(Ai, Ac);
// 
// /* amd */
// 	double Control[AMD_CONTROL], Info[AMD_INFO];
// 	amd_defaults(Control);
// 	amd_order(A->iGetNumCols(), &(Ac[0]), &(Ai[0]), &(invperm[0]), Control, Info);
// 	for (integer i = 0; i < A->iGetNumRows(); i++) {
// 		perm[invperm[i]] = i;
// 	}
// 	
// 
// 
// }
// #endif //HAVE_UMFPACK


//explicit instantiations:
template class NaiveSparsePermSolutionManager<Colamd_ordering>;
#ifdef USE_BOOST
#ifdef HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP
template class NaiveSparsePermSolutionManager<rcmk_ordering>;
#endif /* HAVE_BOOST_GRAPH_CUTHILL_MCKEE_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_KING_ORDERING_HPP
template class NaiveSparsePermSolutionManager<king_ordering>;
#endif /* HAVE_BOOST_GRAPH_KING_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP
template class NaiveSparsePermSolutionManager<md_ordering>;
#endif /* HAVE_BOOST_GRAPH_MINIMUM_DEGREE_ORDERING_HPP */
#ifdef HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP
template class NaiveSparsePermSolutionManager<sloan_ordering>;
#endif /* HAVE_BOOST_GRAPH_SLOAN_ORDERING_HPP */
#endif /* USE_BOOST */
#ifdef USE_METIS
template class NaiveSparsePermSolutionManager<metis_ordering>;
#endif //USE_METIS
// #ifdef HAVE_UMFPACK
// template class NaiveSparsePermSolutionManager<amd_ordering>;
// #endif //HAVE_UMFPACK

