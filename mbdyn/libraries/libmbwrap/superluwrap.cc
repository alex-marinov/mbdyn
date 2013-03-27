/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/*****************************************************************************
 *                                                                           *
 *                          SuperLU C++ WRAPPER                              *
 *                                                                           *
 *****************************************************************************/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */


#ifdef USE_SUPERLU
#ifndef USE_SUPERLU_MT /* SUPERLU and SUPERLU_MT are incompatible */

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

#include "spmh.h"
#include "spmapmh.h"
#include "dirccmh.h"
#include "ccmh.h"


#include "superluwrap.h"

extern "C" {
#include <dsp_defs.h>
#include <util.h>
}

struct SuperLUSolverData {
	SuperMatrix		A,
				AC,
				L,
				U,
				B;
	
	SuperLUStat_t		Gstat;
	std::vector<int>	perm_c, /* column permutation vector */
				perm_r, /* row permutations from partial pivoting */
				etree ; /* elimination tree of A'*A */
				
	superlu_options_t	options;
};


/* SuperLUSolver - begin */

/* Costruttore: si limita ad allocare la memoria */
SuperLUSolver::SuperLUSolver(integer iMatOrd,
		const doublereal &dPivot,
		unsigned ptype)
: Aip(0),
App(0),
Axp(0),
iN(iMatOrd),
iNonZeroes(0),
dPivotFactor(dPivot),
permutation(ptype),
bFirstSol(true),
bRegenerateMatrix(true)
{
	ASSERT(iN > 0);
	ASSERTMSGBREAK(permutation & (SUPERLU_COLAMD|SUPERLU_MMDATA),
		"SuperLU scalar solver unknown permutation strategy");

	SAFENEW(sld, SuperLUSolverData);
	
	/*
	 * This means it's the first run
	 * FIXME: create a dependence on the library's internals
	 */
	sld->A.Store = NULL;

}

/* Distruttore */
SuperLUSolver::~SuperLUSolver(void)
{
}


#ifdef DEBUG
void 
SuperLUSolver::IsValid(void) const
{
	ASSERT(Aip != NULL);
	ASSERT(App != NULL);
	ASSERT(Axp != NULL);
	ASSERT(iN > 0); 
}
#endif /* DEBUG */

/* Fattorizza la matrice */
void
SuperLUSolver::Factor(void)
{
#ifdef DEBUG 
	IsValid();
#endif /* DEBUG */

	ASSERT(iNonZeroes > 0);

	doublereal	drop_tol = 0.0;
	void		*work = NULL;
	int		info = 0, lwork = 0;
	int		panel_size = sp_ienv(1),
			relax = sp_ienv(2);

	if (bRegenerateMatrix) {
		/* NOTE: we could use sld->A.Store == NULL */
		if (bFirstSol) {
			set_default_options(&sld->options);
			sld->options.DiagPivotThresh = dPivotFactor;
			//colperm_t permc_spec = MMD_AT_PLUS_A;
			switch (permutation) {
			case SUPERLU_MMDATA:
				sld->options.ColPerm = MMD_ATA;
				break;
			case SUPERLU_COLAMD: 
			default:
				sld->options.ColPerm = COLAMD;
				break;
			}
			sld->options.Fact = DOFACT;
			sld->options.PrintStat = NO;
			ASSERT(Astore == NULL);

			/* ---------------------------------------------------
			 * Allocate storage and initialize statistics variables. 
			 * ---------------------------------------------------*/
			/* Set up the dense matrix data structure for B. */
			dCreate_Dense_Matrix(&sld->B, iN, 1,
					LinearSolver::pdRhs,
					iN, SLU_DN, SLU_D, SLU_GE);

			/* Set up the sparse matrix data structure for A. */
			dCreate_CompCol_Matrix(&sld->A, iN, iN, iNonZeroes,
					Axp, Aip, App, SLU_NC, SLU_D, SLU_GE);

			StatInit(&sld->Gstat);
			
			sld->perm_c.resize(iN);
			sld->perm_r.resize(iN);
			sld->etree.resize(iN);
			
			bFirstSol = false;	/* never change this again */

		} else {
			sld->options.Fact = DOFACT;
			NCformat *Astore = (NCformat *) sld->A.Store;

			ASSERT(Astore);

			Astore->nnz = iNonZeroes;
			Astore->nzval = Axp;
			Astore->rowind = Aip;
			Astore->colptr = App;
			
			Destroy_CompCol_Permuted(&sld->AC);
			Destroy_SuperNode_Matrix(&sld->L);
			Destroy_CompCol_Matrix(&sld->U);
		}

		int	*pc = &(sld->perm_c[0]);
		get_perm_c(sld->options.ColPerm, &sld->A, pc);

		bRegenerateMatrix = false;
	} else {
		sld->options.Fact = SamePattern;
		Destroy_CompCol_Permuted(&sld->AC);
		Destroy_SuperNode_Matrix(&sld->L);
		Destroy_CompCol_Matrix(&sld->U);
	}


	int	*pr = &(sld->perm_r[0]),
		*pc = &(sld->perm_c[0]),
		*et = &(sld->etree[0]);


	sp_preorder(&sld->options, &sld->A, pc, et, &sld->AC);
	dgstrf(&sld->options, &sld->AC, drop_tol, relax, panel_size, et, work, lwork, pc, pr, 
		&sld->L, &sld->U, &sld->Gstat, &info);

}

/* Risolve */
void
SuperLUSolver::Solve(void) const
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
	
	if (bHasBeenReset) {
      		((SuperLUSolver *)this)->Factor();
      		bHasBeenReset = false;
	}

	/* ------------------------------------------------------------
	 * Solve the system A*X=B, overwriting B with X.
	 * ------------------------------------------------------------*/
	trans_t		trans = NOTRANS;
	int		info = 0;

	int		*pr = &(sld->perm_r[0]),
			*pc = &(sld->perm_c[0]);

	dgstrs(trans, &sld->L, &sld->U, pc, pr,
			&sld->B, &sld->Gstat, &info);
}

/* Index Form */
void
SuperLUSolver::MakeCompactForm(SparseMatrixHandler& mh,
		std::vector<doublereal>& Ax,
		std::vector<integer>& Ar, std::vector<integer>& Ac,
		std::vector<integer>& Ap) const
{
	/* no need to rebuild matrix */
	if (!bHasBeenReset) {
		return;
	}
	
	iNonZeroes = mh.MakeCompressedColumnForm(Ax, Ar, Ap, 0);
	ASSERT(iNonZeroes > 0);

	Axp = &Ax[0];
	Aip = &Ar[0];
	App = &Ap[0];

	/* rebuild matrix ... (CC is broken) */
	bRegenerateMatrix = true;

#if 0
	Destroy_CompCol_Matrix(&sld->A);
#endif
}

/* SuperLUSolver - end */


/* SuperLUSparseSolutionManager - begin: code */

/* Costruttore */
SuperLUSparseSolutionManager::SuperLUSparseSolutionManager(integer iSize,
		const doublereal& dPivotFactor, unsigned ptype)
: iMatSize(iSize), 
Ap(iSize + 1),
xb(iSize),
MH(iSize),
VH(iSize, &xb[0])
{
   	ASSERT(iSize > 0);
   	ASSERT((dPivotFactor >= 0.0) && (dPivotFactor <= 1.0));

   	SAFENEWWITHCONSTRUCTOR(SolutionManager::pLS, 
			       SuperLUSolver,
			       SuperLUSolver(iMatSize, dPivotFactor, ptype));
   
	pLS->pdSetResVec(&(xb[0]));
	pLS->pdSetSolVec(&(xb[0]));
	pLS->SetSolutionManager(this);

#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
}


/* Distruttore; verificare la distruzione degli oggetti piu' complessi */
SuperLUSparseSolutionManager::~SuperLUSparseSolutionManager(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */
   
   	/* Dealloca roba, tra cui i thread */
}

#ifdef DEBUG
/* Test di validita' del manager */
void 
SuperLUSparseSolutionManager::IsValid(void) const
{   
   	ASSERT(iMatSize > 0);
   
#ifdef DEBUG_MEMMANAGER
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(VH.pdGetVec()));
   	ASSERT(defaultMemoryManager.fIsPointerToBlock(pLS));
#endif /* DEBUG_MEMMANAGER */
   
   	ASSERT((VH.IsValid(), 1));
   	ASSERT((pLS->IsValid(), 1));
}
#endif /* DEBUG */

/* Inizializza il gestore delle matrici */
void
SuperLUSparseSolutionManager::MatrReset(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	pLS->Reset();
}

void
SuperLUSparseSolutionManager::MakeCompressedColumnForm(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
	pLS->MakeCompactForm(MH, Ax, Ai, Adummy, Ap);
}

/* Risolve il problema */
void
SuperLUSparseSolutionManager::Solve(void)
{
#ifdef DEBUG
   	IsValid();
#endif /* DEBUG */

	/* FIXME: move this to the matrix handler! */
   	MakeCompressedColumnForm();

#if 0
	std::cerr << "### after MakeIndexForm:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < Ap.size(); i++) {
		std::cerr << i << " " << Ap[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	unsigned c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << Ai[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == Ap[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif

   	pLS->Solve();

#if 0
	std::cerr << "### after Solve:" << std::endl
		<< "{col Ap[col]}={" << std::endl;
	for (unsigned i = 0; i < Ap.size(); i++) {
		std::cerr << i << " " << Ap[i] << std::endl;
	}
	std::cerr << "}" << std::endl;
	
	std::cerr << "{idx Ai[idx] col Ax[idx]}={" << std::endl;
	c = 0;
	for (unsigned i = 0; i < Ax.size(); i++) {
		std::cerr << i << " " << Ai[i] << " " << c << " " << Ax[i] << std::endl;
		if (i == Ap[c]) {
			c++;
		}
	}
	std::cerr << "}" << std::endl;
#endif
}

/* SuperLUSparseSolutionManager - end */

/* SuperLUSparseCCSolutionManager - begin */

template <class CC>
SuperLUSparseCCSolutionManager<CC>::SuperLUSparseCCSolutionManager(integer Dim,
		const doublereal &dPivot, unsigned ptype)
: SuperLUSparseSolutionManager(Dim, dPivot, ptype),
CCReady(false),
Ac(0)
{
	NO_OP;
}

template <class CC>
SuperLUSparseCCSolutionManager<CC>::~SuperLUSparseCCSolutionManager(void) 
{
	if (Ac) {
		SAFEDELETE(Ac);
	}
}

template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MatrReset(void)
{
	pLS->Reset();
}

/* Risolve il sistema  Fattorizzazione + Backward Substitution */
template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MakeCompressedColumnForm(void)
{
	if (!CCReady) {
		pLS->MakeCompactForm(MH, Ax, Ai, Adummy, Ap);

		if (Ac == 0) {
			SAFENEWWITHCONSTRUCTOR(Ac, CC, CC(Ax, Ai, Ap));
		}

		CCReady = true;
	}
}

/* Inizializzatore "speciale" */
template <class CC>
void
SuperLUSparseCCSolutionManager<CC>::MatrInitialize()
{
	CCReady = false;

	MatrReset();
}
	
/* Rende disponibile l'handler per la matrice */
template <class CC>
MatrixHandler*
SuperLUSparseCCSolutionManager<CC>::pMatHdl(void) const
{
	if (!CCReady) {
		return &MH;
	}

	ASSERT(Ac != 0);
	return Ac;
}

template class SuperLUSparseCCSolutionManager<CColMatrixHandler<0> >;
template class SuperLUSparseCCSolutionManager<DirCColMatrixHandler<0> >;

/* SuperLUSparseCCSolutionManager - end */

#endif /* !USE_SUPERLU_MT */
#endif /* USE_SUPERLU */
