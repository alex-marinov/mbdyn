/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2002
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* December 2001 
 * Modified to add a Sparse matrix in row form and to implement methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 1996-2002
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_UMFPACK
#include <umfpackwrap.h>

UmfpackSparseLUSolutionManager::UmfpackSparseLUSolutionManager(integer Dim,
		integer dummy, doublereal dPivot)
: A(Dim),
xVH(0),
bVH(0), 
x(Dim),
b(Dim), 
Symbolic(0),
Numeric(0),
HasBeenReset(true)
{
	SAFENEWWITHCONSTRUCTOR(xVH, MyVectorHandler, 
			MyVectorHandler(Dim, &(x[0])));
	SAFENEWWITHCONSTRUCTOR(bVH, MyVectorHandler, 
			MyVectorHandler(Dim, &(b[0])));
	pdRhs = &(b[0]);
	pdSol = &(x[0]);
	umfpack_defaults(Control);

	if (dPivot != -1. && (dPivot >= 0. && dPivot <= 1.)) {
		/*
		 * 1.0: true partial pivoting
		 * 0.0: treated as 1.0
		 * 
		 * default: 0.1
		 */
		Control[UMFPACK_PIVOT_TOLERANCE] = dPivot;
	}
}

UmfpackSparseLUSolutionManager::~UmfpackSparseLUSolutionManager(void) 
{
	umfpack_free_symbolic(&Symbolic);
	ASSERT(Symbolic == 0);

	umfpack_free_numeric(&Numeric);
	ASSERT(Numeric == 0);
	
	SAFEDELETE(xVH);
	SAFEDELETE(bVH);
}

void
UmfpackSparseLUSolutionManager::MatrInit(const doublereal& d)
{
	A.Reset(d);
	if (Numeric) {
		umfpack_free_numeric(&Numeric);
		ASSERT(Numeric == 0);
	}
	HasBeenReset = true;
}

bool 
UmfpackSparseLUSolutionManager::PrepareSymbolic(void)
{
	const int* const Aip = &(Ai[0]);
	const int* const App = &(Ap[0]);
	int status;

	status = umfpack_symbolic(b.size(), App, Aip, 
			&Symbolic, Control, Info);
	if (status != UMFPACK_OK) {
		umfpack_report_info(Control, Info) ;
		umfpack_report_status(Control, status);
		std::cerr << "umfpack_symbolic failed" << std::endl;

		/* de-allocate memory */
		umfpack_free_symbolic(&Symbolic);
		ASSERT(Symbolic == 0);

		return false;
	}

	return true;
}

/* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
void
UmfpackSparseLUSolutionManager::Solve(const doublereal /* dCoef */)
{
	double t = umfpack_timer() ;
	if (HasBeenReset) {
		A.MakeCompressedColumnForm(Ax, Ai, Ap);
		const double* const Axp = &(Ax[0]);
		const int* const Aip = &(Ai[0]);
		const int* const App = &(Ap[0]);
		int status;
		
		if (Symbolic == 0 && !PrepareSymbolic()) {
			THROW(ErrGeneric());
		}
#ifdef UMFPACK_REPORT
		umfpack_report_symbolic ("Symbolic factorization of A",
				Symbolic, Control) ;
#endif /* UMFPACK_REPORT */

		umfpack_report_info(Control, Info);

#ifdef UMFPACK_REPORT
		double t1 = umfpack_timer() - t;
#endif /* UMFPACK_REPORT */

		status = umfpack_numeric(App, Aip, Axp, Symbolic, 
				&Numeric, Control, Info);
		if (status == UMFPACK_ERROR_different_pattern) {
			umfpack_free_symbolic(&Symbolic);
			if (!PrepareSymbolic()) {
				THROW(ErrGeneric());
			}
			status = umfpack_numeric(App, Aip, Axp, Symbolic, 
					&Numeric, Control, Info);
		}

		if (status != UMFPACK_OK) {
			umfpack_report_info(Control, Info);
			umfpack_report_status(Control, status);
			std::cerr << "umfpack_numeric failed" << std::endl;

			/* de-allocate memory */
			umfpack_free_symbolic(&Symbolic);
			umfpack_free_numeric(&Numeric);
			ASSERT(Numeric == 0);

			THROW(ErrGeneric());
		}
		
#ifdef UMFPACK_REPORT
		umfpack_report_numeric ("Numeric factorization of A",
				Numeric, Control);
#endif /* UMFPACK_REPORT */

		umfpack_report_info(Control, Info);

#ifdef UMFPACK_REPORT
		t1 = umfpack_timer() - t;
#endif /* UMFPACK_REPORT */

		HasBeenReset = false;
	}
	
	BackSub(t);
}

/* Bacward Substitution */
void
UmfpackSparseLUSolutionManager::BackSub(doublereal t_iniz)
{
	const double* const Axp = &(Ax[0]);
	const int* const Aip = &(Ai[0]);
	const int* const App = &(Ap[0]);
	int status;

	ASSERT(HasBeenReset == false);
	
#ifdef UMFPACK_REPORT
	double t = t_iniz;
#endif /* UMFPACK_REPORT */
	Control[UMFPACK_IRSTEP]= 0;
	status = umfpack_solve("Ax=b", App, Aip, Axp, pdSol, pdRhs, 
			Numeric, Control, Info);
	if (status != UMFPACK_OK) {
		umfpack_report_info(Control, Info) ;
		umfpack_report_status(Control, status) ;
		std::cerr << "umfpack_solve failed" << std::endl;
		
		/* de-allocate memory */
		umfpack_free_numeric(&Numeric);
		ASSERT(Numeric == 0);

		THROW(ErrGeneric());
	}
	
	umfpack_report_info(Control, Info);

#ifdef UMFPACK_REPORT
	double t1 = umfpack_timer() - t;
#endif /* UMFPACK_REPORT */
}

/* Rende disponibile l'handler per la matrice */
SpMapMatrixHandler*
UmfpackSparseLUSolutionManager::pMatHdl(void) const
{
	return &A;
}

/* Rende disponibile l'handler per il termine noto */
MyVectorHandler*
UmfpackSparseLUSolutionManager::pResHdl(void) const {
	return bVH;
}

/* Rende disponibile l'handler per la soluzione */
MyVectorHandler*
UmfpackSparseLUSolutionManager::pSolHdl(void) const {
	return xVH;
}

#endif /* USE_UMFPACK */

