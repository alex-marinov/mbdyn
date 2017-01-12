/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * Local headers for lapack library
 *
 * http://netlib.org
 */

#ifndef AC_LAPACK_H
#define AC_LAPACK_H

#ifdef USE_LAPACK

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

#include "ac/f2c.h"

/*
 * Add declarations of lapack routines used by MBDyn
 */

/* Subroutine */ extern int
__FC_DECL__(dggev)(
	const char *JOBVL,
	const char *JOBVR,
	const integer *N,
	doublereal *A,
	const integer *LDA,
	doublereal *B,
	const integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	const integer *LDVL,
	doublereal *VR,
	const integer *LDVR,
	doublereal *WORK,
	const integer *LWORK,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dggevx)(
	const char *BALANC,
	const char *JOBVL,
	const char *JOBVR,
	const char *SENSE,
	const integer *N,
	doublereal *A,
	const integer *LDA,
	doublereal *B,
	const integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	const integer *LDVL,
	doublereal *VR,
	const integer *LDVR,
	integer *ILO,
	integer *IHI,
	doublereal *LSCALE,
	doublereal *RSCALE,
	doublereal *ABNRM,
	doublereal *BBNRM,
	doublereal *RCONDE,
	doublereal *RCONDV,
	doublereal *WORK,
	const integer *LWORK,
	integer *IWORK,
	logical *BWORK,
	integer *INFO);

/* NOTE: according to lapack's documentation, dgegv() is deprecated
 * in favour of dggev()... */
/* Subroutine */ extern int
__FC_DECL__(dgegv)(
	const char *JOBVL,
	const char *JOBVR,
	const integer *N,
	doublereal *A,
	const integer *LDA,
	doublereal *B,
	const integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	const integer *LDVL,
	doublereal *VR,
	const integer *LDVR,
	doublereal *WORK,
	const integer *LWORK,
	integer *INFO);

/* Subroutine */ int
__FC_DECL__(dgeequ)(
	const integer *M,
	const integer *N,
	const doublereal *A,
	const integer *LDA,
	doublereal *R,
	doublereal *C,
	doublereal *ROWCND,
	doublereal *COLCND,
	doublereal *AMAX,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dgelsd)(
	const integer* m,
	const integer* n,
	const integer* nrhs,
	doublereal* a,
	const integer* lda,
	doublereal* b,
	const integer* ldb,
	doublereal* s,
	const doublereal* rcond,
	integer* rank,
	doublereal* work,
	const integer* lwork,
	integer* iwork,
	integer* info);

/* Subroutine */ extern int
__FC_DECL__(dgetrf)(
	const integer *M,
	const integer *N,
	doublereal *A,
	const integer *LDA,
	integer *IPIV,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dgetrs)(
	const char *MODE,
	const integer *N,
	const integer *NRHS,
	const doublereal *A,
	const integer *LDA,
	const integer *IPIV,
	doublereal *B,
	const integer *LDB,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dgecon)(
	const char *NORM,
	const integer *N,
	doublereal *A,
	const integer *LDA,
	const doublereal *ANORM,
	doublereal *RCOND,
	doublereal *WORK,
	integer *IWORK,
	integer *INFO);


/* Function */ extern doublereal
__FC_DECL__(dlamch)(
	const char* CMACH);


#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* USE_LAPACK */

#endif /* AC_LAPACK_H */

