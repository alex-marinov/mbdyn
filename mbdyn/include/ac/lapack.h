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
	char *JOBVL,
	char *JOBVR,
	integer *N,
	doublereal *A,
	integer *LDA,
	doublereal *B,
	integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	integer *LDVL,
	doublereal *VR,
	integer *LDVR,
	doublereal *WORK,
	integer *LWORK,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dggevx)(
	char *BALANC,
	char *JOBVL,
	char *JOBVR,
	char *SENSE,
	integer *N,
	doublereal *A,
	integer *LDA,
	doublereal *B,
	integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	integer *LDVL,
	doublereal *VR,
	integer *LDVR,
	integer *ILO,
	integer *IHI,
	doublereal *LSCALE,
	doublereal *RSCALE,
	doublereal *ABNRM,
	doublereal *BBNRM,
	doublereal *RCONDE,
	doublereal *RCONDV,
	doublereal *WORK,
	integer *LWORK,
	integer *IWORK,
	logical *BWORK,
	integer *INFO);

/* NOTE: according to lapack's documentation, dgegv() is deprecated
 * in favour of dggev()... */
/* Subroutine */ extern int
__FC_DECL__(dgegv)(
	char *JOBVL,
	char *JOBVR,
	integer *N,
	doublereal *A,
	integer *LDA,
	doublereal *B,
	integer *LDB,
	doublereal *ALPHAR,
	doublereal *ALPHAI,
	doublereal *BETA,
	doublereal *VL,
	integer *LDVL,
	doublereal *VR,
	integer *LDVR,
	doublereal *WORK,
	integer *LWORK,
	integer *INFO);

/* Subroutine */ int
__FC_DECL__(dgeequ)(
	integer *M,
	integer *N,
	doublereal *A,
	integer *LDA,
	doublereal *R,
	doublereal *C,
	doublereal *ROWCND,
	doublereal *COLCND,
	doublereal *AMAX,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dgelsd)(
	integer* m,
	integer* n,
	integer* nrhs,
	doublereal* a,
	integer* lda,
	doublereal* b,
	integer* ldb,
	doublereal* s,
	doublereal* rcond,
	integer* rank,
	doublereal* work,
	integer* lwork,
	integer* iwork,
	integer* info);

/* Subroutine */ extern int
__FC_DECL__(dgetrf)(
	integer *N,
	integer *N2,
	doublereal *A,
	integer *LDA,
	integer *IPIV,
	integer *INFO);

/* Subroutine */ extern int
__FC_DECL__(dgetrs)(
	char *MODE,
	integer *N,
	integer *NRHS,
	doublereal *A,
	integer *LDA,
	integer *IPIV,
	doublereal *B,
	integer *LDB,
	integer *INFO);

#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* USE_LAPACK */

#endif /* AC_LAPACK_H */

