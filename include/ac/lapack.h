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

/*
NAME
       DGEMM - one of the matrix-matrix operations   C := alpha*op( A )*op( B ) + beta*C,

SYNOPSIS
       SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)

           DOUBLE                                                         PRECISION ALPHA,BETA

           INTEGER                                                        K,LDA,LDB,LDC,M,N

           CHARACTER                                                      TRANSA,TRANSB

           DOUBLE                                                         PRECISION A(LDA,*),B(LDB,*),C(LDC,*)

PURPOSE
       DGEMM  performs one of the matrix-matrix operations

       where  op( X ) is one of

          op( X ) = X   or   op( X ) = X',

       alpha  and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix,  op( B )  a  k by
       n matrix and  C an m by n matrix.

ARGUMENTS
       TRANSA - CHARACTER*1.  On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplica-
       tion as follows:

       TRANSA = 'N' or 'n',  op( A ) = A.

       TRANSA = 'T' or 't',  op( A ) = A'.

       TRANSA = 'C' or 'c',  op( A ) = A'.

       Unchanged on exit.

       TRANSB - CHARACTER*1.  On entry, TRANSB specifies the form of op( B ) to be used in the matrix multiplica-
       tion as follows:

       TRANSB = 'N' or 'n',  op( B ) = B.

       TRANSB = 'T' or 't',  op( B ) = B'.

       TRANSB = 'C' or 'c',  op( B ) = B'.

       Unchanged on exit.

       M      - INTEGER.
              On entry,  M  specifies  the number  of rows  of the  matrix op( A )  and of  the   matrix   C.   M
              must  be at least  zero.  Unchanged on exit.

       N      - INTEGER.
              On  entry,   N  specifies the number  of columns of the matrix op( B ) and the number of columns of
              the matrix C. N must be at least zero.  Unchanged on exit.

       K      - INTEGER.
              On entry,  K  specifies  the number of columns of the matrix op( A ) and the number of rows of  the
              matrix op( B ). K must be at least  zero.  Unchanged on exit.

       ALPHA  - DOUBLE PRECISION.
              On entry, ALPHA specifies the scalar alpha.  Unchanged on exit.

       A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
              k   when  TRANSA = 'N' or 'n',  and is  m  otherwise.  Before entry with  TRANSA = 'N' or 'n',  the
              leading  m by k part of the array  A  must contain the matrix  A,  otherwise the leading   k  by  m
              part of the array  A  must contain  the matrix A.  Unchanged on exit.

       LDA    - INTEGER.
              On  entry,  LDA  specifies  the first dimension of A as declared in the calling (sub) program. When
              TRANSA = 'N' or 'n' then LDA must be at least  max( 1, m ), otherwise  LDA must be at  least   max(
              1, k ).  Unchanged on exit.

       B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
              n   when  TRANSB = 'N' or 'n',  and is  k  otherwise.  Before entry with  TRANSB = 'N' or 'n',  the
              leading  k by n part of the array  B  must contain the matrix  B,  otherwise the leading   n  by  k
              part of the array  B  must contain  the matrix B.  Unchanged on exit.

       LDB    - INTEGER.
              On  entry,  LDB  specifies  the first dimension of B as declared in the calling (sub) program. When
              TRANSB = 'N' or 'n' then LDB must be at least  max( 1, k ), otherwise  LDB must be at  least   max(
              1, n ).  Unchanged on exit.

       BETA   - DOUBLE PRECISION.
              On  entry,   BETA   specifies the scalar  beta.  When  BETA  is supplied as zero then C need not be
              set on input.  Unchanged on exit.

       C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
              Before entry, the leading  m by n  part of the array  C must contain the matrix   C,   except  when
              beta   is zero, in which case C need not be set on entry.  On exit, the array  C  is overwritten by
              the  m by n  matrix ( alpha*op( A )*op( B ) + beta*C ).

       LDC    - INTEGER.
              On entry, LDC specifies the first dimension of C as declared  in   the   calling   (sub)   program.
              LDC  must  be  at  least max( 1, m ).  Unchanged on exit.

              Level 3 Blas routine.

              --  Written  on 8-February-1989.  Jack Dongarra, Argonne National Laboratory.  Iain Duff, AERE Har-
              well.  Jeremy Du Croz, Numerical Algorithms Group Ltd.  Sven Hammarling, Numerical Algorithms Group
              Ltd.
 */

/* Subroutine */ extern int
__FC_DECL__(dgemm)(const char *const TRANSA,
                   const char *const TRANSB,
                   const integer *const M,
                   const integer *const N,
                   const integer *const K,
                   const doublereal *const ALPHA,
                   const doublereal *const A,
                   const integer *const LDA,
                   const doublereal *const B,
                   const integer *const LDB,
                   const doublereal *const BETA,
                   doublereal *const C,
                   const integer *const LDC);
        
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

/* Subroutine */ extern int
__FC_DECL__(dtrtrs)(const char* UPLO,
                    const char* TRANS,
                    const char* DIAG,
                    const integer* N,
                    const integer* NRHS,
                    const doublereal* A,
                    const integer* LDA,
                    doublereal* B,
                    const integer* LDB,
                    integer* INFO);

#if defined(__cplusplus)
}
#endif /* __cplusplus */

#endif /* USE_LAPACK */

#endif /* AC_LAPACK_H */

