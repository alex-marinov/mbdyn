/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
 * This file contains the subset of the lapack and blas routines
 * that are required to use the dgegv subroutine for the double precision
 * extraction of eigenvalues and eigenvectors of a generalized problem
 *
 *		A * X = lambda * B * X
 *
 * where A and X are real non-symmetric matrices which can be both singular.
 *
 * This package requires the linkning of the f2c libraries:
 * 
 *		cc -o myprog myobj.o libdgegv.c -lf2c -lm
 *
 * as well as the f2c.h header file.
 *
 * The above reported copyright notice refers to the whole package, 
 * while the intellectual property of the lapack and blas libraries
 * remains with their Authors, as indicated in each subroutine's
 * copyright notice, here reported for completeness:
 */

/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999
*/

/*                       
 * 
 * The original package can be found at Netlib:
 *
 *		http://www.netlib.org/lapack
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#if !defined(HAVE_DGEGV) && defined(__HACK_EIG__)

#include <stdlib.h>
#include <stdio.h>

#include <ac/f2c-int.h>
#include <ac/lapack.h>

#ifndef HAVE_POW_DI

static double pow_di(doublereal *d, integer *i)
{
	doublereal ret = *d;

	for ( ; --(*i) > 0;  ) {
		ret *= *d;
	}

	return ret;
}

#endif /* !HAVE_POW_DI */

/* Subroutine */ int dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK driver routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    This routine is deprecated and has been replaced by routine DGGEV.   

    DGEGV computes for a pair of n-by-n real nonsymmetric matrices A and   
    B, the generalized eigenvalues (alphar +/- alphai*i, beta), and   
    optionally, the left and/or right generalized eigenvectors (VL and   
    VR).   

    A generalized eigenvalue for a pair of matrices (A,B) is, roughly   
    speaking, a scalar w or a ratio  alpha/beta = w, such that  A - w*B   
    is singular.  It is usually represented as the pair (alpha,beta),   
    as there is a reasonable interpretation for beta=0, and even for   
    both being zero.  A good beginning reference is the book, "Matrix   
    Computations", by G. Golub & C. van Loan (Johns Hopkins U. Press)   

    A right generalized eigenvector corresponding to a generalized   
    eigenvalue  w  for a pair of matrices (A,B) is a vector  r  such   
    that  (A - w B) r = 0 .  A left generalized eigenvector is a vector   
    l such that l**H * (A - w B) = 0, where l**H is the   
    conjugate-transpose of l.   

    Note: this routine performs "full balancing" on A and B -- see   
    "Further Details", below.   

    Arguments   
    =========   

    JOBVL   (input) CHARACTER*1   
            = 'N':  do not compute the left generalized eigenvectors;   
            = 'V':  compute the left generalized eigenvectors.   

    JOBVR   (input) CHARACTER*1   
            = 'N':  do not compute the right generalized eigenvectors;   
            = 'V':  compute the right generalized eigenvectors.   

    N       (input) INTEGER   
            The order of the matrices A, B, VL, and VR.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the first of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of A on exit, see "Further   
            Details", below.)   

    LDA     (input) INTEGER   
            The leading dimension of A.  LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the second of the pair of matrices whose   
            generalized eigenvalues and (optionally) generalized   
            eigenvectors are to be computed.   
            On exit, the contents will have been destroyed.  (For a   
            description of the contents of B on exit, see "Further   
            Details", below.)   

    LDB     (input) INTEGER   
            The leading dimension of B.  LDB >= max(1,N).   

    ALPHAR  (output) DOUBLE PRECISION array, dimension (N)   
    ALPHAI  (output) DOUBLE PRECISION array, dimension (N)   
    BETA    (output) DOUBLE PRECISION array, dimension (N)   
            On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will   
            be the generalized eigenvalues.  If ALPHAI(j) is zero, then   
            the j-th eigenvalue is real; if positive, then the j-th and   
            (j+1)-st eigenvalues are a complex conjugate pair, with   
            ALPHAI(j+1) negative.   

            Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)   
            may easily over- or underflow, and BETA(j) may even be zero.   
            Thus, the user should avoid naively computing the ratio   
            alpha/beta.  However, ALPHAR and ALPHAI will be always less   
            than and usually comparable with norm(A) in magnitude, and   
            BETA always less than and usually comparable with norm(B).   

    VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)   
            If JOBVL = 'V', the left generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVL = 'N'.   

    LDVL    (input) INTEGER   
            The leading dimension of the matrix VL. LDVL >= 1, and   
            if JOBVL = 'V', LDVL >= N.   

    VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)   
            If JOBVR = 'V', the right generalized eigenvectors.  (See   
            "Purpose", above.)  Real eigenvectors take one column,   
            complex take two columns, the first for the real part and   
            the second for the imaginary part.  Complex eigenvectors   
            correspond to an eigenvalue with positive imaginary part.   
            Each eigenvector will be scaled so the largest component   
            will have abs(real part) + abs(imag. part) = 1, *except*   
            that for eigenvalues with alpha=beta=0, a zero vector will   
            be returned as the corresponding eigenvector.   
            Not referenced if JOBVR = 'N'.   

    LDVR    (input) INTEGER   
            The leading dimension of the matrix VR. LDVR >= 1, and   
            if JOBVR = 'V', LDVR >= N.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,8*N).   
            For good performance, LWORK must generally be larger.   
            To compute the optimal value of LWORK, call ILAENV to get   
            blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:   
            NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR;   
            The optimal LWORK is:   
                2*N + MAX( 6*N, N*(NB+1) ).   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1,...,N:   
                  The QZ iteration failed.  No eigenvectors have been   
                  calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)   
                  should be correct for j=INFO+1,...,N.   
            > N:  errors that usually indicate LAPACK problems:   
                  =N+1: error return from DGGBAL   
                  =N+2: error return from DGEQRF   
                  =N+3: error return from DORMQR   
                  =N+4: error return from DORGQR   
                  =N+5: error return from DGGHRD   
                  =N+6: error return from DHGEQZ (other than failed   
                                                  iteration)   
                  =N+7: error return from DTGEVC   
                  =N+8: error return from DGGBAK (computing VL)   
                  =N+9: error return from DGGBAK (computing VR)   
                  =N+10: error return from DLASCL (various calls)   

    Further Details   
    ===============   

    Balancing   
    ---------   

    This driver calls DGGBAL to both permute and scale rows and columns   
    of A and B.  The permutations PL and PR are chosen so that PL*A*PR   
    and PL*B*R will be upper triangular except for the diagonal blocks   
    A(i:j,i:j) and B(i:j,i:j), with i and j as close together as   
    possible.  The diagonal scaling matrices DL and DR are chosen so   
    that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to   
    one (except for the elements that start out zero.)   

    After the eigenvalues and eigenvectors of the balanced matrices   
    have been computed, DGGBAK transforms the eigenvectors back to what   
    they would have been (in perfect arithmetic) if they had not been   
    balanced.   

    Contents of A and B on Exit   
    -------- -- - --- - -- ----   

    If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or   
    both), then on exit the arrays A and B will contain the real Schur   
    form[*] of the "balanced" versions of A and B.  If no eigenvectors   
    are computed, then only the diagonal blocks will be correct.   

    [*] See DHGEQZ, DGEGS, or read the book "Matrix Computations",   
        by Golub & van Loan, pub. by Johns Hopkins U. Press.   

    =====================================================================   


       Decode the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static doublereal c_b27 = 1.;
    static doublereal c_b38 = 0.;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;
    /* Local variables */
    static doublereal absb, anrm, bnrm;
    static integer itau;
    static doublereal temp;
    static logical ilvl, ilvr;
    static integer lopt;
    static doublereal anrm1, anrm2, bnrm1, bnrm2, absai, scale, absar, sbeta;
    extern logical lsame_(char *, char *);
    static integer ileft, iinfo, icols, iwork, irows, jc;
    extern /* Subroutine */ int dggbak_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *);
    static integer nb;
    extern /* Subroutine */ int dggbal_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer in;
    extern doublereal dlamch_(char *), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *);
    static integer jr;
    static doublereal salfai;
    extern /* Subroutine */ int dgghrd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), dlascl_(char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *);
    static doublereal salfar;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal safmax;
    static char chtemp[1];
    static logical ldumma[1];
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), dtgevc_(char *, char *, 
	    logical *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *), 
	    xerbla_(char *, integer *);
    static integer ijobvl, iright;
    static logical ilimit;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ijobvr;
    extern /* Subroutine */ int dorgqr_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal onepls;
    static integer lwkmin, nb1, nb2, nb3;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static integer lwkopt;
    static logical lquery;
    static integer ihi, ilo;
    static doublereal eps;
    static logical ilv;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define vl_ref(a_1,a_2) vl[(a_2)*vl_dim1 + a_1]
#define vr_ref(a_1,a_2) vr[(a_2)*vr_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1 * 1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1 * 1;
    vr -= vr_offset;
    --work;

    /* Function Body */
    if (lsame_(jobvl, "N")) {
	ijobvl = 1;
	ilvl = FALSE_;
    } else if (lsame_(jobvl, "V")) {
	ijobvl = 2;
	ilvl = TRUE_;
    } else {
	ijobvl = -1;
	ilvl = FALSE_;
    }

    if (lsame_(jobvr, "N")) {
	ijobvr = 1;
	ilvr = FALSE_;
    } else if (lsame_(jobvr, "V")) {
	ijobvr = 2;
	ilvr = TRUE_;
    } else {
	ijobvr = -1;
	ilvr = FALSE_;
    }
    ilv = ilvl || ilvr;

/*     Test the input arguments   

   Computing MAX */
    i__1 = *n << 3;
    lwkmin = max(i__1,1);
    lwkopt = lwkmin;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    *info = 0;
    if (ijobvl <= 0) {
	*info = -1;
    } else if (ijobvr <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    } else if (*ldvl < 1 || ilvl && *ldvl < *n) {
	*info = -12;
    } else if (*ldvr < 1 || ilvr && *ldvr < *n) {
	*info = -14;
    } else if (*lwork < lwkmin && ! lquery) {
	*info = -16;
    }

    if (*info == 0) {
	nb1 = ilaenv_(&c__1, "DGEQRF", " ", n, n, &c_n1, &c_n1, (ftnlen)6, (
		ftnlen)1);
	nb2 = ilaenv_(&c__1, "DORMQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
	nb3 = ilaenv_(&c__1, "DORGQR", " ", n, n, n, &c_n1, (ftnlen)6, (
		ftnlen)1);
/* Computing MAX */
	i__1 = max(nb1,nb2);
	nb = max(i__1,nb3);
/* Computing MAX */
	i__1 = *n * 6, i__2 = *n * (nb + 1);
	lopt = (*n << 1) + max(i__1,i__2);
	work[1] = (doublereal) lopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEGV ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

/*     Get machine constants */

    eps = dlamch_("E") * dlamch_("B");
    safmin = dlamch_("S");
    safmin += safmin;
    safmax = 1. / safmin;
    onepls = eps * 4 + 1.;

/*     Scale A */

    anrm = dlange_("M", n, n, &a[a_offset], lda, &work[1]);
    anrm1 = anrm;
    anrm2 = 1.;
    if (anrm < 1.) {
	if (safmax * anrm < 1.) {
	    anrm1 = safmin;
	    anrm2 = safmax * anrm;
	}
    }

    if (anrm > 0.) {
	dlascl_("G", &c_n1, &c_n1, &anrm, &c_b27, n, n, &a[a_offset], lda, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Scale B */

    bnrm = dlange_("M", n, n, &b[b_offset], ldb, &work[1]);
    bnrm1 = bnrm;
    bnrm2 = 1.;
    if (bnrm < 1.) {
	if (safmax * bnrm < 1.) {
	    bnrm1 = safmin;
	    bnrm2 = safmax * bnrm;
	}
    }

    if (bnrm > 0.) {
	dlascl_("G", &c_n1, &c_n1, &bnrm, &c_b27, n, n, &b[b_offset], ldb, &
		iinfo);
	if (iinfo != 0) {
	    *info = *n + 10;
	    return 0;
	}
    }

/*     Permute the matrix to make it more nearly triangular   
       Workspace layout:  (8*N words -- "work" requires 6*N words)   
          left_permutation, right_permutation, work... */

    ileft = 1;
    iright = *n + 1;
    iwork = iright + *n;
    dggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &work[
	    ileft], &work[iright], &work[iwork], &iinfo);
    if (iinfo != 0) {
	*info = *n + 1;
	goto L120;
    }

/*     Reduce B to triangular form, and initialize VL and/or VR   
       Workspace layout:  ("work..." must have at least N words)   
          left_permutation, right_permutation, tau, work... */

    irows = ihi + 1 - ilo;
    if (ilv) {
	icols = *n + 1 - ilo;
    } else {
	icols = irows;
    }
    itau = iwork;
    iwork = itau + irows;
    i__1 = *lwork + 1 - iwork;
    dgeqrf_(&irows, &icols, &b_ref(ilo, ilo), ldb, &work[itau], &work[iwork], 
	    &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 2;
	goto L120;
    }

    i__1 = *lwork + 1 - iwork;
    dormqr_("L", "T", &irows, &icols, &irows, &b_ref(ilo, ilo), ldb, &work[
	    itau], &a_ref(ilo, ilo), lda, &work[iwork], &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	*info = *n + 3;
	goto L120;
    }

    if (ilvl) {
	dlaset_("Full", n, n, &c_b38, &c_b27, &vl[vl_offset], ldvl)
		;
	i__1 = irows - 1;
	i__2 = irows - 1;
	dlacpy_("L", &i__1, &i__2, &b_ref(ilo + 1, ilo), ldb, &vl_ref(ilo + 1,
		 ilo), ldvl);
	i__1 = *lwork + 1 - iwork;
	dorgqr_(&irows, &irows, &irows, &vl_ref(ilo, ilo), ldvl, &work[itau], 
		&work[iwork], &i__1, &iinfo);
	if (iinfo >= 0) {
/* Computing MAX */
	    i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
	    lwkopt = max(i__1,i__2);
	}
	if (iinfo != 0) {
	    *info = *n + 4;
	    goto L120;
	}
    }

    if (ilvr) {
	dlaset_("Full", n, n, &c_b38, &c_b27, &vr[vr_offset], ldvr)
		;
    }

/*     Reduce to generalized Hessenberg form */

    if (ilv) {

/*        Eigenvectors requested -- work on whole matrix. */

	dgghrd_(jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[b_offset], 
		ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr, &iinfo);
    } else {
	dgghrd_("N", "N", &irows, &c__1, &irows, &a_ref(ilo, ilo), lda, &
		b_ref(ilo, ilo), ldb, &vl[vl_offset], ldvl, &vr[vr_offset], 
		ldvr, &iinfo);
    }
    if (iinfo != 0) {
	*info = *n + 5;
	goto L120;
    }

/*     Perform QZ algorithm   
       Workspace layout:  ("work..." must have at least 1 word)   
          left_permutation, right_permutation, work... */

    iwork = itau;
    if (ilv) {
	*(unsigned char *)chtemp = 'S';
    } else {
	*(unsigned char *)chtemp = 'E';
    }
    i__1 = *lwork + 1 - iwork;
    dhgeqz_(chtemp, jobvl, jobvr, n, &ilo, &ihi, &a[a_offset], lda, &b[
	    b_offset], ldb, &alphar[1], &alphai[1], &beta[1], &vl[vl_offset], 
	    ldvl, &vr[vr_offset], ldvr, &work[iwork], &i__1, &iinfo);
    if (iinfo >= 0) {
/* Computing MAX */
	i__1 = lwkopt, i__2 = (integer) work[iwork] + iwork - 1;
	lwkopt = max(i__1,i__2);
    }
    if (iinfo != 0) {
	if (iinfo > 0 && iinfo <= *n) {
	    *info = iinfo;
	} else if (iinfo > *n && iinfo <= *n << 1) {
	    *info = iinfo - *n;
	} else {
	    *info = *n + 6;
	}
	goto L120;
    }

    if (ilv) {

/*        Compute Eigenvectors  (DTGEVC requires 6*N words of workspace) */

	if (ilvl) {
	    if (ilvr) {
		*(unsigned char *)chtemp = 'B';
	    } else {
		*(unsigned char *)chtemp = 'L';
	    }
	} else {
	    *(unsigned char *)chtemp = 'R';
	}

	dtgevc_(chtemp, "B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
		&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
		iwork], &iinfo);
	if (iinfo != 0) {
	    *info = *n + 7;
	    goto L120;
	}

/*        Undo balancing on VL and VR, rescale */

	if (ilvl) {
	    dggbak_("P", "L", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vl[vl_offset], ldvl, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 8;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		if (alphai[jc] < 0.) {
		    goto L50;
		}
		temp = 0.;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = vl_ref(jr, jc), abs(d__1))
				;
			temp = max(d__2,d__3);
/* L10: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = vl_ref(jr, jc), abs(d__1))
				 + (d__2 = vl_ref(jr, jc + 1), abs(d__2));
			temp = max(d__3,d__4);
/* L20: */
		    }
		}
		if (temp < safmin) {
		    goto L50;
		}
		temp = 1. / temp;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vl_ref(jr, jc) = vl_ref(jr, jc) * temp;
/* L30: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vl_ref(jr, jc) = vl_ref(jr, jc) * temp;
			vl_ref(jr, jc + 1) = vl_ref(jr, jc + 1) * temp;
/* L40: */
		    }
		}
L50:
		;
	    }
	}
	if (ilvr) {
	    dggbak_("P", "R", n, &ilo, &ihi, &work[ileft], &work[iright], n, &
		    vr[vr_offset], ldvr, &iinfo);
	    if (iinfo != 0) {
		*info = *n + 9;
		goto L120;
	    }
	    i__1 = *n;
	    for (jc = 1; jc <= i__1; ++jc) {
		if (alphai[jc] < 0.) {
		    goto L100;
		}
		temp = 0.;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__2 = temp, d__3 = (d__1 = vr_ref(jr, jc), abs(d__1))
				;
			temp = max(d__2,d__3);
/* L60: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
/* Computing MAX */
			d__3 = temp, d__4 = (d__1 = vr_ref(jr, jc), abs(d__1))
				 + (d__2 = vr_ref(jr, jc + 1), abs(d__2));
			temp = max(d__3,d__4);
/* L70: */
		    }
		}
		if (temp < safmin) {
		    goto L100;
		}
		temp = 1. / temp;
		if (alphai[jc] == 0.) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, jc) = vr_ref(jr, jc) * temp;
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, jc) = vr_ref(jr, jc) * temp;
			vr_ref(jr, jc + 1) = vr_ref(jr, jc + 1) * temp;
/* L90: */
		    }
		}
L100:
		;
	    }
	}

/*        End of eigenvector calculation */

    }

/*     Undo scaling in alpha, beta   

       Note: this does not give the alpha and beta for the unscaled   
       problem.   

       Un-scaling is limited to avoid underflow in alpha and beta   
       if they are significant. */

    i__1 = *n;
    for (jc = 1; jc <= i__1; ++jc) {
	absar = (d__1 = alphar[jc], abs(d__1));
	absai = (d__1 = alphai[jc], abs(d__1));
	absb = (d__1 = beta[jc], abs(d__1));
	salfar = anrm * alphar[jc];
	salfai = anrm * alphai[jc];
	sbeta = bnrm * beta[jc];
	ilimit = FALSE_;
	scale = 1.;

/*        Check for significant underflow in ALPHAI   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
	if (abs(salfai) < safmin && absai >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX */
	    d__1 = onepls * safmin, d__2 = anrm2 * absai;
	    scale = onepls * safmin / anrm1 / max(d__1,d__2);

	} else if (salfai == 0.) {

/*           If insignificant underflow in ALPHAI, then make the   
             conjugate eigenvalue real. */

	    if (alphai[jc] < 0. && jc > 1) {
		alphai[jc - 1] = 0.;
	    } else if (alphai[jc] > 0. && jc < *n) {
		alphai[jc + 1] = 0.;
	    }
	}

/*        Check for significant underflow in ALPHAR   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absai, d__1 = max(d__1,d__2), d__2 = eps *
		 absb;
	if (abs(salfar) < safmin && absar >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    d__3 = onepls * safmin, d__4 = anrm2 * absar;
	    d__1 = scale, d__2 = onepls * safmin / anrm1 / max(d__3,d__4);
	    scale = max(d__1,d__2);
	}

/*        Check for significant underflow in BETA   

   Computing MAX */
	d__1 = safmin, d__2 = eps * absar, d__1 = max(d__1,d__2), d__2 = eps *
		 absai;
	if (abs(sbeta) < safmin && absb >= max(d__1,d__2)) {
	    ilimit = TRUE_;
/* Computing MAX   
   Computing MAX */
	    d__3 = onepls * safmin, d__4 = bnrm2 * absb;
	    d__1 = scale, d__2 = onepls * safmin / bnrm1 / max(d__3,d__4);
	    scale = max(d__1,d__2);
	}

/*        Check for possible overflow when limiting scaling */

	if (ilimit) {
/* Computing MAX */
	    d__1 = abs(salfar), d__2 = abs(salfai), d__1 = max(d__1,d__2), 
		    d__2 = abs(sbeta);
	    temp = scale * safmin * max(d__1,d__2);
	    if (temp > 1.) {
		scale /= temp;
	    }
	    if (scale < 1.) {
		ilimit = FALSE_;
	    }
	}

/*        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary. */

	if (ilimit) {
	    salfar = scale * alphar[jc] * anrm;
	    salfai = scale * alphai[jc] * anrm;
	    sbeta = scale * beta[jc] * bnrm;
	}
	alphar[jc] = salfar;
	alphai[jc] = salfai;
	beta[jc] = sbeta;
/* L110: */
    }

L120:
    work[1] = (doublereal) lwkopt;

    return 0;

/*     End of DGEGV */

} /* dgegv_ */

#undef vr_ref
#undef vl_ref
#undef b_ref
#undef a_ref



logical lsame_(char *ca, char *cb)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    LSAME returns .TRUE. if CA is the same letter as CB regardless of   
    case.   

    Arguments   
    =========   

    CA      (input) CHARACTER*1   
    CB      (input) CHARACTER*1   
            CA and CB specify the single characters to be compared.   

   ===================================================================== 
  


       Test if the characters are equal */
    /* System generated locals */
    logical ret_val;
    /* Local variables */
    static integer inta, intb, zcode;


    ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
    if (ret_val) {
	return ret_val;
    }

/*     Now test for equivalence if both characters are alphabetic. */

    zcode = 'Z';

/*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime   
       machines, on which ICHAR returns a value with bit 8 set.   
       ICHAR('A') on Prime machines returns 193 which is the same as   
       ICHAR('A') on an EBCDIC machine. */

    inta = *(unsigned char *)ca;
    intb = *(unsigned char *)cb;

    if (zcode == 90 || zcode == 122) {

/*        ASCII is assumed - ZCODE is the ASCII code of either lower o
r   
          upper case 'Z'. */

	if (inta >= 97 && inta <= 122) {
	    inta += -32;
	}
	if (intb >= 97 && intb <= 122) {
	    intb += -32;
	}

    } else if (zcode == 233 || zcode == 169) {

/*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower
 or   
          upper case 'Z'. */

	if (inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta 
		>= 162 && inta <= 169) {
	    inta += 64;
	}
	if (intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb 
		>= 162 && intb <= 169) {
	    intb += 64;
	}

    } else if (zcode == 218 || zcode == 250) {

/*        ASCII is assumed, on Prime machines - ZCODE is the ASCII cod
e   
          plus 128 of either lower or upper case 'Z'. */

	if (inta >= 225 && inta <= 250) {
	    inta += -32;
	}
	if (intb >= 225 && intb <= 250) {
	    intb += -32;
	}
    }
    ret_val = inta == intb;

/*     RETURN   

       End of LSAME */

    return ret_val;
} /* lsame_ */


integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen 
	opts_len)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    ILAENV is called from the LAPACK routines to choose problem-dependent   
    parameters for the local environment.  See ISPEC for a description of   
    the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    This routine will not function correctly if it is converted to all   
    lower case.  Converting it to all upper case is allowed.   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies the parameter to be returned as the value of   
            ILAENV.   
            = 1: the optimal blocksize; if this value is 1, an unblocked   
                 algorithm will give the best performance.   
            = 2: the minimum block size for which the block routine   
                 should be used; if the usable block size is less than   
                 this value, an unblocked routine should be used.   
            = 3: the crossover point (in a block routine, for N less   
                 than this value, an unblocked routine should be used)   
            = 4: the number of shifts, used in the nonsymmetric   
                 eigenvalue routines   
            = 5: the minimum column dimension for blocking to be used;   
                 rectangular blocks must have dimension at least k by m,   
                 where k is given by ILAENV(2,...) and m by ILAENV(5,...)   
            = 6: the crossover point for the SVD (when reducing an m by n   
                 matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds   
                 this value, a QR factorization is used first to reduce   
                 the matrix to a triangular form.)   
            = 7: the number of processors   
            = 8: the crossover point for the multishift QR and QZ methods   
                 for nonsymmetric eigenvalue problems.   
            = 9: maximum size of the subproblems at the bottom of the   
                 computation tree in the divide-and-conquer algorithm   
                 (used by xGELSD and xGESDD)   
            =10: ieee NaN arithmetic can be trusted not to trap   
            =11: infinity arithmetic can be trusted not to trap   

    NAME    (input) CHARACTER*(*)   
            The name of the calling subroutine, in either upper case or   
            lower case.   

    OPTS    (input) CHARACTER*(*)   
            The character options to the subroutine NAME, concatenated   
            into a single character string.  For example, UPLO = 'U',   
            TRANS = 'T', and DIAG = 'N' for a triangular routine would   
            be specified as OPTS = 'UTN'.   

    N1      (input) INTEGER   
    N2      (input) INTEGER   
    N3      (input) INTEGER   
    N4      (input) INTEGER   
            Problem dimensions for the subroutine NAME; these may not all   
            be required.   

   (ILAENV) (output) INTEGER   
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if ILAENV = -k, the k-th argument had an illegal value.   

    Further Details   
    ===============   

    The following conventions have been used when calling ILAENV from the   
    LAPACK routines:   
    1)  OPTS is a concatenation of all of the character options to   
        subroutine NAME, in the same order that they appear in the   
        argument list for NAME, even if they are not used in determining   
        the value of the parameter specified by ISPEC.   
    2)  The problem dimensions N1, N2, N3, N4 are specified in the order   
        that they appear in the argument list for NAME.  N1 is used   
        first, N2 second, and so on, and unused problem dimensions are   
        passed a value of -1.   
    3)  The parameter value returned by ILAENV is checked for validity in   
        the calling subroutine.  For example, ILAENV is used to retrieve   
        the optimal blocksize for STRTRI as follows:   

        NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )   
        IF( NB.LE.1 ) NB = MAX( 1, N )   

    ===================================================================== */
    /* Table of constant values */
    static integer c__0 = 0;
    static real c_b162 = 0.f;
    static real c_b163 = 1.f;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer ret_val;
    /* Builtin functions   
       Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);
    /* Local variables */
    static integer i__;
    static logical cname, sname;
    static integer nbmin;
    static char c1[1], c2[2], c3[3], c4[2];
    static integer ic, nb;
    extern integer ieeeck_(integer *, real *, real *);
    static integer iz, nx;
    static char subnam[6];




    switch (*ispec) {
	case 1:  goto L100;
	case 2:  goto L100;
	case 3:  goto L100;
	case 4:  goto L400;
	case 5:  goto L500;
	case 6:  goto L600;
	case 7:  goto L700;
	case 8:  goto L800;
	case 9:  goto L900;
	case 10:  goto L1000;
	case 11:  goto L1100;
    }

/*     Invalid value for ISPEC */

    ret_val = -1;
    return ret_val;

L100:

/*     Convert NAME to upper case if the first character is lower case. */

    ret_val = 1;
    s_copy(subnam, name__, (ftnlen)6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

/*        ASCII character set */

	if (ic >= 97 && ic <= 122) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 97 && ic <= 122) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L10: */
	    }
	}

    } else if (iz == 233 || iz == 169) {

/*        EBCDIC character set */

	if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
		ic <= 169) {
	    *(unsigned char *)subnam = (char) (ic + 64);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
			162 && ic <= 169) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
		}
/* L20: */
	    }
	}

    } else if (iz == 218 || iz == 250) {

/*        Prime machines:  ASCII+128 */

	if (ic >= 225 && ic <= 250) {
	    *(unsigned char *)subnam = (char) (ic - 32);
	    for (i__ = 2; i__ <= 6; ++i__) {
		ic = *(unsigned char *)&subnam[i__ - 1];
		if (ic >= 225 && ic <= 250) {
		    *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
		}
/* L30: */
	    }
	}
    }

    *(unsigned char *)c1 = *(unsigned char *)subnam;
    sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
    cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
    if (! (cname || sname)) {
	return ret_val;
    }
    s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
    s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
    s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);

    switch (*ispec) {
	case 1:  goto L110;
	case 2:  goto L200;
	case 3:  goto L300;
    }

L110:

/*     ISPEC = 1:  block size   

       In these examples, separate code is provided for setting NB for   
       real and complex.  We assume that NB will take the same value in   
       single or double precision. */

    nb = 1;

    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
		"RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
		3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
		== 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 32;
	    } else {
		nb = 32;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	} else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 32;
	} else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 64;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nb = 32;
	    }
	}
    } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n4 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    } else {
		if (*n2 <= 64) {
		    nb = 1;
		} else {
		    nb = 32;
		}
	    }
	}
    } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nb = 64;
	    } else {
		nb = 64;
	    }
	}
    } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
	    nb = 1;
	}
    }
    ret_val = nb;
    return ret_val;

L200:

/*     ISPEC = 2:  minimum block size */

    nbmin = 2;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	} else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 2;
	    } else {
		nbmin = 2;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nbmin = 8;
	    } else {
		nbmin = 8;
	    }
	} else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nbmin = 2;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	} else if (*(unsigned char *)c3 == 'M') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nbmin = 2;
	    }
	}
    }
    ret_val = nbmin;
    return ret_val;

L300:

/*     ISPEC = 3:  crossover point */

    nx = 0;
    if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
		ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
		ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
		 {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	} else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
	    if (sname) {
		nx = 128;
	    } else {
		nx = 128;
	    }
	}
    } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
	if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
	if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
	    nx = 32;
	}
    } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
	if (*(unsigned char *)c3 == 'G') {
	    if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
		    (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
		    ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
		     0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
		    c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
		    ftnlen)2, (ftnlen)2) == 0) {
		nx = 128;
	    }
	}
    }
    ret_val = nx;
    return ret_val;

L400:

/*     ISPEC = 4:  number of shifts (used by xHSEQR) */

    ret_val = 6;
    return ret_val;

L500:

/*     ISPEC = 5:  minimum column dimension (not used) */

    ret_val = 2;
    return ret_val;

L600:

/*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

    ret_val = (integer) ((real) min(*n1,*n2) * 1.6f);
    return ret_val;

L700:

/*     ISPEC = 7:  number of processors (not used) */

    ret_val = 1;
    return ret_val;

L800:

/*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    ret_val = 50;
    return ret_val;

L900:

/*     ISPEC = 9:  maximum size of the subproblems at the bottom of the   
                   computation tree in the divide-and-conquer algorithm   
                   (used by xGELSD and xGESDD) */

    ret_val = 25;
    return ret_val;

L1000:

/*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__0, &c_b162, &c_b163);
    }
    return ret_val;

L1100:

/*     ISPEC = 11: infinity arithmetic can be trusted not to trap   

       ILAENV = 0 */
    ret_val = 1;
    if (ret_val == 1) {
	ret_val = ieeeck_(&c__1, &c_b162, &c_b163);
    }
    return ret_val;

/*     End of ILAENV */

} /* ilaenv_ */


/* Subroutine */ int xerbla_(char *srname, integer *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    XERBLA  is an error handler for the LAPACK routines.   
    It is called by an LAPACK routine if an input parameter has an   
    invalid value.  A message is printed and execution stops.   

    Installers may consider modifying the STOP statement in order to   
    call system-specific exception-handling facilities.   

    Arguments   
    =========   

    SRNAME  (input) CHARACTER*6   
            The name of the routine which called XERBLA.   

    INFO    (input) INTEGER   
            The position of the invalid parameter in the parameter list   

            of the calling routine.   

   ===================================================================== 
*/

    printf("** On entry to %6s, parameter number %2li had an illegal value\n",
		srname, *info);

/*     End of XERBLA */

    return 0;
} /* xerbla_ */


doublereal dlamch_(char *cmach)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMCH determines double precision machine parameters.   

    Arguments   
    =========   

    CMACH   (input) CHARACTER*1   
            Specifies the value to be returned by DLAMCH:   
            = 'E' or 'e',   DLAMCH := eps   
            = 'S' or 's ,   DLAMCH := sfmin   
            = 'B' or 'b',   DLAMCH := base   
            = 'P' or 'p',   DLAMCH := eps*base   
            = 'N' or 'n',   DLAMCH := t   
            = 'R' or 'r',   DLAMCH := rnd   
            = 'M' or 'm',   DLAMCH := emin   
            = 'U' or 'u',   DLAMCH := rmin   
            = 'L' or 'l',   DLAMCH := emax   
            = 'O' or 'o',   DLAMCH := rmax   

            where   

            eps   = relative machine precision   
            sfmin = safe minimum, such that 1/sfmin does not overflow   
            base  = base of the machine   
            prec  = eps*base   
            t     = number of (base) digits in the mantissa   
            rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise   
            emin  = minimum exponent before (gradual) underflow   
            rmin  = underflow threshold - base**(emin-1)   
            emax  = largest exponent before overflow   
            rmax  = overflow threshold  - (base**emax)*(1-eps)   

   ===================================================================== 
*/
/* >>Start of File<<   
       Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    doublereal ret_val;
    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    /* Local variables */
    static doublereal base;
    static integer beta;
    static doublereal emin, prec, emax;
    static integer imin, imax;
    static logical lrnd;
    static doublereal rmin, rmax, t, rmach;
    extern logical lsame_(char *, char *);
    static doublereal small, sfmin;
    extern /* Subroutine */ int dlamc2_(integer *, integer *, logical *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *);
    static integer it;
    static doublereal rnd, eps;



    if (first) {
	first = FALSE_;
	dlamc2_(&beta, &it, &lrnd, &eps, &imin, &rmin, &imax, &rmax);
	base = (doublereal) beta;
	t = (doublereal) it;
	if (lrnd) {
	    rnd = 1.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1) / 2;
	} else {
	    rnd = 0.;
	    i__1 = 1 - it;
	    eps = pow_di(&base, &i__1);
	}
	prec = eps * base;
	emin = (doublereal) imin;
	emax = (doublereal) imax;
	sfmin = rmin;
	small = 1. / rmax;
	if (small >= sfmin) {

/*           Use SMALL plus a bit, to avoid the possibility of rou
nding   
             causing overflow when computing  1/sfmin. */

	    sfmin = small * (eps + 1.);
	}
    }

    if (lsame_(cmach, "E")) {
	rmach = eps;
    } else if (lsame_(cmach, "S")) {
	rmach = sfmin;
    } else if (lsame_(cmach, "B")) {
	rmach = base;
    } else if (lsame_(cmach, "P")) {
	rmach = prec;
    } else if (lsame_(cmach, "N")) {
	rmach = t;
    } else if (lsame_(cmach, "R")) {
	rmach = rnd;
    } else if (lsame_(cmach, "M")) {
	rmach = emin;
    } else if (lsame_(cmach, "U")) {
	rmach = rmin;
    } else if (lsame_(cmach, "L")) {
	rmach = emax;
    } else if (lsame_(cmach, "O")) {
	rmach = rmax;
    }

    ret_val = rmach;
    return ret_val;

/*     End of DLAMCH */

} /* dlamch_ */


/* Subroutine */ int dlamc1_(integer *beta, integer *t, logical *rnd, logical 
	*ieee1)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC1 determines the machine parameters given by BETA, T, RND, and   
    IEEE1.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    IEEE1   (output) LOGICAL   
            Specifies whether rounding appears to be done in the IEEE   
            'round to nearest' style.   

    Further Details   
    ===============   

    The routine is based on the routine  ENVRON  by Malcolm and   
    incorporates suggestions by Gentleman and Marovich. See   

       Malcolm M. A. (1972) Algorithms to reveal properties of   
          floating-point arithmetic. Comms. of the ACM, 15, 949-951.   

       Gentleman W. M. and Marovich S. B. (1974) More on algorithms   
          that reveal properties of floating point arithmetic units.   
          Comms. of the ACM, 17, 276-277.   

   ===================================================================== 
*/
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    doublereal d__1, d__2;
    /* Local variables */
    static logical lrnd;
    static doublereal a, b, c, f;
    static integer lbeta;
    static doublereal savec;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;
    static doublereal t1, t2;
    static integer lt;
    static doublereal one, qtr;



    if (first) {
	first = FALSE_;
	one = 1.;

/*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BE
TA,   
          IEEE1, T and RND.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are  stored and not held in registers, 
 or   
          are not affected by optimizers.   

          Compute  a = 2.0**m  with the  smallest positive integer m s
uch   
          that   

             fl( a + 1.0 ) = a. */

	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L10:
	if (c == one) {
	    a *= 2;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L10;
	}
/* +       END WHILE   

          Now compute  b = 2.0**m  with the smallest positive integer 
m   
          such that   

             fl( a + b ) .gt. a. */

	b = 1.;
	c = dlamc3_(&a, &b);

/* +       WHILE( C.EQ.A )LOOP */
L20:
	if (c == a) {
	    b *= 2;
	    c = dlamc3_(&a, &b);
	    goto L20;
	}
/* +       END WHILE   

          Now compute the base.  a and c  are neighbouring floating po
int   
          numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and
 so   
          their difference is beta. Adding 0.25 to c is to ensure that
 it   
          is truncated to beta and not ( beta - 1 ). */

	qtr = one / 4;
	savec = c;
	d__1 = -a;
	c = dlamc3_(&c, &d__1);
	lbeta = (integer) (c + qtr);

/*        Now determine whether rounding or chopping occurs,  by addin
g a   
          bit  less  than  beta/2  and a  bit  more  than  beta/2  to 
 a. */

	b = (doublereal) lbeta;
	d__1 = b / 2;
	d__2 = -b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (c == a) {
	    lrnd = TRUE_;
	} else {
	    lrnd = FALSE_;
	}
	d__1 = b / 2;
	d__2 = b / 100;
	f = dlamc3_(&d__1, &d__2);
	c = dlamc3_(&f, &a);
	if (lrnd && c == a) {
	    lrnd = FALSE_;
	}

/*        Try and decide whether rounding is done in the  IEEE  'round
 to   
          nearest' style. B/2 is half a unit in the last place of the 
two   
          numbers A and SAVEC. Furthermore, A is even, i.e. has last  
bit   
          zero, and SAVEC is odd. Thus adding B/2 to A should not  cha
nge   
          A, but adding B/2 to SAVEC should change SAVEC. */

	d__1 = b / 2;
	t1 = dlamc3_(&d__1, &a);
	d__1 = b / 2;
	t2 = dlamc3_(&d__1, &savec);
	lieee1 = t1 == a && t2 > savec && lrnd;

/*        Now find  the  mantissa, t.  It should  be the  integer part
 of   
          log to the base beta of a,  however it is safer to determine
  t   
          by powering.  So we find t as the smallest positive integer 
for   
          which   

             fl( beta**t + 1.0 ) = 1.0. */

	lt = 0;
	a = 1.;
	c = 1.;

/* +       WHILE( C.EQ.ONE )LOOP */
L30:
	if (c == one) {
	    ++lt;
	    a *= lbeta;
	    c = dlamc3_(&a, &one);
	    d__1 = -a;
	    c = dlamc3_(&c, &d__1);
	    goto L30;
	}
/* +       END WHILE */

    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *ieee1 = lieee1;
    return 0;

/*     End of DLAMC1 */

} /* dlamc1_ */


/* Subroutine */ int dlamc2_(integer *beta, integer *t, logical *rnd, 
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax, 
	doublereal *rmax)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC2 determines the machine parameters specified in its argument   
    list.   

    Arguments   
    =========   

    BETA    (output) INTEGER   
            The base of the machine.   

    T       (output) INTEGER   
            The number of ( BETA ) digits in the mantissa.   

    RND     (output) LOGICAL   
            Specifies whether proper rounding  ( RND = .TRUE. )  or   
            chopping  ( RND = .FALSE. )  occurs in addition. This may not 
  
            be a reliable guide to the way in which the machine performs 
  
            its arithmetic.   

    EPS     (output) DOUBLE PRECISION   
            The smallest positive number such that   

               fl( 1.0 - EPS ) .LT. 1.0,   

            where fl denotes the computed value.   

    EMIN    (output) INTEGER   
            The minimum exponent before (gradual) underflow occurs.   

    RMIN    (output) DOUBLE PRECISION   
            The smallest normalized number for the machine, given by   
            BASE**( EMIN - 1 ), where  BASE  is the floating point value 
  
            of BETA.   

    EMAX    (output) INTEGER   
            The maximum exponent before overflow occurs.   

    RMAX    (output) DOUBLE PRECISION   
            The largest positive number for the machine, given by   
            BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point 
  
            value of BETA.   

    Further Details   
    ===============   

    The computation of  EPS  is based on a routine PARANOIA by   
    W. Kahan of the University of California at Berkeley.   

   ===================================================================== 
*/
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* Initialized data */
    static logical first = TRUE_;
    static logical iwarn = FALSE_;
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5;
    /* Builtin functions */
    double pow_di(doublereal *, integer *);
    /* Local variables */
    static logical ieee;
    static doublereal half;
    static logical lrnd;
    static doublereal leps, zero, a, b, c;
    static integer i, lbeta;
    static doublereal rbase;
    static integer lemin, lemax, gnmin;
    static doublereal small;
    static integer gpmin;
    static doublereal third, lrmin, lrmax, sixth;
    extern /* Subroutine */ int dlamc1_(integer *, integer *, logical *, 
	    logical *);
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static logical lieee1;
    extern /* Subroutine */ int dlamc4_(integer *, doublereal *, integer *), 
	    dlamc5_(integer *, integer *, integer *, logical *, integer *, 
	    doublereal *);
    static integer lt, ngnmin, ngpmin;
    static doublereal one, two;



    if (first) {
	first = FALSE_;
	zero = 0.;
	one = 1.;
	two = 2.;

/*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values
 of   
          BETA, T, RND, EPS, EMIN and RMIN.   

          Throughout this routine  we use the function  DLAMC3  to ens
ure   
          that relevant values are stored  and not held in registers, 
 or   
          are not affected by optimizers.   

          DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. 
*/

	dlamc1_(&lbeta, &lt, &lrnd, &lieee1);

/*        Start to find EPS. */

	b = (doublereal) lbeta;
	i__1 = -lt;
	a = pow_di(&b, &i__1);
	leps = a;

/*        Try some tricks to see whether or not this is the correct  E
PS. */

	b = two / 3;
	half = one / 2;
	d__1 = -half;
	sixth = dlamc3_(&b, &d__1);
	third = dlamc3_(&sixth, &sixth);
	d__1 = -half;
	b = dlamc3_(&third, &d__1);
	b = dlamc3_(&b, &sixth);
	b = abs(b);
	if (b < leps) {
	    b = leps;
	}

	leps = 1.;

/* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
L10:
	if (leps > b && b > zero) {
	    leps = b;
	    d__1 = half * leps;
/* Computing 5th power */
	    d__3 = two, d__4 = d__3, d__3 *= d__3;
/* Computing 2nd power */
	    d__5 = leps;
	    d__2 = d__4 * (d__3 * d__3) * (d__5 * d__5);
	    c = dlamc3_(&d__1, &d__2);
	    d__1 = -c;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    d__1 = -b;
	    c = dlamc3_(&half, &d__1);
	    b = dlamc3_(&half, &c);
	    goto L10;
	}
/* +       END WHILE */

	if (a < leps) {
	    leps = a;
	}

/*        Computation of EPS complete.   

          Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3
)).   
          Keep dividing  A by BETA until (gradual) underflow occurs. T
his   
          is detected when we cannot recover the previous A. */

	rbase = one / lbeta;
	small = one;
	for (i = 1; i <= 3; ++i) {
	    d__1 = small * rbase;
	    small = dlamc3_(&d__1, &zero);
/* L20: */
	}
	a = dlamc3_(&one, &small);
	dlamc4_(&ngpmin, &one, &lbeta);
	d__1 = -one;
	dlamc4_(&ngnmin, &d__1, &lbeta);
	dlamc4_(&gpmin, &a, &lbeta);
	d__1 = -a;
	dlamc4_(&gnmin, &d__1, &lbeta);
	ieee = FALSE_;

	if (ngpmin == ngnmin && gpmin == gnmin) {
	    if (ngpmin == gpmin) {
		lemin = ngpmin;
/*            ( Non twos-complement machines, no gradual under
flow;   
                e.g.,  VAX ) */
	    } else if (gpmin - ngpmin == 3) {
		lemin = ngpmin - 1 + lt;
		ieee = TRUE_;
/*            ( Non twos-complement machines, with gradual und
erflow;   
                e.g., IEEE standard followers ) */
	    } else {
		lemin = min(ngpmin,gpmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if (ngpmin == gpmin && ngnmin == gnmin) {
	    if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1) {
		lemin = max(ngpmin,ngnmin);
/*            ( Twos-complement machines, no gradual underflow
;   
                e.g., CYBER 205 ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else if ((i__1 = ngpmin - ngnmin, abs(i__1)) == 1 && gpmin == gnmin)
		 {
	    if (gpmin - min(ngpmin,ngnmin) == 3) {
		lemin = max(ngpmin,ngnmin) - 1 + lt;
/*            ( Twos-complement machines with gradual underflo
w;   
                no known machine ) */
	    } else {
		lemin = min(ngpmin,ngnmin);
/*            ( A guess; no known machine ) */
		iwarn = TRUE_;
	    }

	} else {
/* Computing MIN */
	    i__1 = min(ngpmin,ngnmin), i__1 = min(i__1,gpmin);
	    lemin = min(i__1,gnmin);
/*         ( A guess; no known machine ) */
	    iwarn = TRUE_;
	}
/* **   
   Comment out this if block if EMIN is ok */
	if (iwarn) {
	    first = TRUE_;
	    printf("\n\n WARNING. The value EMIN may be incorrect:- ");
	    printf("EMIN = %8li\n",lemin);
	    printf("If, after inspection, the value EMIN looks acceptable");
            printf("please comment out \n the IF block as marked within the"); 
            printf("code of routine DLAMC2, \n otherwise supply EMIN"); 
            printf("explicitly.\n");
	}
/* **   

          Assume IEEE arithmetic if we found denormalised  numbers abo
ve,   
          or if arithmetic seems to round in the  IEEE style,  determi
ned   
          in routine DLAMC1. A true IEEE machine should have both  thi
ngs   
          true; however, faulty machines may have one or the other. */

	ieee = ieee || lieee1;

/*        Compute  RMIN by successive division by  BETA. We could comp
ute   
          RMIN as BASE**( EMIN - 1 ),  but some machines underflow dur
ing   
          this computation. */

	lrmin = 1.;
	i__1 = 1 - lemin;
	for (i = 1; i <= 1-lemin; ++i) {
	    d__1 = lrmin * rbase;
	    lrmin = dlamc3_(&d__1, &zero);
/* L30: */
	}

/*        Finally, call DLAMC5 to compute EMAX and RMAX. */

	dlamc5_(&lbeta, &lt, &lemin, &ieee, &lemax, &lrmax);
    }

    *beta = lbeta;
    *t = lt;
    *rnd = lrnd;
    *eps = leps;
    *emin = lemin;
    *rmin = lrmin;
    *emax = lemax;
    *rmax = lrmax;

    return 0;


/*     End of DLAMC2 */

} /* dlamc2_ */


doublereal dlamc3_(doublereal *a, doublereal *b)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC3  is intended to force  A  and  B  to be stored prior to doing 
  
    the addition of  A  and  B ,  for use in situations where optimizers 
  
    might hold one of these in a register.   

    Arguments   
    =========   

    A, B    (input) DOUBLE PRECISION   
            The values A and B.   

   ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    doublereal ret_val;



    ret_val = *a + *b;

    return ret_val;

/*     End of DLAMC3 */

} /* dlamc3_ */


/* Subroutine */ int dlamc4_(integer *emin, doublereal *start, integer *base)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC4 is a service routine for DLAMC2.   

    Arguments   
    =========   

    EMIN    (output) EMIN   
            The minimum exponent before (gradual) underflow, computed by 
  
            setting A = START and dividing by BASE until the previous A   
            can not be recovered.   

    START   (input) DOUBLE PRECISION   
            The starting point for determining EMIN.   

    BASE    (input) INTEGER   
            The base of the machine.   

   ===================================================================== 
*/
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Local variables */
    static doublereal zero, a;
    static integer i;
    static doublereal rbase, b1, b2, c1, c2, d1, d2;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static doublereal one;



    a = *start;
    one = 1.;
    rbase = one / *base;
    zero = 0.;
    *emin = 1;
    d__1 = a * rbase;
    b1 = dlamc3_(&d__1, &zero);
    c1 = a;
    c2 = a;
    d1 = a;
    d2 = a;
/* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.   
      $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
L10:
    if (c1 == a && c2 == a && d1 == a && d2 == a) {
	--(*emin);
	a = b1;
	d__1 = a / *base;
	b1 = dlamc3_(&d__1, &zero);
	d__1 = b1 * *base;
	c1 = dlamc3_(&d__1, &zero);
	d1 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d1 += b1;
/* L20: */
	}
	d__1 = a * rbase;
	b2 = dlamc3_(&d__1, &zero);
	d__1 = b2 / rbase;
	c2 = dlamc3_(&d__1, &zero);
	d2 = zero;
	i__1 = *base;
	for (i = 1; i <= *base; ++i) {
	    d2 += b2;
/* L30: */
	}
	goto L10;
    }
/* +    END WHILE */

    return 0;

/*     End of DLAMC4 */

} /* dlamc4_ */


/* Subroutine */ int dlamc5_(integer *beta, integer *p, integer *emin, 
	logical *ieee, integer *emax, doublereal *rmax)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAMC5 attempts to compute RMAX, the largest machine floating-point   
    number, without overflow.  It assumes that EMAX + abs(EMIN) sum   
    approximately to a power of 2.  It will fail on machines where this   
    assumption does not hold, for example, the Cyber 205 (EMIN = -28625, 
  
    EMAX = 28718).  It will also fail if the value supplied for EMIN is   
    too large (i.e. too close to zero), probably with overflow.   

    Arguments   
    =========   

    BETA    (input) INTEGER   
            The base of floating-point arithmetic.   

    P       (input) INTEGER   
            The number of base BETA digits in the mantissa of a   
            floating-point value.   

    EMIN    (input) INTEGER   
            The minimum exponent before (gradual) underflow.   

    IEEE    (input) LOGICAL   
            A logical flag specifying whether or not the arithmetic   
            system is thought to comply with the IEEE standard.   

    EMAX    (output) INTEGER   
            The largest exponent before overflow   

    RMAX    (output) DOUBLE PRECISION   
            The largest machine floating-point number.   

   ===================================================================== 
  


       First compute LEXP and UEXP, two powers of 2 that bound   
       abs(EMIN). We then assume that EMAX + abs(EMIN) will sum   
       approximately to the bound that is closest to abs(EMIN).   
       (EMAX is the exponent of the required number RMAX). */
    /* Table of constant values */
    static doublereal c_b5 = 0.;
    
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Local variables */
    static integer lexp;
    static doublereal oldy;
    static integer uexp, i;
    static doublereal y, z;
    static integer nbits;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    static doublereal recbas;
    static integer exbits, expsum, try__;



    lexp = 1;
    exbits = 1;
L10:
    try__ = lexp << 1;
    if (try__ <= -(*emin)) {
	lexp = try__;
	++exbits;
	goto L10;
    }
    if (lexp == -(*emin)) {
	uexp = lexp;
    } else {
	uexp = try__;
	++exbits;
    }

/*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater   
       than or equal to EMIN. EXBITS is the number of bits needed to   
       store the exponent. */

    if (uexp + *emin > -lexp - *emin) {
	expsum = lexp << 1;
    } else {
	expsum = uexp << 1;
    }

/*     EXPSUM is the exponent range, approximately equal to   
       EMAX - EMIN + 1 . */

    *emax = expsum + *emin - 1;
    nbits = exbits + 1 + *p;

/*     NBITS is the total number of bits needed to store a   
       floating-point number. */

    if (nbits % 2 == 1 && *beta == 2) {

/*        Either there are an odd number of bits used to store a   
          floating-point number, which is unlikely, or some bits are 
  
          not used in the representation of numbers, which is possible
,   
          (e.g. Cray machines) or the mantissa has an implicit bit,   
          (e.g. IEEE machines, Dec Vax machines), which is perhaps the
   
          most likely. We have to assume the last alternative.   
          If this is true, then we need to reduce EMAX by one because 
  
          there must be some way of representing zero in an implicit-b
it   
          system. On machines like Cray, we are reducing EMAX by one 
  
          unnecessarily. */

	--(*emax);
    }

    if (*ieee) {

/*        Assume we are on an IEEE machine which reserves one exponent
   
          for infinity and NaN. */

	--(*emax);
    }

/*     Now create RMAX, the largest machine number, which should   
       be equal to (1.0 - BETA**(-P)) * BETA**EMAX .   

       First compute 1.0 - BETA**(-P), being careful that the   
       result is less than 1.0 . */

    recbas = 1. / *beta;
    z = *beta - 1.;
    y = 0.;
    i__1 = *p;
    for (i = 1; i <= *p; ++i) {
	z *= recbas;
	if (y < 1.) {
	    oldy = y;
	}
	y = dlamc3_(&y, &z);
/* L20: */
    }
    if (y >= 1.) {
	y = oldy;
    }

/*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = *emax;
    for (i = 1; i <= *emax; ++i) {
	d__1 = y * *beta;
	y = dlamc3_(&d__1, &c_b5);
/* L30: */
    }

    *rmax = y;
    return 0;

/*     End of DLAMC5 */

} /* dlamc5_ */


doublereal dlange_(char *norm, integer *m, integer *n, doublereal *a, integer 
	*lda, doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANGE  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    real matrix A.   

    Description   
    ===========   

    DLANGE returns the value   

       DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANGE as described   
            above.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.  When M = 0,   
            DLANGE is set to zero.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.  When N = 0,   
            DLANGE is set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(M,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= M when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer i__, j;
    static doublereal scale;
    extern logical lsame_(char *, char *);
    static doublereal value;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal sum;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (min(*m,*n) == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = a_ref(i__, j), abs(d__1));
		value = max(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum += (d__1 = a_ref(i__, j), abs(d__1));
/* L30: */
	    }
	    value = max(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[i__] += (d__1 = a_ref(i__, j), abs(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *m;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = value, d__2 = work[i__];
	    value = max(d__1,d__2);
/* L80: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    dlassq_(m, &a_ref(1, j), &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANGE */

} /* dlange_ */

#undef a_ref



/* Subroutine */ int dlascl_(char *type__, integer *kl, integer *ku, 
	doublereal *cfrom, doublereal *cto, integer *m, integer *n, 
	doublereal *a, integer *lda, integer *info)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLASCL multiplies the M by N real matrix A by the real scalar   
    CTO/CFROM.  This is done without over/underflow as long as the final   
    result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that   
    A may be full, upper triangular, lower triangular, upper Hessenberg,   
    or banded.   

    Arguments   
    =========   

    TYPE    (input) CHARACTER*1   
            TYPE indices the storage type of the input matrix.   
            = 'G':  A is a full matrix.   
            = 'L':  A is a lower triangular matrix.   
            = 'U':  A is an upper triangular matrix.   
            = 'H':  A is an upper Hessenberg matrix.   
            = 'B':  A is a symmetric band matrix with lower bandwidth KL   
                    and upper bandwidth KU and with the only the lower   
                    half stored.   
            = 'Q':  A is a symmetric band matrix with lower bandwidth KL   
                    and upper bandwidth KU and with the only the upper   
                    half stored.   
            = 'Z':  A is a band matrix with lower bandwidth KL and upper   
                    bandwidth KU.   

    KL      (input) INTEGER   
            The lower bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    KU      (input) INTEGER   
            The upper bandwidth of A.  Referenced only if TYPE = 'B',   
            'Q' or 'Z'.   

    CFROM   (input) DOUBLE PRECISION   
    CTO     (input) DOUBLE PRECISION   
            The matrix A is multiplied by CTO/CFROM. A(I,J) is computed   
            without over/underflow if the final result CTO*A(I,J)/CFROM   
            can be represented without over/underflow.  CFROM must be   
            nonzero.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)   
            The matrix to be multiplied by CTO/CFROM.  See TYPE for the   
            storage type.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    INFO    (output) INTEGER   
            0  - successful exit   
            <0 - if INFO = -i, the i-th argument had an illegal value.   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    /* Local variables */
    static logical done;
    static doublereal ctoc;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer itype, k1, k2, k3, k4;
    static doublereal cfrom1;
    extern doublereal dlamch_(char *);
    static doublereal cfromc;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal bignum, smlnum, mul, cto1;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    *info = 0;

    if (lsame_(type__, "G")) {
	itype = 0;
    } else if (lsame_(type__, "L")) {
	itype = 1;
    } else if (lsame_(type__, "U")) {
	itype = 2;
    } else if (lsame_(type__, "H")) {
	itype = 3;
    } else if (lsame_(type__, "B")) {
	itype = 4;
    } else if (lsame_(type__, "Q")) {
	itype = 5;
    } else if (lsame_(type__, "Z")) {
	itype = 6;
    } else {
	itype = -1;
    }

    if (itype == -1) {
	*info = -1;
    } else if (*cfrom == 0.) {
	*info = -4;
    } else if (*m < 0) {
	*info = -6;
    } else if (*n < 0 || itype == 4 && *n != *m || itype == 5 && *n != *m) {
	*info = -7;
    } else if (itype <= 3 && *lda < max(1,*m)) {
	*info = -9;
    } else if (itype >= 4) {
/* Computing MAX */
	i__1 = *m - 1;
	if (*kl < 0 || *kl > max(i__1,0)) {
	    *info = -2;
	} else /* if(complicated condition) */ {
/* Computing MAX */
	    i__1 = *n - 1;
	    if (*ku < 0 || *ku > max(i__1,0) || (itype == 4 || itype == 5) && 
		    *kl != *ku) {
		*info = -3;
	    } else if (itype == 4 && *lda < *kl + 1 || itype == 5 && *lda < *
		    ku + 1 || itype == 6 && *lda < (*kl << 1) + *ku + 1) {
		*info = -9;
	    }
	}
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DLASCL", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *m == 0) {
	return 0;
    }

/*     Get machine parameters */

    smlnum = dlamch_("S");
    bignum = 1. / smlnum;

    cfromc = *cfrom;
    ctoc = *cto;

L10:
    cfrom1 = cfromc * smlnum;
    cto1 = ctoc / bignum;
    if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
	mul = smlnum;
	done = FALSE_;
	cfromc = cfrom1;
    } else if (abs(cto1) > abs(cfromc)) {
	mul = bignum;
	done = FALSE_;
	ctoc = cto1;
    } else {
	mul = ctoc / cfromc;
	done = TRUE_;
    }

    if (itype == 0) {

/*        Full matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L20: */
	    }
/* L30: */
	}

    } else if (itype == 1) {

/*        Lower triangular matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L40: */
	    }
/* L50: */
	}

    } else if (itype == 2) {

/*        Upper triangular matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = min(j,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L60: */
	    }
/* L70: */
	}

    } else if (itype == 3) {

/*        Upper Hessenberg matrix */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = j + 1;
	    i__2 = min(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L80: */
	    }
/* L90: */
	}

    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

	k3 = *kl + 1;
	k4 = *n + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = k3, i__4 = k4 - j;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L100: */
	    }
/* L110: */
	}

    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

	k1 = *ku + 2;
	k3 = *ku + 1;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__2 = k1 - j;
	    i__3 = k3;
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L120: */
	    }
/* L130: */
	}

    } else if (itype == 6) {

/*        Band matrix */

	k1 = *kl + *ku + 2;
	k2 = *kl + 1;
	k3 = (*kl << 1) + *ku + 1;
	k4 = *kl + *ku + 1 + *m;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	    i__3 = k1 - j;
/* Computing MIN */
	    i__4 = k3, i__5 = k4 - j;
	    i__2 = min(i__4,i__5);
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
		a_ref(i__, j) = a_ref(i__, j) * mul;
/* L140: */
	    }
/* L150: */
	}

    }

    if (! done) {
	goto L10;
    }

    return 0;

/*     End of DLASCL */

} /* dlascl_ */

#undef a_ref



/* Subroutine */ int dggbal_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi, 
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAL balances a pair of general real matrices (A,B).  This   
    involves, first, permuting A and B by similarity transformations to   
    isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N   
    elements on the diagonal; and second, applying a diagonal similarity   
    transformation to rows and columns ILO to IHI to make the rows   
    and columns as close in norm as possible. Both steps are optional.   

    Balancing may reduce the 1-norm of the matrices, and improve the   
    accuracy of the computed eigenvalues and/or eigenvectors in the   
    generalized eigenvalue problem A*x = lambda*B*x.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the operations to be performed on A and B:   
            = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0   
                    and RSCALE(I) = 1.0 for i = 1,...,N.   
            = 'P':  permute only;   
            = 'S':  scale only;   
            = 'B':  both permute and scale.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the input matrix A.   
            On exit,  A is overwritten by the balanced matrix.   
            If JOB = 'N', A is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the input matrix B.   
            On exit,  B is overwritten by the balanced matrix.   
            If JOB = 'N', B is not referenced.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= max(1,N).   

    ILO     (output) INTEGER   
    IHI     (output) INTEGER   
            ILO and IHI are set to integers such that on exit   
            A(i,j) = 0 and B(i,j) = 0 if i > j and   
            j = 1,...,ILO-1 or i = IHI+1,...,N.   
            If JOB = 'N' or 'S', ILO = 1 and IHI = N.   

    LSCALE  (output) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the left side of A and B.  If P(j) is the index of the   
            row interchanged with row j, and D(j)   
            is the scaling factor applied to row j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    RSCALE  (output) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and scaling factors applied   
            to the right side of A and B.  If P(j) is the index of the   
            column interchanged with column j, and D(j)   
            is the scaling factor applied to column j, then   
              LSCALE(j) = P(j)    for J = 1,...,ILO-1   
                        = D(j)    for J = ILO,...,IHI   
                        = P(j)    for J = IHI+1,...,N.   
            The order in which the interchanges are made is N to IHI+1,   
            then 1 to ILO-1.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. WARD, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b34 = 10.;
    static doublereal c_b70 = .5;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3;
    /* Builtin functions */
    double d_lg10(doublereal *), d_sign(doublereal *, doublereal *), pow_di(
	    doublereal *, integer *);
    /* Local variables */
    static integer lcab;
    static doublereal beta, coef;
    static integer irab, lrab;
    static doublereal basl, cmax;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal coef2, coef5;
    static integer i__, j, k, l, m;
    static doublereal gamma, t, alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    static doublereal sfmin, sfmax;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer iflow;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer kount, jc;
    static doublereal ta, tb, tc;
    extern doublereal dlamch_(char *);
    static integer ir, it;
    static doublereal ew;
    static integer nr;
    static doublereal pgamma;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static integer lsfmin, lsfmax, ip1, jp1, lm1;
    static doublereal cab, rab, ewc, cor, sum;
    static integer nrp2, icab;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --lscale;
    --rscale;
    --work;

    /* Function Body */
    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(job, "S") 
	    && ! lsame_(job, "B")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGBAL", &i__1);
	return 0;
    }

    k = 1;
    l = *n;

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }

    if (lsame_(job, "N")) {
	*ilo = 1;
	*ihi = *n;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lscale[i__] = 1.;
	    rscale[i__] = 1.;
/* L10: */
	}
	return 0;
    }

    if (k == l) {
	*ilo = 1;
	*ihi = 1;
	lscale[1] = 1.;
	rscale[1] = 1.;
	return 0;
    }

    if (lsame_(job, "S")) {
	goto L190;
    }

    goto L30;

/*     Permute the matrices A and B to isolate the eigenvalues.   

       Find row with one nonzero in columns 1 through L */

L20:
    l = lm1;
    if (l != 1) {
	goto L30;
    }

    rscale[1] = 1.;
    lscale[1] = 1.;
    goto L190;

L30:
    lm1 = l - 1;
    for (i__ = l; i__ >= 1; --i__) {
	i__1 = lm1;
	for (j = 1; j <= i__1; ++j) {
	    jp1 = j + 1;
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L50;
	    }
/* L40: */
	}
	j = l;
	goto L70;

L50:
	i__1 = l;
	for (j = jp1; j <= i__1; ++j) {
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L80;
	    }
/* L60: */
	}
	j = jp1 - 1;

L70:
	m = l;
	iflow = 1;
	goto L160;
L80:
	;
    }
    goto L100;

/*     Find column with one nonzero in rows K through N */

L90:
    ++k;

L100:
    i__1 = l;
    for (j = k; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = k; i__ <= i__2; ++i__) {
	    ip1 = i__ + 1;
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L120;
	    }
/* L110: */
	}
	i__ = l;
	goto L140;
L120:
	i__2 = l;
	for (i__ = ip1; i__ <= i__2; ++i__) {
	    if (a_ref(i__, j) != 0. || b_ref(i__, j) != 0.) {
		goto L150;
	    }
/* L130: */
	}
	i__ = ip1 - 1;
L140:
	m = k;
	iflow = 2;
	goto L160;
L150:
	;
    }
    goto L190;

/*     Permute rows M and I */

L160:
    lscale[m] = (doublereal) i__;
    if (i__ == m) {
	goto L170;
    }
    i__1 = *n - k + 1;
    dswap_(&i__1, &a_ref(i__, k), lda, &a_ref(m, k), lda);
    i__1 = *n - k + 1;
    dswap_(&i__1, &b_ref(i__, k), ldb, &b_ref(m, k), ldb);

/*     Permute columns M and J */

L170:
    rscale[m] = (doublereal) j;
    if (j == m) {
	goto L180;
    }
    dswap_(&l, &a_ref(1, j), &c__1, &a_ref(1, m), &c__1);
    dswap_(&l, &b_ref(1, j), &c__1, &b_ref(1, m), &c__1);

L180:
    switch (iflow) {
	case 1:  goto L20;
	case 2:  goto L90;
    }

L190:
    *ilo = k;
    *ihi = l;

    if (*ilo == *ihi) {
	return 0;
    }

    if (lsame_(job, "P")) {
	return 0;
    }

/*     Balance the submatrix in rows ILO to IHI. */

    nr = *ihi - *ilo + 1;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	rscale[i__] = 0.;
	lscale[i__] = 0.;

	work[i__] = 0.;
	work[i__ + *n] = 0.;
	work[i__ + (*n << 1)] = 0.;
	work[i__ + *n * 3] = 0.;
	work[i__ + (*n << 2)] = 0.;
	work[i__ + *n * 5] = 0.;
/* L200: */
    }

/*     Compute right side vector in resulting linear equations */

    basl = d_lg10(&c_b34);
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *ihi;
	for (j = *ilo; j <= i__2; ++j) {
	    tb = b_ref(i__, j);
	    ta = a_ref(i__, j);
	    if (ta == 0.) {
		goto L210;
	    }
	    d__1 = abs(ta);
	    ta = d_lg10(&d__1) / basl;
L210:
	    if (tb == 0.) {
		goto L220;
	    }
	    d__1 = abs(tb);
	    tb = d_lg10(&d__1) / basl;
L220:
	    work[i__ + (*n << 2)] = work[i__ + (*n << 2)] - ta - tb;
	    work[j + *n * 5] = work[j + *n * 5] - ta - tb;
/* L230: */
	}
/* L240: */
    }

    coef = 1. / (doublereal) (nr << 1);
    coef2 = coef * coef;
    coef5 = coef2 * .5;
    nrp2 = nr + 2;
    beta = 0.;
    it = 1;

/*     Start generalized conjugate gradient iteration */

L250:

    gamma = ddot_(&nr, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1) + ddot_(&nr, &work[*ilo + *n * 5], &c__1, &work[*ilo + *
	    n * 5], &c__1);

    ew = 0.;
    ewc = 0.;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	ew += work[i__ + (*n << 2)];
	ewc += work[i__ + *n * 5];
/* L260: */
    }

/* Computing 2nd power */
    d__1 = ew;
/* Computing 2nd power */
    d__2 = ewc;
/* Computing 2nd power */
    d__3 = ew - ewc;
    gamma = coef * gamma - coef2 * (d__1 * d__1 + d__2 * d__2) - coef5 * (
	    d__3 * d__3);
    if (gamma == 0.) {
	goto L350;
    }
    if (it != 1) {
	beta = gamma / pgamma;
    }
    t = coef5 * (ewc - ew * 3.);
    tc = coef5 * (ew - ewc * 3.);

    dscal_(&nr, &beta, &work[*ilo], &c__1);
    dscal_(&nr, &beta, &work[*ilo + *n], &c__1);

    daxpy_(&nr, &coef, &work[*ilo + (*n << 2)], &c__1, &work[*ilo + *n], &
	    c__1);
    daxpy_(&nr, &coef, &work[*ilo + *n * 5], &c__1, &work[*ilo], &c__1);

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	work[i__] += tc;
	work[i__ + *n] += t;
/* L270: */
    }

/*     Apply matrix to vector */

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (j = *ilo; j <= i__2; ++j) {
	    if (a_ref(i__, j) == 0.) {
		goto L280;
	    }
	    ++kount;
	    sum += work[j];
L280:
	    if (b_ref(i__, j) == 0.) {
		goto L290;
	    }
	    ++kount;
	    sum += work[j];
L290:
	    ;
	}
	work[i__ + (*n << 1)] = (doublereal) kount * work[i__ + *n] + sum;
/* L300: */
    }

    i__1 = *ihi;
    for (j = *ilo; j <= i__1; ++j) {
	kount = 0;
	sum = 0.;
	i__2 = *ihi;
	for (i__ = *ilo; i__ <= i__2; ++i__) {
	    if (a_ref(i__, j) == 0.) {
		goto L310;
	    }
	    ++kount;
	    sum += work[i__ + *n];
L310:
	    if (b_ref(i__, j) == 0.) {
		goto L320;
	    }
	    ++kount;
	    sum += work[i__ + *n];
L320:
	    ;
	}
	work[j + *n * 3] = (doublereal) kount * work[j] + sum;
/* L330: */
    }

    sum = ddot_(&nr, &work[*ilo + *n], &c__1, &work[*ilo + (*n << 1)], &c__1) 
	    + ddot_(&nr, &work[*ilo], &c__1, &work[*ilo + *n * 3], &c__1);
    alpha = gamma / sum;

/*     Determine correction to current iteration */

    cmax = 0.;
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	cor = alpha * work[i__ + *n];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	lscale[i__] += cor;
	cor = alpha * work[i__];
	if (abs(cor) > cmax) {
	    cmax = abs(cor);
	}
	rscale[i__] += cor;
/* L340: */
    }
    if (cmax < .5) {
	goto L350;
    }

    d__1 = -alpha;
    daxpy_(&nr, &d__1, &work[*ilo + (*n << 1)], &c__1, &work[*ilo + (*n << 2)]
	    , &c__1);
    d__1 = -alpha;
    daxpy_(&nr, &d__1, &work[*ilo + *n * 3], &c__1, &work[*ilo + *n * 5], &
	    c__1);

    pgamma = gamma;
    ++it;
    if (it <= nrp2) {
	goto L250;
    }

/*     End generalized conjugate gradient iteration */

L350:
    sfmin = dlamch_("S");
    sfmax = 1. / sfmin;
    lsfmin = (integer) (d_lg10(&sfmin) / basl + 1.);
    lsfmax = (integer) (d_lg10(&sfmax) / basl);
    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *n - *ilo + 1;
	irab = idamax_(&i__2, &a_ref(i__, *ilo), lda);
	rab = (d__1 = a_ref(i__, irab + *ilo - 1), abs(d__1));
	i__2 = *n - *ilo + 1;
	irab = idamax_(&i__2, &b_ref(i__, *ilo), lda);
/* Computing MAX */
	d__2 = rab, d__3 = (d__1 = b_ref(i__, irab + *ilo - 1), abs(d__1));
	rab = max(d__2,d__3);
	d__1 = rab + sfmin;
	lrab = (integer) (d_lg10(&d__1) / basl + 1.);
	ir = (integer) (lscale[i__] + d_sign(&c_b70, &lscale[i__]));
/* Computing MIN */
	i__2 = max(ir,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lrab;
	ir = min(i__2,i__3);
	lscale[i__] = pow_di(&c_b34, &ir);
	icab = idamax_(ihi, &a_ref(1, i__), &c__1);
	cab = (d__1 = a_ref(icab, i__), abs(d__1));
	icab = idamax_(ihi, &b_ref(1, i__), &c__1);
/* Computing MAX */
	d__2 = cab, d__3 = (d__1 = b_ref(icab, i__), abs(d__1));
	cab = max(d__2,d__3);
	d__1 = cab + sfmin;
	lcab = (integer) (d_lg10(&d__1) / basl + 1.);
	jc = (integer) (rscale[i__] + d_sign(&c_b70, &rscale[i__]));
/* Computing MIN */
	i__2 = max(jc,lsfmin), i__2 = min(i__2,lsfmax), i__3 = lsfmax - lcab;
	jc = min(i__2,i__3);
	rscale[i__] = pow_di(&c_b34, &jc);
/* L360: */
    }

/*     Row scaling of matrices A and B */

    i__1 = *ihi;
    for (i__ = *ilo; i__ <= i__1; ++i__) {
	i__2 = *n - *ilo + 1;
	dscal_(&i__2, &lscale[i__], &a_ref(i__, *ilo), lda);
	i__2 = *n - *ilo + 1;
	dscal_(&i__2, &lscale[i__], &b_ref(i__, *ilo), ldb);
/* L370: */
    }

/*     Column scaling of matrices A and B */

    i__1 = *ihi;
    for (j = *ilo; j <= i__1; ++j) {
	dscal_(ihi, &rscale[j], &a_ref(1, j), &c__1);
	dscal_(ihi, &rscale[j], &b_ref(1, j), &c__1);
/* L380: */
    }

    return 0;

/*     End of DGGBAL */

} /* dggbal_ */

#undef b_ref
#undef a_ref



/* Subroutine */ int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DGEQRF computes a QR factorization of a real M-by-N matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and above the diagonal of the array   
            contain the min(M,N)-by-N upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal,   
            with the array TAU, represent the orthogonal matrix Q as a   
            product of min(m,n) elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is   
            the optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),   
    and tau in TAU(i).   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    /* Local variables */
    static integer i__, k, nbmin, iinfo;
    extern /* Subroutine */ int dgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    static integer ib, nb;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nx;
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "DGEQRF", " ", m, n, &c_n1, &c_n1, (ftnlen)6, (ftnlen)
	    1);
    lwkopt = *n * nb;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQRF", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    k = min(*m,*n);
    if (k == 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGEQRF", " ", m, n, &c_n1, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", m, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__3 = k - i__ + 1;
	    ib = min(i__3,nb);

/*           Compute the QR factorization of the current block   
             A(i:m,i:i+ib-1) */

	    i__3 = *m - i__ + 1;
	    dgeqr2_(&i__3, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
		    iinfo);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *m - i__ + 1;
		dlarft_("Forward", "Columnwise", &i__3, &ib, &a_ref(i__, i__),
			 lda, &tau[i__], &work[1], &ldwork);

/*              Apply H' to A(i:m,i+ib:n) from the left */

		i__3 = *m - i__ + 1;
		i__4 = *n - i__ - ib + 1;
		dlarfb_("Left", "Transpose", "Forward", "Columnwise", &i__3, &
			i__4, &ib, &a_ref(i__, i__), lda, &work[1], &ldwork, &
			a_ref(i__, i__ + ib), lda, &work[ib + 1], &ldwork);
	    }
/* L10: */
	}
    } else {
	i__ = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i__ <= k) {
	i__2 = *m - i__ + 1;
	i__1 = *n - i__ + 1;
	dgeqr2_(&i__2, &i__1, &a_ref(i__, i__), lda, &tau[i__], &work[1], &
		iinfo);
    }

    work[1] = (doublereal) iws;
    return 0;

/*     End of DGEQRF */

} /* dgeqrf_ */

#undef a_ref



/* Subroutine */ int dormqr_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DORMQR overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(1) H(2) . . . H(k)   

    as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by   
            DGEQRF in the first k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= max(1,M);   
            if SIDE = 'R', LDA >= max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= max(1,N);   
            if SIDE = 'R', LWORK >= max(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__2 = 2;
    static integer c__65 = 65;
    
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2, i__3[2], i__4, 
	    i__5;
    char ch__1[2];
    /* Builtin functions   
       Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    /* Local variables */
    static logical left;
    static integer i__;
    static doublereal t[4160]	/* was [65][64] */;
    extern logical lsame_(char *, char *);
    static integer nbmin, iinfo, i1, i2, i3;
    extern /* Subroutine */ int dorm2r_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer ib, ic, jc, nb, mi, ni;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nq, nw;
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical notran;
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    lquery = *lwork == -1;

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! left && ! lsame_(side, "R")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < max(1,nq)) {
	*info = -7;
    } else if (*ldc < max(1,*m)) {
	*info = -10;
    } else if (*lwork < max(1,nw) && ! lquery) {
	*info = -12;
    }

    if (*info == 0) {

/*        Determine the block size.  NB may be at most NBMAX, where NBMAX   
          is used to define the local array T.   

   Computing MIN   
   Writing concatenation */
	i__3[0] = 1, a__1[0] = side;
	i__3[1] = 1, a__1[1] = trans;
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	i__1 = 64, i__2 = ilaenv_(&c__1, "DORMQR", ch__1, m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)2);
	nb = min(i__1,i__2);
	lwkopt = max(1,nw) * nb;
	work[1] = (doublereal) lwkopt;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMQR", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX   
   Writing concatenation */
	    i__3[0] = 1, a__1[0] = side;
	    i__3[1] = 1, a__1[1] = trans;
	    s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQR", ch__1, m, n, k, &c_n1, (
		    ftnlen)6, (ftnlen)2);
	    nbmin = max(i__1,i__2);
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */

	dorm2r_(side, trans, m, n, k, &a[a_offset], lda, &tau[1], &c__[
		c_offset], ldc, &work[1], &iinfo);
    } else {

/*        Use blocked code */

	if (left && ! notran || ! left && notran) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__1 = i2;
	i__2 = i3;
	for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
	    i__4 = nb, i__5 = *k - i__ + 1;
	    ib = min(i__4,i__5);

/*           Form the triangular factor of the block reflector   
             H = H(i) H(i+1) . . . H(i+ib-1) */

	    i__4 = nq - i__ + 1;
	    dlarft_("Forward", "Columnwise", &i__4, &ib, &a_ref(i__, i__), 
		    lda, &tau[i__], t, &c__65);
	    if (left) {

/*              H or H' is applied to C(i:m,1:n) */

		mi = *m - i__ + 1;
		ic = i__;
	    } else {

/*              H or H' is applied to C(1:m,i:n) */

		ni = *n - i__ + 1;
		jc = i__;
	    }

/*           Apply H or H' */

	    dlarfb_(side, trans, "Forward", "Columnwise", &mi, &ni, &ib, &
		    a_ref(i__, i__), lda, t, &c__65, &c___ref(ic, jc), ldc, &
		    work[1], &ldwork);
/* L10: */
	}
    }
    work[1] = (doublereal) lwkopt;
    return 0;

/*     End of DORMQR */

} /* dormqr_ */

#undef c___ref
#undef a_ref



/* Subroutine */ int dlaset_(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower   
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper   
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) DOUBLE PRECISION   
            The constant to which the offdiagonal elements are to be set.   

    BETA    (input) DOUBLE PRECISION   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows:   

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    if (lsame_(uplo, "U")) {

/*        Set the strictly upper triangular or trapezoidal part of the   
          array to ALPHA. */

	i__1 = *n;
	for (j = 2; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = min(i__3,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (lsame_(uplo, "L")) {

/*        Set the strictly lower triangular or trapezoidal part of the   
          array to ALPHA. */

	i__1 = min(*m,*n);
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first min(M,N) diagonal elements to BETA. */

    i__1 = min(*m,*n);
    for (i__ = 1; i__ <= i__1; ++i__) {
	a_ref(i__, i__) = *beta;
/* L70: */
    }

    return 0;

/*     End of DLASET */

} /* dlaset_ */

#undef a_ref



/* Subroutine */ int dlacpy_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLACPY copies all or part of a two-dimensional matrix A to another   
    matrix B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be copied to B.   
            = 'U':      Upper triangular part   
            = 'L':      Lower triangular part   
            Otherwise:  All of the matrix A   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The m by n matrix A.  If UPLO = 'U', only the upper triangle   
            or trapezoid is accessed; if UPLO = 'L', only the lower   
            triangle or trapezoid is accessed.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    B       (output) DOUBLE PRECISION array, dimension (LDB,N)   
            On exit, B = A in the locations specified by UPLO.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,M).   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    if (lsame_(uplo, "U")) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = min(j,*m);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(uplo, "L")) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = j; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L30: */
	    }
/* L40: */
	}
    } else {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = a_ref(i__, j);
/* L50: */
	    }
/* L60: */
	}
    }
    return 0;

/*     End of DLACPY */

} /* dlacpy_ */

#undef b_ref
#undef a_ref



/* Subroutine */ int dorgqr_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DORGQR generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the first N columns of a product of K elementary   
    reflectors of order M   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= max(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static integer c_n1 = -1;
    static integer c__3 = 3;
    static integer c__2 = 2;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, j, l, nbmin, iinfo;
    extern /* Subroutine */ int dorg2r_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *);
    static integer ib, nb, ki, kk;
    extern /* Subroutine */ int dlarfb_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer nx;
    extern /* Subroutine */ int dlarft_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *), xerbla_(char *, integer *);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static integer ldwork, lwkopt;
    static logical lquery;
    static integer iws;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    nb = ilaenv_(&c__1, "DORGQR", " ", m, n, k, &c_n1, (ftnlen)6, (ftnlen)1);
    lwkopt = max(1,*n) * nb;
    work[1] = (doublereal) lwkopt;
    lquery = *lwork == -1;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGQR", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1] = 1.;
	return 0;
    }

    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQR", " ", m, n, k, &c_n1, (
		ftnlen)6, (ftnlen)1);
	nx = max(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduce NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", m, n, k, &c_n1,
			 (ftnlen)6, (ftnlen)1);
		nbmin = max(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk columns are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = min(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

	i__1 = *n;
	for (j = kk + 1; j <= i__1; ++j) {
	    i__2 = kk;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		a_ref(i__, j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;
	dorg2r_(&i__1, &i__2, &i__3, &a_ref(kk + 1, kk + 1), lda, &tau[kk + 1]
		, &work[1], &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i__ = ki + 1; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i__ + 1;
	    ib = min(i__2,i__3);
	    if (i__ + ib <= *n) {

/*              Form the triangular factor of the block reflector   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *m - i__ + 1;
		dlarft_("Forward", "Columnwise", &i__2, &ib, &a_ref(i__, i__),
			 lda, &tau[i__], &work[1], &ldwork);

/*              Apply H to A(i:m,i+ib:n) from the left */

		i__2 = *m - i__ + 1;
		i__3 = *n - i__ - ib + 1;
		dlarfb_("Left", "No transpose", "Forward", "Columnwise", &
			i__2, &i__3, &ib, &a_ref(i__, i__), lda, &work[1], &
			ldwork, &a_ref(i__, i__ + ib), lda, &work[ib + 1], &
			ldwork);
	    }

/*           Apply H to rows i:m of current block */

	    i__2 = *m - i__ + 1;
	    dorg2r_(&i__2, &ib, &ib, &a_ref(i__, i__), lda, &tau[i__], &work[
		    1], &iinfo);

/*           Set rows 1:i-1 of current block to zero */

	    i__2 = i__ + ib - 1;
	    for (j = i__; j <= i__2; ++j) {
		i__3 = i__ - 1;
		for (l = 1; l <= i__3; ++l) {
		    a_ref(l, j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    work[1] = (doublereal) iws;
    return 0;

/*     End of DORGQR */

} /* dorgqr_ */

#undef a_ref



/* Subroutine */ int dgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGHRD reduces a pair of real matrices (A,B) to generalized upper   
    Hessenberg form using orthogonal transformations, where A is a   
    general matrix and B is upper triangular:  Q' * A * Z = H and   
    Q' * B * Z = T, where H is upper Hessenberg, T is upper triangular,   
    and Q and Z are orthogonal, and ' means transpose.   

    The orthogonal matrices Q and Z are determined as products of Givens   
    rotations.  They may either be formed explicitly, or they may be   
    postmultiplied into input matrices Q1 and Z1, so that   

         Q1 * A * Z1' = (Q1*Q) * H * (Z1*Z)'   
         Q1 * B * Z1' = (Q1*Q) * T * (Z1*Z)'   

    Arguments   
    =========   

    COMPQ   (input) CHARACTER*1   
            = 'N': do not compute Q;   
            = 'I': Q is initialized to the unit matrix, and the   
                   orthogonal matrix Q is returned;   
            = 'V': Q must contain an orthogonal matrix Q1 on entry,   
                   and the product Q1*Q is returned.   

    COMPZ   (input) CHARACTER*1   
            = 'N': do not compute Z;   
            = 'I': Z is initialized to the unit matrix, and the   
                   orthogonal matrix Z is returned;   
            = 'V': Z must contain an orthogonal matrix Z1 on entry,   
                   and the product Z1*Z is returned.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows and   
            columns 1:ILO-1 and IHI+1:N.  ILO and IHI are normally set   
            by a previous call to DGGBAL; otherwise they should be set   
            to 1 and N respectively.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the N-by-N general matrix to be reduced.   
            On exit, the upper triangle and the first subdiagonal of A   
            are overwritten with the upper Hessenberg matrix H, and the   
            rest is set to zero.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,N).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.   
            On exit, the upper triangular matrix T = Q' B Z.  The   
            elements below the diagonal are set to zero.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max(1,N).   

    Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)   
            If COMPQ='N':  Q is not referenced.   
            If COMPQ='I':  on entry, Q need not be set, and on exit it   
                           contains the orthogonal matrix Q, where Q'   
                           is the product of the Givens transformations   
                           which are applied to A and B on the left.   
            If COMPQ='V':  on entry, Q must contain an orthogonal matrix   
                           Q1, and on exit this is overwritten by Q1*Q.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.   
            LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.   

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   
            If COMPZ='N':  Z is not referenced.   
            If COMPZ='I':  on entry, Z need not be set, and on exit it   
                           contains the orthogonal matrix Z, which is   
                           the product of the Givens transformations   
                           which are applied to A and B on the right.   
            If COMPZ='V':  on entry, Z must contain an orthogonal matrix   
                           Z1, and on exit this is overwritten by Z1*Z.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.   
            LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    This routine reduces A to Hessenberg and B to triangular form by   
    an unblocked reduction, as described in _Matrix_Computations_,   
    by Golub and Van Loan (Johns Hopkins Press.)   

    =====================================================================   


       Decode COMPQ   

       Parameter adjustments */
    /* Table of constant values */
    static doublereal c_b10 = 0.;
    static doublereal c_b11 = 1.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer jcol;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jrow;
    static doublereal c__, s;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *);
    static integer icompq, icompz;
    static logical ilq, ilz;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define q_ref(a_1,a_2) q[(a_2)*q_dim1 + a_1]
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;

    /* Function Body */
    if (lsame_(compq, "N")) {
	ilq = FALSE_;
	icompq = 1;
    } else if (lsame_(compq, "V")) {
	ilq = TRUE_;
	icompq = 2;
    } else if (lsame_(compq, "I")) {
	ilq = TRUE_;
	icompq = 3;
    } else {
	icompq = 0;
    }

/*     Decode COMPZ */

    if (lsame_(compz, "N")) {
	ilz = FALSE_;
	icompz = 1;
    } else if (lsame_(compz, "V")) {
	ilz = TRUE_;
	icompz = 2;
    } else if (lsame_(compz, "I")) {
	ilz = TRUE_;
	icompz = 3;
    } else {
	icompz = 0;
    }

/*     Test the input parameters. */

    *info = 0;
    if (icompq <= 0) {
	*info = -1;
    } else if (icompz <= 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -5;
    } else if (*lda < max(1,*n)) {
	*info = -7;
    } else if (*ldb < max(1,*n)) {
	*info = -9;
    } else if (ilq && *ldq < *n || *ldq < 1) {
	*info = -11;
    } else if (ilz && *ldz < *n || *ldz < 1) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGHRD", &i__1);
	return 0;
    }

/*     Initialize Q and Z if desired. */

    if (icompq == 3) {
	dlaset_("Full", n, n, &c_b10, &c_b11, &q[q_offset], ldq);
    }
    if (icompz == 3) {
	dlaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz);
    }

/*     Quick return if possible */

    if (*n <= 1) {
	return 0;
    }

/*     Zero out lower triangle of B */

    i__1 = *n - 1;
    for (jcol = 1; jcol <= i__1; ++jcol) {
	i__2 = *n;
	for (jrow = jcol + 1; jrow <= i__2; ++jrow) {
	    b_ref(jrow, jcol) = 0.;
/* L10: */
	}
/* L20: */
    }

/*     Reduce A and B */

    i__1 = *ihi - 2;
    for (jcol = *ilo; jcol <= i__1; ++jcol) {

	i__2 = jcol + 2;
	for (jrow = *ihi; jrow >= i__2; --jrow) {

/*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL) */

	    temp = a_ref(jrow - 1, jcol);
	    dlartg_(&temp, &a_ref(jrow, jcol), &c__, &s, &a_ref(jrow - 1, 
		    jcol));
	    a_ref(jrow, jcol) = 0.;
	    i__3 = *n - jcol;
	    drot_(&i__3, &a_ref(jrow - 1, jcol + 1), lda, &a_ref(jrow, jcol + 
		    1), lda, &c__, &s);
	    i__3 = *n + 2 - jrow;
	    drot_(&i__3, &b_ref(jrow - 1, jrow - 1), ldb, &b_ref(jrow, jrow - 
		    1), ldb, &c__, &s);
	    if (ilq) {
		drot_(n, &q_ref(1, jrow - 1), &c__1, &q_ref(1, jrow), &c__1, &
			c__, &s);
	    }

/*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1) */

	    temp = b_ref(jrow, jrow);
	    dlartg_(&temp, &b_ref(jrow, jrow - 1), &c__, &s, &b_ref(jrow, 
		    jrow));
	    b_ref(jrow, jrow - 1) = 0.;
	    drot_(ihi, &a_ref(1, jrow), &c__1, &a_ref(1, jrow - 1), &c__1, &
		    c__, &s);
	    i__3 = jrow - 1;
	    drot_(&i__3, &b_ref(1, jrow), &c__1, &b_ref(1, jrow - 1), &c__1, &
		    c__, &s);
	    if (ilz) {
		drot_(n, &z___ref(1, jrow), &c__1, &z___ref(1, jrow - 1), &
			c__1, &c__, &s);
	    }
/* L30: */
	}
/* L40: */
    }

    return 0;

/*     End of DGGHRD */

} /* dgghrd_ */

#undef z___ref
#undef q_ref
#undef b_ref
#undef a_ref



/* Subroutine */ int dhgeqz_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz, 
	doublereal *work, integer *lwork, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DHGEQZ implements a single-/double-shift version of the QZ method for   
    finding the generalized eigenvalues   

    w(j)=(ALPHAR(j) + i*ALPHAI(j))/BETAR(j)   of the equation   

         det( A - w(i) B ) = 0   

    In addition, the pair A,B may be reduced to generalized Schur form:   
    B is upper triangular, and A is block upper triangular, where the   
    diagonal blocks are either 1-by-1 or 2-by-2, the 2-by-2 blocks having   
    complex generalized eigenvalues (see the description of the argument   
    JOB.)   

    If JOB='S', then the pair (A,B) is simultaneously reduced to Schur   
    form by applying one orthogonal tranformation (usually called Q) on   
    the left and another (usually called Z) on the right.  The 2-by-2   
    upper-triangular diagonal blocks of B corresponding to 2-by-2 blocks   
    of A will be reduced to positive diagonal matrices.  (I.e.,   
    if A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) and   
    B(j+1,j+1) will be positive.)   

    If JOB='E', then at each iteration, the same transformations   
    are computed, but they are only applied to those parts of A and B   
    which are needed to compute ALPHAR, ALPHAI, and BETAR.   

    If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the orthogonal   
    transformations used to reduce (A,B) are accumulated into the arrays   
    Q and Z s.t.:   

         Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*   
         Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*   

    Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix   
         Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),   
         pp. 241--256.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            = 'E': compute only ALPHAR, ALPHAI, and BETA.  A and B will   
                   not necessarily be put into generalized Schur form.   
            = 'S': put A and B into generalized Schur form, as well   
                   as computing ALPHAR, ALPHAI, and BETA.   

    COMPQ   (input) CHARACTER*1   
            = 'N': do not modify Q.   
            = 'V': multiply the array Q on the right by the transpose of   
                   the orthogonal tranformation that is applied to the   
                   left side of A and B to reduce them to Schur form.   
            = 'I': like COMPQ='V', except that Q will be initialized to   
                   the identity first.   

    COMPZ   (input) CHARACTER*1   
            = 'N': do not modify Z.   
            = 'V': multiply the array Z on the right by the orthogonal   
                   tranformation that is applied to the right side of   
                   A and B to reduce them to Schur form.   
            = 'I': like COMPZ='V', except that Z will be initialized to   
                   the identity first.   

    N       (input) INTEGER   
            The order of the matrices A, B, Q, and Z.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            It is assumed that A is already upper triangular in rows and   
            columns 1:ILO-1 and IHI+1:N.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)   
            On entry, the N-by-N upper Hessenberg matrix A.  Elements   
            below the subdiagonal must be zero.   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to generalized Schur form.   
            If JOB='E', then on exit A will have been destroyed.   
               The diagonal blocks will be correct, but the off-diagonal   
               portion will be meaningless.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max( 1, N ).   

    B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)   
            On entry, the N-by-N upper triangular matrix B.  Elements   
            below the diagonal must be zero.  2-by-2 blocks in B   
            corresponding to 2-by-2 blocks in A will be reduced to   
            positive diagonal form.  (I.e., if A(j+1,j) is non-zero,   
            then B(j+1,j)=B(j,j+1)=0 and B(j,j) and B(j+1,j+1) will be   
            positive.)   
            If JOB='S', then on exit A and B will have been   
               simultaneously reduced to Schur form.   
            If JOB='E', then on exit B will have been destroyed.   
               Elements corresponding to diagonal blocks of A will be   
               correct, but the off-diagonal portion will be meaningless.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= max( 1, N ).   

    ALPHAR  (output) DOUBLE PRECISION array, dimension (N)   
            ALPHAR(1:N) will be set to real parts of the diagonal   
            elements of A that would result from reducing A and B to   
            Schur form and then further reducing them both to triangular   
            form using unitary transformations s.t. the diagonal of B   
            was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block   
            (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=A(j,j).   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   

    ALPHAI  (output) DOUBLE PRECISION array, dimension (N)   
            ALPHAI(1:N) will be set to imaginary parts of the diagonal   
            elements of A that would result from reducing A and B to   
            Schur form and then further reducing them both to triangular   
            form using unitary transformations s.t. the diagonal of B   
            was non-negative real.  Thus, if A(j,j) is in a 1-by-1 block   
            (i.e., A(j+1,j)=A(j,j+1)=0), then ALPHAR(j)=0.   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   

    BETA    (output) DOUBLE PRECISION array, dimension (N)   
            BETA(1:N) will be set to the (real) diagonal elements of B   
            that would result from reducing A and B to Schur form and   
            then further reducing them both to triangular form using   
            unitary transformations s.t. the diagonal of B was   
            non-negative real.  Thus, if A(j,j) is in a 1-by-1 block   
            (i.e., A(j+1,j)=A(j,j+1)=0), then BETA(j)=B(j,j).   
            Note that the (real or complex) values   
            (ALPHAR(j) + i*ALPHAI(j))/BETA(j), j=1,...,N, are the   
            generalized eigenvalues of the matrix pencil A - wB.   
            (Note that BETA(1:N) will always be non-negative, and no   
            BETAI is necessary.)   

    Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)   
            If COMPQ='N', then Q will not be referenced.   
            If COMPQ='V' or 'I', then the transpose of the orthogonal   
               transformations which are applied to A and B on the left   
               will be applied to the array Q on the right.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q.  LDQ >= 1.   
            If COMPQ='V' or 'I', then LDQ >= N.   

    Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   
            If COMPZ='N', then Z will not be referenced.   
            If COMPZ='V' or 'I', then the orthogonal transformations   
               which are applied to A and B on the right will be applied   
               to the array Z on the right.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1.   
            If COMPZ='V' or 'I', then LDZ >= N.   

    WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)   
            On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= max(1,N).   

            If LWORK = -1, then a workspace query is assumed; the routine   
            only calculates the optimal size of the WORK array, returns   
            this value as the first entry of the WORK array, and no error   
            message related to LWORK is issued by XERBLA.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            = 1,...,N: the QZ iteration did not converge.  (A,B) is not   
                       in Schur form, but ALPHAR(i), ALPHAI(i), and   
                       BETA(i), i=INFO+1,...,N should be correct.   
            = N+1,...,2*N: the shift calculation failed.  (A,B) is not   
                       in Schur form, but ALPHAR(i), ALPHAI(i), and   
                       BETA(i), i=INFO-N+1,...,N should be correct.   
            > 2*N:     various "impossible" errors.   

    Further Details   
    ===============   

    Iteration counters:   

    JITER  -- counts iterations.   
    IITER  -- counts iterations run since ILAST was last   
              changed.  This is therefore reset only when a 1-by-1 or   
              2-by-2 block deflates off the bottom.   

    =====================================================================   

      $                     SAFETY = 1.0E+0 )   

       Decode JOB, COMPQ, COMPZ   

       Parameter adjustments */
    /* Table of constant values */
    static doublereal c_b12 = 0.;
    static doublereal c_b13 = 1.;
    static integer c__1 = 1;
    static integer c__3 = 3;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublereal ad11l, ad12l, ad21l, ad22l, ad32l, wabs, atol, btol, 
	    temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlag2_(
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal temp2, s1inv, c__;
    static integer j;
    static doublereal s, t, v[3], scale;
    extern logical lsame_(char *, char *);
    static integer iiter, ilast, jiter;
    static doublereal anorm, bnorm;
    static integer maxit;
    static doublereal tempi, tempr, s1, s2, u1, u2;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static logical ilazr2;
    static doublereal a11, a12, a21, a22, b11, b22, c12, c21;
    static integer jc;
    static doublereal an, bn, cl, cq, cr;
    static integer in;
    static doublereal ascale, bscale, u12, w11;
    static integer jr;
    static doublereal cz, sl, w12, w21, w22, wi;
    extern doublereal dlamch_(char *);
    static doublereal sr;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal vs, wr;
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal safmin;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal safmax;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static doublereal eshift;
    static logical ilschr;
    static doublereal b1a, b2a;
    static integer icompq, ilastm;
    static doublereal a1i;
    static integer ischur;
    static doublereal a2i, b1i;
    static logical ilazro;
    static integer icompz, ifirst;
    static doublereal b2i;
    static integer ifrstm;
    static doublereal a1r;
    static integer istart;
    static logical ilpivt;
    static doublereal a2r, b1r, b2r;
    static logical lquery;
    static doublereal wr2, ad11, ad12, ad21, ad22, c11i, c22i;
    static integer jch;
    static doublereal c11r, c22r, u12l;
    static logical ilq;
    static doublereal tau, sqi;
    static logical ilz;
    static doublereal ulp, sqr, szi, szr;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define q_ref(a_1,a_2) q[(a_2)*q_dim1 + a_1]
#define z___ref(a_1,a_2) z__[(a_2)*z_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    --alphar;
    --alphai;
    --beta;
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1 * 1;
    q -= q_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1 * 1;
    z__ -= z_offset;
    --work;

    /* Function Body */
    if (lsame_(job, "E")) {
	ilschr = FALSE_;
	ischur = 1;
    } else if (lsame_(job, "S")) {
	ilschr = TRUE_;
	ischur = 2;
    } else {
	ischur = 0;
    }

    if (lsame_(compq, "N")) {
	ilq = FALSE_;
	icompq = 1;
    } else if (lsame_(compq, "V")) {
	ilq = TRUE_;
	icompq = 2;
    } else if (lsame_(compq, "I")) {
	ilq = TRUE_;
	icompq = 3;
    } else {
	icompq = 0;
    }

    if (lsame_(compz, "N")) {
	ilz = FALSE_;
	icompz = 1;
    } else if (lsame_(compz, "V")) {
	ilz = TRUE_;
	icompz = 2;
    } else if (lsame_(compz, "I")) {
	ilz = TRUE_;
	icompz = 3;
    } else {
	icompz = 0;
    }

/*     Check Argument Values */

    *info = 0;
    work[1] = (doublereal) max(1,*n);
    lquery = *lwork == -1;
    if (ischur == 0) {
	*info = -1;
    } else if (icompq == 0) {
	*info = -2;
    } else if (icompz == 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1) {
	*info = -5;
    } else if (*ihi > *n || *ihi < *ilo - 1) {
	*info = -6;
    } else if (*lda < *n) {
	*info = -8;
    } else if (*ldb < *n) {
	*info = -10;
    } else if (*ldq < 1 || ilq && *ldq < *n) {
	*info = -15;
    } else if (*ldz < 1 || ilz && *ldz < *n) {
	*info = -17;
    } else if (*lwork < max(1,*n) && ! lquery) {
	*info = -19;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DHGEQZ", &i__1);
	return 0;
    } else if (lquery) {
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	work[1] = 1.;
	return 0;
    }

/*     Initialize Q and Z */

    if (icompq == 3) {
	dlaset_("Full", n, n, &c_b12, &c_b13, &q[q_offset], ldq);
    }
    if (icompz == 3) {
	dlaset_("Full", n, n, &c_b12, &c_b13, &z__[z_offset], ldz);
    }

/*     Machine Constants */

    in = *ihi + 1 - *ilo;
    safmin = dlamch_("S");
    safmax = 1. / safmin;
    ulp = dlamch_("E") * dlamch_("B");
    anorm = dlanhs_("F", &in, &a_ref(*ilo, *ilo), lda, &work[1]);
    bnorm = dlanhs_("F", &in, &b_ref(*ilo, *ilo), ldb, &work[1]);
/* Computing MAX */
    d__1 = safmin, d__2 = ulp * anorm;
    atol = max(d__1,d__2);
/* Computing MAX */
    d__1 = safmin, d__2 = ulp * bnorm;
    btol = max(d__1,d__2);
    ascale = 1. / max(safmin,anorm);
    bscale = 1. / max(safmin,bnorm);

/*     Set Eigenvalues IHI+1:N */

    i__1 = *n;
    for (j = *ihi + 1; j <= i__1; ++j) {
	if (b_ref(j, j) < 0.) {
	    if (ilschr) {
		i__2 = j;
		for (jr = 1; jr <= i__2; ++jr) {
		    a_ref(jr, j) = -a_ref(jr, j);
		    b_ref(jr, j) = -b_ref(jr, j);
/* L10: */
		}
	    } else {
		a_ref(j, j) = -a_ref(j, j);
		b_ref(j, j) = -b_ref(j, j);
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    z___ref(jr, j) = -z___ref(jr, j);
/* L20: */
		}
	    }
	}
	alphar[j] = a_ref(j, j);
	alphai[j] = 0.;
	beta[j] = b_ref(j, j);
/* L30: */
    }

/*     If IHI < ILO, skip QZ steps */

    if (*ihi < *ilo) {
	goto L380;
    }

/*     MAIN QZ ITERATION LOOP   

       Initialize dynamic indices   

       Eigenvalues ILAST+1:N have been found.   
          Column operations modify rows IFRSTM:whatever.   
          Row operations modify columns whatever:ILASTM.   

       If only eigenvalues are being computed, then   
          IFRSTM is the row of the last splitting row above row ILAST;   
          this is always at least ILO.   
       IITER counts iterations since the last eigenvalue was found,   
          to tell when to use an extraordinary shift.   
       MAXIT is the maximum number of QZ sweeps allowed. */

    ilast = *ihi;
    if (ilschr) {
	ifrstm = 1;
	ilastm = *n;
    } else {
	ifrstm = *ilo;
	ilastm = *ihi;
    }
    iiter = 0;
    eshift = 0.;
    maxit = (*ihi - *ilo + 1) * 30;

    i__1 = maxit;
    for (jiter = 1; jiter <= i__1; ++jiter) {

/*        Split the matrix if possible.   

          Two tests:   
             1: A(j,j-1)=0  or  j=ILO   
             2: B(j,j)=0 */

	if (ilast == *ilo) {

/*           Special case: j=ILAST */

	    goto L80;
	} else {
	    if ((d__1 = a_ref(ilast, ilast - 1), abs(d__1)) <= atol) {
		a_ref(ilast, ilast - 1) = 0.;
		goto L80;
	    }
	}

	if ((d__1 = b_ref(ilast, ilast), abs(d__1)) <= btol) {
	    b_ref(ilast, ilast) = 0.;
	    goto L70;
	}

/*        General case: j<ILAST */

	i__2 = *ilo;
	for (j = ilast - 1; j >= i__2; --j) {

/*           Test 1: for A(j,j-1)=0 or j=ILO */

	    if (j == *ilo) {
		ilazro = TRUE_;
	    } else {
		if ((d__1 = a_ref(j, j - 1), abs(d__1)) <= atol) {
		    a_ref(j, j - 1) = 0.;
		    ilazro = TRUE_;
		} else {
		    ilazro = FALSE_;
		}
	    }

/*           Test 2: for B(j,j)=0 */

	    if ((d__1 = b_ref(j, j), abs(d__1)) < btol) {
		b_ref(j, j) = 0.;

/*              Test 1a: Check for 2 consecutive small subdiagonals in A */

		ilazr2 = FALSE_;
		if (! ilazro) {
		    temp = (d__1 = a_ref(j, j - 1), abs(d__1));
		    temp2 = (d__1 = a_ref(j, j), abs(d__1));
		    tempr = max(temp,temp2);
		    if (tempr < 1. && tempr != 0.) {
			temp /= tempr;
			temp2 /= tempr;
		    }
		    if (temp * (ascale * (d__1 = a_ref(j + 1, j), abs(d__1))) 
			    <= temp2 * (ascale * atol)) {
			ilazr2 = TRUE_;
		    }
		}

/*              If both tests pass (1 & 2), i.e., the leading diagonal   
                element of B in the block is zero, split a 1x1 block off   
                at the top. (I.e., at the J-th row/column) The leading   
                diagonal element of the remainder can also be zero, so   
                this may have to be done repeatedly. */

		if (ilazro || ilazr2) {
		    i__3 = ilast - 1;
		    for (jch = j; jch <= i__3; ++jch) {
			temp = a_ref(jch, jch);
			dlartg_(&temp, &a_ref(jch + 1, jch), &c__, &s, &a_ref(
				jch, jch));
			a_ref(jch + 1, jch) = 0.;
			i__4 = ilastm - jch;
			drot_(&i__4, &a_ref(jch, jch + 1), lda, &a_ref(jch + 
				1, jch + 1), lda, &c__, &s);
			i__4 = ilastm - jch;
			drot_(&i__4, &b_ref(jch, jch + 1), ldb, &b_ref(jch + 
				1, jch + 1), ldb, &c__, &s);
			if (ilq) {
			    drot_(n, &q_ref(1, jch), &c__1, &q_ref(1, jch + 1)
				    , &c__1, &c__, &s);
			}
			if (ilazr2) {
			    a_ref(jch, jch - 1) = a_ref(jch, jch - 1) * c__;
			}
			ilazr2 = FALSE_;
			if ((d__1 = b_ref(jch + 1, jch + 1), abs(d__1)) >= 
				btol) {
			    if (jch + 1 >= ilast) {
				goto L80;
			    } else {
				ifirst = jch + 1;
				goto L110;
			    }
			}
			b_ref(jch + 1, jch + 1) = 0.;
/* L40: */
		    }
		    goto L70;
		} else {

/*                 Only test 2 passed -- chase the zero to B(ILAST,ILAST)   
                   Then process as in the case B(ILAST,ILAST)=0 */

		    i__3 = ilast - 1;
		    for (jch = j; jch <= i__3; ++jch) {
			temp = b_ref(jch, jch + 1);
			dlartg_(&temp, &b_ref(jch + 1, jch + 1), &c__, &s, &
				b_ref(jch, jch + 1));
			b_ref(jch + 1, jch + 1) = 0.;
			if (jch < ilastm - 1) {
			    i__4 = ilastm - jch - 1;
			    drot_(&i__4, &b_ref(jch, jch + 2), ldb, &b_ref(
				    jch + 1, jch + 2), ldb, &c__, &s);
			}
			i__4 = ilastm - jch + 2;
			drot_(&i__4, &a_ref(jch, jch - 1), lda, &a_ref(jch + 
				1, jch - 1), lda, &c__, &s);
			if (ilq) {
			    drot_(n, &q_ref(1, jch), &c__1, &q_ref(1, jch + 1)
				    , &c__1, &c__, &s);
			}
			temp = a_ref(jch + 1, jch);
			dlartg_(&temp, &a_ref(jch + 1, jch - 1), &c__, &s, &
				a_ref(jch + 1, jch));
			a_ref(jch + 1, jch - 1) = 0.;
			i__4 = jch + 1 - ifrstm;
			drot_(&i__4, &a_ref(ifrstm, jch), &c__1, &a_ref(
				ifrstm, jch - 1), &c__1, &c__, &s);
			i__4 = jch - ifrstm;
			drot_(&i__4, &b_ref(ifrstm, jch), &c__1, &b_ref(
				ifrstm, jch - 1), &c__1, &c__, &s);
			if (ilz) {
			    drot_(n, &z___ref(1, jch), &c__1, &z___ref(1, jch 
				    - 1), &c__1, &c__, &s);
			}
/* L50: */
		    }
		    goto L70;
		}
	    } else if (ilazro) {

/*              Only test 1 passed -- work on J:ILAST */

		ifirst = j;
		goto L110;
	    }

/*           Neither test passed -- try next J   

   L60: */
	}

/*        (Drop-through is "impossible") */

	*info = *n + 1;
	goto L420;

/*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a   
          1x1 block. */

L70:
	temp = a_ref(ilast, ilast);
	dlartg_(&temp, &a_ref(ilast, ilast - 1), &c__, &s, &a_ref(ilast, 
		ilast));
	a_ref(ilast, ilast - 1) = 0.;
	i__2 = ilast - ifrstm;
	drot_(&i__2, &a_ref(ifrstm, ilast), &c__1, &a_ref(ifrstm, ilast - 1), 
		&c__1, &c__, &s);
	i__2 = ilast - ifrstm;
	drot_(&i__2, &b_ref(ifrstm, ilast), &c__1, &b_ref(ifrstm, ilast - 1), 
		&c__1, &c__, &s);
	if (ilz) {
	    drot_(n, &z___ref(1, ilast), &c__1, &z___ref(1, ilast - 1), &c__1,
		     &c__, &s);
	}

/*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,   
                                and BETA */

L80:
	if (b_ref(ilast, ilast) < 0.) {
	    if (ilschr) {
		i__2 = ilast;
		for (j = ifrstm; j <= i__2; ++j) {
		    a_ref(j, ilast) = -a_ref(j, ilast);
		    b_ref(j, ilast) = -b_ref(j, ilast);
/* L90: */
		}
	    } else {
		a_ref(ilast, ilast) = -a_ref(ilast, ilast);
		b_ref(ilast, ilast) = -b_ref(ilast, ilast);
	    }
	    if (ilz) {
		i__2 = *n;
		for (j = 1; j <= i__2; ++j) {
		    z___ref(j, ilast) = -z___ref(j, ilast);
/* L100: */
		}
	    }
	}
	alphar[ilast] = a_ref(ilast, ilast);
	alphai[ilast] = 0.;
	beta[ilast] = b_ref(ilast, ilast);

/*        Go to next block -- exit if finished. */

	--ilast;
	if (ilast < *ilo) {
	    goto L380;
	}

/*        Reset counters */

	iiter = 0;
	eshift = 0.;
	if (! ilschr) {
	    ilastm = ilast;
	    if (ifrstm > ilast) {
		ifrstm = *ilo;
	    }
	}
	goto L350;

/*        QZ step   

          This iteration only involves rows/columns IFIRST:ILAST. We   
          assume IFIRST < ILAST, and that the diagonal of B is non-zero. */

L110:
	++iiter;
	if (! ilschr) {
	    ifrstm = ifirst;
	}

/*        Compute single shifts.   

          At this point, IFIRST < ILAST, and the diagonal elements of   
          B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in   
          magnitude) */

	if (iiter / 10 * 10 == iiter) {

/*           Exceptional shift.  Chosen for no particularly good reason.   
             (Single shift only.) */

	    if ((doublereal) maxit * safmin * (d__1 = a_ref(ilast - 1, ilast),
		     abs(d__1)) < (d__2 = b_ref(ilast - 1, ilast - 1), abs(
		    d__2))) {
		eshift += a_ref(ilast - 1, ilast) / b_ref(ilast - 1, ilast - 
			1);
	    } else {
		eshift += 1. / (safmin * (doublereal) maxit);
	    }
	    s1 = 1.;
	    wr = eshift;

	} else {

/*           Shifts based on the generalized eigenvalues of the   
             bottom-right 2x2 block of A and B. The first eigenvalue   
             returned by DLAG2 is the Wilkinson shift (AEP p.512), */

	    d__1 = safmin * 100.;
	    dlag2_(&a_ref(ilast - 1, ilast - 1), lda, &b_ref(ilast - 1, ilast 
		    - 1), ldb, &d__1, &s1, &s2, &wr, &wr2, &wi);

/* Computing MAX   
   Computing MAX */
	    d__3 = 1., d__4 = abs(wr), d__3 = max(d__3,d__4), d__4 = abs(wi);
	    d__1 = s1, d__2 = safmin * max(d__3,d__4);
	    temp = max(d__1,d__2);
	    if (wi != 0.) {
		goto L200;
	    }
	}

/*        Fiddle with shift to avoid overflow */

	temp = min(ascale,1.) * (safmax * .5);
	if (s1 > temp) {
	    scale = temp / s1;
	} else {
	    scale = 1.;
	}

	temp = min(bscale,1.) * (safmax * .5);
	if (abs(wr) > temp) {
/* Computing MIN */
	    d__1 = scale, d__2 = temp / abs(wr);
	    scale = min(d__1,d__2);
	}
	s1 = scale * s1;
	wr = scale * wr;

/*        Now check for two consecutive small subdiagonals. */

	i__2 = ifirst + 1;
	for (j = ilast - 1; j >= i__2; --j) {
	    istart = j;
	    temp = (d__1 = s1 * a_ref(j, j - 1), abs(d__1));
	    temp2 = (d__1 = s1 * a_ref(j, j) - wr * b_ref(j, j), abs(d__1));
	    tempr = max(temp,temp2);
	    if (tempr < 1. && tempr != 0.) {
		temp /= tempr;
		temp2 /= tempr;
	    }
	    if ((d__1 = ascale * a_ref(j + 1, j) * temp, abs(d__1)) <= ascale 
		    * atol * temp2) {
		goto L130;
	    }
/* L120: */
	}

	istart = ifirst;
L130:

/*        Do an implicit single-shift QZ sweep.   

          Initial Q */

	temp = s1 * a_ref(istart, istart) - wr * b_ref(istart, istart);
	temp2 = s1 * a_ref(istart + 1, istart);
	dlartg_(&temp, &temp2, &c__, &s, &tempr);

/*        Sweep */

	i__2 = ilast - 1;
	for (j = istart; j <= i__2; ++j) {
	    if (j > istart) {
		temp = a_ref(j, j - 1);
		dlartg_(&temp, &a_ref(j + 1, j - 1), &c__, &s, &a_ref(j, j - 
			1));
		a_ref(j + 1, j - 1) = 0.;
	    }

	    i__3 = ilastm;
	    for (jc = j; jc <= i__3; ++jc) {
		temp = c__ * a_ref(j, jc) + s * a_ref(j + 1, jc);
		a_ref(j + 1, jc) = -s * a_ref(j, jc) + c__ * a_ref(j + 1, jc);
		a_ref(j, jc) = temp;
		temp2 = c__ * b_ref(j, jc) + s * b_ref(j + 1, jc);
		b_ref(j + 1, jc) = -s * b_ref(j, jc) + c__ * b_ref(j + 1, jc);
		b_ref(j, jc) = temp2;
/* L140: */
	    }
	    if (ilq) {
		i__3 = *n;
		for (jr = 1; jr <= i__3; ++jr) {
		    temp = c__ * q_ref(jr, j) + s * q_ref(jr, j + 1);
		    q_ref(jr, j + 1) = -s * q_ref(jr, j) + c__ * q_ref(jr, j 
			    + 1);
		    q_ref(jr, j) = temp;
/* L150: */
		}
	    }

	    temp = b_ref(j + 1, j + 1);
	    dlartg_(&temp, &b_ref(j + 1, j), &c__, &s, &b_ref(j + 1, j + 1));
	    b_ref(j + 1, j) = 0.;

/* Computing MIN */
	    i__4 = j + 2;
	    i__3 = min(i__4,ilast);
	    for (jr = ifrstm; jr <= i__3; ++jr) {
		temp = c__ * a_ref(jr, j + 1) + s * a_ref(jr, j);
		a_ref(jr, j) = -s * a_ref(jr, j + 1) + c__ * a_ref(jr, j);
		a_ref(jr, j + 1) = temp;
/* L160: */
	    }
	    i__3 = j;
	    for (jr = ifrstm; jr <= i__3; ++jr) {
		temp = c__ * b_ref(jr, j + 1) + s * b_ref(jr, j);
		b_ref(jr, j) = -s * b_ref(jr, j + 1) + c__ * b_ref(jr, j);
		b_ref(jr, j + 1) = temp;
/* L170: */
	    }
	    if (ilz) {
		i__3 = *n;
		for (jr = 1; jr <= i__3; ++jr) {
		    temp = c__ * z___ref(jr, j + 1) + s * z___ref(jr, j);
		    z___ref(jr, j) = -s * z___ref(jr, j + 1) + c__ * z___ref(
			    jr, j);
		    z___ref(jr, j + 1) = temp;
/* L180: */
		}
	    }
/* L190: */
	}

	goto L350;

/*        Use Francis double-shift   

          Note: the Francis double-shift should work with real shifts,   
                but only if the block is at least 3x3.   
                This code may break if this point is reached with   
                a 2x2 block with real eigenvalues. */

L200:
	if (ifirst + 1 == ilast) {

/*           Special case -- 2x2 block with complex eigenvectors   

             Step 1: Standardize, that is, rotate so that   

                         ( B11  0  )   
                     B = (         )  with B11 non-negative.   
                         (  0  B22 ) */

	    dlasv2_(&b_ref(ilast - 1, ilast - 1), &b_ref(ilast - 1, ilast), &
		    b_ref(ilast, ilast), &b22, &b11, &sr, &cr, &sl, &cl);

	    if (b11 < 0.) {
		cr = -cr;
		sr = -sr;
		b11 = -b11;
		b22 = -b22;
	    }

	    i__2 = ilastm + 1 - ifirst;
	    drot_(&i__2, &a_ref(ilast - 1, ilast - 1), lda, &a_ref(ilast, 
		    ilast - 1), lda, &cl, &sl);
	    i__2 = ilast + 1 - ifrstm;
	    drot_(&i__2, &a_ref(ifrstm, ilast - 1), &c__1, &a_ref(ifrstm, 
		    ilast), &c__1, &cr, &sr);

	    if (ilast < ilastm) {
		i__2 = ilastm - ilast;
		drot_(&i__2, &b_ref(ilast - 1, ilast + 1), ldb, &b_ref(ilast, 
			ilast + 1), lda, &cl, &sl);
	    }
	    if (ifrstm < ilast - 1) {
		i__2 = ifirst - ifrstm;
		drot_(&i__2, &b_ref(ifrstm, ilast - 1), &c__1, &b_ref(ifrstm, 
			ilast), &c__1, &cr, &sr);
	    }

	    if (ilq) {
		drot_(n, &q_ref(1, ilast - 1), &c__1, &q_ref(1, ilast), &c__1,
			 &cl, &sl);
	    }
	    if (ilz) {
		drot_(n, &z___ref(1, ilast - 1), &c__1, &z___ref(1, ilast), &
			c__1, &cr, &sr);
	    }

	    b_ref(ilast - 1, ilast - 1) = b11;
	    b_ref(ilast - 1, ilast) = 0.;
	    b_ref(ilast, ilast - 1) = 0.;
	    b_ref(ilast, ilast) = b22;

/*           If B22 is negative, negate column ILAST */

	    if (b22 < 0.) {
		i__2 = ilast;
		for (j = ifrstm; j <= i__2; ++j) {
		    a_ref(j, ilast) = -a_ref(j, ilast);
		    b_ref(j, ilast) = -b_ref(j, ilast);
/* L210: */
		}

		if (ilz) {
		    i__2 = *n;
		    for (j = 1; j <= i__2; ++j) {
			z___ref(j, ilast) = -z___ref(j, ilast);
/* L220: */
		    }
		}
	    }

/*           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)   

             Recompute shift */

	    d__1 = safmin * 100.;
	    dlag2_(&a_ref(ilast - 1, ilast - 1), lda, &b_ref(ilast - 1, ilast 
		    - 1), ldb, &d__1, &s1, &temp, &wr, &temp2, &wi);

/*           If standardization has perturbed the shift onto real line,   
             do another (real single-shift) QR step. */

	    if (wi == 0.) {
		goto L350;
	    }
	    s1inv = 1. / s1;

/*           Do EISPACK (QZVAL) computation of alpha and beta */

	    a11 = a_ref(ilast - 1, ilast - 1);
	    a21 = a_ref(ilast, ilast - 1);
	    a12 = a_ref(ilast - 1, ilast);
	    a22 = a_ref(ilast, ilast);

/*           Compute complex Givens rotation on right   
             (Assume some element of C = (sA - wB) > unfl )   
                              __   
             (sA - wB) ( CZ   -SZ )   
                       ( SZ    CZ ) */

	    c11r = s1 * a11 - wr * b11;
	    c11i = -wi * b11;
	    c12 = s1 * a12;
	    c21 = s1 * a21;
	    c22r = s1 * a22 - wr * b22;
	    c22i = -wi * b22;

	    if (abs(c11r) + abs(c11i) + abs(c12) > abs(c21) + abs(c22r) + abs(
		    c22i)) {
		t = dlapy3_(&c12, &c11r, &c11i);
		cz = c12 / t;
		szr = -c11r / t;
		szi = -c11i / t;
	    } else {
		cz = dlapy2_(&c22r, &c22i);
		if (cz <= safmin) {
		    cz = 0.;
		    szr = 1.;
		    szi = 0.;
		} else {
		    tempr = c22r / cz;
		    tempi = c22i / cz;
		    t = dlapy2_(&cz, &c21);
		    cz /= t;
		    szr = -c21 * tempr / t;
		    szi = c21 * tempi / t;
		}
	    }

/*           Compute Givens rotation on left   

             (  CQ   SQ )   
             (  __      )  A or B   
             ( -SQ   CQ ) */

	    an = abs(a11) + abs(a12) + abs(a21) + abs(a22);
	    bn = abs(b11) + abs(b22);
	    wabs = abs(wr) + abs(wi);
	    if (s1 * an > wabs * bn) {
		cq = cz * b11;
		sqr = szr * b22;
		sqi = -szi * b22;
	    } else {
		a1r = cz * a11 + szr * a12;
		a1i = szi * a12;
		a2r = cz * a21 + szr * a22;
		a2i = szi * a22;
		cq = dlapy2_(&a1r, &a1i);
		if (cq <= safmin) {
		    cq = 0.;
		    sqr = 1.;
		    sqi = 0.;
		} else {
		    tempr = a1r / cq;
		    tempi = a1i / cq;
		    sqr = tempr * a2r + tempi * a2i;
		    sqi = tempi * a2r - tempr * a2i;
		}
	    }
	    t = dlapy3_(&cq, &sqr, &sqi);
	    cq /= t;
	    sqr /= t;
	    sqi /= t;

/*           Compute diagonal elements of QBZ */

	    tempr = sqr * szr - sqi * szi;
	    tempi = sqr * szi + sqi * szr;
	    b1r = cq * cz * b11 + tempr * b22;
	    b1i = tempi * b22;
	    b1a = dlapy2_(&b1r, &b1i);
	    b2r = cq * cz * b22 + tempr * b11;
	    b2i = -tempi * b11;
	    b2a = dlapy2_(&b2r, &b2i);

/*           Normalize so beta > 0, and Im( alpha1 ) > 0 */

	    beta[ilast - 1] = b1a;
	    beta[ilast] = b2a;
	    alphar[ilast - 1] = wr * b1a * s1inv;
	    alphai[ilast - 1] = wi * b1a * s1inv;
	    alphar[ilast] = wr * b2a * s1inv;
	    alphai[ilast] = -(wi * b2a) * s1inv;

/*           Step 3: Go to next block -- exit if finished. */

	    ilast = ifirst - 1;
	    if (ilast < *ilo) {
		goto L380;
	    }

/*           Reset counters */

	    iiter = 0;
	    eshift = 0.;
	    if (! ilschr) {
		ilastm = ilast;
		if (ifrstm > ilast) {
		    ifrstm = *ilo;
		}
	    }
	    goto L350;
	} else {

/*           Usual case: 3x3 or larger block, using Francis implicit   
                         double-shift   

                                      2   
             Eigenvalue equation is  w  - c w + d = 0,   

                                           -1 2        -1   
             so compute 1st column of  (A B  )  - c A B   + d   
             using the formula in QZIT (from EISPACK)   

             We assume that the block is at least 3x3 */

	    ad11 = ascale * a_ref(ilast - 1, ilast - 1) / (bscale * b_ref(
		    ilast - 1, ilast - 1));
	    ad21 = ascale * a_ref(ilast, ilast - 1) / (bscale * b_ref(ilast - 
		    1, ilast - 1));
	    ad12 = ascale * a_ref(ilast - 1, ilast) / (bscale * b_ref(ilast, 
		    ilast));
	    ad22 = ascale * a_ref(ilast, ilast) / (bscale * b_ref(ilast, 
		    ilast));
	    u12 = b_ref(ilast - 1, ilast) / b_ref(ilast, ilast);
	    ad11l = ascale * a_ref(ifirst, ifirst) / (bscale * b_ref(ifirst, 
		    ifirst));
	    ad21l = ascale * a_ref(ifirst + 1, ifirst) / (bscale * b_ref(
		    ifirst, ifirst));
	    ad12l = ascale * a_ref(ifirst, ifirst + 1) / (bscale * b_ref(
		    ifirst + 1, ifirst + 1));
	    ad22l = ascale * a_ref(ifirst + 1, ifirst + 1) / (bscale * b_ref(
		    ifirst + 1, ifirst + 1));
	    ad32l = ascale * a_ref(ifirst + 2, ifirst + 1) / (bscale * b_ref(
		    ifirst + 1, ifirst + 1));
	    u12l = b_ref(ifirst, ifirst + 1) / b_ref(ifirst + 1, ifirst + 1);

	    v[0] = (ad11 - ad11l) * (ad22 - ad11l) - ad12 * ad21 + ad21 * u12 
		    * ad11l + (ad12l - ad11l * u12l) * ad21l;
	    v[1] = (ad22l - ad11l - ad21l * u12l - (ad11 - ad11l) - (ad22 - 
		    ad11l) + ad21 * u12) * ad21l;
	    v[2] = ad32l * ad21l;

	    istart = ifirst;

	    dlarfg_(&c__3, v, &v[1], &c__1, &tau);
	    v[0] = 1.;

/*           Sweep */

	    i__2 = ilast - 2;
	    for (j = istart; j <= i__2; ++j) {

/*              All but last elements: use 3x3 Householder transforms.   

                Zero (j-1)st column of A */

		if (j > istart) {
		    v[0] = a_ref(j, j - 1);
		    v[1] = a_ref(j + 1, j - 1);
		    v[2] = a_ref(j + 2, j - 1);

		    dlarfg_(&c__3, &a_ref(j, j - 1), &v[1], &c__1, &tau);
		    v[0] = 1.;
		    a_ref(j + 1, j - 1) = 0.;
		    a_ref(j + 2, j - 1) = 0.;
		}

		i__3 = ilastm;
		for (jc = j; jc <= i__3; ++jc) {
		    temp = tau * (a_ref(j, jc) + v[1] * a_ref(j + 1, jc) + v[
			    2] * a_ref(j + 2, jc));
		    a_ref(j, jc) = a_ref(j, jc) - temp;
		    a_ref(j + 1, jc) = a_ref(j + 1, jc) - temp * v[1];
		    a_ref(j + 2, jc) = a_ref(j + 2, jc) - temp * v[2];
		    temp2 = tau * (b_ref(j, jc) + v[1] * b_ref(j + 1, jc) + v[
			    2] * b_ref(j + 2, jc));
		    b_ref(j, jc) = b_ref(j, jc) - temp2;
		    b_ref(j + 1, jc) = b_ref(j + 1, jc) - temp2 * v[1];
		    b_ref(j + 2, jc) = b_ref(j + 2, jc) - temp2 * v[2];
/* L230: */
		}
		if (ilq) {
		    i__3 = *n;
		    for (jr = 1; jr <= i__3; ++jr) {
			temp = tau * (q_ref(jr, j) + v[1] * q_ref(jr, j + 1) 
				+ v[2] * q_ref(jr, j + 2));
			q_ref(jr, j) = q_ref(jr, j) - temp;
			q_ref(jr, j + 1) = q_ref(jr, j + 1) - temp * v[1];
			q_ref(jr, j + 2) = q_ref(jr, j + 2) - temp * v[2];
/* L240: */
		    }
		}

/*              Zero j-th column of B (see DLAGBC for details)   

                Swap rows to pivot */

		ilpivt = FALSE_;
/* Computing MAX */
		d__3 = (d__1 = b_ref(j + 1, j + 1), abs(d__1)), d__4 = (d__2 =
			 b_ref(j + 1, j + 2), abs(d__2));
		temp = max(d__3,d__4);
/* Computing MAX */
		d__3 = (d__1 = b_ref(j + 2, j + 1), abs(d__1)), d__4 = (d__2 =
			 b_ref(j + 2, j + 2), abs(d__2));
		temp2 = max(d__3,d__4);
		if (max(temp,temp2) < safmin) {
		    scale = 0.;
		    u1 = 1.;
		    u2 = 0.;
		    goto L250;
		} else if (temp >= temp2) {
		    w11 = b_ref(j + 1, j + 1);
		    w21 = b_ref(j + 2, j + 1);
		    w12 = b_ref(j + 1, j + 2);
		    w22 = b_ref(j + 2, j + 2);
		    u1 = b_ref(j + 1, j);
		    u2 = b_ref(j + 2, j);
		} else {
		    w21 = b_ref(j + 1, j + 1);
		    w11 = b_ref(j + 2, j + 1);
		    w22 = b_ref(j + 1, j + 2);
		    w12 = b_ref(j + 2, j + 2);
		    u2 = b_ref(j + 1, j);
		    u1 = b_ref(j + 2, j);
		}

/*              Swap columns if nec. */

		if (abs(w12) > abs(w11)) {
		    ilpivt = TRUE_;
		    temp = w12;
		    temp2 = w22;
		    w12 = w11;
		    w22 = w21;
		    w11 = temp;
		    w21 = temp2;
		}

/*              LU-factor */

		temp = w21 / w11;
		u2 -= temp * u1;
		w22 -= temp * w12;
		w21 = 0.;

/*              Compute SCALE */

		scale = 1.;
		if (abs(w22) < safmin) {
		    scale = 0.;
		    u2 = 1.;
		    u1 = -w12 / w11;
		    goto L250;
		}
		if (abs(w22) < abs(u2)) {
		    scale = (d__1 = w22 / u2, abs(d__1));
		}
		if (abs(w11) < abs(u1)) {
/* Computing MIN */
		    d__2 = scale, d__3 = (d__1 = w11 / u1, abs(d__1));
		    scale = min(d__2,d__3);
		}

/*              Solve */

		u2 = scale * u2 / w22;
		u1 = (scale * u1 - w12 * u2) / w11;

L250:
		if (ilpivt) {
		    temp = u2;
		    u2 = u1;
		    u1 = temp;
		}

/*              Compute Householder Vector   

   Computing 2nd power */
		d__1 = scale;
/* Computing 2nd power */
		d__2 = u1;
/* Computing 2nd power */
		d__3 = u2;
		t = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		tau = scale / t + 1.;
		vs = -1. / (scale + t);
		v[0] = 1.;
		v[1] = vs * u1;
		v[2] = vs * u2;

/*              Apply transformations from the right.   

   Computing MIN */
		i__4 = j + 3;
		i__3 = min(i__4,ilast);
		for (jr = ifrstm; jr <= i__3; ++jr) {
		    temp = tau * (a_ref(jr, j) + v[1] * a_ref(jr, j + 1) + v[
			    2] * a_ref(jr, j + 2));
		    a_ref(jr, j) = a_ref(jr, j) - temp;
		    a_ref(jr, j + 1) = a_ref(jr, j + 1) - temp * v[1];
		    a_ref(jr, j + 2) = a_ref(jr, j + 2) - temp * v[2];
/* L260: */
		}
		i__3 = j + 2;
		for (jr = ifrstm; jr <= i__3; ++jr) {
		    temp = tau * (b_ref(jr, j) + v[1] * b_ref(jr, j + 1) + v[
			    2] * b_ref(jr, j + 2));
		    b_ref(jr, j) = b_ref(jr, j) - temp;
		    b_ref(jr, j + 1) = b_ref(jr, j + 1) - temp * v[1];
		    b_ref(jr, j + 2) = b_ref(jr, j + 2) - temp * v[2];
/* L270: */
		}
		if (ilz) {
		    i__3 = *n;
		    for (jr = 1; jr <= i__3; ++jr) {
			temp = tau * (z___ref(jr, j) + v[1] * z___ref(jr, j + 
				1) + v[2] * z___ref(jr, j + 2));
			z___ref(jr, j) = z___ref(jr, j) - temp;
			z___ref(jr, j + 1) = z___ref(jr, j + 1) - temp * v[1];
			z___ref(jr, j + 2) = z___ref(jr, j + 2) - temp * v[2];
/* L280: */
		    }
		}
		b_ref(j + 1, j) = 0.;
		b_ref(j + 2, j) = 0.;
/* L290: */
	    }

/*           Last elements: Use Givens rotations   

             Rotations from the left */

	    j = ilast - 1;
	    temp = a_ref(j, j - 1);
	    dlartg_(&temp, &a_ref(j + 1, j - 1), &c__, &s, &a_ref(j, j - 1));
	    a_ref(j + 1, j - 1) = 0.;

	    i__2 = ilastm;
	    for (jc = j; jc <= i__2; ++jc) {
		temp = c__ * a_ref(j, jc) + s * a_ref(j + 1, jc);
		a_ref(j + 1, jc) = -s * a_ref(j, jc) + c__ * a_ref(j + 1, jc);
		a_ref(j, jc) = temp;
		temp2 = c__ * b_ref(j, jc) + s * b_ref(j + 1, jc);
		b_ref(j + 1, jc) = -s * b_ref(j, jc) + c__ * b_ref(j + 1, jc);
		b_ref(j, jc) = temp2;
/* L300: */
	    }
	    if (ilq) {
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    temp = c__ * q_ref(jr, j) + s * q_ref(jr, j + 1);
		    q_ref(jr, j + 1) = -s * q_ref(jr, j) + c__ * q_ref(jr, j 
			    + 1);
		    q_ref(jr, j) = temp;
/* L310: */
		}
	    }

/*           Rotations from the right. */

	    temp = b_ref(j + 1, j + 1);
	    dlartg_(&temp, &b_ref(j + 1, j), &c__, &s, &b_ref(j + 1, j + 1));
	    b_ref(j + 1, j) = 0.;

	    i__2 = ilast;
	    for (jr = ifrstm; jr <= i__2; ++jr) {
		temp = c__ * a_ref(jr, j + 1) + s * a_ref(jr, j);
		a_ref(jr, j) = -s * a_ref(jr, j + 1) + c__ * a_ref(jr, j);
		a_ref(jr, j + 1) = temp;
/* L320: */
	    }
	    i__2 = ilast - 1;
	    for (jr = ifrstm; jr <= i__2; ++jr) {
		temp = c__ * b_ref(jr, j + 1) + s * b_ref(jr, j);
		b_ref(jr, j) = -s * b_ref(jr, j + 1) + c__ * b_ref(jr, j);
		b_ref(jr, j + 1) = temp;
/* L330: */
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    temp = c__ * z___ref(jr, j + 1) + s * z___ref(jr, j);
		    z___ref(jr, j) = -s * z___ref(jr, j + 1) + c__ * z___ref(
			    jr, j);
		    z___ref(jr, j + 1) = temp;
/* L340: */
		}
	    }

/*           End of Double-Shift code */

	}

	goto L350;

/*        End of iteration loop */

L350:
/* L360: */
	;
    }

/*     Drop-through = non-convergence   

   L370: */
    *info = ilast;
    goto L420;

/*     Successful completion of all QZ steps */

L380:

/*     Set Eigenvalues 1:ILO-1 */

    i__1 = *ilo - 1;
    for (j = 1; j <= i__1; ++j) {
	if (b_ref(j, j) < 0.) {
	    if (ilschr) {
		i__2 = j;
		for (jr = 1; jr <= i__2; ++jr) {
		    a_ref(jr, j) = -a_ref(jr, j);
		    b_ref(jr, j) = -b_ref(jr, j);
/* L390: */
		}
	    } else {
		a_ref(j, j) = -a_ref(j, j);
		b_ref(j, j) = -b_ref(j, j);
	    }
	    if (ilz) {
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    z___ref(jr, j) = -z___ref(jr, j);
/* L400: */
		}
	    }
	}
	alphar[j] = a_ref(j, j);
	alphai[j] = 0.;
	beta[j] = b_ref(j, j);
/* L410: */
    }

/*     Normal Termination */

    *info = 0;

/*     Exit (other than argument error) -- return optimal workspace size */

L420:
    work[1] = (doublereal) (*n);
    return 0;

/*     End of DHGEQZ */

} /* dhgeqz_ */

#undef z___ref
#undef q_ref
#undef b_ref
#undef a_ref



/* Subroutine */ int dtgevc_(char *side, char *howmny, logical *select, 
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer 
	*mm, integer *m, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   



    Purpose   
    =======   

    DTGEVC computes some or all of the right and/or left generalized   
    eigenvectors of a pair of real upper triangular matrices (A,B).   

    The right generalized eigenvector x and the left generalized   
    eigenvector y of (A,B) corresponding to a generalized eigenvalue   
    w are defined by:   

            (A - wB) * x = 0  and  y**H * (A - wB) = 0   

    where y**H denotes the conjugate tranpose of y.   

    If an eigenvalue w is determined by zero diagonal elements of both A   
    and B, a unit vector is returned as the corresponding eigenvector.   

    If all eigenvectors are requested, the routine may either return   
    the matrices X and/or Y of right or left eigenvectors of (A,B), or   
    the products Z*X and/or Q*Y, where Z and Q are input orthogonal   
    matrices.  If (A,B) was obtained from the generalized real-Schur   
    factorization of an original pair of matrices   
       (A0,B0) = (Q*A*Z**H,Q*B*Z**H),   
    then Z*X and Q*Y are the matrices of right or left eigenvectors of   
    A.   

    A must be block upper triangular, with 1-by-1 and 2-by-2 diagonal   
    blocks.  Corresponding to each 2-by-2 diagonal block is a complex   
    conjugate pair of eigenvalues and eigenvectors; only one   
    eigenvector of the pair is computed, namely the one corresponding   
    to the eigenvalue with positive imaginary part.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'R': compute right eigenvectors only;   
            = 'L': compute left eigenvectors only;   
            = 'B': compute both right and left eigenvectors.   

    HOWMNY  (input) CHARACTER*1   
            = 'A': compute all right and/or left eigenvectors;   
            = 'B': compute all right and/or left eigenvectors, and   
                   backtransform them using the input matrices supplied   
                   in VR and/or VL;   
            = 'S': compute selected right and/or left eigenvectors,   
                   specified by the logical array SELECT.   

    SELECT  (input) LOGICAL array, dimension (N)   
            If HOWMNY='S', SELECT specifies the eigenvectors to be   
            computed.   
            If HOWMNY='A' or 'B', SELECT is not referenced.   
            To select the real eigenvector corresponding to the real   
            eigenvalue w(j), SELECT(j) must be set to .TRUE.  To select   
            the complex eigenvector corresponding to a complex conjugate   
            pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must   
            be set to .TRUE..   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The upper quasi-triangular matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of array A.  LDA >= max(1, N).   

    B       (input) DOUBLE PRECISION array, dimension (LDB,N)   
            The upper triangular matrix B.  If A has a 2-by-2 diagonal   
            block, then the corresponding 2-by-2 block of B must be   
            diagonal with positive elements.   

    LDB     (input) INTEGER   
            The leading dimension of array B.  LDB >= max(1,N).   

    VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)   
            On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must   
            contain an N-by-N matrix Q (usually the orthogonal matrix Q   
            of left Schur vectors returned by DHGEQZ).   
            On exit, if SIDE = 'L' or 'B', VL contains:   
            if HOWMNY = 'A', the matrix Y of left eigenvectors of (A,B);   
            if HOWMNY = 'B', the matrix Q*Y;   
            if HOWMNY = 'S', the left eigenvectors of (A,B) specified by   
                        SELECT, stored consecutively in the columns of   
                        VL, in the same order as their eigenvalues.   
            If SIDE = 'R', VL is not referenced.   

            A complex eigenvector corresponding to a complex eigenvalue   
            is stored in two consecutive columns, the first holding the   
            real part, and the second the imaginary part.   

    LDVL    (input) INTEGER   
            The leading dimension of array VL.   
            LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.   

    VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)   
            On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must   
            contain an N-by-N matrix Q (usually the orthogonal matrix Z   
            of right Schur vectors returned by DHGEQZ).   
            On exit, if SIDE = 'R' or 'B', VR contains:   
            if HOWMNY = 'A', the matrix X of right eigenvectors of (A,B);   
            if HOWMNY = 'B', the matrix Z*X;   
            if HOWMNY = 'S', the right eigenvectors of (A,B) specified by   
                        SELECT, stored consecutively in the columns of   
                        VR, in the same order as their eigenvalues.   
            If SIDE = 'L', VR is not referenced.   

            A complex eigenvector corresponding to a complex eigenvalue   
            is stored in two consecutive columns, the first holding the   
            real part and the second the imaginary part.   

    LDVR    (input) INTEGER   
            The leading dimension of the array VR.   
            LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.   

    MM      (input) INTEGER   
            The number of columns in the arrays VL and/or VR. MM >= M.   

    M       (output) INTEGER   
            The number of columns in the arrays VL and/or VR actually   
            used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M   
            is set to N.  Each selected real eigenvector occupies one   
            column and each selected complex eigenvector occupies two   
            columns.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex   
                  eigenvalue.   

    Further Details   
    ===============   

    Allocation of workspace:   
    ---------- -- ---------   

       WORK( j ) = 1-norm of j-th column of A, above the diagonal   
       WORK( N+j ) = 1-norm of j-th column of B, above the diagonal   
       WORK( 2*N+1:3*N ) = real part of eigenvector   
       WORK( 3*N+1:4*N ) = imaginary part of eigenvector   
       WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector   
       WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector   

    Rowwise vs. columnwise solution methods:   
    ------- --  ---------- -------- -------   

    Finding a generalized eigenvector consists basically of solving the   
    singular triangular system   

     (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)   

    Consider finding the i-th right eigenvector (assume all eigenvalues   
    are real). The equation to be solved is:   
         n                   i   
    0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1   
        k=j                 k=j   

    where  C = (A - w B)  (The components v(i+1:n) are 0.)   

    The "rowwise" method is:   

    (1)  v(i) := 1   
    for j = i-1,. . .,1:   
                            i   
        (2) compute  s = - sum C(j,k) v(k)   and   
                          k=j+1   

        (3) v(j) := s / C(j,j)   

    Step 2 is sometimes called the "dot product" step, since it is an   
    inner product between the j-th row and the portion of the eigenvector   
    that has been computed so far.   

    The "columnwise" method consists basically in doing the sums   
    for all the rows in parallel.  As each v(j) is computed, the   
    contribution of v(j) times the j-th column of C is added to the   
    partial sums.  Since FORTRAN arrays are stored columnwise, this has   
    the advantage that at each step, the elements of C that are accessed   
    are adjacent to one another, whereas with the rowwise method, the   
    elements accessed at a step are spaced LDA (and LDB) words apart.   

    When finding left eigenvectors, the matrix in question is the   
    transpose of the one in storage, so the rowwise method then   
    actually accesses columns of A and B at each step, and so is the   
    preferred method.   

    =====================================================================   


       Decode and Test the input parameters   

       Parameter adjustments */
    /* Table of constant values */
    static logical c_true = TRUE_;
    static integer c__2 = 2;
    static doublereal c_b35 = 1.;
    static integer c__1 = 1;
    static doublereal c_b37 = 0.;
    static logical c_false = FALSE_;
    
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, vl_dim1, vl_offset, vr_dim1, 
	    vr_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    /* Local variables */
    static integer ibeg, ieig, iend;
    static doublereal dmin__, temp, suma[4]	/* was [2][2] */, sumb[4]	
	    /* was [2][2] */, xmax;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal cim2a, cim2b, cre2a, cre2b, temp2, bdiag[2];
    static integer i__, j;
    static doublereal acoef, scale;
    static logical ilall;
    static integer iside;
    static doublereal sbeta;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);
    static logical il2by2;
    static integer iinfo;
    static doublereal small;
    static logical compl;
    static doublereal anorm, bnorm;
    static logical compr;
    extern /* Subroutine */ int dlaln2_(logical *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, doublereal *
	    , doublereal *, integer *, doublereal *, doublereal *, integer *);
    static doublereal temp2i;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    static doublereal temp2r;
    static integer ja;
    static logical ilabad, ilbbad;
    static integer jc, je, na;
    static doublereal acoefa, bcoefa, cimaga, cimagb;
    static logical ilback;
    static integer im;
    static doublereal bcoefi, ascale, bscale, creala;
    static integer jr;
    static doublereal crealb;
    extern doublereal dlamch_(char *);
    static doublereal bcoefr;
    static integer jw, nw;
    static doublereal salfar, safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal xscale, bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical ilcomp, ilcplx;
    static integer ihwmny;
    static doublereal big;
    static logical lsa, lsb;
    static doublereal ulp, sum[4]	/* was [2][2] */;
#define suma_ref(a_1,a_2) suma[(a_2)*2 + a_1 - 3]
#define sumb_ref(a_1,a_2) sumb[(a_2)*2 + a_1 - 3]
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define vl_ref(a_1,a_2) vl[(a_2)*vl_dim1 + a_1]
#define vr_ref(a_1,a_2) vr[(a_2)*vr_dim1 + a_1]
#define sum_ref(a_1,a_2) sum[(a_2)*2 + a_1 - 3]


    --select;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    vl_dim1 = *ldvl;
    vl_offset = 1 + vl_dim1 * 1;
    vl -= vl_offset;
    vr_dim1 = *ldvr;
    vr_offset = 1 + vr_dim1 * 1;
    vr -= vr_offset;
    --work;

    /* Function Body */
    if (lsame_(howmny, "A")) {
	ihwmny = 1;
	ilall = TRUE_;
	ilback = FALSE_;
    } else if (lsame_(howmny, "S")) {
	ihwmny = 2;
	ilall = FALSE_;
	ilback = FALSE_;
    } else if (lsame_(howmny, "B") || lsame_(howmny, 
	    "T")) {
	ihwmny = 3;
	ilall = TRUE_;
	ilback = TRUE_;
    } else {
	ihwmny = -1;
	ilall = TRUE_;
    }

    if (lsame_(side, "R")) {
	iside = 1;
	compl = FALSE_;
	compr = TRUE_;
    } else if (lsame_(side, "L")) {
	iside = 2;
	compl = TRUE_;
	compr = FALSE_;
    } else if (lsame_(side, "B")) {
	iside = 3;
	compl = TRUE_;
	compr = TRUE_;
    } else {
	iside = -1;
    }

    *info = 0;
    if (iside < 0) {
	*info = -1;
    } else if (ihwmny < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -4;
    } else if (*lda < max(1,*n)) {
	*info = -6;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTGEVC", &i__1);
	return 0;
    }

/*     Count the number of eigenvectors to be computed */

    if (! ilall) {
	im = 0;
	ilcplx = FALSE_;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L10;
	    }
	    if (j < *n) {
		if (a_ref(j + 1, j) != 0.) {
		    ilcplx = TRUE_;
		}
	    }
	    if (ilcplx) {
		if (select[j] || select[j + 1]) {
		    im += 2;
		}
	    } else {
		if (select[j]) {
		    ++im;
		}
	    }
L10:
	    ;
	}
    } else {
	im = *n;
    }

/*     Check 2-by-2 diagonal blocks of A, B */

    ilabad = FALSE_;
    ilbbad = FALSE_;
    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j) {
	if (a_ref(j + 1, j) != 0.) {
	    if (b_ref(j, j) == 0. || b_ref(j + 1, j + 1) == 0. || b_ref(j, j 
		    + 1) != 0.) {
		ilbbad = TRUE_;
	    }
	    if (j < *n - 1) {
		if (a_ref(j + 2, j + 1) != 0.) {
		    ilabad = TRUE_;
		}
	    }
	}
/* L20: */
    }

    if (ilabad) {
	*info = -5;
    } else if (ilbbad) {
	*info = -7;
    } else if (compl && *ldvl < *n || *ldvl < 1) {
	*info = -10;
    } else if (compr && *ldvr < *n || *ldvr < 1) {
	*info = -12;
    } else if (*mm < im) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTGEVC", &i__1);
	return 0;
    }

/*     Quick return if possible */

    *m = im;
    if (*n == 0) {
	return 0;
    }

/*     Machine Constants */

    safmin = dlamch_("Safe minimum");
    big = 1. / safmin;
    dlabad_(&safmin, &big);
    ulp = dlamch_("Epsilon") * dlamch_("Base");
    small = safmin * *n / ulp;
    big = 1. / small;
    bignum = 1. / (safmin * *n);

/*     Compute the 1-norm of each column of the strictly upper triangular   
       part (i.e., excluding all elements belonging to the diagonal   
       blocks) of A and B to check for possible overflow in the   
       triangular solver. */

    anorm = (d__1 = a_ref(1, 1), abs(d__1));
    if (*n > 1) {
	anorm += (d__1 = a_ref(2, 1), abs(d__1));
    }
    bnorm = (d__1 = b_ref(1, 1), abs(d__1));
    work[1] = 0.;
    work[*n + 1] = 0.;

    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = 0.;
	temp2 = 0.;
	if (a_ref(j, j - 1) == 0.) {
	    iend = j - 1;
	} else {
	    iend = j - 2;
	}
	i__2 = iend;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    temp += (d__1 = a_ref(i__, j), abs(d__1));
	    temp2 += (d__1 = b_ref(i__, j), abs(d__1));
/* L30: */
	}
	work[j] = temp;
	work[*n + j] = temp2;
/* Computing MIN */
	i__3 = j + 1;
	i__2 = min(i__3,*n);
	for (i__ = iend + 1; i__ <= i__2; ++i__) {
	    temp += (d__1 = a_ref(i__, j), abs(d__1));
	    temp2 += (d__1 = b_ref(i__, j), abs(d__1));
/* L40: */
	}
	anorm = max(anorm,temp);
	bnorm = max(bnorm,temp2);
/* L50: */
    }

    ascale = 1. / max(anorm,safmin);
    bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

    if (compl) {
	ieig = 0;

/*        Main loop over eigenvalues */

	ilcplx = FALSE_;
	i__1 = *n;
	for (je = 1; je <= i__1; ++je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or   
             (b) this would be the second of a complex pair.   
             Check for complex eigenvalue, so as to be sure of which   
             entry(-ies) of SELECT to look at. */

	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L220;
	    }
	    nw = 1;
	    if (je < *n) {
		if (a_ref(je + 1, je) != 0.) {
		    ilcplx = TRUE_;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = TRUE_;
	    } else if (ilcplx) {
		ilcomp = select[je] || select[je + 1];
	    } else {
		ilcomp = select[je];
	    }
	    if (! ilcomp) {
		goto L220;
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or   
             (c) complex eigenvalue. */

	    if (! ilcplx) {
		if ((d__1 = a_ref(je, je), abs(d__1)) <= safmin && (d__2 = 
			b_ref(je, je), abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

		    ++ieig;
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vl_ref(jr, ieig) = 0.;
/* L60: */
		    }
		    vl_ref(ieig, ieig) = 1.;
		    goto L220;
		}
	    }

/*           Clear vector */

	    i__2 = nw * *n;
	    for (jr = 1; jr <= i__2; ++jr) {
		work[(*n << 1) + jr] = 0.;
/* L70: */
	    }
/*                                                 T   
             Compute coefficients in  ( a A - b B )  y = 0   
                a  is  ACOEF   
                b  is  BCOEFR + i*BCOEFI */

	    if (! ilcplx) {

/*              Real eigenvalue   

   Computing MAX */
		d__3 = (d__1 = a_ref(je, je), abs(d__1)) * ascale, d__4 = (
			d__2 = b_ref(je, je), abs(d__2)) * bscale, d__3 = max(
			d__3,d__4);
		temp = 1. / max(d__3,safmin);
		salfar = temp * a_ref(je, je) * ascale;
		sbeta = temp * b_ref(je, je) * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = 0.;

/*              Scale to avoid underflow */

		scale = 1.;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
		    scale = max(d__1,d__2);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
		    scale = min(d__1,d__2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);

/*              First component is 1 */

		work[(*n << 1) + je] = 1.;
		xmax = 1.;
	    } else {

/*              Complex eigenvalue */

		d__1 = safmin * 100.;
		dlag2_(&a_ref(je, je), lda, &b_ref(je, je), ldb, &d__1, &
			acoef, &temp, &bcoefr, &temp2, &bcoefi);
		bcoefi = -bcoefi;
		if (bcoefi == 0.) {
		    *info = je;
		    return 0;
		}

/*              Scale to avoid over/underflow */

		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = 1.;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
		    scale = max(d__1,d__2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
		    scale = min(d__1,d__2);
		}
		if (scale != 1.) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}

/*              Compute first two components of eigenvector */

		temp = acoef * a_ref(je + 1, je);
		temp2r = acoef * a_ref(je, je) - bcoefr * b_ref(je, je);
		temp2i = -bcoefi * b_ref(je, je);
		if (abs(temp) > abs(temp2r) + abs(temp2i)) {
		    work[(*n << 1) + je] = 1.;
		    work[*n * 3 + je] = 0.;
		    work[(*n << 1) + je + 1] = -temp2r / temp;
		    work[*n * 3 + je + 1] = -temp2i / temp;
		} else {
		    work[(*n << 1) + je + 1] = 1.;
		    work[*n * 3 + je + 1] = 0.;
		    temp = acoef * a_ref(je, je + 1);
		    work[(*n << 1) + je] = (bcoefr * b_ref(je + 1, je + 1) - 
			    acoef * a_ref(je + 1, je + 1)) / temp;
		    work[*n * 3 + je] = bcoefi * b_ref(je + 1, je + 1) / temp;
		}
/* Computing MAX */
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je + 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je + 1], abs(d__4));
		xmax = max(d__5,d__6);
	    }

/* Computing MAX */
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
	    dmin__ = max(d__1,safmin);

/*                                           T   
             Triangular solve of  (a A - b B)  y = 0   

                                     T   
             (rowwise in  (a A - b B) , or columnwise in (a A - b B) ) */

	    il2by2 = FALSE_;

	    i__2 = *n;
	    for (j = je + nw; j <= i__2; ++j) {
		if (il2by2) {
		    il2by2 = FALSE_;
		    goto L160;
		}

		na = 1;
		bdiag[0] = b_ref(j, j);
		if (j < *n) {
		    if (a_ref(j + 1, j) != 0.) {
			il2by2 = TRUE_;
			bdiag[1] = b_ref(j + 1, j + 1);
			na = 2;
		    }
		}

/*              Check whether scaling is necessary for dot products */

		xscale = 1. / max(1.,xmax);
/* Computing MAX */
		d__1 = work[j], d__2 = work[*n + j], d__1 = max(d__1,d__2), 
			d__2 = acoefa * work[j] + bcoefa * work[*n + j];
		temp = max(d__1,d__2);
		if (il2by2) {
/* Computing MAX */
		    d__1 = temp, d__2 = work[j + 1], d__1 = max(d__1,d__2), 
			    d__2 = work[*n + j + 1], d__1 = max(d__1,d__2), 
			    d__2 = acoefa * work[j + 1] + bcoefa * work[*n + 
			    j + 1];
		    temp = max(d__1,d__2);
		}
		if (temp > bignum * xscale) {
		    i__3 = nw - 1;
		    for (jw = 0; jw <= i__3; ++jw) {
			i__4 = j - 1;
			for (jr = je; jr <= i__4; ++jr) {
			    work[(jw + 2) * *n + jr] = xscale * work[(jw + 2) 
				    * *n + jr];
/* L80: */
			}
/* L90: */
		    }
		    xmax *= xscale;
		}

/*              Compute dot products   

                      j-1   
                SUM = sum  conjg( a*A(k,j) - b*B(k,j) )*x(k)   
                      k=je   

                To reduce the op count, this is done as   

                _        j-1                  _        j-1   
                a*conjg( sum  A(k,j)*x(k) ) - b*conjg( sum  B(k,j)*x(k) )   
                         k=je                          k=je   

                which may cause underflow problems if A or B are close   
                to underflow.  (E.g., less than SMALL.)   


                A series of compiler directives to defeat vectorization   
                for the next loop   

   $PL$ CMCHAR=' '   
   DIR$          NEXTSCALAR   
   $DIR          SCALAR   
   DIR$          NEXT SCALAR   
   VD$L          NOVECTOR   
   DEC$          NOVECTOR   
   VD$           NOVECTOR   
   VDIR          NOVECTOR   
   VOCL          LOOP,SCALAR   
   IBM           PREFER SCALAR   
   $PL$ CMCHAR='*' */

		i__3 = nw;
		for (jw = 1; jw <= i__3; ++jw) {

/* $PL$ CMCHAR=' '   
   DIR$             NEXTSCALAR   
   $DIR             SCALAR   
   DIR$             NEXT SCALAR   
   VD$L             NOVECTOR   
   DEC$             NOVECTOR   
   VD$              NOVECTOR   
   VDIR             NOVECTOR   
   VOCL             LOOP,SCALAR   
   IBM              PREFER SCALAR   
   $PL$ CMCHAR='*' */

		    i__4 = na;
		    for (ja = 1; ja <= i__4; ++ja) {
			suma_ref(ja, jw) = 0.;
			sumb_ref(ja, jw) = 0.;

			i__5 = j - 1;
			for (jr = je; jr <= i__5; ++jr) {
			    suma_ref(ja, jw) = suma_ref(ja, jw) + a_ref(jr, j 
				    + ja - 1) * work[(jw + 1) * *n + jr];
			    sumb_ref(ja, jw) = sumb_ref(ja, jw) + b_ref(jr, j 
				    + ja - 1) * work[(jw + 1) * *n + jr];
/* L100: */
			}
/* L110: */
		    }
/* L120: */
		}

/* $PL$ CMCHAR=' '   
   DIR$          NEXTSCALAR   
   $DIR          SCALAR   
   DIR$          NEXT SCALAR   
   VD$L          NOVECTOR   
   DEC$          NOVECTOR   
   VD$           NOVECTOR   
   VDIR          NOVECTOR   
   VOCL          LOOP,SCALAR   
   IBM           PREFER SCALAR   
   $PL$ CMCHAR='*' */

		i__3 = na;
		for (ja = 1; ja <= i__3; ++ja) {
		    if (ilcplx) {
			sum_ref(ja, 1) = -acoef * suma_ref(ja, 1) + bcoefr * 
				sumb_ref(ja, 1) - bcoefi * sumb_ref(ja, 2);
			sum_ref(ja, 2) = -acoef * suma_ref(ja, 2) + bcoefr * 
				sumb_ref(ja, 2) + bcoefi * sumb_ref(ja, 1);
		    } else {
			sum_ref(ja, 1) = -acoef * suma_ref(ja, 1) + bcoefr * 
				sumb_ref(ja, 1);
		    }
/* L130: */
		}

/*                                  T   
                Solve  ( a A - b B )  y = SUM(,)   
                with scaling and perturbation of the denominator */

		dlaln2_(&c_true, &na, &nw, &dmin__, &acoef, &a_ref(j, j), lda,
			 bdiag, &bdiag[1], sum, &c__2, &bcoefr, &bcoefi, &
			work[(*n << 1) + j], n, &scale, &temp, &iinfo);
		if (scale < 1.) {
		    i__3 = nw - 1;
		    for (jw = 0; jw <= i__3; ++jw) {
			i__4 = j - 1;
			for (jr = je; jr <= i__4; ++jr) {
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
/* L140: */
			}
/* L150: */
		    }
		    xmax = scale * xmax;
		}
		xmax = max(xmax,temp);
L160:
		;
	    }

/*           Copy eigenvector to VL, back transforming if   
             HOWMNY='B'. */

	    ++ieig;
	    if (ilback) {
		i__2 = nw - 1;
		for (jw = 0; jw <= i__2; ++jw) {
		    i__3 = *n + 1 - je;
		    dgemv_("N", n, &i__3, &c_b35, &vl_ref(1, je), ldvl, &work[
			    (jw + 2) * *n + je], &c__1, &c_b37, &work[(jw + 4)
			     * *n + 1], &c__1);
/* L170: */
		}
		dlacpy_(" ", n, &nw, &work[(*n << 2) + 1], n, &vl_ref(1, je), 
			ldvl);
		ibeg = 1;
	    } else {
		dlacpy_(" ", n, &nw, &work[(*n << 1) + 1], n, &vl_ref(1, ieig)
			, ldvl);
		ibeg = je;
	    }

/*           Scale eigenvector */

	    xmax = 0.;
	    if (ilcplx) {
		i__2 = *n;
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
		    d__3 = xmax, d__4 = (d__1 = vl_ref(j, ieig), abs(d__1)) + 
			    (d__2 = vl_ref(j, ieig + 1), abs(d__2));
		    xmax = max(d__3,d__4);
/* L180: */
		}
	    } else {
		i__2 = *n;
		for (j = ibeg; j <= i__2; ++j) {
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = vl_ref(j, ieig), abs(d__1));
		    xmax = max(d__2,d__3);
/* L190: */
		}
	    }

	    if (xmax > safmin) {
		xscale = 1. / xmax;

		i__2 = nw - 1;
		for (jw = 0; jw <= i__2; ++jw) {
		    i__3 = *n;
		    for (jr = ibeg; jr <= i__3; ++jr) {
			vl_ref(jr, ieig + jw) = xscale * vl_ref(jr, ieig + jw)
				;
/* L200: */
		    }
/* L210: */
		}
	    }
	    ieig = ieig + nw - 1;

L220:
	    ;
	}
    }

/*     Right eigenvectors */

    if (compr) {
	ieig = im + 1;

/*        Main loop over eigenvalues */

	ilcplx = FALSE_;
	for (je = *n; je >= 1; --je) {

/*           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or   
             (b) this would be the second of a complex pair.   
             Check for complex eigenvalue, so as to be sure of which   
             entry(-ies) of SELECT to look at -- if complex, SELECT(JE)   
             or SELECT(JE-1).   
             If this is a complex pair, the 2-by-2 diagonal block   
             corresponding to the eigenvalue is in rows/columns JE-1:JE */

	    if (ilcplx) {
		ilcplx = FALSE_;
		goto L500;
	    }
	    nw = 1;
	    if (je > 1) {
		if (a_ref(je, je - 1) != 0.) {
		    ilcplx = TRUE_;
		    nw = 2;
		}
	    }
	    if (ilall) {
		ilcomp = TRUE_;
	    } else if (ilcplx) {
		ilcomp = select[je] || select[je - 1];
	    } else {
		ilcomp = select[je];
	    }
	    if (! ilcomp) {
		goto L500;
	    }

/*           Decide if (a) singular pencil, (b) real eigenvalue, or   
             (c) complex eigenvalue. */

	    if (! ilcplx) {
		if ((d__1 = a_ref(je, je), abs(d__1)) <= safmin && (d__2 = 
			b_ref(je, je), abs(d__2)) <= safmin) {

/*                 Singular matrix pencil -- unit eigenvector */

		    --ieig;
		    i__1 = *n;
		    for (jr = 1; jr <= i__1; ++jr) {
			vr_ref(jr, ieig) = 0.;
/* L230: */
		    }
		    vr_ref(ieig, ieig) = 1.;
		    goto L500;
		}
	    }

/*           Clear vector */

	    i__1 = nw - 1;
	    for (jw = 0; jw <= i__1; ++jw) {
		i__2 = *n;
		for (jr = 1; jr <= i__2; ++jr) {
		    work[(jw + 2) * *n + jr] = 0.;
/* L240: */
		}
/* L250: */
	    }

/*           Compute coefficients in  ( a A - b B ) x = 0   
                a  is  ACOEF   
                b  is  BCOEFR + i*BCOEFI */

	    if (! ilcplx) {

/*              Real eigenvalue   

   Computing MAX */
		d__3 = (d__1 = a_ref(je, je), abs(d__1)) * ascale, d__4 = (
			d__2 = b_ref(je, je), abs(d__2)) * bscale, d__3 = max(
			d__3,d__4);
		temp = 1. / max(d__3,safmin);
		salfar = temp * a_ref(je, je) * ascale;
		sbeta = temp * b_ref(je, je) * bscale;
		acoef = sbeta * ascale;
		bcoefr = salfar * bscale;
		bcoefi = 0.;

/*              Scale to avoid underflow */

		scale = 1.;
		lsa = abs(sbeta) >= safmin && abs(acoef) < small;
		lsb = abs(salfar) >= safmin && abs(bcoefr) < small;
		if (lsa) {
		    scale = small / abs(sbeta) * min(anorm,big);
		}
		if (lsb) {
/* Computing MAX */
		    d__1 = scale, d__2 = small / abs(salfar) * min(bnorm,big);
		    scale = max(d__1,d__2);
		}
		if (lsa || lsb) {
/* Computing MIN   
   Computing MAX */
		    d__3 = 1., d__4 = abs(acoef), d__3 = max(d__3,d__4), d__4 
			    = abs(bcoefr);
		    d__1 = scale, d__2 = 1. / (safmin * max(d__3,d__4));
		    scale = min(d__1,d__2);
		    if (lsa) {
			acoef = ascale * (scale * sbeta);
		    } else {
			acoef = scale * acoef;
		    }
		    if (lsb) {
			bcoefr = bscale * (scale * salfar);
		    } else {
			bcoefr = scale * bcoefr;
		    }
		}
		acoefa = abs(acoef);
		bcoefa = abs(bcoefr);

/*              First component is 1 */

		work[(*n << 1) + je] = 1.;
		xmax = 1.;

/*              Compute contribution from column JE of A and B to sum   
                (See "Further Details", above.) */

		i__1 = je - 1;
		for (jr = 1; jr <= i__1; ++jr) {
		    work[(*n << 1) + jr] = bcoefr * b_ref(jr, je) - acoef * 
			    a_ref(jr, je);
/* L260: */
		}
	    } else {

/*              Complex eigenvalue */

		d__1 = safmin * 100.;
		dlag2_(&a_ref(je - 1, je - 1), lda, &b_ref(je - 1, je - 1), 
			ldb, &d__1, &acoef, &temp, &bcoefr, &temp2, &bcoefi);
		if (bcoefi == 0.) {
		    *info = je - 1;
		    return 0;
		}

/*              Scale to avoid over/underflow */

		acoefa = abs(acoef);
		bcoefa = abs(bcoefr) + abs(bcoefi);
		scale = 1.;
		if (acoefa * ulp < safmin && acoefa >= safmin) {
		    scale = safmin / ulp / acoefa;
		}
		if (bcoefa * ulp < safmin && bcoefa >= safmin) {
/* Computing MAX */
		    d__1 = scale, d__2 = safmin / ulp / bcoefa;
		    scale = max(d__1,d__2);
		}
		if (safmin * acoefa > ascale) {
		    scale = ascale / (safmin * acoefa);
		}
		if (safmin * bcoefa > bscale) {
/* Computing MIN */
		    d__1 = scale, d__2 = bscale / (safmin * bcoefa);
		    scale = min(d__1,d__2);
		}
		if (scale != 1.) {
		    acoef = scale * acoef;
		    acoefa = abs(acoef);
		    bcoefr = scale * bcoefr;
		    bcoefi = scale * bcoefi;
		    bcoefa = abs(bcoefr) + abs(bcoefi);
		}

/*              Compute first two components of eigenvector   
                and contribution to sums */

		temp = acoef * a_ref(je, je - 1);
		temp2r = acoef * a_ref(je, je) - bcoefr * b_ref(je, je);
		temp2i = -bcoefi * b_ref(je, je);
		if (abs(temp) >= abs(temp2r) + abs(temp2i)) {
		    work[(*n << 1) + je] = 1.;
		    work[*n * 3 + je] = 0.;
		    work[(*n << 1) + je - 1] = -temp2r / temp;
		    work[*n * 3 + je - 1] = -temp2i / temp;
		} else {
		    work[(*n << 1) + je - 1] = 1.;
		    work[*n * 3 + je - 1] = 0.;
		    temp = acoef * a_ref(je - 1, je);
		    work[(*n << 1) + je] = (bcoefr * b_ref(je - 1, je - 1) - 
			    acoef * a_ref(je - 1, je - 1)) / temp;
		    work[*n * 3 + je] = bcoefi * b_ref(je - 1, je - 1) / temp;
		}

/* Computing MAX */
		d__5 = (d__1 = work[(*n << 1) + je], abs(d__1)) + (d__2 = 
			work[*n * 3 + je], abs(d__2)), d__6 = (d__3 = work[(*
			n << 1) + je - 1], abs(d__3)) + (d__4 = work[*n * 3 + 
			je - 1], abs(d__4));
		xmax = max(d__5,d__6);

/*              Compute contribution from columns JE and JE-1   
                of A and B to the sums. */

		creala = acoef * work[(*n << 1) + je - 1];
		cimaga = acoef * work[*n * 3 + je - 1];
		crealb = bcoefr * work[(*n << 1) + je - 1] - bcoefi * work[*n 
			* 3 + je - 1];
		cimagb = bcoefi * work[(*n << 1) + je - 1] + bcoefr * work[*n 
			* 3 + je - 1];
		cre2a = acoef * work[(*n << 1) + je];
		cim2a = acoef * work[*n * 3 + je];
		cre2b = bcoefr * work[(*n << 1) + je] - bcoefi * work[*n * 3 
			+ je];
		cim2b = bcoefi * work[(*n << 1) + je] + bcoefr * work[*n * 3 
			+ je];
		i__1 = je - 2;
		for (jr = 1; jr <= i__1; ++jr) {
		    work[(*n << 1) + jr] = -creala * a_ref(jr, je - 1) + 
			    crealb * b_ref(jr, je - 1) - cre2a * a_ref(jr, je)
			     + cre2b * b_ref(jr, je);
		    work[*n * 3 + jr] = -cimaga * a_ref(jr, je - 1) + cimagb *
			     b_ref(jr, je - 1) - cim2a * a_ref(jr, je) + 
			    cim2b * b_ref(jr, je);
/* L270: */
		}
	    }

/* Computing MAX */
	    d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, d__1 = 
		    max(d__1,d__2);
	    dmin__ = max(d__1,safmin);

/*           Columnwise triangular solve of  (a A - b B)  x = 0 */

	    il2by2 = FALSE_;
	    for (j = je - nw; j >= 1; --j) {

/*              If a 2-by-2 block, is in position j-1:j, wait until   
                next iteration to process it (when it will be j:j+1) */

		if (! il2by2 && j > 1) {
		    if (a_ref(j, j - 1) != 0.) {
			il2by2 = TRUE_;
			goto L370;
		    }
		}
		bdiag[0] = b_ref(j, j);
		if (il2by2) {
		    na = 2;
		    bdiag[1] = b_ref(j + 1, j + 1);
		} else {
		    na = 1;
		}

/*              Compute x(j) (and x(j+1), if 2-by-2 block) */

		dlaln2_(&c_false, &na, &nw, &dmin__, &acoef, &a_ref(j, j), 
			lda, bdiag, &bdiag[1], &work[(*n << 1) + j], n, &
			bcoefr, &bcoefi, sum, &c__2, &scale, &temp, &iinfo);
		if (scale < 1.) {

		    i__1 = nw - 1;
		    for (jw = 0; jw <= i__1; ++jw) {
			i__2 = je;
			for (jr = 1; jr <= i__2; ++jr) {
			    work[(jw + 2) * *n + jr] = scale * work[(jw + 2) *
				     *n + jr];
/* L280: */
			}
/* L290: */
		    }
		}
/* Computing MAX */
		d__1 = scale * xmax;
		xmax = max(d__1,temp);

		i__1 = nw;
		for (jw = 1; jw <= i__1; ++jw) {
		    i__2 = na;
		    for (ja = 1; ja <= i__2; ++ja) {
			work[(jw + 1) * *n + j + ja - 1] = sum_ref(ja, jw);
/* L300: */
		    }
/* L310: */
		}

/*              w = w + x(j)*(a A(*,j) - b B(*,j) ) with scaling */

		if (j > 1) {

/*                 Check whether scaling is necessary for sum. */

		    xscale = 1. / max(1.,xmax);
		    temp = acoefa * work[j] + bcoefa * work[*n + j];
		    if (il2by2) {
/* Computing MAX */
			d__1 = temp, d__2 = acoefa * work[j + 1] + bcoefa * 
				work[*n + j + 1];
			temp = max(d__1,d__2);
		    }
/* Computing MAX */
		    d__1 = max(temp,acoefa);
		    temp = max(d__1,bcoefa);
		    if (temp > bignum * xscale) {

			i__1 = nw - 1;
			for (jw = 0; jw <= i__1; ++jw) {
			    i__2 = je;
			    for (jr = 1; jr <= i__2; ++jr) {
				work[(jw + 2) * *n + jr] = xscale * work[(jw 
					+ 2) * *n + jr];
/* L320: */
			    }
/* L330: */
			}
			xmax *= xscale;
		    }

/*                 Compute the contributions of the off-diagonals of   
                   column j (and j+1, if 2-by-2 block) of A and B to the   
                   sums. */


		    i__1 = na;
		    for (ja = 1; ja <= i__1; ++ja) {
			if (ilcplx) {
			    creala = acoef * work[(*n << 1) + j + ja - 1];
			    cimaga = acoef * work[*n * 3 + j + ja - 1];
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1] - 
				    bcoefi * work[*n * 3 + j + ja - 1];
			    cimagb = bcoefi * work[(*n << 1) + j + ja - 1] + 
				    bcoefr * work[*n * 3 + j + ja - 1];
			    i__2 = j - 1;
			    for (jr = 1; jr <= i__2; ++jr) {
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * a_ref(jr, j + ja - 1) + 
					crealb * b_ref(jr, j + ja - 1);
				work[*n * 3 + jr] = work[*n * 3 + jr] - 
					cimaga * a_ref(jr, j + ja - 1) + 
					cimagb * b_ref(jr, j + ja - 1);
/* L340: */
			    }
			} else {
			    creala = acoef * work[(*n << 1) + j + ja - 1];
			    crealb = bcoefr * work[(*n << 1) + j + ja - 1];
			    i__2 = j - 1;
			    for (jr = 1; jr <= i__2; ++jr) {
				work[(*n << 1) + jr] = work[(*n << 1) + jr] - 
					creala * a_ref(jr, j + ja - 1) + 
					crealb * b_ref(jr, j + ja - 1);
/* L350: */
			    }
			}
/* L360: */
		    }
		}

		il2by2 = FALSE_;
L370:
		;
	    }

/*           Copy eigenvector to VR, back transforming if   
             HOWMNY='B'. */

	    ieig -= nw;
	    if (ilback) {

		i__1 = nw - 1;
		for (jw = 0; jw <= i__1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			work[(jw + 4) * *n + jr] = work[(jw + 2) * *n + 1] * 
				vr_ref(jr, 1);
/* L380: */
		    }

/*                 A series of compiler directives to defeat   
                   vectorization for the next loop */


		    i__2 = je;
		    for (jc = 2; jc <= i__2; ++jc) {
			i__3 = *n;
			for (jr = 1; jr <= i__3; ++jr) {
			    work[(jw + 4) * *n + jr] += work[(jw + 2) * *n + 
				    jc] * vr_ref(jr, jc);
/* L390: */
			}
/* L400: */
		    }
/* L410: */
		}

		i__1 = nw - 1;
		for (jw = 0; jw <= i__1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, ieig + jw) = work[(jw + 4) * *n + jr];
/* L420: */
		    }
/* L430: */
		}

		iend = *n;
	    } else {
		i__1 = nw - 1;
		for (jw = 0; jw <= i__1; ++jw) {
		    i__2 = *n;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, ieig + jw) = work[(jw + 2) * *n + jr];
/* L440: */
		    }
/* L450: */
		}

		iend = je;
	    }

/*           Scale eigenvector */

	    xmax = 0.;
	    if (ilcplx) {
		i__1 = iend;
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		    d__3 = xmax, d__4 = (d__1 = vr_ref(j, ieig), abs(d__1)) + 
			    (d__2 = vr_ref(j, ieig + 1), abs(d__2));
		    xmax = max(d__3,d__4);
/* L460: */
		}
	    } else {
		i__1 = iend;
		for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		    d__2 = xmax, d__3 = (d__1 = vr_ref(j, ieig), abs(d__1));
		    xmax = max(d__2,d__3);
/* L470: */
		}
	    }

	    if (xmax > safmin) {
		xscale = 1. / xmax;
		i__1 = nw - 1;
		for (jw = 0; jw <= i__1; ++jw) {
		    i__2 = iend;
		    for (jr = 1; jr <= i__2; ++jr) {
			vr_ref(jr, ieig + jw) = xscale * vr_ref(jr, ieig + jw)
				;
/* L480: */
		    }
/* L490: */
		}
	    }
L500:
	    ;
	}
    }

    return 0;

/*     End of DTGEVC */

} /* dtgevc_ */

#undef sum_ref
#undef vr_ref
#undef vl_ref
#undef b_ref
#undef a_ref
#undef sumb_ref
#undef suma_ref



/* Subroutine */ int dggbak_(char *job, char *side, integer *n, integer *ilo, 
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m, 
	doublereal *v, integer *ldv, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGBAK forms the right or left eigenvectors of a real generalized   
    eigenvalue problem A*x = lambda*B*x, by backward transformation on   
    the computed eigenvectors of the balanced pair of matrices output by   
    DGGBAL.   

    Arguments   
    =========   

    JOB     (input) CHARACTER*1   
            Specifies the type of backward transformation required:   
            = 'N':  do nothing, return immediately;   
            = 'P':  do backward transformation for permutation only;   
            = 'S':  do backward transformation for scaling only;   
            = 'B':  do backward transformations for both permutation and   
                    scaling.   
            JOB must be the same as the argument JOB supplied to DGGBAL.   

    SIDE    (input) CHARACTER*1   
            = 'R':  V contains right eigenvectors;   
            = 'L':  V contains left eigenvectors.   

    N       (input) INTEGER   
            The number of rows of the matrix V.  N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            The integers ILO and IHI determined by DGGBAL.   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    LSCALE  (input) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the left side of A and B, as returned by DGGBAL.   

    RSCALE  (input) DOUBLE PRECISION array, dimension (N)   
            Details of the permutations and/or scaling factors applied   
            to the right side of A and B, as returned by DGGBAL.   

    M       (input) INTEGER   
            The number of columns of the matrix V.  M >= 0.   

    V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)   
            On entry, the matrix of right or left eigenvectors to be   
            transformed, as returned by DTGEVC.   
            On exit, V is overwritten by the transformed eigenvectors.   

    LDV     (input) INTEGER   
            The leading dimension of the matrix V. LDV >= max(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    See R.C. Ward, Balancing the generalized eigenvalue problem,   
                   SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.   

    =====================================================================   


       Test the input parameters   

       Parameter adjustments */
    /* System generated locals */
    integer v_dim1, v_offset, i__1;
    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical leftv;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical rightv;
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]

    --lscale;
    --rscale;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;

    /* Function Body */
    rightv = lsame_(side, "R");
    leftv = lsame_(side, "L");

    *info = 0;
    if (! lsame_(job, "N") && ! lsame_(job, "P") && ! lsame_(job, "S") 
	    && ! lsame_(job, "B")) {
	*info = -1;
    } else if (! rightv && ! leftv) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*ilo < 1) {
	*info = -4;
    } else if (*ihi < *ilo || *ihi > max(1,*n)) {
	*info = -5;
    } else if (*m < 0) {
	*info = -6;
    } else if (*ldv < max(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGBAK", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return 0;
    }
    if (*m == 0) {
	return 0;
    }
    if (lsame_(job, "N")) {
	return 0;
    }

    if (*ilo == *ihi) {
	goto L30;
    }

/*     Backward balance */

    if (lsame_(job, "S") || lsame_(job, "B")) {

/*        Backward transformation on right eigenvectors */

	if (rightv) {
	    i__1 = *ihi;
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
		dscal_(m, &rscale[i__], &v_ref(i__, 1), ldv);
/* L10: */
	    }
	}

/*        Backward transformation on left eigenvectors */

	if (leftv) {
	    i__1 = *ihi;
	    for (i__ = *ilo; i__ <= i__1; ++i__) {
		dscal_(m, &lscale[i__], &v_ref(i__, 1), ldv);
/* L20: */
	    }
	}
    }

/*     Backward permutation */

L30:
    if (lsame_(job, "P") || lsame_(job, "B")) {

/*        Backward permutation on right eigenvectors */

	if (rightv) {
	    if (*ilo == 1) {
		goto L50;
	    }

	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
		k = (integer) rscale[i__];
		if (k == i__) {
		    goto L40;
		}
		dswap_(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L40:
		;
	    }

L50:
	    if (*ihi == *n) {
		goto L70;
	    }
	    i__1 = *n;
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
		k = (integer) rscale[i__];
		if (k == i__) {
		    goto L60;
		}
		dswap_(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L60:
		;
	    }
	}

/*        Backward permutation on left eigenvectors */

L70:
	if (leftv) {
	    if (*ilo == 1) {
		goto L90;
	    }
	    for (i__ = *ilo - 1; i__ >= 1; --i__) {
		k = (integer) lscale[i__];
		if (k == i__) {
		    goto L80;
		}
		dswap_(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L80:
		;
	    }

L90:
	    if (*ihi == *n) {
		goto L110;
	    }
	    i__1 = *n;
	    for (i__ = *ihi + 1; i__ <= i__1; ++i__) {
		k = (integer) lscale[i__];
		if (k == i__) {
		    goto L100;
		}
		dswap_(m, &v_ref(i__, 1), ldv, &v_ref(k, 1), ldv);
L100:
		;
	    }
	}
    }

L110:

    return 0;

/*     End of DGGBAK */

} /* dggbak_ */

#undef v_ref



integer ieeeck_(integer *ispec, real *zero, real *one)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1998   


    Purpose   
    =======   

    IEEECK is called from the ILAENV to verify that Infinity and   
    possibly NaN arithmetic is safe (i.e. will not trap).   

    Arguments   
    =========   

    ISPEC   (input) INTEGER   
            Specifies whether to test just for inifinity arithmetic   
            or whether to test for infinity and NaN arithmetic.   
            = 0: Verify infinity arithmetic only.   
            = 1: Verify infinity and NaN arithmetic.   

    ZERO    (input) REAL   
            Must contain the value 0.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    ONE     (input) REAL   
            Must contain the value 1.0   
            This is passed to prevent the compiler from optimizing   
            away this code.   

    RETURN VALUE:  INTEGER   
            = 0:  Arithmetic failed to produce the correct answers   
            = 1:  Arithmetic produced the correct answers */
    /* System generated locals */
    integer ret_val;
    /* Local variables */
    static real neginf, posinf, negzro, newzro, nan1, nan2, nan3, nan4, nan5, 
	    nan6;


    ret_val = 1;

    posinf = *one / *zero;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf = -(*one) / *zero;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    negzro = *one / (neginf + *one);
    if (negzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    neginf = *one / negzro;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    newzro = negzro + *zero;
    if (newzro != *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf = *one / newzro;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }

    neginf *= posinf;
    if (neginf >= *zero) {
	ret_val = 0;
	return ret_val;
    }

    posinf *= posinf;
    if (posinf <= *one) {
	ret_val = 0;
	return ret_val;
    }




/*     Return if we were only asked to check infinity arithmetic */

    if (*ispec == 0) {
	return ret_val;
    }

    nan1 = posinf + neginf;

    nan2 = posinf / neginf;

    nan3 = posinf / posinf;

    nan4 = posinf * *zero;

    nan5 = neginf * negzro;

    nan6 = nan5 * 0.f;

    if (nan1 == nan1) {
	ret_val = 0;
	return ret_val;
    }

    if (nan2 == nan2) {
	ret_val = 0;
	return ret_val;
    }

    if (nan3 == nan3) {
	ret_val = 0;
	return ret_val;
    }

    if (nan4 == nan4) {
	ret_val = 0;
	return ret_val;
    }

    if (nan5 == nan5) {
	ret_val = 0;
	return ret_val;
    }

    if (nan6 == nan6) {
	ret_val = 0;
	return ret_val;
    }

    return ret_val;
} /* ieeeck_ */


/* Subroutine */ int dlassq_(integer *n, doublereal *x, integer *incx, 
	doublereal *scale, doublereal *sumsq)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1999   


    Purpose   
    =======   

    DLASSQ  returns the values  scl  and  smsq  such that   

       ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,   

    where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is   
    assumed to be non-negative and  scl  returns the value   

       scl = max( scale, abs( x( i ) ) ).   

    scale and sumsq must be supplied in SCALE and SUMSQ and   
    scl and smsq are overwritten on SCALE and SUMSQ respectively.   

    The routine makes only one pass through the vector x.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of elements to be used from the vector X.   

    X       (input) DOUBLE PRECISION array, dimension (N)   
            The vector for which a scaled sum of squares is computed.   
               x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.   

    INCX    (input) INTEGER   
            The increment between successive values of the vector X.   
            INCX > 0.   

    SCALE   (input/output) DOUBLE PRECISION   
            On entry, the value  scale  in the equation above.   
            On exit, SCALE is overwritten with  scl , the scaling factor   
            for the sum of squares.   

    SUMSQ   (input/output) DOUBLE PRECISION   
            On entry, the value  sumsq  in the equation above.   
            On exit, SUMSQ is overwritten with  smsq , the basic sum of   
            squares from which  scl  has been factored out.   

   =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    /* Local variables */
    static doublereal absxi;
    static integer ix;

    --x;

    /* Function Body */
    if (*n > 0) {
	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
		absxi = (d__1 = x[ix], abs(d__1));
		if (*scale < absxi) {
/* Computing 2nd power */
		    d__1 = *scale / absxi;
		    *sumsq = *sumsq * (d__1 * d__1) + 1;
		    *scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / *scale;
		    *sumsq += d__1 * d__1;
		}
	    }
/* L10: */
	}
    }
    return 0;

/*     End of DLASSQ */

} /* dlassq_ */


/* Subroutine */ int dgeqr2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGEQR2 computes a QR factorization of a real m by n matrix A:   
    A = Q * R.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and above the diagonal of the array   
            contain the min(m,n) by n upper trapezoidal matrix R (R is   
            upper triangular if m >= n); the elements below the diagonal,   
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(1,M).   

    TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))   
            The scalar factors of the elementary reflectors (see Further   
            Details).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = min(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),   
    and tau in TAU(i).   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *), dlarfg_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *), xerbla_(char *, integer *);
    static doublereal aii;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < max(1,*m)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQR2", &i__1);
	return 0;
    }

    k = min(*m,*n);

    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)   

   Computing MIN */
	i__2 = i__ + 1;
	i__3 = *m - i__ + 1;
	dlarfg_(&i__3, &a_ref(i__, i__), &a_ref(min(i__2,*m), i__), &c__1, &
		tau[i__]);
	if (i__ < *n) {

/*           Apply H(i) to A(i:m,i+1:n) from the left */

	    aii = a_ref(i__, i__);
	    a_ref(i__, i__) = 1.;
	    i__2 = *m - i__ + 1;
	    i__3 = *n - i__;
	    dlarf_("Left", &i__2, &i__3, &a_ref(i__, i__), &c__1, &tau[i__], &
		    a_ref(i__, i__ + 1), lda, &work[1]);
	    a_ref(i__, i__) = aii;
	}
/* L10: */
    }
    return 0;

/*     End of DGEQR2 */

} /* dgeqr2_ */

#undef a_ref



/* Subroutine */ int dlarft_(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t, 
	integer *ldt)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFT forms the triangular factor T of a real block reflector H   
    of order n, which is defined as a product of k elementary reflectors.   

    If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;   

    If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.   

    If STOREV = 'C', the vector which defines the elementary reflector   
    H(i) is stored in the i-th column of the array V, and   

       H  =  I - V * T * V'   

    If STOREV = 'R', the vector which defines the elementary reflector   
    H(i) is stored in the i-th row of the array V, and   

       H  =  I - V' * T * V   

    Arguments   
    =========   

    DIRECT  (input) CHARACTER*1   
            Specifies the order in which the elementary reflectors are   
            multiplied to form the block reflector:   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Specifies how the vectors which define the elementary   
            reflectors are stored (see also Further Details):   
            = 'C': columnwise   
            = 'R': rowwise   

    N       (input) INTEGER   
            The order of the block reflector H. N >= 0.   

    K       (input) INTEGER   
            The order of the triangular factor T (= the number of   
            elementary reflectors). K >= 1.   

    V       (input/output) DOUBLE PRECISION array, dimension   
                                 (LDV,K) if STOREV = 'C'   
                                 (LDV,N) if STOREV = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i).   

    T       (output) DOUBLE PRECISION array, dimension (LDT,K)   
            The k by k triangular factor T of the block reflector.   
            If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is   
            lower triangular. The rest of the array is not used.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    Further Details   
    ===============   

    The shape of the matrix V and the storage of the vectors which define   
    the H(i) is best illustrated by the following example with n = 5 and   
    k = 3. The elements equal to 1 are not stored; the corresponding   
    array elements are modified but restored on exit. The rest of the   
    array is not used.   

    DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':   

                 V = (  1       )                 V = (  1 v1 v1 v1 v1 )   
                     ( v1  1    )                     (     1 v2 v2 v2 )   
                     ( v1 v2  1 )                     (        1 v3 v3 )   
                     ( v1 v2 v3 )   
                     ( v1 v2 v3 )   

    DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':   

                 V = ( v1 v2 v3 )                 V = ( v1 v1  1       )   
                     ( v1 v2 v3 )                     ( v2 v2 v2  1    )   
                     (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )   
                     (     1 v3 )   
                     (        1 )   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b8 = 0.;
    
    /* System generated locals */
    integer t_dim1, t_offset, v_dim1, v_offset, i__1, i__2, i__3;
    doublereal d__1;
    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dtrmv_(char *, 
	    char *, char *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal vii;
#define t_ref(a_1,a_2) t[(a_2)*t_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]


    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    --tau;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;

    /* Function Body */
    if (*n == 0) {
	return 0;
    }

    if (lsame_(direct, "F")) {
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (tau[i__] == 0.) {

/*              H(i)  =  I */

		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    t_ref(j, i__) = 0.;
/* L10: */
		}
	    } else {

/*              general case */

		vii = v_ref(i__, i__);
		v_ref(i__, i__) = 1.;
		if (lsame_(storev, "C")) {

/*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i) */

		    i__2 = *n - i__ + 1;
		    i__3 = i__ - 1;
		    d__1 = -tau[i__];
		    dgemv_("Transpose", &i__2, &i__3, &d__1, &v_ref(i__, 1), 
			    ldv, &v_ref(i__, i__), &c__1, &c_b8, &t_ref(1, 
			    i__), &c__1);
		} else {

/*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)' */

		    i__2 = i__ - 1;
		    i__3 = *n - i__ + 1;
		    d__1 = -tau[i__];
		    dgemv_("No transpose", &i__2, &i__3, &d__1, &v_ref(1, i__)
			    , ldv, &v_ref(i__, i__), ldv, &c_b8, &t_ref(1, 
			    i__), &c__1);
		}
		v_ref(i__, i__) = vii;

/*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i) */

		i__2 = i__ - 1;
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[
			t_offset], ldt, &t_ref(1, i__), &c__1);
		t_ref(i__, i__) = tau[i__];
	    }
/* L20: */
	}
    } else {
	for (i__ = *k; i__ >= 1; --i__) {
	    if (tau[i__] == 0.) {

/*              H(i)  =  I */

		i__1 = *k;
		for (j = i__; j <= i__1; ++j) {
		    t_ref(j, i__) = 0.;
/* L30: */
		}
	    } else {

/*              general case */

		if (i__ < *k) {
		    if (lsame_(storev, "C")) {
			vii = v_ref(*n - *k + i__, i__);
			v_ref(*n - *k + i__, i__) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i) */

			i__1 = *n - *k + i__;
			i__2 = *k - i__;
			d__1 = -tau[i__];
			dgemv_("Transpose", &i__1, &i__2, &d__1, &v_ref(1, 
				i__ + 1), ldv, &v_ref(1, i__), &c__1, &c_b8, &
				t_ref(i__ + 1, i__), &c__1);
			v_ref(*n - *k + i__, i__) = vii;
		    } else {
			vii = v_ref(i__, *n - *k + i__);
			v_ref(i__, *n - *k + i__) = 1.;

/*                    T(i+1:k,i) :=   
                              - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)' */

			i__1 = *k - i__;
			i__2 = *n - *k + i__;
			d__1 = -tau[i__];
			dgemv_("No transpose", &i__1, &i__2, &d__1, &v_ref(
				i__ + 1, 1), ldv, &v_ref(i__, 1), ldv, &c_b8, 
				&t_ref(i__ + 1, i__), &c__1);
			v_ref(i__, *n - *k + i__) = vii;
		    }

/*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i) */

		    i__1 = *k - i__;
		    dtrmv_("Lower", "No transpose", "Non-unit", &i__1, &t_ref(
			    i__ + 1, i__ + 1), ldt, &t_ref(i__ + 1, i__), &
			    c__1);
		}
		t_ref(i__, i__) = tau[i__];
	    }
/* L40: */
	}
    }
    return 0;

/*     End of DLARFT */

} /* dlarft_ */

#undef v_ref
#undef t_ref



/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc, 
	doublereal *work, integer *ldwork)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARFB applies a real block reflector H or its transpose H' to a   
    real m by n matrix C, from either the left or the right.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply H or H' from the Left   
            = 'R': apply H or H' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply H (No transpose)   
            = 'T': apply H' (Transpose)   

    DIRECT  (input) CHARACTER*1   
            Indicates how H is formed from a product of elementary   
            reflectors   
            = 'F': H = H(1) H(2) . . . H(k) (Forward)   
            = 'B': H = H(k) . . . H(2) H(1) (Backward)   

    STOREV  (input) CHARACTER*1   
            Indicates how the vectors which define the elementary   
            reflectors are stored:   
            = 'C': Columnwise   
            = 'R': Rowwise   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    K       (input) INTEGER   
            The order of the matrix T (= the number of elementary   
            reflectors whose product defines the block reflector).   

    V       (input) DOUBLE PRECISION array, dimension   
                                  (LDV,K) if STOREV = 'C'   
                                  (LDV,M) if STOREV = 'R' and SIDE = 'L'   
                                  (LDV,N) if STOREV = 'R' and SIDE = 'R'   
            The matrix V. See further details.   

    LDV     (input) INTEGER   
            The leading dimension of the array V.   
            If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);   
            if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);   
            if STOREV = 'R', LDV >= K.   

    T       (input) DOUBLE PRECISION array, dimension (LDT,K)   
            The triangular k by k matrix T in the representation of the   
            block reflector.   

    LDT     (input) INTEGER   
            The leading dimension of the array T. LDT >= K.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by H*C or H'*C or C*H or C*H'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDA >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)   

    LDWORK  (input) INTEGER   
            The leading dimension of the array WORK.   
            If SIDE = 'L', LDWORK >= max(1,N);   
            if SIDE = 'R', LDWORK >= max(1,M).   

    =====================================================================   


       Quick return if possible   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    static doublereal c_b14 = 1.;
    static doublereal c_b25 = -1.;
    
    /* System generated locals */
    integer c_dim1, c_offset, t_dim1, t_offset, v_dim1, v_offset, work_dim1, 
	    work_offset, i__1, i__2;
    /* Local variables */
    static integer i__, j;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char transt[1];
#define work_ref(a_1,a_2) work[(a_2)*work_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
#define v_ref(a_1,a_2) v[(a_2)*v_dim1 + a_1]


    v_dim1 = *ldv;
    v_offset = 1 + v_dim1 * 1;
    v -= v_offset;
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    work_dim1 = *ldwork;
    work_offset = 1 + work_dim1 * 1;
    work -= work_offset;

    /* Function Body */
    if (*m <= 0 || *n <= 0) {
	return 0;
    }

    if (lsame_(trans, "N")) {
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

    if (lsame_(storev, "C")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1 )    (first K rows)   
                       ( V2 )   
             where  V1  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c___ref(j, 1), ldc, &work_ref(1, j), &c__1);
/* L10: */
		}

/*              W := W * V1 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2 */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c___ref(*k + 1, 1), ldc, &v_ref(*k + 1, 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2 * W' */

		    i__1 = *m - *k;
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v_ref(*k + 1, 1), ldv, &work[work_offset], ldwork,
			     &c_b14, &c___ref(*k + 1, 1), ldc);
		}

/*              W := W * V1' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(j, i__) = c___ref(j, i__) - work_ref(i__, j);
/* L20: */
		    }
/* L30: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c___ref(1, j), &c__1, &work_ref(1, j), &c__1);
/* L40: */
		}

/*              W := W * V1 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c___ref(1, *k + 1), ldc, &v_ref(*k + 1, 1)
			    , ldv, &c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C2 := C2 - W * V2' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v_ref(*k + 1, 1), ldv,
			     &c_b14, &c___ref(1, *k + 1), ldc);
		}

/*              W := W * V1' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = c___ref(i__, j) - work_ref(i__, j);
/* L50: */
		    }
/* L60: */
		}
	    }

	} else {

/*           Let  V =  ( V1 )   
                       ( V2 )    (last K rows)   
             where  V2  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c___ref(*m - *k + j, 1), ldc, &work_ref(1, j), 
			    &c__1);
/* L70: */
		}

/*              W := W * V2 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v_ref(*m - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1 */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "No transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1 * W' */

		    i__1 = *m - *k;
		    dgemm_("No transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v[v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc)
			    ;
		}

/*              W := W * V2' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v_ref(*m - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(*m - *k + j, i__) = c___ref(*m - *k + j, i__) 
				- work_ref(i__, j);
/* L80: */
		    }
/* L90: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c___ref(1, *n - *k + j), &c__1, &work_ref(1, j)
			    , &c__1);
/* L100: */
		}

/*              W := W * V2 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v_ref(*n - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, k, &i__1, &
			    c_b14, &c__[c_offset], ldc, &v[v_offset], ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V' */

		if (*n > *k) {

/*                 C1 := C1 - W * V1' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, &i__1, k, &c_b25, &
			    work[work_offset], ldwork, &v[v_offset], ldv, &
			    c_b14, &c__[c_offset], ldc)
			    ;
		}

/*              W := W * V2' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v_ref(*n - *k + 1, 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, *n - *k + j) = c___ref(i__, *n - *k + j) 
				- work_ref(i__, j);
/* L110: */
		    }
/* L120: */
		}
	    }
	}

    } else if (lsame_(storev, "R")) {

	if (lsame_(direct, "F")) {

/*           Let  V =  ( V1  V2 )    (V1: first K columns)   
             where  V1  is unit upper triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

                W := C1' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c___ref(j, 1), ldc, &work_ref(1, j), &c__1);
/* L130: */
		}

/*              W := W * V1' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*m > *k) {

/*                 W := W + C2'*V2' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c___ref(*k + 1, 1), ldc, &v_ref(1, *k + 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Upper", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C2 := C2 - V2' * W' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &
			    v_ref(1, *k + 1), ldv, &work[work_offset], ldwork,
			     &c_b14, &c___ref(*k + 1, 1), ldc);
		}

/*              W := W * V1 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", n, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(j, i__) = c___ref(j, i__) - work_ref(i__, j);
/* L140: */
		    }
/* L150: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

                W := C1 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c___ref(1, j), &c__1, &work_ref(1, j), &c__1);
/* L160: */
		}

/*              W := W * V1' */

		dtrmm_("Right", "Upper", "Transpose", "Unit", m, k, &c_b14, &
			v[v_offset], ldv, &work[work_offset], ldwork);
		if (*n > *k) {

/*                 W := W + C2 * V2' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c___ref(1, *k + 1), ldc, &v_ref(1, *k + 1), ldv, &
			    c_b14, &work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Upper", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C2 := C2 - W * V2 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v_ref(1, *k + 
			    1), ldv, &c_b14, &c___ref(1, *k + 1), ldc);
		}

/*              W := W * V1 */

		dtrmm_("Right", "Upper", "No transpose", "Unit", m, k, &c_b14,
			 &v[v_offset], ldv, &work[work_offset], ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = c___ref(i__, j) - work_ref(i__, j);
/* L170: */
		    }
/* L180: */
		}

	    }

	} else {

/*           Let  V =  ( V1  V2 )    (V2: last K columns)   
             where  V2  is unit lower triangular. */

	    if (lsame_(side, "L")) {

/*              Form  H * C  or  H' * C  where  C = ( C1 )   
                                                    ( C2 )   

                W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)   

                W := C2' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(n, &c___ref(*m - *k + j, 1), ldc, &work_ref(1, j), 
			    &c__1);
/* L190: */
		}

/*              W := W * V2' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b14, &
			v_ref(1, *m - *k + 1), ldv, &work[work_offset], 
			ldwork);
		if (*m > *k) {

/*                 W := W + C1'*V1' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T'  or  W * T */

		dtrmm_("Right", "Lower", transt, "Non-unit", n, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - V' * W' */

		if (*m > *k) {

/*                 C1 := C1 - V1' * W' */

		    i__1 = *m - *k;
		    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b25, &v[
			    v_offset], ldv, &work[work_offset], ldwork, &
			    c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b14,
			 &v_ref(1, *m - *k + 1), ldv, &work[work_offset], 
			ldwork);

/*              C2 := C2 - W' */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *n;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(*m - *k + j, i__) = c___ref(*m - *k + j, i__) 
				- work_ref(i__, j);
/* L200: */
		    }
/* L210: */
		}

	    } else if (lsame_(side, "R")) {

/*              Form  C * H  or  C * H'  where  C = ( C1  C2 )   

                W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)   

                W := C2 */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    dcopy_(m, &c___ref(1, *n - *k + j), &c__1, &work_ref(1, j)
			    , &c__1);
/* L220: */
		}

/*              W := W * V2' */

		dtrmm_("Right", "Lower", "Transpose", "Unit", m, k, &c_b14, &
			v_ref(1, *n - *k + 1), ldv, &work[work_offset], 
			ldwork);
		if (*n > *k) {

/*                 W := W + C1 * V1' */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "Transpose", m, k, &i__1, &c_b14, &
			    c__[c_offset], ldc, &v[v_offset], ldv, &c_b14, &
			    work[work_offset], ldwork);
		}

/*              W := W * T  or  W * T' */

		dtrmm_("Right", "Lower", trans, "Non-unit", m, k, &c_b14, &t[
			t_offset], ldt, &work[work_offset], ldwork);

/*              C := C - W * V */

		if (*n > *k) {

/*                 C1 := C1 - W * V1 */

		    i__1 = *n - *k;
		    dgemm_("No transpose", "No transpose", m, &i__1, k, &
			    c_b25, &work[work_offset], ldwork, &v[v_offset], 
			    ldv, &c_b14, &c__[c_offset], ldc);
		}

/*              W := W * V2 */

		dtrmm_("Right", "Lower", "No transpose", "Unit", m, k, &c_b14,
			 &v_ref(1, *n - *k + 1), ldv, &work[work_offset], 
			ldwork);

/*              C1 := C1 - W */

		i__1 = *k;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, *n - *k + j) = c___ref(i__, *n - *k + j) 
				- work_ref(i__, j);
/* L230: */
		    }
/* L240: */
		}

	    }

	}
    }

    return 0;

/*     End of DLARFB */

} /* dlarfb_ */

#undef v_ref
#undef c___ref
#undef work_ref



/* Subroutine */ int dorm2r_(char *side, char *trans, integer *m, integer *n, 
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORM2R overwrites the general real m by n matrix C with   

          Q * C  if SIDE = 'L' and TRANS = 'N', or   

          Q'* C  if SIDE = 'L' and TRANS = 'T', or   

          C * Q  if SIDE = 'R' and TRANS = 'N', or   

          C * Q' if SIDE = 'R' and TRANS = 'T',   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(1) H(2) . . . H(k)   

    as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q' from the Left   
            = 'R': apply Q or Q' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply Q  (No transpose)   
            = 'T': apply Q' (Transpose)   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by   
            DGEQRF in the first k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= max(1,M);   
            if SIDE = 'R', LDA >= max(1,N).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                                     (N) if SIDE = 'L',   
                                     (M) if SIDE = 'R'   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;
    /* Local variables */
    static logical left;
    static integer i__;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *);
    extern logical lsame_(char *, char *);
    static integer i1, i2, i3, ic, jc, mi, ni, nq;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical notran;
    static doublereal aii;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! lsame_(side, "R")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < max(1,nq)) {
	*info = -7;
    } else if (*ldc < max(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORM2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return 0;
    }

    if (left && ! notran || ! left && notran) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
    }

    if (left) {
	ni = *n;
	jc = 1;
    } else {
	mi = *m;
	ic = 1;
    }

    i__1 = i2;
    i__2 = i3;
    for (i__ = i1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	if (left) {

/*           H(i) is applied to C(i:m,1:n) */

	    mi = *m - i__ + 1;
	    ic = i__;
	} else {

/*           H(i) is applied to C(1:m,i:n) */

	    ni = *n - i__ + 1;
	    jc = i__;
	}

/*        Apply H(i) */

	aii = a_ref(i__, i__);
	a_ref(i__, i__) = 1.;
	dlarf_(side, &mi, &ni, &a_ref(i__, i__), &c__1, &tau[i__], &c___ref(
		ic, jc), ldc, &work[1]);
	a_ref(i__, i__) = aii;
/* L10: */
    }
    return 0;

/*     End of DORM2R */

} /* dorm2r_ */

#undef c___ref
#undef a_ref



/* Subroutine */ int dorg2r_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info)
{
/*  -- LAPACK routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORG2R generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary   
    reflectors of order m   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the   
            matrix Q. N >= K >= 0.   

    A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as   
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= max(1,M).   

    TAU     (input) DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace) DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    =====================================================================   


       Test the input arguments   

       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    /* Local variables */
    static integer i__, j, l;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *), xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --tau;
    --work;

    /* Function Body */
    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < max(1,*m)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORG2R", &i__1);
	return 0;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return 0;
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

    i__1 = *n;
    for (j = *k + 1; j <= i__1; ++j) {
	i__2 = *m;
	for (l = 1; l <= i__2; ++l) {
	    a_ref(l, j) = 0.;
/* L10: */
	}
	a_ref(j, j) = 1.;
/* L20: */
    }

    for (i__ = *k; i__ >= 1; --i__) {

/*        Apply H(i) to A(i:m,i:n) from the left */

	if (i__ < *n) {
	    a_ref(i__, i__) = 1.;
	    i__1 = *m - i__ + 1;
	    i__2 = *n - i__;
	    dlarf_("Left", &i__1, &i__2, &a_ref(i__, i__), &c__1, &tau[i__], &
		    a_ref(i__, i__ + 1), lda, &work[1]);
	}
	if (i__ < *m) {
	    i__1 = *m - i__;
	    d__1 = -tau[i__];
	    dscal_(&i__1, &d__1, &a_ref(i__ + 1, i__), &c__1);
	}
	a_ref(i__, i__) = 1. - tau[i__];

/*        Set A(1:i-1,i) to zero */

	i__1 = i__ - 1;
	for (l = 1; l <= i__1; ++l) {
	    a_ref(l, i__) = 0.;
/* L30: */
	}
/* L40: */
    }
    return 0;

/*     End of DORG2R */

} /* dorg2r_ */

#undef a_ref



/* Subroutine */ int dlartg_(doublereal *f, doublereal *g, doublereal *cs, 
	doublereal *sn, doublereal *r__)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARTG generate a plane rotation so that   

       [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.   
       [ -SN  CS  ]     [ G ]     [ 0 ]   

    This is a slower, more accurate version of the BLAS1 routine DROTG,   
    with the following other differences:   
       F and G are unchanged on return.   
       If G=0, then CS=1 and SN=0.   
       If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any   
          floating point operations (saves work in DBDSQR when   
          there are zeros on the diagonal).   

    If F exceeds G in magnitude, CS will be positive.   

    Arguments   
    =========   

    F       (input) DOUBLE PRECISION   
            The first component of vector to be rotated.   

    G       (input) DOUBLE PRECISION   
            The second component of vector to be rotated.   

    CS      (output) DOUBLE PRECISION   
            The cosine of the rotation.   

    SN      (output) DOUBLE PRECISION   
            The sine of the rotation.   

    R       (output) DOUBLE PRECISION   
            The nonzero component of the rotated vector.   

    ===================================================================== */
    /* Initialized data */
    static logical first = TRUE_;
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;
    /* Builtin functions */
    double log(doublereal), pow_di(doublereal *, integer *), sqrt(doublereal);
    /* Local variables */
    static integer i__;
    static doublereal scale;
    static integer count;
    static doublereal f1, g1, safmn2, safmx2;
    extern doublereal dlamch_(char *);
    static doublereal safmin, eps;



    if (first) {
	first = FALSE_;
	safmin = dlamch_("S");
	eps = dlamch_("E");
	d__1 = dlamch_("B");
	i__1 = (integer) (log(safmin / eps) / log(dlamch_("B")) / 
		2.);
	safmn2 = pow_di(&d__1, &i__1);
	safmx2 = 1. / safmn2;
    }
    if (*g == 0.) {
	*cs = 1.;
	*sn = 0.;
	*r__ = *f;
    } else if (*f == 0.) {
	*cs = 0.;
	*sn = 1.;
	*r__ = *g;
    } else {
	f1 = *f;
	g1 = *g;
/* Computing MAX */
	d__1 = abs(f1), d__2 = abs(g1);
	scale = max(d__1,d__2);
	if (scale >= safmx2) {
	    count = 0;
L10:
	    ++count;
	    f1 *= safmn2;
	    g1 *= safmn2;
/* Computing MAX */
	    d__1 = abs(f1), d__2 = abs(g1);
	    scale = max(d__1,d__2);
	    if (scale >= safmx2) {
		goto L10;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	    i__1 = count;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		*r__ *= safmx2;
/* L20: */
	    }
	} else if (scale <= safmn2) {
	    count = 0;
L30:
	    ++count;
	    f1 *= safmx2;
	    g1 *= safmx2;
/* Computing MAX */
	    d__1 = abs(f1), d__2 = abs(g1);
	    scale = max(d__1,d__2);
	    if (scale <= safmn2) {
		goto L30;
	    }
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	    i__1 = count;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		*r__ *= safmn2;
/* L40: */
	    }
	} else {
/* Computing 2nd power */
	    d__1 = f1;
/* Computing 2nd power */
	    d__2 = g1;
	    *r__ = sqrt(d__1 * d__1 + d__2 * d__2);
	    *cs = f1 / *r__;
	    *sn = g1 / *r__;
	}
	if (abs(*f) > abs(*g) && *cs < 0.) {
	    *cs = -(*cs);
	    *sn = -(*sn);
	    *r__ = -(*r__);
	}
    }
    return 0;

/*     End of DLARTG */

} /* dlartg_ */


doublereal dlanhs_(char *norm, integer *n, doublereal *a, integer *lda, 
	doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANHS  returns the value of the one norm,  or the Frobenius norm, or   
    the  infinity norm,  or the  element of  largest absolute value  of a   
    Hessenberg matrix A.   

    Description   
    ===========   

    DLANHS returns the value   

       DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum),   
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and   
    normF  denotes the  Frobenius norm of a matrix (square root of sum of   
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANHS as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANHS is   
            set to zero.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
            The n by n upper Hessenberg matrix A; the part of A below the   
            first sub-diagonal is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= max(N,1).   

    WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static integer c__1 = 1;
    
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static integer i__, j;
    static doublereal scale;
    extern logical lsame_(char *, char *);
    static doublereal value;
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);
    static doublereal sum;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]


    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --work;

    /* Function Body */
    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find max(abs(A(i,j))). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = a_ref(i__, j), abs(d__1));
		value = max(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)
	    norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    sum = 0.;
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sum += (d__1 = a_ref(i__, j), abs(d__1));
/* L30: */
	    }
	    value = max(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    work[i__] = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    for (i__ = 1; i__ <= i__2; ++i__) {
		work[i__] += (d__1 = a_ref(i__, j), abs(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	    d__1 = value, d__2 = work[i__];
	    value = max(d__1,d__2);
/* L80: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = min(i__3,i__4);
	    dlassq_(&i__2, &a_ref(1, j), &c__1, &scale, &sum);
/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANHS */

} /* dlanhs_ */

#undef a_ref



/* Subroutine */ int dlag2_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	scale2, doublereal *wr1, doublereal *wr2, doublereal *wi)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue   
    problem  A - w B, with scaling as necessary to avoid over-/underflow.   

    The scaling factor "s" results in a modified eigenvalue equation   

        s A - w B   

    where  s  is a non-negative scaling factor chosen so that  w,  w B,   
    and  s A  do not overflow and, if possible, do not underflow, either.   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION array, dimension (LDA, 2)   
            On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm   
            is less than 1/SAFMIN.  Entries less than   
            sqrt(SAFMIN)*norm(A) are subject to being treated as zero.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= 2.   

    B       (input) DOUBLE PRECISION array, dimension (LDB, 2)   
            On entry, the 2 x 2 upper triangular matrix B.  It is   
            assumed that the one-norm of B is less than 1/SAFMIN.  The   
            diagonals should be at least sqrt(SAFMIN) times the largest   
            element of B (in absolute value); if a diagonal is smaller   
            than that, then  +/- sqrt(SAFMIN) will be used instead of   
            that diagonal.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= 2.   

    SAFMIN  (input) DOUBLE PRECISION   
            The smallest positive number s.t. 1/SAFMIN does not   
            overflow.  (This should always be DLAMCH('S') -- it is an   
            argument in order to avoid having to call DLAMCH frequently.)   

    SCALE1  (output) DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the first eigenvalue.  If   
            the eigenvalues are complex, then the eigenvalues are   
            ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the   
            exponent range of the machine), SCALE1=SCALE2, and SCALE1   
            will always be positive.  If the eigenvalues are real, then   
            the first (real) eigenvalue is  WR1 / SCALE1 , but this may   
            overflow or underflow, and in fact, SCALE1 may be zero or   
            less than the underflow threshhold if the exact eigenvalue   
            is sufficiently large.   

    SCALE2  (output) DOUBLE PRECISION   
            A scaling factor used to avoid over-/underflow in the   
            eigenvalue equation which defines the second eigenvalue.  If   
            the eigenvalues are complex, then SCALE2=SCALE1.  If the   
            eigenvalues are real, then the second (real) eigenvalue is   
            WR2 / SCALE2 , but this may overflow or underflow, and in   
            fact, SCALE2 may be zero or less than the underflow   
            threshhold if the exact eigenvalue is sufficiently large.   

    WR1     (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WR1 is SCALE1 times the   
            eigenvalue closest to the (2,2) element of A B**(-1).  If the   
            eigenvalue is complex, then WR1=WR2 is SCALE1 times the real   
            part of the eigenvalues.   

    WR2     (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WR2 is SCALE2 times the   
            other eigenvalue.  If the eigenvalue is complex, then   
            WR1=WR2 is SCALE1 times the real part of the eigenvalues.   

    WI      (output) DOUBLE PRECISION   
            If the eigenvalue is real, then WI is zero.  If the   
            eigenvalue is complex, then WI is SCALE1 times the imaginary   
            part of the eigenvalues.  WI will always be non-negative.   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    static doublereal diff, bmin, wbig, wabs, wdet, r__, binv11, binv22, 
	    discr, anorm, bnorm, bsize, shift, c1, c2, c3, c4, c5, rtmin, 
	    rtmax, wsize, s1, s2, a11, a12, a21, a22, b11, b12, b22, ascale, 
	    bscale, pp, qq, ss, wscale, safmax, wsmall, as11, as12, as22, sum,
	     abi22;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;

    /* Function Body */
    rtmin = sqrt(*safmin);
    rtmax = 1. / rtmin;
    safmax = 1. / *safmin;

/*     Scale A   

   Computing MAX */
    d__5 = (d__1 = a_ref(1, 1), abs(d__1)) + (d__2 = a_ref(2, 1), abs(d__2)), 
	    d__6 = (d__3 = a_ref(1, 2), abs(d__3)) + (d__4 = a_ref(2, 2), abs(
	    d__4)), d__5 = max(d__5,d__6);
    anorm = max(d__5,*safmin);
    ascale = 1. / anorm;
    a11 = ascale * a_ref(1, 1);
    a21 = ascale * a_ref(2, 1);
    a12 = ascale * a_ref(1, 2);
    a22 = ascale * a_ref(2, 2);

/*     Perturb B if necessary to insure non-singularity */

    b11 = b_ref(1, 1);
    b12 = b_ref(1, 2);
    b22 = b_ref(2, 2);
/* Computing MAX */
    d__1 = abs(b11), d__2 = abs(b12), d__1 = max(d__1,d__2), d__2 = abs(b22), 
	    d__1 = max(d__1,d__2);
    bmin = rtmin * max(d__1,rtmin);
    if (abs(b11) < bmin) {
	b11 = d_sign(&bmin, &b11);
    }
    if (abs(b22) < bmin) {
	b22 = d_sign(&bmin, &b22);
    }

/*     Scale B   

   Computing MAX */
    d__1 = abs(b11), d__2 = abs(b12) + abs(b22), d__1 = max(d__1,d__2);
    bnorm = max(d__1,*safmin);
/* Computing MAX */
    d__1 = abs(b11), d__2 = abs(b22);
    bsize = max(d__1,d__2);
    bscale = 1. / bsize;
    b11 *= bscale;
    b12 *= bscale;
    b22 *= bscale;

/*     Compute larger eigenvalue by method described by C. van Loan   

       ( AS is A shifted by -SHIFT*B ) */

    binv11 = 1. / b11;
    binv22 = 1. / b22;
    s1 = a11 * binv11;
    s2 = a22 * binv22;
    if (abs(s1) <= abs(s2)) {
	as12 = a12 - s1 * b12;
	as22 = a22 - s1 * b22;
	ss = a21 * (binv11 * binv22);
	abi22 = as22 * binv22 - ss * b12;
	pp = abi22 * .5;
	shift = s1;
    } else {
	as12 = a12 - s2 * b12;
	as11 = a11 - s2 * b11;
	ss = a21 * (binv11 * binv22);
	abi22 = -ss * b12;
	pp = (as11 * binv11 + abi22) * .5;
	shift = s2;
    }
    qq = ss * as12;
    if ((d__1 = pp * rtmin, abs(d__1)) >= 1.) {
/* Computing 2nd power */
	d__1 = rtmin * pp;
	discr = d__1 * d__1 + qq * *safmin;
	r__ = sqrt((abs(discr))) * rtmax;
    } else {
/* Computing 2nd power */
	d__1 = pp;
	if (d__1 * d__1 + abs(qq) <= *safmin) {
/* Computing 2nd power */
	    d__1 = rtmax * pp;
	    discr = d__1 * d__1 + qq * safmax;
	    r__ = sqrt((abs(discr))) * rtmin;
	} else {
/* Computing 2nd power */
	    d__1 = pp;
	    discr = d__1 * d__1 + qq;
	    r__ = sqrt((abs(discr)));
	}
    }

/*     Note: the test of R in the following IF is to cover the case when   
             DISCR is small and negative and is flushed to zero during   
             the calculation of R.  On machines which have a consistent   
             flush-to-zero threshhold and handle numbers above that   
             threshhold correctly, it would not be necessary. */

    if (discr >= 0. || r__ == 0.) {
	sum = pp + d_sign(&r__, &pp);
	diff = pp - d_sign(&r__, &pp);
	wbig = shift + sum;

/*        Compute smaller eigenvalue */

	wsmall = shift + diff;
/* Computing MAX */
	d__1 = abs(wsmall);
	if (abs(wbig) * .5 > max(d__1,*safmin)) {
	    wdet = (a11 * a22 - a12 * a21) * (binv11 * binv22);
	    wsmall = wdet / wbig;
	}

/*        Choose (real) eigenvalue closest to 2,2 element of A*B**(-1)   
          for WR1. */

	if (pp > abi22) {
	    *wr1 = min(wbig,wsmall);
	    *wr2 = max(wbig,wsmall);
	} else {
	    *wr1 = max(wbig,wsmall);
	    *wr2 = min(wbig,wsmall);
	}
	*wi = 0.;
    } else {

/*        Complex eigenvalues */

	*wr1 = shift + pp;
	*wr2 = *wr1;
	*wi = r__;
    }

/*     Further scaling to avoid underflow and overflow in computing   
       SCALE1 and overflow in computing w*B.   

       This scale factor (WSCALE) is bounded from above using C1 and C2,   
       and from below using C3 and C4.   
          C1 implements the condition  s A  must never overflow.   
          C2 implements the condition  w B  must never overflow.   
          C3, with C2,   
             implement the condition that s A - w B must never overflow.   
          C4 implements the condition  s    should not underflow.   
          C5 implements the condition  max(s,|w|) should be at least 2. */

    c1 = bsize * (*safmin * max(1.,ascale));
    c2 = *safmin * max(1.,bnorm);
    c3 = bsize * *safmin;
    if (ascale <= 1. && bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale / *safmin * bsize;
	c4 = min(d__1,d__2);
    } else {
	c4 = 1.;
    }
    if (ascale <= 1. || bsize <= 1.) {
/* Computing MIN */
	d__1 = 1., d__2 = ascale * bsize;
	c5 = min(d__1,d__2);
    } else {
	c5 = 1.;
    }

/*     Scale first eigenvalue */

    wabs = abs(*wr1) + abs(*wi);
/* Computing MAX   
   Computing MIN */
    d__3 = c4, d__4 = max(wabs,c5) * .5;
    d__1 = max(*safmin,c1), d__2 = (wabs * c2 + c3) * 1.0000100000000001, 
	    d__1 = max(d__1,d__2), d__2 = min(d__3,d__4);
    wsize = max(d__1,d__2);
    if (wsize != 1.) {
	wscale = 1. / wsize;
	if (wsize > 1.) {
	    *scale1 = max(ascale,bsize) * wscale * min(ascale,bsize);
	} else {
	    *scale1 = min(ascale,bsize) * wscale * max(ascale,bsize);
	}
	*wr1 *= wscale;
	if (*wi != 0.) {
	    *wi *= wscale;
	    *wr2 = *wr1;
	    *scale2 = *scale1;
	}
    } else {
	*scale1 = ascale * bsize;
	*scale2 = *scale1;
    }

/*     Scale second eigenvalue (if real) */

    if (*wi == 0.) {
/* Computing MAX   
   Computing MIN   
   Computing MAX */
	d__5 = abs(*wr2);
	d__3 = c4, d__4 = max(d__5,c5) * .5;
	d__1 = max(*safmin,c1), d__2 = (abs(*wr2) * c2 + c3) * 
		1.0000100000000001, d__1 = max(d__1,d__2), d__2 = min(d__3,
		d__4);
	wsize = max(d__1,d__2);
	if (wsize != 1.) {
	    wscale = 1. / wsize;
	    if (wsize > 1.) {
		*scale2 = max(ascale,bsize) * wscale * min(ascale,bsize);
	    } else {
		*scale2 = min(ascale,bsize) * wscale * max(ascale,bsize);
	    }
	    *wr2 *= wscale;
	} else {
	    *scale2 = ascale * bsize;
	}
    }

/*     End of DLAG2 */

    return 0;
} /* dlag2_ */

#undef b_ref
#undef a_ref



/* Subroutine */ int dlasv2_(doublereal *f, doublereal *g, doublereal *h__, 
	doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	csr, doublereal *snl, doublereal *csl)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASV2 computes the singular value decomposition of a 2-by-2   
    triangular matrix   
       [  F   G  ]   
       [  0   H  ].   
    On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the   
    smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and   
    right singular vectors for abs(SSMAX), giving the decomposition   

       [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]   
       [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].   

    Arguments   
    =========   

    F       (input) DOUBLE PRECISION   
            The (1,1) element of the 2-by-2 matrix.   

    G       (input) DOUBLE PRECISION   
            The (1,2) element of the 2-by-2 matrix.   

    H       (input) DOUBLE PRECISION   
            The (2,2) element of the 2-by-2 matrix.   

    SSMIN   (output) DOUBLE PRECISION   
            abs(SSMIN) is the smaller singular value.   

    SSMAX   (output) DOUBLE PRECISION   
            abs(SSMAX) is the larger singular value.   

    SNL     (output) DOUBLE PRECISION   
    CSL     (output) DOUBLE PRECISION   
            The vector (CSL, SNL) is a unit left singular vector for the   
            singular value abs(SSMAX).   

    SNR     (output) DOUBLE PRECISION   
    CSR     (output) DOUBLE PRECISION   
            The vector (CSR, SNR) is a unit right singular vector for the   
            singular value abs(SSMAX).   

    Further Details   
    ===============   

    Any input parameter may be aliased with any output parameter.   

    Barring over/underflow and assuming a guard digit in subtraction, all   
    output quantities are correct to within a few units in the last   
    place (ulps).   

    In IEEE arithmetic, the code works correctly if one matrix element is   
    infinite.   

    Overflow will not occur unless the largest singular value itself   
    overflows or is within a few ulps of overflow. (On machines with   
    partial overflow, like the Cray, overflow may occur if the largest   
    singular value is within a factor of 2 of overflow.)   

    Underflow is harmless if underflow is gradual. Otherwise, results   
    may correspond to a matrix modified by perturbations of size near   
    the underflow threshold.   

   ===================================================================== */
    /* Table of constant values */
    static doublereal c_b3 = 2.;
    static doublereal c_b4 = 1.;
    
    /* System generated locals */
    doublereal d__1;
    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);
    /* Local variables */
    static integer pmax;
    static doublereal temp;
    static logical swap;
    static doublereal a, d__, l, m, r__, s, t, tsign, fa, ga, ha;
    extern doublereal dlamch_(char *);
    static doublereal ft, gt, ht, mm;
    static logical gasmal;
    static doublereal tt, clt, crt, slt, srt;




    ft = *f;
    fa = abs(ft);
    ht = *h__;
    ha = abs(*h__);

/*     PMAX points to the maximum absolute element of matrix   
         PMAX = 1 if F largest in absolute values   
         PMAX = 2 if G largest in absolute values   
         PMAX = 3 if H largest in absolute values */

    pmax = 1;
    swap = ha > fa;
    if (swap) {
	pmax = 3;
	temp = ft;
	ft = ht;
	ht = temp;
	temp = fa;
	fa = ha;
	ha = temp;

/*        Now FA .ge. HA */

    }
    gt = *g;
    ga = abs(gt);
    if (ga == 0.) {

/*        Diagonal matrix */

	*ssmin = ha;
	*ssmax = fa;
	clt = 1.;
	crt = 1.;
	slt = 0.;
	srt = 0.;
    } else {
	gasmal = TRUE_;
	if (ga > fa) {
	    pmax = 2;
	    if (fa / ga < dlamch_("EPS")) {

/*              Case of very large GA */

		gasmal = FALSE_;
		*ssmax = ga;
		if (ha > 1.) {
		    *ssmin = fa / (ga / ha);
		} else {
		    *ssmin = fa / ga * ha;
		}
		clt = 1.;
		slt = ht / gt;
		srt = 1.;
		crt = ft / gt;
	    }
	}
	if (gasmal) {

/*           Normal case */

	    d__ = fa - ha;
	    if (d__ == fa) {

/*              Copes with infinite F or H */

		l = 1.;
	    } else {
		l = d__ / fa;
	    }

/*           Note that 0 .le. L .le. 1 */

	    m = gt / ft;

/*           Note that abs(M) .le. 1/macheps */

	    t = 2. - l;

/*           Note that T .ge. 1 */

	    mm = m * m;
	    tt = t * t;
	    s = sqrt(tt + mm);

/*           Note that 1 .le. S .le. 1 + 1/macheps */

	    if (l == 0.) {
		r__ = abs(m);
	    } else {
		r__ = sqrt(l * l + mm);
	    }

/*           Note that 0 .le. R .le. 1 + 1/macheps */

	    a = (s + r__) * .5;

/*           Note that 1 .le. A .le. 1 + abs(M) */

	    *ssmin = ha / a;
	    *ssmax = fa * a;
	    if (mm == 0.) {

/*              Note that M is very tiny */

		if (l == 0.) {
		    t = d_sign(&c_b3, &ft) * d_sign(&c_b4, &gt);
		} else {
		    t = gt / d_sign(&d__, &ft) + m / t;
		}
	    } else {
		t = (m / (s + t) + m / (r__ + l)) * (a + 1.);
	    }
	    l = sqrt(t * t + 4.);
	    crt = 2. / l;
	    srt = t / l;
	    clt = (crt + srt * m) / a;
	    slt = ht / ft * srt / a;
	}
    }
    if (swap) {
	*csl = srt;
	*snl = crt;
	*csr = slt;
	*snr = clt;
    } else {
	*csl = clt;
	*snl = slt;
	*csr = crt;
	*snr = srt;
    }

/*     Correct signs of SSMAX and SSMIN */

    if (pmax == 1) {
	tsign = d_sign(&c_b4, csr) * d_sign(&c_b4, csl) * d_sign(&c_b4, f);
    }
    if (pmax == 2) {
	tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, csl) * d_sign(&c_b4, g);
    }
    if (pmax == 3) {
	tsign = d_sign(&c_b4, snr) * d_sign(&c_b4, snl) * d_sign(&c_b4, h__);
    }
    *ssmax = d_sign(ssmax, &tsign);
    d__1 = tsign * d_sign(&c_b4, f) * d_sign(&c_b4, h__);
    *ssmin = d_sign(ssmin, &d__1);
    return 0;

/*     End of DLASV2 */

} /* dlasv2_ */


doublereal dlapy3_(doublereal *x, doublereal *y, doublereal *z__)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause   
    unnecessary overflow.   

    Arguments   
    =========   

    X       (input) DOUBLE PRECISION   
    Y       (input) DOUBLE PRECISION   
    Z       (input) DOUBLE PRECISION   
            X, Y and Z specify the values x, y and z.   

    ===================================================================== */
    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublereal xabs, yabs, zabs, w;



    xabs = abs(*x);
    yabs = abs(*y);
    zabs = abs(*z__);
/* Computing MAX */
    d__1 = max(xabs,yabs);
    w = max(d__1,zabs);
    if (w == 0.) {
	ret_val = 0.;
    } else {
/* Computing 2nd power */
	d__1 = xabs / w;
/* Computing 2nd power */
	d__2 = yabs / w;
/* Computing 2nd power */
	d__3 = zabs / w;
	ret_val = w * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    return ret_val;

/*     End of DLAPY3 */

} /* dlapy3_ */


doublereal dlapy2_(doublereal *x, doublereal *y)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary   
    overflow.   

    Arguments   
    =========   

    X       (input) DOUBLE PRECISION   
    Y       (input) DOUBLE PRECISION   
            X and Y specify the values x and y.   

    ===================================================================== */
    /* System generated locals */
    doublereal ret_val, d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublereal xabs, yabs, w, z__;



    xabs = abs(*x);
    yabs = abs(*y);
    w = max(xabs,yabs);
    z__ = min(xabs,yabs);
    if (z__ == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z__ / w;
	ret_val = w * sqrt(d__1 * d__1 + 1.);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */


/* Subroutine */ int dlarfg_(integer *n, doublereal *alpha, doublereal *x, 
	integer *incx, doublereal *tau)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DLARFG generates a real elementary reflector H of order n, such   
    that   

          H * ( alpha ) = ( beta ),   H' * H = I.   
              (   x   )   (   0  )   

    where alpha and beta are scalars, and x is an (n-1)-element real   
    vector. H is represented in the form   

          H = I - tau * ( 1 ) * ( 1 v' ) ,   
                        ( v )   

    where tau is a real scalar and v is a real (n-1)-element   
    vector.   

    If the elements of x are all zero, then tau = 0 and H is taken to be   
    the unit matrix.   

    Otherwise  1 <= tau <= 2.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the elementary reflector.   

    ALPHA   (input/output) DOUBLE PRECISION   
            On entry, the value alpha.   
            On exit, it is overwritten with the value beta.   

    X       (input/output) DOUBLE PRECISION array, dimension   
                           (1+(N-2)*abs(INCX))   
            On entry, the vector x.   
            On exit, it is overwritten with the vector v.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    TAU     (output) DOUBLE PRECISION   
            The value tau.   

    =====================================================================   


       Parameter adjustments */
    /* System generated locals */
    integer i__1;
    doublereal d__1;
    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);
    /* Local variables */
    static doublereal beta;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer j;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal xnorm;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *);
    static doublereal safmin, rsafmn;
    static integer knt;

    --x;

    /* Function Body */
    if (*n <= 1) {
	*tau = 0.;
	return 0;
    }

    i__1 = *n - 1;
    xnorm = dnrm2_(&i__1, &x[1], incx);

    if (xnorm == 0.) {

/*        H  =  I */

	*tau = 0.;
    } else {

/*        general case */

	d__1 = dlapy2_(alpha, &xnorm);
	beta = -d_sign(&d__1, alpha);
	safmin = dlamch_("S") / dlamch_("E");
	if (abs(beta) < safmin) {

/*           XNORM, BETA may be inaccurate; scale X and recompute them */

	    rsafmn = 1. / safmin;
	    knt = 0;
L10:
	    ++knt;
	    i__1 = *n - 1;
	    dscal_(&i__1, &rsafmn, &x[1], incx);
	    beta *= rsafmn;
	    *alpha *= rsafmn;
	    if (abs(beta) < safmin) {
		goto L10;
	    }

/*           New BETA is at most 1, at least SAFMIN */

	    i__1 = *n - 1;
	    xnorm = dnrm2_(&i__1, &x[1], incx);
	    d__1 = dlapy2_(alpha, &xnorm);
	    beta = -d_sign(&d__1, alpha);
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &x[1], incx);

/*           If ALPHA is subnormal, it may lose relative accuracy */

	    *alpha = beta;
	    i__1 = knt;
	    for (j = 1; j <= i__1; ++j) {
		*alpha *= safmin;
/* L20: */
	    }
	} else {
	    *tau = (beta - *alpha) / beta;
	    i__1 = *n - 1;
	    d__1 = 1. / (*alpha - beta);
	    dscal_(&i__1, &d__1, &x[1], incx);
	    *alpha = beta;
	}
    }

    return 0;

/*     End of DLARFG */

} /* dlarfg_ */


/* Subroutine */ int dlabad_(doublereal *small, doublereal *large)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLABAD takes as input the values computed by DLAMCH for underflow and   
    overflow, and returns the square root of each of these values if the   
    log of LARGE is sufficiently large.  This subroutine is intended to   
    identify machines with a large exponent range, such as the Crays, and   
    redefine the underflow and overflow limits to be the square roots of   
    the values computed by DLAMCH.  This subroutine is needed because   
    DLAMCH does not compensate for poor arithmetic in the upper half of   
    the exponent range, as is found on a Cray.   

    Arguments   
    =========   

    SMALL   (input/output) DOUBLE PRECISION   
            On entry, the underflow threshold as computed by DLAMCH.   
            On exit, if LOG10(LARGE) is sufficiently large, the square   
            root of SMALL, otherwise unchanged.   

    LARGE   (input/output) DOUBLE PRECISION   
            On entry, the overflow threshold as computed by DLAMCH.   
            On exit, if LOG10(LARGE) is sufficiently large, the square   
            root of LARGE, otherwise unchanged.   

    =====================================================================   


       If it looks like we're on a Cray, take the square root of   
       SMALL and LARGE to avoid overflow and underflow problems. */
    /* Builtin functions */
    double d_lg10(doublereal *), sqrt(doublereal);


    if (d_lg10(large) > 2e3) {
	*small = sqrt(*small);
	*large = sqrt(*large);
    }

    return 0;

/*     End of DLABAD */

} /* dlabad_ */


/* Subroutine */ int dlaln2_(logical *ltrans, integer *na, integer *nw, 
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda, 
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb, 
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx, 
	doublereal *scale, doublereal *xnorm, integer *info)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLALN2 solves a system of the form  (ca A - w D ) X = s B   
    or (ca A' - w D) X = s B   with possible scaling ("s") and   
    perturbation of A.  (A' means A-transpose.)   

    A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA   
    real diagonal matrix, w is a real or complex value, and X and B are   
    NA x 1 matrices -- real if w is real, complex if w is complex.  NA   
    may be 1 or 2.   

    If w is complex, X and B are represented as NA x 2 matrices,   
    the first column of each being the real part and the second   
    being the imaginary part.   

    "s" is a scaling factor (.LE. 1), computed by DLALN2, which is   
    so chosen that X can be computed without overflow.  X is further   
    scaled if necessary to assure that norm(ca A - w D)*norm(X) is less   
    than overflow.   

    If both singular values of (ca A - w D) are less than SMIN,   
    SMIN*identity will be used instead of (ca A - w D).  If only one   
    singular value is less than SMIN, one element of (ca A - w D) will be   
    perturbed enough to make the smallest singular value roughly SMIN.   
    If both singular values are at least SMIN, (ca A - w D) will not be   
    perturbed.  In any case, the perturbation will be at most some small   
    multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values   
    are computed by infinity-norm approximations, and thus will only be   
    correct to a factor of 2 or so.   

    Note: all input quantities are assumed to be smaller than overflow   
    by a reasonable factor.  (See BIGNUM.)   

    Arguments   
    ==========   

    LTRANS  (input) LOGICAL   
            =.TRUE.:  A-transpose will be used.   
            =.FALSE.: A will be used (not transposed.)   

    NA      (input) INTEGER   
            The size of the matrix A.  It may (only) be 1 or 2.   

    NW      (input) INTEGER   
            1 if "w" is real, 2 if "w" is complex.  It may only be 1   
            or 2.   

    SMIN    (input) DOUBLE PRECISION   
            The desired lower bound on the singular values of A.  This   
            should be a safe distance away from underflow or overflow,   
            say, between (underflow/machine precision) and  (machine   
            precision * overflow ).  (See BIGNUM and ULP.)   

    CA      (input) DOUBLE PRECISION   
            The coefficient c, which A is multiplied by.   

    A       (input) DOUBLE PRECISION array, dimension (LDA,NA)   
            The NA x NA matrix A.   

    LDA     (input) INTEGER   
            The leading dimension of A.  It must be at least NA.   

    D1      (input) DOUBLE PRECISION   
            The 1,1 element in the diagonal matrix D.   

    D2      (input) DOUBLE PRECISION   
            The 2,2 element in the diagonal matrix D.  Not used if NW=1.   

    B       (input) DOUBLE PRECISION array, dimension (LDB,NW)   
            The NA x NW matrix B (right-hand side).  If NW=2 ("w" is   
            complex), column 1 contains the real part of B and column 2   
            contains the imaginary part.   

    LDB     (input) INTEGER   
            The leading dimension of B.  It must be at least NA.   

    WR      (input) DOUBLE PRECISION   
            The real part of the scalar "w".   

    WI      (input) DOUBLE PRECISION   
            The imaginary part of the scalar "w".  Not used if NW=1.   

    X       (output) DOUBLE PRECISION array, dimension (LDX,NW)   
            The NA x NW matrix X (unknowns), as computed by DLALN2.   
            If NW=2 ("w" is complex), on exit, column 1 will contain   
            the real part of X and column 2 will contain the imaginary   
            part.   

    LDX     (input) INTEGER   
            The leading dimension of X.  It must be at least NA.   

    SCALE   (output) DOUBLE PRECISION   
            The scale factor that B must be multiplied by to insure   
            that overflow does not occur when computing X.  Thus,   
            (ca A - w D) X  will be SCALE*B, not B (ignoring   
            perturbations of A.)  It will be at most 1.   

    XNORM   (output) DOUBLE PRECISION   
            The infinity-norm of X, when X is regarded as an NA x NW   
            real matrix.   

    INFO    (output) INTEGER   
            An error flag.  It will be set to zero if no error occurs,   
            a negative number if an argument is in error, or a positive   
            number if  ca A - w D  had to be perturbed.   
            The possible values are:   
            = 0: No error occurred, and (ca A - w D) did not have to be   
                   perturbed.   
            = 1: (ca A - w D) had to be perturbed to make its smallest   
                 (or only) singular value greater than SMIN.   
            NOTE: In the interests of speed, this routine does not   
                  check the inputs for errors.   

   =====================================================================   

       Parameter adjustments */
    /* Initialized data */
    static logical zswap[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
    static logical rswap[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };
    static integer ipivot[16]	/* was [4][4] */ = { 1,2,3,4,2,1,4,3,3,4,1,2,
	    4,3,2,1 };
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;
    static doublereal equiv_0[4], equiv_1[4];
    /* Local variables */
    static doublereal bbnd, cmax, ui11r, ui12s, temp, ur11r, ur12s;
    static integer j;
    static doublereal u22abs;
    static integer icmax;
    static doublereal bnorm, cnorm, smini;
#define ci (equiv_0)
#define cr (equiv_1)
    extern doublereal dlamch_(char *);
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum, bi1, bi2, br1, br2, smlnum, xi1, xi2, xr1, xr2, 
	    ci21, ci22, cr21, cr22, li21, csi, ui11, lr21, ui12, ui22;
#define civ (equiv_0)
    static doublereal csr, ur11, ur12, ur22;
#define crv (equiv_1)
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define x_ref(a_1,a_2) x[(a_2)*x_dim1 + a_1]
#define ci_ref(a_1,a_2) ci[(a_2)*2 + a_1 - 3]
#define cr_ref(a_1,a_2) cr[(a_2)*2 + a_1 - 3]
#define ipivot_ref(a_1,a_2) ipivot[(a_2)*4 + a_1 - 5]

    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1 * 1;
    x -= x_offset;

    /* Function Body   

       Compute BIGNUM */

    smlnum = 2. * dlamch_("Safe minimum");
    bignum = 1. / smlnum;
    smini = max(*smin,smlnum);

/*     Don't check for input errors */

    *info = 0;

/*     Standard Initializations */

    *scale = 1.;

    if (*na == 1) {

/*        1 x 1  (i.e., scalar) system   C X = B */

	if (*nw == 1) {

/*           Real 1x1 system.   

             C = ca A - w D */

	    csr = *ca * a_ref(1, 1) - *wr * *d1;
	    cnorm = abs(csr);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = b_ref(1, 1), abs(d__1));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    x_ref(1, 1) = b_ref(1, 1) * *scale / csr;
	    *xnorm = (d__1 = x_ref(1, 1), abs(d__1));
	} else {

/*           Complex 1x1 system (w is complex)   

             C = ca A - w D */

	    csr = *ca * a_ref(1, 1) - *wr * *d1;
	    csi = -(*wi) * *d1;
	    cnorm = abs(csr) + abs(csi);

/*           If | C | < SMINI, use C = SMINI */

	    if (cnorm < smini) {
		csr = smini;
		csi = 0.;
		cnorm = smini;
		*info = 1;
	    }

/*           Check scaling for  X = B / C */

	    bnorm = (d__1 = b_ref(1, 1), abs(d__1)) + (d__2 = b_ref(1, 2), 
		    abs(d__2));
	    if (cnorm < 1. && bnorm > 1.) {
		if (bnorm > bignum * cnorm) {
		    *scale = 1. / bnorm;
		}
	    }

/*           Compute X */

	    d__1 = *scale * b_ref(1, 1);
	    d__2 = *scale * b_ref(1, 2);
	    dladiv_(&d__1, &d__2, &csr, &csi, &x_ref(1, 1), &x_ref(1, 2));
	    *xnorm = (d__1 = x_ref(1, 1), abs(d__1)) + (d__2 = x_ref(1, 2), 
		    abs(d__2));
	}

    } else {

/*        2x2 System   

          Compute the real part of  C = ca A - w D  (or  ca A' - w D ) */

	cr_ref(1, 1) = *ca * a_ref(1, 1) - *wr * *d1;
	cr_ref(2, 2) = *ca * a_ref(2, 2) - *wr * *d2;
	if (*ltrans) {
	    cr_ref(1, 2) = *ca * a_ref(2, 1);
	    cr_ref(2, 1) = *ca * a_ref(1, 2);
	} else {
	    cr_ref(2, 1) = *ca * a_ref(2, 1);
	    cr_ref(1, 2) = *ca * a_ref(1, 2);
	}

	if (*nw == 1) {

/*           Real 2x2 system  (w is real)   

             Find the largest element in C */

	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], abs(d__1)) > cmax) {
		    cmax = (d__1 = crv[j - 1], abs(d__1));
		    icmax = j;
		}
/* L10: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__3 = (d__1 = b_ref(1, 1), abs(d__1)), d__4 = (d__2 = b_ref(
			2, 1), abs(d__2));
		bnorm = max(d__3,d__4);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		x_ref(1, 1) = temp * b_ref(1, 1);
		x_ref(2, 1) = temp * b_ref(2, 1);
		*xnorm = temp * bnorm;
		*info = 1;
		return 0;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    cr21 = crv[ipivot_ref(2, icmax) - 1];
	    ur12 = crv[ipivot_ref(3, icmax) - 1];
	    cr22 = crv[ipivot_ref(4, icmax) - 1];
	    ur11r = 1. / ur11;
	    lr21 = ur11r * cr21;
	    ur22 = cr22 - ur12 * lr21;

/*           If smaller pivot < SMINI, use SMINI */

	    if (abs(ur22) < smini) {
		ur22 = smini;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br1 = b_ref(2, 1);
		br2 = b_ref(1, 1);
	    } else {
		br1 = b_ref(1, 1);
		br2 = b_ref(2, 1);
	    }
	    br2 -= lr21 * br1;
/* Computing MAX */
	    d__2 = (d__1 = br1 * (ur22 * ur11r), abs(d__1)), d__3 = abs(br2);
	    bbnd = max(d__2,d__3);
	    if (bbnd > 1. && abs(ur22) < 1.) {
		if (bbnd >= bignum * abs(ur22)) {
		    *scale = 1. / bbnd;
		}
	    }

	    xr2 = br2 * *scale / ur22;
	    xr1 = *scale * br1 * ur11r - xr2 * (ur11r * ur12);
	    if (zswap[icmax - 1]) {
		x_ref(1, 1) = xr2;
		x_ref(2, 1) = xr1;
	    } else {
		x_ref(1, 1) = xr1;
		x_ref(2, 1) = xr2;
	    }
/* Computing MAX */
	    d__1 = abs(xr1), d__2 = abs(xr2);
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x_ref(1, 1) = temp * x_ref(1, 1);
		    x_ref(2, 1) = temp * x_ref(2, 1);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	} else {

/*           Complex 2x2 system  (w is complex)   

             Find the largest element in C */

	    ci_ref(1, 1) = -(*wi) * *d1;
	    ci_ref(2, 1) = 0.;
	    ci_ref(1, 2) = 0.;
	    ci_ref(2, 2) = -(*wi) * *d2;
	    cmax = 0.;
	    icmax = 0;

	    for (j = 1; j <= 4; ++j) {
		if ((d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1], abs(
			d__2)) > cmax) {
		    cmax = (d__1 = crv[j - 1], abs(d__1)) + (d__2 = civ[j - 1]
			    , abs(d__2));
		    icmax = j;
		}
/* L20: */
	    }

/*           If norm(C) < SMINI, use SMINI*identity. */

	    if (cmax < smini) {
/* Computing MAX */
		d__5 = (d__1 = b_ref(1, 1), abs(d__1)) + (d__2 = b_ref(1, 2), 
			abs(d__2)), d__6 = (d__3 = b_ref(2, 1), abs(d__3)) + (
			d__4 = b_ref(2, 2), abs(d__4));
		bnorm = max(d__5,d__6);
		if (smini < 1. && bnorm > 1.) {
		    if (bnorm > bignum * smini) {
			*scale = 1. / bnorm;
		    }
		}
		temp = *scale / smini;
		x_ref(1, 1) = temp * b_ref(1, 1);
		x_ref(2, 1) = temp * b_ref(2, 1);
		x_ref(1, 2) = temp * b_ref(1, 2);
		x_ref(2, 2) = temp * b_ref(2, 2);
		*xnorm = temp * bnorm;
		*info = 1;
		return 0;
	    }

/*           Gaussian elimination with complete pivoting. */

	    ur11 = crv[icmax - 1];
	    ui11 = civ[icmax - 1];
	    cr21 = crv[ipivot_ref(2, icmax) - 1];
	    ci21 = civ[ipivot_ref(2, icmax) - 1];
	    ur12 = crv[ipivot_ref(3, icmax) - 1];
	    ui12 = civ[ipivot_ref(3, icmax) - 1];
	    cr22 = crv[ipivot_ref(4, icmax) - 1];
	    ci22 = civ[ipivot_ref(4, icmax) - 1];
	    if (icmax == 1 || icmax == 4) {

/*              Code when off-diagonals of pivoted C are real */

		if (abs(ur11) > abs(ui11)) {
		    temp = ui11 / ur11;
/* Computing 2nd power */
		    d__1 = temp;
		    ur11r = 1. / (ur11 * (d__1 * d__1 + 1.));
		    ui11r = -temp * ur11r;
		} else {
		    temp = ur11 / ui11;
/* Computing 2nd power */
		    d__1 = temp;
		    ui11r = -1. / (ui11 * (d__1 * d__1 + 1.));
		    ur11r = -temp * ui11r;
		}
		lr21 = cr21 * ur11r;
		li21 = cr21 * ui11r;
		ur12s = ur12 * ur11r;
		ui12s = ur12 * ui11r;
		ur22 = cr22 - ur12 * lr21;
		ui22 = ci22 - ur12 * li21;
	    } else {

/*              Code when diagonals of pivoted C are real */

		ur11r = 1. / ur11;
		ui11r = 0.;
		lr21 = cr21 * ur11r;
		li21 = ci21 * ur11r;
		ur12s = ur12 * ur11r;
		ui12s = ui12 * ur11r;
		ur22 = cr22 - ur12 * lr21 + ui12 * li21;
		ui22 = -ur12 * li21 - ui12 * lr21;
	    }
	    u22abs = abs(ur22) + abs(ui22);

/*           If smaller pivot < SMINI, use SMINI */

	    if (u22abs < smini) {
		ur22 = smini;
		ui22 = 0.;
		*info = 1;
	    }
	    if (rswap[icmax - 1]) {
		br2 = b_ref(1, 1);
		br1 = b_ref(2, 1);
		bi2 = b_ref(1, 2);
		bi1 = b_ref(2, 2);
	    } else {
		br1 = b_ref(1, 1);
		br2 = b_ref(2, 1);
		bi1 = b_ref(1, 2);
		bi2 = b_ref(2, 2);
	    }
	    br2 = br2 - lr21 * br1 + li21 * bi1;
	    bi2 = bi2 - li21 * br1 - lr21 * bi1;
/* Computing MAX */
	    d__1 = (abs(br1) + abs(bi1)) * (u22abs * (abs(ur11r) + abs(ui11r))
		    ), d__2 = abs(br2) + abs(bi2);
	    bbnd = max(d__1,d__2);
	    if (bbnd > 1. && u22abs < 1.) {
		if (bbnd >= bignum * u22abs) {
		    *scale = 1. / bbnd;
		    br1 = *scale * br1;
		    bi1 = *scale * bi1;
		    br2 = *scale * br2;
		    bi2 = *scale * bi2;
		}
	    }

	    dladiv_(&br2, &bi2, &ur22, &ui22, &xr2, &xi2);
	    xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
	    xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
	    if (zswap[icmax - 1]) {
		x_ref(1, 1) = xr2;
		x_ref(2, 1) = xr1;
		x_ref(1, 2) = xi2;
		x_ref(2, 2) = xi1;
	    } else {
		x_ref(1, 1) = xr1;
		x_ref(2, 1) = xr2;
		x_ref(1, 2) = xi1;
		x_ref(2, 2) = xi2;
	    }
/* Computing MAX */
	    d__1 = abs(xr1) + abs(xi1), d__2 = abs(xr2) + abs(xi2);
	    *xnorm = max(d__1,d__2);

/*           Further scaling if  norm(A) norm(X) > overflow */

	    if (*xnorm > 1. && cmax > 1.) {
		if (*xnorm > bignum / cmax) {
		    temp = cmax / bignum;
		    x_ref(1, 1) = temp * x_ref(1, 1);
		    x_ref(2, 1) = temp * x_ref(2, 1);
		    x_ref(1, 2) = temp * x_ref(1, 2);
		    x_ref(2, 2) = temp * x_ref(2, 2);
		    *xnorm = temp * *xnorm;
		    *scale = temp * *scale;
		}
	    }
	}
    }

    return 0;

/*     End of DLALN2 */

} /* dlaln2_ */

#undef ipivot_ref
#undef cr_ref
#undef ci_ref
#undef x_ref
#undef b_ref
#undef a_ref
#undef crv
#undef civ
#undef cr
#undef ci



/* Subroutine */ int dlarf_(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc, 
	doublereal *work)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARF applies a real elementary reflector H to a real m by n matrix   
    C, from either the left or the right. H is represented in the form   

          H = I - tau * v * v'   

    where tau is a real scalar and v is a real vector.   

    If tau = 0, then H is taken to be the unit matrix.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': form  H * C   
            = 'R': form  C * H   

    M       (input) INTEGER   
            The number of rows of the matrix C.   

    N       (input) INTEGER   
            The number of columns of the matrix C.   

    V       (input) DOUBLE PRECISION array, dimension   
                       (1 + (M-1)*abs(INCV)) if SIDE = 'L'   
                    or (1 + (N-1)*abs(INCV)) if SIDE = 'R'   
            The vector v in the representation of H. V is not used if   
            TAU = 0.   

    INCV    (input) INTEGER   
            The increment between elements of v. INCV <> 0.   

    TAU     (input) DOUBLE PRECISION   
            The value tau in the representation of H.   

    C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by the matrix H * C if SIDE = 'L',   
            or C * H if SIDE = 'R'.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= max(1,M).   

    WORK    (workspace) DOUBLE PRECISION array, dimension   
                           (N) if SIDE = 'L'   
                        or (M) if SIDE = 'R'   

    =====================================================================   


       Parameter adjustments */
    /* Table of constant values */
    static doublereal c_b4 = 1.;
    static doublereal c_b5 = 0.;
    static integer c__1 = 1;
    
    /* System generated locals */
    integer c_dim1, c_offset;
    doublereal d__1;
    /* Local variables */
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *);


    --v;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    --work;

    /* Function Body */
    if (lsame_(side, "L")) {

/*        Form  H * C */

	if (*tau != 0.) {

/*           w := C' * v */

	    dgemv_("Transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], incv,
		     &c_b5, &work[1], &c__1);

/*           C := C - v * w' */

	    d__1 = -(*tau);
	    dger_(m, n, &d__1, &v[1], incv, &work[1], &c__1, &c__[c_offset], 
		    ldc);
	}
    } else {

/*        Form  C * H */

	if (*tau != 0.) {

/*           w := C * v */

	    dgemv_("No transpose", m, n, &c_b4, &c__[c_offset], ldc, &v[1], 
		    incv, &c_b5, &work[1], &c__1);

/*           C := C - w * v' */

	    d__1 = -(*tau);
	    dger_(m, n, &d__1, &work[1], &c__1, &v[1], incv, &c__[c_offset], 
		    ldc);
	}
    }
    return 0;

/*     End of DLARF */

} /* dlarf_ */


/* Subroutine */ int dladiv_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *p, doublereal *q)
{
/*  -- LAPACK auxiliary routine (version 3.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLADIV performs complex division in  real arithmetic   

                          a + i*b   
               p + i*q = ---------   
                          c + i*d   

    The algorithm is due to Robert L. Smith and can be found   
    in D. Knuth, The art of Computer Programming, Vol.2, p.195   

    Arguments   
    =========   

    A       (input) DOUBLE PRECISION   
    B       (input) DOUBLE PRECISION   
    C       (input) DOUBLE PRECISION   
    D       (input) DOUBLE PRECISION   
            The scalars a, b, c, and d in the above expression.   

    P       (output) DOUBLE PRECISION   
    Q       (output) DOUBLE PRECISION   
            The scalars p and q in the above expression.   

    ===================================================================== */
    static doublereal e, f;



    if (abs(*d__) < abs(*c__)) {
	e = *d__ / *c__;
	f = *c__ + *d__ * e;
	*p = (*a + *b * e) / f;
	*q = (*b - *a * e) / f;
    } else {
	e = *c__ / *d__;
	f = *d__ + *c__ * e;
	*p = (*b + *a * e) / f;
	*q = (-(*a) + *b * e) / f;
    }

    return 0;

/*     End of DLADIV */

} /* dladiv_ */


/* Subroutine */ int dswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;
/*     interchanges two vectors.   
       uses unrolled loops for increments equal one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dy;
    --dx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*       code for unequal increments or equal increments not equal   
           to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*       code for both increments equal to 1   
         clean-up loop */
L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
	dtemp = dx[i__ + 1];
	dx[i__ + 1] = dy[i__ + 1];
	dy[i__ + 1] = dtemp;
	dtemp = dx[i__ + 2];
	dx[i__ + 2] = dy[i__ + 2];
	dy[i__ + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* dswap_ */


/* Subroutine */ int dger_(integer *m, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *y, integer *incy, 
	doublereal *a, integer *lda)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, ix, jy, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    DGER   performs the rank 1 operation   
       A := alpha*x*y' + A,   
    where alpha is a scalar, x is an m element vector, y is an n element   
    vector and A is an m by n matrix.   
    Parameters   
    ==========   
    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( m - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the m   
             element vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Y      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients. On exit, A is   
             overwritten by the updated matrix.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    --x;
    --y;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    /* Function Body */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < max(1,*m)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DGER  ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0.) {
	return 0;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */
    if (*incy > 0) {
	jy = 1;
    } else {
	jy = 1 - (*n - 1) * *incy;
    }
    if (*incx == 1) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a_ref(i__, j) = a_ref(i__, j) + x[i__] * temp;
/* L10: */
		}
	    }
	    jy += *incy;
/* L20: */
	}
    } else {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*m - 1) * *incx;
	}
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    if (y[jy] != 0.) {
		temp = *alpha * y[jy];
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    a_ref(i__, j) = a_ref(i__, j) + x[ix] * temp;
		    ix += *incx;
/* L30: */
		}
	    }
	    jy += *incy;
/* L40: */
	}
    }
    return 0;
/*     End of DGER  . */
} /* dger_ */
#undef a_ref


/* Subroutine */ int dgemv_(char *trans, integer *m, integer *n, doublereal *
	alpha, doublereal *a, integer *lda, doublereal *x, integer *incx, 
	doublereal *beta, doublereal *y, integer *incy)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer lenx, leny, i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    DGEMV  performs one of the matrix-vector operations   
       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n matrix.   
    Parameters   
    ==========   
    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   
                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   
                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading m by n part of the array A must   
             contain the matrix of coefficients.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             max( 1, m ).   
             Unchanged on exit.   
    X      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    BETA   - DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   
    Y      - DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the   
             updated vector y.   
    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    --y;
    /* Function Body */
    info = 0;
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! lsame_(trans, "C")
	    ) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*lda < max(1,*m)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DGEMV ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
	return 0;
    }
/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set   
       up the start points in  X  and  Y. */
    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A.   
       First form  y := beta*y. */
    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[i__] = *beta * y[i__];
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    y[iy] = *beta * y[iy];
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return 0;
    }
    if (lsame_(trans, "N")) {
/*        Form  y := alpha*A*x + y. */
	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a_ref(i__, j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != 0.) {
		    temp = *alpha * x[jx];
		    iy = ky;
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a_ref(i__, j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
/* L80: */
	    }
	}
    } else {
/*        Form  y := alpha*A'*x + y. */
	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a_ref(i__, j) * x[i__];
/* L90: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		temp = 0.;
		ix = kx;
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a_ref(i__, j) * x[ix];
		    ix += *incx;
/* L110: */
		}
		y[jy] += *alpha * temp;
		jy += *incy;
/* L120: */
	    }
	}
    }
    return 0;
/*     End of DGEMV . */
} /* dgemv_ */
#undef a_ref


/* Subroutine */ int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    /* Local variables */
    static integer i__, m, nincx, mp1;
/*     scales a vector by a constant.   
       uses unrolled loops for increment equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dx;
    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }
/*        code for increment not equal to 1 */
    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;
/*        code for increment equal to 1   
          clean-up loop */
L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */


doublereal dnrm2_(integer *n, doublereal *x, integer *incx)
{
/*        The following loop is equivalent to this call to the LAPACK   
          auxiliary routine:   
          CALL DLASSQ( N, X, INCX, SCALE, SSQ ) */
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;
    /* Builtin functions */
    double sqrt(doublereal);
    /* Local variables */
    static doublereal norm, scale, absxi;
    static integer ix;
    static doublereal ssq;
/*  DNRM2 returns the euclidean norm of a vector via the function   
    name, so that   
       DNRM2 := sqrt( x'*x )   
    -- This version written on 25-October-1982.   
       Modified on 14-October-1993 to inline the call to DLASSQ.   
       Sven Hammarling, Nag Ltd.   
       Parameter adjustments */
    --x;
    /* Function Body */
    if (*n < 1 || *incx < 1) {
	norm = 0.;
    } else if (*n == 1) {
	norm = abs(x[1]);
    } else {
	scale = 0.;
	ssq = 1.;


	i__1 = (*n - 1) * *incx + 1;
	i__2 = *incx;
	for (ix = 1; i__2 < 0 ? ix >= i__1 : ix <= i__1; ix += i__2) {
	    if (x[ix] != 0.) {
		absxi = (d__1 = x[ix], abs(d__1));
		if (scale < absxi) {
/* Computing 2nd power */
		    d__1 = scale / absxi;
		    ssq = ssq * (d__1 * d__1) + 1.;
		    scale = absxi;
		} else {
/* Computing 2nd power */
		    d__1 = absxi / scale;
		    ssq += d__1 * d__1;
		}
	    }
/* L10: */
	}
	norm = scale * sqrt(ssq);
    }

    ret_val = norm;
    return ret_val;

/*     End of DNRM2. */

} /* dnrm2_ */


/* Subroutine */ int dtrmm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3;
    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j, k;
    static logical lside;
    extern logical lsame_(char *, char *);
    static integer nrowa;
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
/*  Purpose   
    =======   
    DTRMM  performs one of the matrix-matrix operations   
       B := alpha*op( A )*B,   or   B := alpha*B*op( A ),   
    where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or   
    non-unit,  upper or lower triangular matrix  and  op( A )  is one  of   
       op( A ) = A   or   op( A ) = A'.   
    Parameters   
    ==========   
    SIDE   - CHARACTER*1.   
             On entry,  SIDE specifies whether  op( A ) multiplies B from   
             the left or right as follows:   
                SIDE = 'L' or 'l'   B := alpha*op( A )*B.   
                SIDE = 'R' or 'r'   B := alpha*B*op( A ).   
             Unchanged on exit.   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix A is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n'   op( A ) = A.   
                TRANSA = 'T' or 't'   op( A ) = A'.   
                TRANSA = 'C' or 'c'   op( A ) = A'.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit triangular   
             as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry, M specifies the number of rows of B. M must be at   
             least zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the number of columns of B.  N must be   
             at least zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry,  ALPHA specifies the scalar  alpha. When  alpha is   
             zero then  A is not referenced and  B need not be set before   
             entry.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m   
             when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.   
             Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k   
             upper triangular part of the array  A must contain the upper   
             triangular matrix  and the strictly lower triangular part of   
             A is not referenced.   
             Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k   
             lower triangular part of the array  A must contain the lower   
             triangular matrix  and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u',  the diagonal elements of   
             A  are not referenced either,  but are assumed to be  unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program.  When  SIDE = 'L' or 'l'  then   
             LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'   
             then LDA must be at least max( 1, n ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).   
             Before entry,  the leading  m by n part of the array  B must   
             contain the matrix  B,  and  on exit  is overwritten  by the   
             transformed matrix.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in  the  calling  (sub)  program.   LDB  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    /* Function Body */
    lside = lsame_(side, "L");
    if (lside) {
	nrowa = *m;
    } else {
	nrowa = *n;
    }
    nounit = lsame_(diag, "N");
    upper = lsame_(uplo, "U");
    info = 0;
    if (! lside && ! lsame_(side, "R")) {
	info = 1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	info = 2;
    } else if (! lsame_(transa, "N") && ! lsame_(transa,
	     "T") && ! lsame_(transa, "C")) {
	info = 3;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, 
	    "N")) {
	info = 4;
    } else if (*m < 0) {
	info = 5;
    } else if (*n < 0) {
	info = 6;
    } else if (*lda < max(1,nrowa)) {
	info = 9;
    } else if (*ldb < max(1,*m)) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DTRMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
/*     And when  alpha.eq.zero. */
    if (*alpha == 0.) {
	i__1 = *n;
	for (j = 1; j <= i__1; ++j) {
	    i__2 = *m;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		b_ref(i__, j) = 0.;
/* L10: */
	    }
/* L20: */
	}
	return 0;
    }
/*     Start the operations. */
    if (lside) {
	if (lsame_(transa, "N")) {
/*           Form  B := alpha*A*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (k = 1; k <= i__2; ++k) {
			if (b_ref(k, j) != 0.) {
			    temp = *alpha * b_ref(k, j);
			    i__3 = k - 1;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L30: */
			    }
			    if (nounit) {
				temp *= a_ref(k, k);
			    }
			    b_ref(k, j) = temp;
			}
/* L40: */
		    }
/* L50: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (k = *m; k >= 1; --k) {
			if (b_ref(k, j) != 0.) {
			    temp = *alpha * b_ref(k, j);
			    b_ref(k, j) = temp;
			    if (nounit) {
				b_ref(k, j) = b_ref(k, j) * a_ref(k, k);
			    }
			    i__2 = *m;
			    for (i__ = k + 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * a_ref(
					i__, k);
/* L60: */
			    }
			}
/* L70: */
		    }
/* L80: */
		}
	    }
	} else {
/*           Form  B := alpha*A'*B. */
	    if (upper) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    for (i__ = *m; i__ >= 1; --i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__2 = i__ - 1;
			for (k = 1; k <= i__2; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L90: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L100: */
		    }
/* L110: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			temp = b_ref(i__, j);
			if (nounit) {
			    temp *= a_ref(i__, i__);
			}
			i__3 = *m;
			for (k = i__ + 1; k <= i__3; ++k) {
			    temp += a_ref(k, i__) * b_ref(k, j);
/* L120: */
			}
			b_ref(i__, j) = *alpha * temp;
/* L130: */
		    }
/* L140: */
		}
	    }
	}
    } else {
	if (lsame_(transa, "N")) {
/*           Form  B := alpha*B*A. */
	    if (upper) {
		for (j = *n; j >= 1; --j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__1 = *m;
		    for (i__ = 1; i__ <= i__1; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L150: */
		    }
		    i__1 = j - 1;
		    for (k = 1; k <= i__1; ++k) {
			if (a_ref(k, j) != 0.) {
			    temp = *alpha * a_ref(k, j);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L160: */
			    }
			}
/* L170: */
		    }
/* L180: */
		}
	    } else {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			b_ref(i__, j) = temp * b_ref(i__, j);
/* L190: */
		    }
		    i__2 = *n;
		    for (k = j + 1; k <= i__2; ++k) {
			if (a_ref(k, j) != 0.) {
			    temp = *alpha * a_ref(k, j);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L200: */
			    }
			}
/* L210: */
		    }
/* L220: */
		}
	    }
	} else {
/*           Form  B := alpha*B*A'. */
	    if (upper) {
		i__1 = *n;
		for (k = 1; k <= i__1; ++k) {
		    i__2 = k - 1;
		    for (j = 1; j <= i__2; ++j) {
			if (a_ref(j, k) != 0.) {
			    temp = *alpha * a_ref(j, k);
			    i__3 = *m;
			    for (i__ = 1; i__ <= i__3; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L230: */
			    }
			}
/* L240: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.) {
			i__2 = *m;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L250: */
			}
		    }
/* L260: */
		}
	    } else {
		for (k = *n; k >= 1; --k) {
		    i__1 = *n;
		    for (j = k + 1; j <= i__1; ++j) {
			if (a_ref(j, k) != 0.) {
			    temp = *alpha * a_ref(j, k);
			    i__2 = *m;
			    for (i__ = 1; i__ <= i__2; ++i__) {
				b_ref(i__, j) = b_ref(i__, j) + temp * b_ref(
					i__, k);
/* L270: */
			    }
			}
/* L280: */
		    }
		    temp = *alpha;
		    if (nounit) {
			temp *= a_ref(k, k);
		    }
		    if (temp != 1.) {
			i__1 = *m;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    b_ref(i__, k) = temp * b_ref(i__, k);
/* L290: */
			}
		    }
/* L300: */
		}
	    }
	}
    }
    return 0;
/*     End of DTRMM . */
} /* dtrmm_ */
#undef b_ref
#undef a_ref


/* Subroutine */ int dgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
	integer *ldc)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    /* Local variables */
    static integer info;
    static logical nota, notb;
    static doublereal temp;
    static integer i__, j, l, ncola;
    extern logical lsame_(char *, char *);
    static integer nrowa, nrowb;
    extern /* Subroutine */ int xerbla_(char *, integer *);
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
#define b_ref(a_1,a_2) b[(a_2)*b_dim1 + a_1]
#define c___ref(a_1,a_2) c__[(a_2)*c_dim1 + a_1]
/*  Purpose   
    =======   
    DGEMM  performs one of the matrix-matrix operations   
       C := alpha*op( A )*op( B ) + beta*C,   
    where  op( X ) is one of   
       op( X ) = X   or   op( X ) = X',   
    alpha and beta are scalars, and A, B and C are matrices, with op( A )   
    an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.   
    Parameters   
    ==========   
    TRANSA - CHARACTER*1.   
             On entry, TRANSA specifies the form of op( A ) to be used in   
             the matrix multiplication as follows:   
                TRANSA = 'N' or 'n',  op( A ) = A.   
                TRANSA = 'T' or 't',  op( A ) = A'.   
                TRANSA = 'C' or 'c',  op( A ) = A'.   
             Unchanged on exit.   
    TRANSB - CHARACTER*1.   
             On entry, TRANSB specifies the form of op( B ) to be used in   
             the matrix multiplication as follows:   
                TRANSB = 'N' or 'n',  op( B ) = B.   
                TRANSB = 'T' or 't',  op( B ) = B'.   
                TRANSB = 'C' or 'c',  op( B ) = B'.   
             Unchanged on exit.   
    M      - INTEGER.   
             On entry,  M  specifies  the number  of rows  of the  matrix   
             op( A )  and of the  matrix  C.  M  must  be at least  zero.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry,  N  specifies the number  of columns of the matrix   
             op( B ) and the number of columns of the matrix C. N must be   
             at least zero.   
             Unchanged on exit.   
    K      - INTEGER.   
             On entry,  K  specifies  the number of columns of the matrix   
             op( A ) and the number of rows of the matrix op( B ). K must   
             be at least  zero.   
             Unchanged on exit.   
    ALPHA  - DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is   
             k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.   
             Before entry with  TRANSA = 'N' or 'n',  the leading  m by k   
             part of the array  A  must contain the matrix  A,  otherwise   
             the leading  k by m  part of the array  A  must contain  the   
             matrix A.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. When  TRANSA = 'N' or 'n' then   
             LDA must be at least  max( 1, m ), otherwise  LDA must be at   
             least  max( 1, k ).   
             Unchanged on exit.   
    B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is   
             n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.   
             Before entry with  TRANSB = 'N' or 'n',  the leading  k by n   
             part of the array  B  must contain the matrix  B,  otherwise   
             the leading  n by k  part of the array  B  must contain  the   
             matrix B.   
             Unchanged on exit.   
    LDB    - INTEGER.   
             On entry, LDB specifies the first dimension of B as declared   
             in the calling (sub) program. When  TRANSB = 'N' or 'n' then   
             LDB must be at least  max( 1, k ), otherwise  LDB must be at   
             least  max( 1, n ).   
             Unchanged on exit.   
    BETA   - DOUBLE PRECISION.   
             On entry,  BETA  specifies the scalar  beta.  When  BETA  is   
             supplied as zero then C need not be set on input.   
             Unchanged on exit.   
    C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).   
             Before entry, the leading  m by n  part of the array  C must   
             contain the matrix  C,  except when  beta  is zero, in which   
             case C need not be set on entry.   
             On exit, the array  C  is overwritten by the  m by n  matrix   
             ( alpha*op( A )*op( B ) + beta*C ).   
    LDC    - INTEGER.   
             On entry, LDC specifies the first dimension of C as declared   
             in  the  calling  (sub)  program.   LDC  must  be  at  least   
             max( 1, m ).   
             Unchanged on exit.   
    Level 3 Blas routine.   
    -- Written on 8-February-1989.   
       Jack Dongarra, Argonne National Laboratory.   
       Iain Duff, AERE Harwell.   
       Jeremy Du Croz, Numerical Algorithms Group Ltd.   
       Sven Hammarling, Numerical Algorithms Group Ltd.   
       Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not   
       transposed and set  NROWA, NCOLA and  NROWB  as the number of rows   
       and  columns of  A  and the  number of  rows  of  B  respectively.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1 * 1;
    b -= b_offset;
    c_dim1 = *ldc;
    c_offset = 1 + c_dim1 * 1;
    c__ -= c_offset;
    /* Function Body */
    nota = lsame_(transa, "N");
    notb = lsame_(transb, "N");
    if (nota) {
	nrowa = *m;
	ncola = *k;
    } else {
	nrowa = *k;
	ncola = *m;
    }
    if (notb) {
	nrowb = *k;
    } else {
	nrowb = *n;
    }
/*     Test the input parameters. */
    info = 0;
    if (! nota && ! lsame_(transa, "C") && ! lsame_(
	    transa, "T")) {
	info = 1;
    } else if (! notb && ! lsame_(transb, "C") && ! 
	    lsame_(transb, "T")) {
	info = 2;
    } else if (*m < 0) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*k < 0) {
	info = 5;
    } else if (*lda < max(1,nrowa)) {
	info = 8;
    } else if (*ldb < max(1,nrowb)) {
	info = 10;
    } else if (*ldc < max(1,*m)) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("DGEMM ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*m == 0 || *n == 0 || (*alpha == 0. || *k == 0) && *beta == 1.) {
	return 0;
    }
/*     And if  alpha.eq.zero. */
    if (*alpha == 0.) {
	if (*beta == 0.) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = 0.;
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    c___ref(i__, j) = *beta * c___ref(i__, j);
/* L30: */
		}
/* L40: */
	    }
	}
	return 0;
    }
/*     Start the operations. */
    if (notb) {
	if (nota) {
/*           Form  C := alpha*A*B + beta*C. */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L50: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L60: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(l, j) != 0.) {
			temp = *alpha * b_ref(l, j);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L70: */
			}
		    }
/* L80: */
		}
/* L90: */
	    }
	} else {
/*           Form  C := alpha*A'*B + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(l, j);
/* L100: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L110: */
		}
/* L120: */
	    }
	}
    } else {
	if (nota) {
/*           Form  C := alpha*A*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		if (*beta == 0.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = 0.;
/* L130: */
		    }
		} else if (*beta != 1.) {
		    i__2 = *m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			c___ref(i__, j) = *beta * c___ref(i__, j);
/* L140: */
		    }
		}
		i__2 = *k;
		for (l = 1; l <= i__2; ++l) {
		    if (b_ref(j, l) != 0.) {
			temp = *alpha * b_ref(j, l);
			i__3 = *m;
			for (i__ = 1; i__ <= i__3; ++i__) {
			    c___ref(i__, j) = c___ref(i__, j) + temp * a_ref(
				    i__, l);
/* L150: */
			}
		    }
/* L160: */
		}
/* L170: */
	    }
	} else {
/*           Form  C := alpha*A'*B' + beta*C */
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp = 0.;
		    i__3 = *k;
		    for (l = 1; l <= i__3; ++l) {
			temp += a_ref(l, i__) * b_ref(j, l);
/* L180: */
		    }
		    if (*beta == 0.) {
			c___ref(i__, j) = *alpha * temp;
		    } else {
			c___ref(i__, j) = *alpha * temp + *beta * c___ref(i__,
				 j);
		    }
/* L190: */
		}
/* L200: */
	    }
	}
    }
    return 0;
/*     End of DGEMM . */
} /* dgemm_ */
#undef c___ref
#undef b_ref
#undef a_ref


/* Subroutine */ int dcopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m, ix, iy, mp1;
/*     copies a vector, x, to a vector, y.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dy;
    --dx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */


/* Subroutine */ int dtrmv_(char *uplo, char *trans, char *diag, integer *n, 
	doublereal *a, integer *lda, doublereal *x, integer *incx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    /* Local variables */
    static integer info;
    static doublereal temp;
    static integer i__, j;
    extern logical lsame_(char *, char *);
    static integer ix, jx, kx;
    extern /* Subroutine */ int xerbla_(char *, integer *);
    static logical nounit;
#define a_ref(a_1,a_2) a[(a_2)*a_dim1 + a_1]
/*  Purpose   
    =======   
    DTRMV  performs one of the matrix-vector operations   
       x := A*x,   or   x := A'*x,   
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
    Parameters   
    ==========   
    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   
                UPLO = 'U' or 'u'   A is an upper triangular matrix.   
                UPLO = 'L' or 'l'   A is a lower triangular matrix.   
             Unchanged on exit.   
    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   x := A*x.   
                TRANS = 'T' or 't'   x := A'*x.   
                TRANS = 'C' or 'c'   x := A'*x.   
             Unchanged on exit.   
    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   
                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   
                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   
             Unchanged on exit.   
    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   
    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper   
             triangular matrix and the strictly lower triangular part of   
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower   
             triangular matrix and the strictly upper triangular part of   
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of   
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   
    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared   
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   
    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x. On exit, X is overwritten with the   
             tranformed vector x.   
    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   
    Level 2 Blas routine.   
    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   
       Test the input parameters.   
       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;
    --x;
    /* Function Body */
    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, 
	    "T") && ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, 
	    "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRMV ", &info);
	return 0;
    }
/*     Quick return if possible. */
    if (*n == 0) {
	return 0;
    }
    nounit = lsame_(diag, "N");
/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */
    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }
/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */
    if (lsame_(trans, "N")) {
/*        Form  x := A*x. */
	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[i__] += temp * a_ref(i__, j);
/* L10: */
			}
			if (nounit) {
			    x[j] *= a_ref(j, j);
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__2 = j - 1;
			for (i__ = 1; i__ <= i__2; ++i__) {
			    x[ix] += temp * a_ref(i__, j);
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    x[jx] *= a_ref(j, j);
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (x[j] != 0.) {
			temp = x[j];
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[i__] += temp * a_ref(i__, j);
/* L50: */
			}
			if (nounit) {
			    x[j] *= a_ref(j, j);
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (x[jx] != 0.) {
			temp = x[jx];
			ix = kx;
			i__1 = j + 1;
			for (i__ = *n; i__ >= i__1; --i__) {
			    x[ix] += temp * a_ref(i__, j);
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    x[jx] *= a_ref(j, j);
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {
/*        Form  x := A'*x. */
	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			temp += a_ref(i__, j) * x[i__];
/* L90: */
		    }
		    x[j] = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    for (i__ = j - 1; i__ >= 1; --i__) {
			ix -= *incx;
			temp += a_ref(i__, j) * x[ix];
/* L110: */
		    }
		    x[jx] = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[j];
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			temp += a_ref(i__, j) * x[i__];
/* L130: */
		    }
		    x[j] = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= i__1; ++j) {
		    temp = x[jx];
		    ix = jx;
		    if (nounit) {
			temp *= a_ref(j, j);
		    }
		    i__2 = *n;
		    for (i__ = j + 1; i__ <= i__2; ++i__) {
			ix += *incx;
			temp += a_ref(i__, j) * x[ix];
/* L150: */
		    }
		    x[jx] = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }
    return 0;
/*     End of DTRMV . */
} /* dtrmv_ */
#undef a_ref


/* Subroutine */ int drot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__;
    static doublereal dtemp;
    static integer ix, iy;
/*     applies a plane rotation.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dy;
    --dx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*       code for unequal increments or equal increments not equal   
           to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[ix] + *s * dy[iy];
	dy[iy] = *c__ * dy[iy] - *s * dx[ix];
	dx[ix] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*       code for both increments equal to 1 */
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[i__] + *s * dy[i__];
	dy[i__] = *c__ * dy[i__] - *s * dx[i__];
	dx[i__] = dtemp;
/* L30: */
    }
    return 0;
} /* drot_ */


integer idamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;
    /* Local variables */
    static doublereal dmax__;
    static integer i__, ix;
/*     finds the index of element having max. absolute value.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dx;
    /* Function Body */
    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }
/*        code for increment not equal to 1 */
    ix = 1;
    dmax__ = abs(dx[1]);
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[ix], abs(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;
/*        code for increment equal to 1 */
L20:
    dmax__ = abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[i__], abs(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[i__], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */


/* Subroutine */ int daxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;
    /* Local variables */
    static integer i__, m, ix, iy, mp1;
/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dy;
    --dx;
    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */


doublereal ddot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;
    /* Local variables */
    static integer i__, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;
/*     forms the dot product of two vectors.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   
       Parameter adjustments */
    --dy;
    --dx;
    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }
/*        code for unequal increments or equal increments   
            not equal to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;
/*        code for both increments equal to 1   
          clean-up loop */
L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

#endif /* HAVE_DGEGV && __HACK_EIG__ */

