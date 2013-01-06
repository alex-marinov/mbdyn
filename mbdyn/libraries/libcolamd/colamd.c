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

/*
 * Colamd is distributed with MBDyn because its license seems to allow it.
 * The original Copyright is preserved and reported below; credit of course
 * goes to the original Authors.
 *
 * Colamd is used by the naive solver.  The functions have been renamed
 * in the mbdyn_* namespace for compatibility with other solvers that may
 * require linking the original colamd functions with incompatible types.
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */

#include "ac/f2c.h"

#else /* !HAVE_CONFIG_H */
/* to ease compilation outside of MBDyn...
 * replace long and double with the preferred types */
#include <math.h>
typedef long int integer;
typedef double doublereal;
#endif /* !HAVE_CONFIG_H */

/* ========================================================================== */
/* === colamd/symamd - a sparse matrix column ordering algorithm ============ */
/* ========================================================================== */

/*
    colamd:  an approximate minimum degree column ordering algorithm,
    	for LU factorization of symmetric or unsymmetric matrices,
	QR factorization, least squares, interior point methods for
	linear programming problems, and other related problems.

    symamd:  an approximate minimum degree ordering algorithm for Cholesky
    	factorization of symmetric matrices.

    Purpose:

	Colamd computes a permutation Q such that the Cholesky factorization of
	(AQ)'(AQ) has less fill-in and requires fewer floating point operations
	than A'A.  This also provides a good ordering for sparse partial
	pivoting methods, P(AQ) = LU, where Q is computed prior to numerical
	factorization, and P is computed during numerical factorization via
	conventional partial pivoting with row interchanges.  Colamd is the
	column ordering method used in SuperLU, part of the ScaLAPACK library.
	It is also available as built-in function in Matlab Version 6,
	available from MathWorks, Inc. (http://www.mathworks.com).  This
	routine can be used in place of colmmd in Matlab.  By default, the \
	and / operators in Matlab perform a column ordering (using colmmd
	or colamd) prior to LU factorization using sparse partial pivoting,
	in the built-in Matlab lu(A) routine.

    	Symamd computes a permutation P of a symmetric matrix A such that the
	Cholesky factorization of PAP' has less fill-in and requires fewer
	floating point operations than A.  Symamd constructs a matrix M such
	that M'M has the same nonzero pattern of A, and then orders the columns
	of M using colmmd.  The column ordering of M is then returned as the
	row and column ordering P of A. 

    Authors:

	The authors of the code itself are Stefan I. Larimore and Timothy A.
	Davis (davis@cise.ufl.edu), University of Florida.  The algorithm was
	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
	Ng, Oak Ridge National Laboratory.

    Date:

	May 4, 2001.  Version 2.1.

    Acknowledgements:

	This work was supported by the National Science Foundation, under
	grants DMS-9504974 and DMS-9803599.

    Notice:

	Copyright (c) 1998-2001 by the University of Florida.
	All Rights Reserved.

	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.

	Permission is hereby granted to use or copy this program for any
	purpose, provided the above notices are retained on all copies.
	User documentation of any code that uses this code must cite the
	Authors, the Copyright, and "Used by permission."  If this code is
	accessible from within Matlab, then typing "help colamd" and "help
	symamd" must cite the Authors.  Permission to modify the code and to
	distribute modified code is granted, provided the above notices are
	retained, and a notice that the code was modified is included with the
	above copyright notice.  You must also retain the Availability
	information below, of the original version.

	This software is provided free of charge.

    Availability:

	The colamd/symamd library is available at

	    http://www.cise.ufl.edu/research/sparse/colamd/

	This is the http://www.cise.ufl.edu/research/sparse/colamd/colamd.c
	file.  It requires the colamd.h file.  It is required by the colamdmex.c
	and symamdmex.c files, for the Matlab interface to colamd and symamd.

    Changes to the colamd library since Version 1.0 and 1.1:

	No bugs were found in version 1.1.  These changes merely add new
	functionality.

    	* added the COLAMD_RECOMMENDED (nnz, n_row, n_col) macro.

	* moved the output statistics, from A, to a separate output argument.
		The arguments changed for the C-callable routines.

	* added colamd_report and symamd_report.

	* added a C-callable symamd routine.  Formerly, symamd was only
		available as a mexFunction from Matlab.

	* added error-checking to symamd.  Formerly, it assumed its input
		was error-free.

	* added the optional stats and knobs arguments to the symamd mexFunction

	* deleted colamd_help.  A help message is still available from
		"help colamd" and "help symamd" in Matlab.

	* deleted colamdtree.m and symamdtree.m.  Now, colamd.m and symamd.m
		also do the elimination tree post-ordering.  The Version 1.1
		colamd and symamd mexFunctions, which do not do the post-
		ordering, are now visible as colamdmex and symamdmex from
		Matlab.  Essentialy, the post-ordering is now the default
		behavior of colamd.m and symamd.m, to match the behavior of
		colmmd and symmmd.  The post-ordering is only available in the
		Matlab interface, not the C-callable interface.

	* made a slight change to the dense row/column detection in symamd,
		to match the stated specifications.

    Changes from Version 2.0 to 2.1:

	* TRUE and FALSE are predefined on some systems, so they are defined
		here only if not already defined.
	
	* web site changed

	* UNIX Makefile modified, to handle the case if "." is not in your path.

*/

/* ========================================================================== */
/* === Description of user-callable routines ================================ */
/* ========================================================================== */

/*
    ----------------------------------------------------------------------------
    colamd_recommended:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    integer colamd_recommended (integer nnz, integer n_row, integer n_col) ;

	    or as a C macro

	    #include "colamd.h"
	    Alen = COLAMD_RECOMMENDED (integer nnz, integer n_row, integer n_col) ;

	Purpose:

	    Returns recommended value of Alen for use by colamd.  Returns -1
	    if any input argument is negative.  The use of this routine
	    or macro is optional.  Note that the macro uses its arguments
	    more than once, so be careful for side effects, if you pass
	    expressions as arguments to COLAMD_RECOMMENDED.  Not needed for
	    symamd, which dynamically allocates its own memory.

	Arguments (all input arguments):

	    integer nnz ;		Number of nonzeros in the matrix A.  This must
				be the same value as p [n_col] in the call to
				colamd - otherwise you will get a wrong value
				of the recommended memory to use.

	    integer n_row ;		Number of rows in the matrix A.

	    integer n_col ;		Number of columns in the matrix A.

    ----------------------------------------------------------------------------
    colamd_set_defaults:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    colamd_set_defaults (integer knobs [COLAMD_KNOBS]) ;

	Purpose:

	    Sets the default parameters.  The use of this routine is optional.

	Arguments:

	    double knobs [COLAMD_KNOBS] ;	Output only.

		Colamd: rows with more than (knobs [COLAMD_DENSE_ROW] * n_col)
		entries are removed prior to ordering.  Columns with more than
		(knobs [COLAMD_DENSE_COL] * n_row) entries are removed prior to
		ordering, and placed last in the output column ordering. 

		Symamd: uses only knobs [COLAMD_DENSE_ROW], which is knobs [0].
		Rows and columns with more than (knobs [COLAMD_DENSE_ROW] * n)
		entries are removed prior to ordering, and placed last in the
		output ordering.

		COLAMD_DENSE_ROW and COLAMD_DENSE_COL are defined as 0 and 1,
		respectively, in colamd.h.  Default values of these two knobs
		are both 0.5.  Currently, only knobs [0] and knobs [1] are
		used, but future versions may use more knobs.  If so, they will
		be properly set to their defaults by the future version of
		colamd_set_defaults, so that the code that calls colamd will
		not need to change, assuming that you either use
		colamd_set_defaults, or pass a (double *) NULL pointer as the
		knobs array to colamd or symamd.

    ----------------------------------------------------------------------------
    colamd:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    integer colamd (integer n_row, integer n_col, integer Alen, integer *A, integer *p,
	    	double knobs [COLAMD_KNOBS], integer stats [COLAMD_STATS]) ;

	Purpose:

	    Computes a column ordering (Q) of A such that P(AQ)=LU or
	    (AQ)'AQ=LL' have less fill-in and require fewer floating point
	    operations than factorizing the unpermuted matrix A or A'A,
	    respectively.
	    
	Returns:

	    TRUE (1) if successful, FALSE (0) otherwise.

	Arguments:

	    integer n_row ;		Input argument.

		Number of rows in the matrix A.
		Restriction:  n_row >= 0.
		Colamd returns FALSE if n_row is negative.

	    integer n_col ;		Input argument.

		Number of columns in the matrix A.
		Restriction:  n_col >= 0.
		Colamd returns FALSE if n_col is negative.

	    integer Alen ;		Input argument.

		Restriction (see note):
		Alen >= 2*nnz + 6*(n_col+1) + 4*(n_row+1) + n_col
		Colamd returns FALSE if these conditions are not met.

		Note:  this restriction makes an modest assumption regarding
		the size of the two typedef's structures in colamd.h.
		We do, however, guarantee that

			Alen >= colamd_recommended (nnz, n_row, n_col)
		
		or equivalently as a C preprocessor macro: 

			Alen >= COLAMD_RECOMMENDED (nnz, n_row, n_col)

		will be sufficient.

	    integer A [Alen] ;	Input argument, undefined on output.

		A is an integer array of size Alen.  Alen must be at least as
		large as the bare minimum value given above, but this is very
		low, and can result in excessive run time.  For best
		performance, we recommend that Alen be greater than or equal to
		colamd_recommended (nnz, n_row, n_col), which adds
		nnz/5 to the bare minimum value given above.

		On input, the row indices of the entries in column c of the
		matrix are held in A [(p [c]) ... (p [c+1]-1)].  The row indices
		in a given column c need not be in ascending order, and
		duplicate row indices may be be present.  However, colamd will
		work a little faster if both of these conditions are met
		(Colamd puts the matrix into this format, if it finds that the
		the conditions are not met).

		The matrix is 0-based.  That is, rows are in the range 0 to
		n_row-1, and columns are in the range 0 to n_col-1.  Colamd
		returns FALSE if any row index is out of range.

		The contents of A are modified during ordering, and are
		undefined on output.

	    integer p [n_col+1] ;	Both input and output argument.

		p is an integer array of size n_col+1.  On input, it holds the
		"pointers" for the column form of the matrix A.  Column c of
		the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
		entry, p [0], must be zero, and p [c] <= p [c+1] must hold
		for all c in the range 0 to n_col-1.  The value p [n_col] is
		thus the total number of entries in the pattern of the matrix A.
		Colamd returns FALSE if these conditions are not met.

		On output, if colamd returns TRUE, the array p holds the column
		permutation (Q, for P(AQ)=LU or (AQ)'(AQ)=LL'), where p [0] is
		the first column index in the new ordering, and p [n_col-1] is
		the last.  That is, p [k] = j means that column j of A is the
		kth pivot column, in AQ, where k is in the range 0 to n_col-1
		(p [0] = j means that column j of A is the first column in AQ).

		If colamd returns FALSE, then no permutation is returned, and
		p is undefined on output.

	    double knobs [COLAMD_KNOBS] ;	Input argument.

		See colamd_set_defaults for a description.

	    integer stats [COLAMD_STATS] ;		Output argument.

		Statistics on the ordering, and error status.
		See colamd.h for related definitions.
		Colamd returns FALSE if stats is not present.

		stats [0]:  number of dense or empty rows ignored.

		stats [1]:  number of dense or empty columns ignored (and
				ordered last in the output permutation p)
				Note that a row can become "empty" if it
				contains only "dense" and/or "empty" columns,
				and similarly a column can become "empty" if it
				only contains "dense" and/or "empty" rows.

		stats [2]:  number of garbage collections performed.
				This can be excessively high if Alen is close
				to the minimum required value.

		stats [3]:  status code.  < 0 is an error code.
			    > 1 is a warning or notice.

			0	OK.  Each column of the input matrix contained
				row indices in increasing order, with no
				duplicates.

			1	OK, but columns of input matrix were jumbled
				(unsorted columns or duplicate entries).  Colamd
				had to do some extra work to sort the matrix
				first and remove duplicate entries, but it
				still was able to return a valid permutation
				(return value of colamd was TRUE).

					stats [4]: highest numbered column that
						is unsorted or has duplicate
						entries.
					stats [5]: last seen duplicate or
						unsorted row index.
					stats [6]: number of duplicate or
						unsorted row indices.

			-1	A is a null pointer

			-2	p is a null pointer

			-3 	n_row is negative

					stats [4]: n_row

			-4	n_col is negative

					stats [4]: n_col

			-5	number of nonzeros in matrix is negative

					stats [4]: number of nonzeros, p [n_col]

			-6	p [0] is nonzero

					stats [4]: p [0]

			-7	A is too small

					stats [4]: required size
					stats [5]: actual size (Alen)

			-8	a column has a negative number of entries

					stats [4]: column with < 0 entries
					stats [5]: number of entries in col

			-9	a row index is out of bounds

					stats [4]: column with bad row index
					stats [5]: bad row index
					stats [6]: n_row, # of rows of matrx

			-10	(unused; see symamd.c)

			-999	(unused; see symamd.c)

		Future versions may return more statistics in the stats array.

	Example:
	
	    See http://www.cise.ufl.edu/research/sparse/colamd/example.c
	    for a complete example.

	    To order the columns of a 5-by-4 matrix with 11 nonzero entries in
	    the following nonzero pattern

	    	x 0 x 0
		x 0 x x
		0 x x 0
		0 0 x x
		x x 0 0

	    with default knobs and no output statistics, do the following:

		#include "colamd.h"
		#define ALEN COLAMD_RECOMMENDED (11, 5, 4)
		integer A [ALEN] = {1, 2, 5, 3, 5, 1, 2, 3, 4, 2, 4} ;
		integer p [ ] = {0, 3, 5, 9, 11} ;
		integer stats [COLAMD_STATS] ;
		colamd (5, 4, ALEN, A, p, (double *) NULL, stats) ;

	    The permutation is returned in the array p, and A is destroyed.

    ----------------------------------------------------------------------------
    symamd:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    integer symamd (integer n, integer *A, integer *p, integer *perm,
	    	integer knobs [COLAMD_KNOBS], integer stats [COLAMD_STATS],
		void (*allocate) (size_t, size_t), void (*release) (void *)) ;

	Purpose:

    	    The symamd routine computes an ordering P of a symmetric sparse
	    matrix A such that the Cholesky factorization PAP' = LL' remains
	    sparse.  It is based on a column ordering of a matrix M constructed
	    so that the nonzero pattern of M'M is the same as A.  The matrix A
	    is assumed to be symmetric; only the strictly lower triangular part
	    is accessed.  You must pass your selected memory allocator (usually
	    calloc/free or mxCalloc/mxFree) to symamd, for it to allocate
	    memory for the temporary matrix M.

	Returns:

	    TRUE (1) if successful, FALSE (0) otherwise.

	Arguments:

	    integer n ;		Input argument.

	    	Number of rows and columns in the symmetrix matrix A.
		Restriction:  n >= 0.
		Symamd returns FALSE if n is negative.

	    integer A [nnz] ;	Input argument.

	    	A is an integer array of size nnz, where nnz = p [n].
		
		The row indices of the entries in column c of the matrix are
		held in A [(p [c]) ... (p [c+1]-1)].  The row indices in a
		given column c need not be in ascending order, and duplicate
		row indices may be present.  However, symamd will run faster
		if the columns are in sorted order with no duplicate entries. 

		The matrix is 0-based.  That is, rows are in the range 0 to
		n-1, and columns are in the range 0 to n-1.  Symamd
		returns FALSE if any row index is out of range.

		The contents of A are not modified.

	    integer p [n+1] ;   	Input argument.

		p is an integer array of size n+1.  On input, it holds the
		"pointers" for the column form of the matrix A.  Column c of
		the matrix A is held in A [(p [c]) ... (p [c+1]-1)].  The first
		entry, p [0], must be zero, and p [c] <= p [c+1] must hold
		for all c in the range 0 to n-1.  The value p [n] is
		thus the total number of entries in the pattern of the matrix A.
		Symamd returns FALSE if these conditions are not met.

		The contents of p are not modified.

	    integer perm [n+1] ;   	Output argument.

		On output, if symamd returns TRUE, the array perm holds the
		permutation P, where perm [0] is the first index in the new
		ordering, and perm [n-1] is the last.  That is, perm [k] = j
		means that row and column j of A is the kth column in PAP',
		where k is in the range 0 to n-1 (perm [0] = j means
		that row and column j of A are the first row and column in
		PAP').  The array is used as a workspace during the ordering,
		which is why it must be of length n+1, not just n.

	    double knobs [COLAMD_KNOBS] ;	Input argument.

		See colamd_set_defaults for a description.

	    integer stats [COLAMD_STATS] ;		Output argument.

		Statistics on the ordering, and error status.
		See colamd.h for related definitions.
		Symamd returns FALSE if stats is not present.

		stats [0]:  number of dense or empty row and columns ignored
				(and ordered last in the output permutation 
				perm).  Note that a row/column can become
				"empty" if it contains only "dense" and/or
				"empty" columns/rows.

		stats [1]:  (same as stats [0])

		stats [2]:  number of garbage collections performed.

		stats [3]:  status code.  < 0 is an error code.
			    > 1 is a warning or notice.

			0	OK.  Each column of the input matrix contained
				row indices in increasing order, with no
				duplicates.

			1	OK, but columns of input matrix were jumbled
				(unsorted columns or duplicate entries).  Symamd
				had to do some extra work to sort the matrix
				first and remove duplicate entries, but it
				still was able to return a valid permutation
				(return value of symamd was TRUE).

					stats [4]: highest numbered column that
						is unsorted or has duplicate
						entries.
					stats [5]: last seen duplicate or
						unsorted row index.
					stats [6]: number of duplicate or
						unsorted row indices.

			-1	A is a null pointer

			-2	p is a null pointer

			-3	(unused, see colamd.c)

			-4 	n is negative

					stats [4]: n

			-5	number of nonzeros in matrix is negative

					stats [4]: # of nonzeros (p [n]).

			-6	p [0] is nonzero

					stats [4]: p [0]

			-7	(unused)

			-8	a column has a negative number of entries

					stats [4]: column with < 0 entries
					stats [5]: number of entries in col

			-9	a row index is out of bounds

					stats [4]: column with bad row index
					stats [5]: bad row index
					stats [6]: n_row, # of rows of matrx

			-10	out of memory (unable to allocate temporary
				workspace for M or count arrays using the
				"allocate" routine passed into symamd).

			-999	internal error.  colamd failed to order the
				matrix M, when it should have succeeded.  This
				indicates a bug.  If this (and *only* this)
				error code occurs, please contact the authors.
				Don't contact the authors if you get any other
				error code.

		Future versions may return more statistics in the stats array.

	    void * (*allocate) (size_t, size_t)

	    	A pointer to a function providing memory allocation.  The
		allocated memory must be returned initialized to zero.  For a
		C application, this argument should normally be a pointer to
		calloc.  For a Matlab mexFunction, the routine mxCalloc is
		passed instead.

	    void (*release) (size_t, size_t)

	    	A pointer to a function that frees memory allocated by the
		memory allocation routine above.  For a C application, this
		argument should normally be a pointer to free.  For a Matlab
		mexFunction, the routine mxFree is passed instead.


    ----------------------------------------------------------------------------
    colamd_report:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    colamd_report (integer stats [COLAMD_STATS]) ;

	Purpose:

	    Prints the error status and statistics recorded in the stats
	    array on the standard error output (for a standard C routine)
	    or on the Matlab output (for a mexFunction).

	Arguments:

	    integer stats [COLAMD_STATS] ;	Input only.  Statistics from colamd.


    ----------------------------------------------------------------------------
    symamd_report:
    ----------------------------------------------------------------------------

	C syntax:

	    #include "colamd.h"
	    symamd_report (integer stats [COLAMD_STATS]) ;

	Purpose:

	    Prints the error status and statistics recorded in the stats
	    array on the standard error output (for a standard C routine)
	    or on the Matlab output (for a mexFunction).

	Arguments:

	    integer stats [COLAMD_STATS] ;	Input only.  Statistics from symamd.


*/

/* ========================================================================== */
/* === Scaffolding code definitions  ======================================== */
/* ========================================================================== */

/* Ensure that debugging is turned off: */
#ifndef NDEBUG
#define NDEBUG
#endif /* NDEBUG */

/*
   Our "scaffolding code" philosophy:  In our opinion, well-written library
   code should keep its "debugging" code, and just normally have it turned off
   by the compiler so as not to interfere with performance.  This serves
   several purposes:

   (1) assertions act as comments to the reader, telling you what the code
	expects at that point.  All assertions will always be true (unless
	there really is a bug, of course).

   (2) leaving in the scaffolding code assists anyone who would like to modify
	the code, or understand the algorithm (by reading the debugging output,
	one can get a glimpse into what the code is doing).

   (3) (gasp!) for actually finding bugs.  This code has been heavily tested
	and "should" be fully functional and bug-free ... but you never know...

    To enable debugging, comment out the "#define NDEBUG" above.  For a Matlab
    mexFunction, you will also need to modify mexopts.sh to remove the -DNDEBUG
    definition.  The code will become outrageously slow when debugging is
    enabled.  To control the level of debugging output, set an environment
    variable D to 0 (little), 1 (some), 2, 3, or 4 (lots).  When debugging,
    you should see the following message on the standard output:

    	colamd: debug version, D = 1 (THIS WILL BE SLOW!)

    or a similar message for symamd.  If you don't, then debugging has not
    been enabled.

*/

/* ========================================================================== */
/* === Include files ======================================================== */
/* ========================================================================== */

#include "colamd.h"
#include <limits.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#else
#include <stdio.h>
#include <assert.h>
#endif /* MATLAB_MEX_FILE */

/* ========================================================================== */
/* === Definitions ========================================================== */
/* ========================================================================== */

/* Routines are either PUBLIC (user-callable) or PRIVATE (not user-callable) */
#define PUBLIC
#define PRIVATE static

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define ONES_COMPLEMENT(r) (-(r)-1)

/* -------------------------------------------------------------------------- */
/* Change for version 2.1:  define TRUE and FALSE only if not yet defined */  
/* -------------------------------------------------------------------------- */

#ifndef TRUE
#define TRUE (1)
#endif

#ifndef FALSE
#define FALSE (0)
#endif

/* -------------------------------------------------------------------------- */

#define EMPTY	(-1)

/* Row and column status */
#define ALIVE	(0)
#define DEAD	(-1)

/* Column status */
#define DEAD_PRINCIPAL		(-1)
#define DEAD_NON_PRINCIPAL	(-2)

/* Macros for row and column status update and checking. */
#define ROW_IS_DEAD(r)			ROW_IS_MARKED_DEAD (Row[r].shared2.mark)
#define ROW_IS_MARKED_DEAD(row_mark)	(row_mark < ALIVE)
#define ROW_IS_ALIVE(r)			(Row [r].shared2.mark >= ALIVE)
#define COL_IS_DEAD(c)			(Col [c].start < ALIVE)
#define COL_IS_ALIVE(c)			(Col [c].start >= ALIVE)
#define COL_IS_DEAD_PRINCIPAL(c)	(Col [c].start == DEAD_PRINCIPAL)
#define KILL_ROW(r)			{ Row [r].shared2.mark = DEAD ; }
#define KILL_PRINCIPAL_COL(c)		{ Col [c].start = DEAD_PRINCIPAL ; }
#define KILL_NON_PRINCIPAL_COL(c)	{ Col [c].start = DEAD_NON_PRINCIPAL ; }

/* ========================================================================== */
/* === Colamd reporting mechanism =========================================== */
/* ========================================================================== */

#ifdef MATLAB_MEX_FILE

/* use mexPrintf in a Matlab mexFunction, for debugging and statistics output */
#define PRINTF mexPrintf

/* In Matlab, matrices are 1-based to the user, but 0-based internally */
#define INDEX(i) ((i)+1)

#else

/* Use printf in standard C environment, for debugging and statistics output. */
/* Output is generated only if debugging is enabled at compile time, or if */
/* the caller explicitly calls colamd_report or symamd_report. */
#define PRINTF printf

/* In C, matrices are 0-based and indices are reported as such in *_report */
#define INDEX(i) (i)

#endif /* MATLAB_MEX_FILE */

/* ========================================================================== */
/* === Prototypes of PRIVATE routines ======================================= */
/* ========================================================================== */

PRIVATE integer init_rows_cols
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer p [],
    integer stats [COLAMD_STATS]
) ;

PRIVATE void init_scoring
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer head [],
    double knobs [COLAMD_KNOBS],
    integer *p_n_row2,
    integer *p_n_col2,
    integer *p_max_deg
) ;

PRIVATE integer find_ordering
(
    integer n_row,
    integer n_col,
    integer Alen,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer head [],
    integer n_col2,
    integer max_deg,
    integer pfree
) ;

PRIVATE void order_children
(
    integer n_col,
    mbdyn_Colamd_Col Col [],
    integer p []
) ;

PRIVATE void detect_super_cols
(

#ifndef NDEBUG
    integer n_col,
    mbdyn_Colamd_Row Row [],
#endif /* NDEBUG */

    mbdyn_Colamd_Col Col [],
    integer A [],
    integer head [],
    integer row_start,
    integer row_length
) ;

PRIVATE integer garbage_collection
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer *pfree
) ;

PRIVATE integer clear_mark
(
    integer n_row,
    mbdyn_Colamd_Row Row []
) ;

PRIVATE void print_report
(
    char *method,
    integer stats [COLAMD_STATS]
) ;

/* ========================================================================== */
/* === Debugging prototypes and definitions ================================= */
/* ========================================================================== */

#ifndef NDEBUG

/* colamd_debug is the *ONLY* global variable, and is only */
/* present when debugging */

PRIVATE integer colamd_debug ;	/* debug print level */

#define DEBUG0(params) { (void) PRINTF params ; }
#define DEBUG1(params) { if (colamd_debug >= 1) (void) PRINTF params ; }
#define DEBUG2(params) { if (colamd_debug >= 2) (void) PRINTF params ; }
#define DEBUG3(params) { if (colamd_debug >= 3) (void) PRINTF params ; }
#define DEBUG4(params) { if (colamd_debug >= 4) (void) PRINTF params ; }

#ifdef MATLAB_MEX_FILE
#define ASSERT(expression) (mxAssert ((expression), ""))
#else
#define ASSERT(expression) (assert (expression))
#endif /* MATLAB_MEX_FILE */

PRIVATE void colamd_get_debug	/* gets the debug print level from getenv */
(
    char *method
) ;

PRIVATE void debug_deg_lists
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer head [],
    integer min_score,
    integer should,
    integer max_deg
) ;

PRIVATE void debug_mark
(
    integer n_row,
    mbdyn_Colamd_Row Row [],
    integer tag_mark,
    integer max_mark
) ;

PRIVATE void debug_matrix
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A []
) ;

PRIVATE void debug_structures
(
    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer n_col2
) ;

#else /* NDEBUG */

/* === No debugging ========================================================= */

#define DEBUG0(params) ;
#define DEBUG1(params) ;
#define DEBUG2(params) ;
#define DEBUG3(params) ;
#define DEBUG4(params) ;

#define ASSERT(expression) ((void) 0)

#endif /* NDEBUG */

/* ========================================================================== */



/* ========================================================================== */
/* === USER-CALLABLE ROUTINES: ============================================== */
/* ========================================================================== */


/* ========================================================================== */
/* === colamd_recommended =================================================== */
/* ========================================================================== */

/*
    The colamd_recommended routine returns the suggested size for Alen.  This
    value has been determined to provide good balance between the number of
    garbage collections and the memory requirements for colamd.  If any
    argument is negative, a -1 is returned as an error condition.  This
    function is also available as a macro defined in colamd.h, so that you
    can use it for a statically-allocated array size.
*/

PUBLIC integer mbdyn_colamd_recommended	/* returns recommended value of Alen. */
(
    /* === Parameters ======================================================= */

    integer nnz,			/* number of nonzeros in A */
    integer n_row,			/* number of rows in A */
    integer n_col			/* number of columns in A */
)
{
    return (COLAMD_RECOMMENDED (nnz, n_row, n_col)) ; 
}


/* ========================================================================== */
/* === colamd_set_defaults ================================================== */
/* ========================================================================== */

/*
    The colamd_set_defaults routine sets the default values of the user-
    controllable parameters for colamd:

	knobs [0]	rows with knobs[0]*n_col entries or more are removed
			prior to ordering in colamd.  Rows and columns with
			knobs[0]*n_col entries or more are removed prior to
			ordering in symamd and placed last in the output
			ordering.

	knobs [1]	columns with knobs[1]*n_row entries or more are removed
			prior to ordering in colamd, and placed last in the
			column permutation.  Symamd ignores this knob.

	knobs [2..19]	unused, but future versions might use this
*/

PUBLIC void mbdyn_colamd_set_defaults
(
    /* === Parameters ======================================================= */

    double knobs [COLAMD_KNOBS]		/* knob array */
)
{
    /* === Local variables ================================================== */

    integer i ;

    if (!knobs)
    {
	return ;			/* no knobs to initialize */
    }
    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	knobs [i] = 0 ;
    }
    knobs [COLAMD_DENSE_ROW] = 0.5 ;	/* ignore rows over 50% dense */
    knobs [COLAMD_DENSE_COL] = 0.5 ;	/* ignore columns over 50% dense */
}


/* ========================================================================== */
/* === symamd =============================================================== */
/* ========================================================================== */

PUBLIC integer mbdyn_symamd			/* return TRUE if OK, FALSE otherwise */
(
    /* === Parameters ======================================================= */

    integer n,				/* number of rows and columns of A */
    integer A [],				/* row indices of A */
    integer p [],				/* column pointers of A */
    integer perm [],			/* output permutation, size n+1 */
    double knobs [COLAMD_KNOBS],	/* parameters (uses defaults if NULL) */
    integer stats [COLAMD_STATS],		/* output statistics and error codes */
    void * (*allocate) (size_t, size_t),
    					/* pointer to calloc (ANSI C) or */
					/* mxCalloc (for Matlab mexFunction) */
    void (*release) (void *)
    					/* pointer to free (ANSI C) or */
    					/* mxFree (for Matlab mexFunction) */
)
{
    /* === Local variables ================================================== */

    integer *count ;		/* length of each column of M, and col pointer*/
    integer *mark ;			/* mark array for finding duplicate entries */
    integer *M ;			/* row indices of matrix M */
    integer Mlen ;			/* length of M */
    integer n_row ;			/* number of rows in M */
    integer nnz ;			/* number of entries in A */
    integer i ;			/* row index of A */
    integer j ;			/* column index of A */
    integer k ;			/* row index of M */ 
    integer mnz ;			/* number of nonzeros in M */
    integer pp ;			/* index into a column of A */
    integer last_row ;		/* last row seen in the current column */
    integer length ;		/* number of nonzeros in a column */

    double cknobs [COLAMD_KNOBS] ;		/* knobs for colamd */
    double default_knobs [COLAMD_KNOBS] ;	/* default knobs for colamd */
    integer cstats [COLAMD_STATS] ;			/* colamd stats */

#ifndef NDEBUG
    colamd_get_debug ("symamd") ;
#endif /* NDEBUG */

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("symamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    if (!A)
    {
    	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	DEBUG0 (("symamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	DEBUG0 (("symamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n < 0)		/* n must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n ;
	DEBUG0 (("symamd: n negative %d\n", n)) ;
    	return (FALSE) ;
    }

    nnz = p [n] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	DEBUG0 (("symamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero ;
	stats [COLAMD_INFO1] = p [0] ;
	DEBUG0 (("symamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	mbdyn_colamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    /* === Allocate count and mark ========================================== */

    count = (integer *) ((*allocate) (n+1, sizeof (int))) ;
    if (!count)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	DEBUG0 (("symamd: allocate count (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }

    mark = (integer *) ((*allocate) (n+1, sizeof (int))) ;
    if (!mark)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	DEBUG0 (("symamd: allocate mark (size %d) failed\n", n+1)) ;
	return (FALSE) ;
    }

    /* === Compute column counts of M, check if A is valid ================== */

    stats [COLAMD_INFO3] = 0 ;  /* number of duplicate or unsorted row indices*/

    for (i = 0 ; i < n ; i++)
    {
    	mark [i] = -1 ;
    }

    for (j = 0 ; j < n ; j++)
    {
	last_row = -1 ;

	length = p [j+1] - p [j] ;
	if (length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = j ;
	    stats [COLAMD_INFO2] = length ;
	    (*release) ((void *) count) ;
	    (*release) ((void *) mark) ;
	    DEBUG0 (("symamd: col %d negative length %d\n", j, length)) ;
	    return (FALSE) ;
	}

	for (pp = p [j] ; pp < p [j+1] ; pp++)
	{
	    i = A [pp] ;
	    if (i < 0 || i >= n)
	    {
		/* row index i, in column j, is out of bounds */
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		stats [COLAMD_INFO3] = n ;
		(*release) ((void *) count) ;
		(*release) ((void *) mark) ;
		DEBUG0 (("symamd: row %d col %d out of bounds\n", i, j)) ;
		return (FALSE) ;
	    }

	    if (i <= last_row || mark [i] == j)
	    {
		/* row index is unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = j ;
		stats [COLAMD_INFO2] = i ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("symamd: row %d col %d unsorted/duplicate\n", i, j)) ;
	    }

	    if (i > j && mark [i] != j)
	    {
		/* row k of M will contain column indices i and j */
		count [i]++ ;
		count [j]++ ;
	    }

	    /* mark the row as having been seen in this column */
	    mark [i] = j ;

	    last_row = i ;
	}
    }

    if (stats [COLAMD_STATUS] == COLAMD_OK)
    {
	/* if there are no duplicate entries, then mark is no longer needed */
	(*release) ((void *) mark) ;
    }

    /* === Compute column pointers of M ===================================== */

    /* use output permutation, perm, for column pointers of M */
    perm [0] = 0 ;
    for (j = 1 ; j <= n ; j++)
    {
	perm [j] = perm [j-1] + count [j-1] ;
    }
    for (j = 0 ; j < n ; j++)
    {
	count [j] = perm [j] ;
    }

    /* === Construct M ====================================================== */

    mnz = perm [n] ;
    n_row = mnz / 2 ;
    Mlen = mbdyn_colamd_recommended (mnz, n_row, n) ;
    M = (integer *) ((*allocate) (Mlen, sizeof (int))) ;
    DEBUG0 (("symamd: M is %d-by-%d with %d entries, Mlen = %d\n",
    	n_row, n, mnz, Mlen)) ;

    if (!M)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_out_of_memory ;
	(*release) ((void *) count) ;
	(*release) ((void *) mark) ;
	DEBUG0 (("symamd: allocate M (size %d) failed\n", Mlen)) ;
	return (FALSE) ;
    }

    k = 0 ;

    if (stats [COLAMD_STATUS] == COLAMD_OK)
    {
	/* Matrix is OK */
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j)
		{
		    /* row k of M contains column indices i and j */
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		}
	    }
	}
    }
    else
    {
	/* Matrix is jumbled.  Do not add duplicates to M.  Unsorted cols OK. */
	DEBUG0 (("symamd: Duplicates in A.\n")) ;
	for (i = 0 ; i < n ; i++)
	{
	    mark [i] = -1 ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    ASSERT (p [j+1] - p [j] >= 0) ;
	    for (pp = p [j] ; pp < p [j+1] ; pp++)
	    {
		i = A [pp] ;
		ASSERT (i >= 0 && i < n) ;
		if (i > j && mark [i] != j)
		{
		    /* row k of M contains column indices i and j */
		    M [count [i]++] = k ;
		    M [count [j]++] = k ;
		    k++ ;
		    mark [i] = j ;
		}
	    }
	}
	(*release) ((void *) mark) ;
    }

    /* count and mark no longer needed */
    (*release) ((void *) count) ;
    ASSERT (k == n_row) ;

    /* === Adjust the knobs for M =========================================== */

    for (i = 0 ; i < COLAMD_KNOBS ; i++)
    {
	cknobs [i] = knobs [i] ;
    }

    /* there are no dense rows in M */
    cknobs [COLAMD_DENSE_ROW] = 1.0 ;

    if (n_row != 0 && n < n_row)
    {
	/* On input, the knob is a fraction of 1..n, the number of rows of A. */
	/* Convert it to a fraction of 1..n_row, of the number of rows of M. */
    	cknobs [COLAMD_DENSE_COL] = (knobs [COLAMD_DENSE_ROW] * n) / n_row ;
    }
    else
    {
	/* no dense columns in M */
    	cknobs [COLAMD_DENSE_COL] = 1.0 ;
    }

    DEBUG0 (("symamd: dense col knob for M: %g\n", cknobs [COLAMD_DENSE_COL])) ;

    /* === Order the columns of M =========================================== */

    if (!mbdyn_colamd (n_row, n, Mlen, M, perm, cknobs, cstats))
    {
	/* This "cannot" happen, unless there is a bug in the code. */
	stats [COLAMD_STATUS] = COLAMD_ERROR_internal_error ;
	(*release) ((void *) M) ;
	DEBUG0 (("symamd: internal error!\n")) ;
	return (FALSE) ;
    }

    /* Note that the output permutation is now in perm */

    /* === get the statistics for symamd from colamd ======================== */

    /* note that a dense column in colamd means a dense row and col in symamd */
    stats [COLAMD_DENSE_ROW]    = cstats [COLAMD_DENSE_COL] ;
    stats [COLAMD_DENSE_COL]    = cstats [COLAMD_DENSE_COL] ;
    stats [COLAMD_DEFRAG_COUNT] = cstats [COLAMD_DEFRAG_COUNT] ;

    /* === Free M =========================================================== */

    (*release) ((void *) M) ;
    DEBUG0 (("symamd: done.\n")) ;
    return (TRUE) ;

}

/* ========================================================================== */
/* === colamd =============================================================== */
/* ========================================================================== */

/*
    The colamd routine computes a column ordering Q of a sparse matrix
    A such that the LU factorization P(AQ) = LU remains sparse, where P is
    selected via partial pivoting.   The routine can also be viewed as
    providing a permutation Q such that the Cholesky factorization
    (AQ)'(AQ) = LL' remains sparse.
*/

PUBLIC integer mbdyn_colamd		/* returns TRUE if successful, FALSE otherwise*/
(
    /* === Parameters ======================================================= */

    integer n_row,			/* number of rows in A */
    integer n_col,			/* number of columns in A */
    integer Alen,			/* length of A */
    integer A [],			/* row indices of A */
    integer p [],			/* pointers to columns in A */
    double knobs [COLAMD_KNOBS],/* parameters (uses defaults if NULL) */
    integer stats [COLAMD_STATS]	/* output statistics and error codes */
)
{
    /* === Local variables ================================================== */

    integer i ;			/* loop index */
    integer nnz ;			/* nonzeros in A */
    integer Row_size ;		/* size of Row [], in integers */
    integer Col_size ;		/* size of Col [], in integers */
    integer need ;			/* minimum required length of A */
    mbdyn_Colamd_Row *Row ;		/* pointer into A of Row [0..n_row] array */
    mbdyn_Colamd_Col *Col ;		/* pointer into A of Col [0..n_col] array */
    integer n_col2 ;		/* number of non-dense, non-empty columns */
    integer n_row2 ;		/* number of non-dense, non-empty rows */
    integer ngarbage ;		/* number of garbage collections performed */
    integer max_deg ;		/* maximum row degree */
    double default_knobs [COLAMD_KNOBS] ;	/* default knobs array */

#ifndef NDEBUG
    colamd_get_debug ("colamd") ;
#endif /* NDEBUG */

    /* === Check the input arguments ======================================== */

    if (!stats)
    {
	DEBUG0 (("colamd: stats not present\n")) ;
	return (FALSE) ;
    }
    for (i = 0 ; i < COLAMD_STATS ; i++)
    {
	stats [i] = 0 ;
    }
    stats [COLAMD_STATUS] = COLAMD_OK ;
    stats [COLAMD_INFO1] = -1 ;
    stats [COLAMD_INFO2] = -1 ;

    if (!A)		/* A is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_not_present ;
	DEBUG0 (("colamd: A not present\n")) ;
	return (FALSE) ;
    }

    if (!p)		/* p is not present */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p_not_present ;
	DEBUG0 (("colamd: p not present\n")) ;
    	return (FALSE) ;
    }

    if (n_row < 0)	/* n_row must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nrow_negative ;
	stats [COLAMD_INFO1] = n_row ;
	DEBUG0 (("colamd: nrow negative %d\n", n_row)) ;
    	return (FALSE) ;
    }

    if (n_col < 0)	/* n_col must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_ncol_negative ;
	stats [COLAMD_INFO1] = n_col ;
	DEBUG0 (("colamd: ncol negative %d\n", n_col)) ;
    	return (FALSE) ;
    }

    nnz = p [n_col] ;
    if (nnz < 0)	/* nnz must be >= 0 */
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_nnz_negative ;
	stats [COLAMD_INFO1] = nnz ;
	DEBUG0 (("colamd: number of entries negative %d\n", nnz)) ;
	return (FALSE) ;
    }

    if (p [0] != 0)
    {
	stats [COLAMD_STATUS] = COLAMD_ERROR_p0_nonzero	;
	stats [COLAMD_INFO1] = p [0] ;
	DEBUG0 (("colamd: p[0] not zero %d\n", p [0])) ;
	return (FALSE) ;
    }

    /* === If no knobs, set default knobs =================================== */

    if (!knobs)
    {
	mbdyn_colamd_set_defaults (default_knobs) ;
	knobs = default_knobs ;
    }

    /* === Allocate the Row and Col arrays from array A ===================== */

    Col_size = COLAMD_C (n_col) ;
    Row_size = COLAMD_R (n_row) ;
    need = 2*nnz + n_col + Col_size + Row_size ;

    if (need > Alen)
    {
	/* not enough space in array A to perform the ordering */
	stats [COLAMD_STATUS] = COLAMD_ERROR_A_too_small ;
	stats [COLAMD_INFO1] = need ;
	stats [COLAMD_INFO2] = Alen ;
	DEBUG0 (("colamd: Need Alen >= %d, given only Alen = %d\n", need,Alen));
	return (FALSE) ;
    }

    Alen -= Col_size + Row_size ;
    Col = (mbdyn_Colamd_Col *) &A [Alen] ;
    Row = (mbdyn_Colamd_Row *) &A [Alen + Col_size] ;

    /* === Construct the row and column data structures ===================== */

    if (!init_rows_cols (n_row, n_col, Row, Col, A, p, stats))
    {
	/* input matrix is invalid */
	DEBUG0 (("colamd: Matrix invalid\n")) ;
	return (FALSE) ;
    }

    /* === Initialize scores, kill dense rows/columns ======================= */

    init_scoring (n_row, n_col, Row, Col, A, p, knobs,
	&n_row2, &n_col2, &max_deg) ;

    /* === Order the supercolumns =========================================== */

    ngarbage = find_ordering (n_row, n_col, Alen, Row, Col, A, p,
	n_col2, max_deg, 2*nnz) ;

    /* === Order the non-principal columns ================================== */

    order_children (n_col, Col, p) ;

    /* === Return statistics in stats ======================================= */

    stats [COLAMD_DENSE_ROW] = n_row - n_row2 ;
    stats [COLAMD_DENSE_COL] = n_col - n_col2 ;
    stats [COLAMD_DEFRAG_COUNT] = ngarbage ;
    DEBUG0 (("colamd: done.\n")) ; 
    return (TRUE) ;
}


/* ========================================================================== */
/* === colamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void mbdyn_colamd_report
(
    integer stats [COLAMD_STATS]
)
{
    print_report ("colamd", stats) ;
}


/* ========================================================================== */
/* === symamd_report ======================================================== */
/* ========================================================================== */

PUBLIC void mbdyn_symamd_report
(
    integer stats [COLAMD_STATS]
)
{
    print_report ("symamd", stats) ;
}



/* ========================================================================== */
/* === NON-USER-CALLABLE ROUTINES: ========================================== */
/* ========================================================================== */

/* There are no user-callable routines beyond this point in the file */


/* ========================================================================== */
/* === init_rows_cols ======================================================= */
/* ========================================================================== */

/*
    Takes the column form of the matrix in A and creates the row form of the
    matrix.  Also, row and column attributes are stored in the Col and Row
    structs.  If the columns are un-sorted or contain duplicate row indices,
    this routine will also sort and remove duplicate row indices from the
    column form of the matrix.  Returns FALSE if the matrix is invalid,
    TRUE otherwise.  Not user-callable.
*/

PRIVATE integer init_rows_cols	/* returns TRUE if OK, or FALSE otherwise */
(
    /* === Parameters ======================================================= */

    integer n_row,			/* number of rows of A */
    integer n_col,			/* number of columns of A */
    mbdyn_Colamd_Row Row [],		/* of size n_row+1 */
    mbdyn_Colamd_Col Col [],		/* of size n_col+1 */
    integer A [],			/* row indices of A, of size Alen */
    integer p [],			/* pointers to columns in A, of size n_col+1 */
    integer stats [COLAMD_STATS]	/* colamd statistics */ 
)
{
    /* === Local variables ================================================== */

    integer col ;			/* a column index */
    integer row ;			/* a row index */
    integer *cp ;			/* a column pointer */
    integer *cp_end ;		/* a pointer to the end of a column */
    integer *rp ;			/* a row pointer */
    integer *rp_end ;		/* a pointer to the end of a row */
    integer last_row ;		/* previous row */

    /* === Initialize columns, and check column pointers ==================== */

    for (col = 0 ; col < n_col ; col++)
    {
	Col [col].start = p [col] ;
	Col [col].length = p [col+1] - p [col] ;

	if (Col [col].length < 0)
	{
	    /* column pointers must be non-decreasing */
	    stats [COLAMD_STATUS] = COLAMD_ERROR_col_length_negative ;
	    stats [COLAMD_INFO1] = col ;
	    stats [COLAMD_INFO2] = Col [col].length ;
	    DEBUG0 (("colamd: col %d length %d < 0\n", col, Col [col].length)) ;
	    return (FALSE) ;
	}

	Col [col].shared1.thickness = 1 ;
	Col [col].shared2.score = 0 ;
	Col [col].shared3.prev = EMPTY ;
	Col [col].shared4.degree_next = EMPTY ;
    }

    /* p [0..n_col] no longer needed, used as "head" in subsequent routines */

    /* === Scan columns, compute row degrees, and check row indices ========= */

    stats [COLAMD_INFO3] = 0 ;	/* number of duplicate or unsorted row indices*/

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].length = 0 ;
	Row [row].shared2.mark = -1 ;
    }

    for (col = 0 ; col < n_col ; col++)
    {
	last_row = -1 ;

	cp = &A [p [col]] ;
	cp_end = &A [p [col+1]] ;

	while (cp < cp_end)
	{
	    row = *cp++ ;

	    /* make sure row indices within range */
	    if (row < 0 || row >= n_row)
	    {
		stats [COLAMD_STATUS] = COLAMD_ERROR_row_index_out_of_bounds ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		stats [COLAMD_INFO3] = n_row ;
		DEBUG0 (("colamd: row %d col %d out of bounds\n", row, col)) ;
		return (FALSE) ;
	    }

	    if (row <= last_row || Row [row].shared2.mark == col)
	    {
		/* row index are unsorted or repeated (or both), thus col */
		/* is jumbled.  This is a notice, not an error condition. */
		stats [COLAMD_STATUS] = COLAMD_OK_BUT_JUMBLED ;
		stats [COLAMD_INFO1] = col ;
		stats [COLAMD_INFO2] = row ;
		(stats [COLAMD_INFO3]) ++ ;
		DEBUG1 (("colamd: row %d col %d unsorted/duplicate\n",row,col));
	    }

	    if (Row [row].shared2.mark != col)
	    {
		Row [row].length++ ;
	    }
	    else
	    {
		/* this is a repeated entry in the column, */
		/* it will be removed */
		Col [col].length-- ;
	    }

	    /* mark the row as having been seen in this column */
	    Row [row].shared2.mark = col ;

	    last_row = row ;
	}
    }

    /* === Compute row pointers ============================================= */

    /* row form of the matrix starts directly after the column */
    /* form of matrix in A */
    Row [0].start = p [n_col] ;
    Row [0].shared1.p = Row [0].start ;
    Row [0].shared2.mark = -1 ;
    for (row = 1 ; row < n_row ; row++)
    {
	Row [row].start = Row [row-1].start + Row [row-1].length ;
	Row [row].shared1.p = Row [row].start ;
	Row [row].shared2.mark = -1 ;
    }

    /* === Create row form ================================================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
	/* if cols jumbled, watch for repeated row indices */
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		row = *cp++ ;
		if (Row [row].shared2.mark != col)
		{
		    A [(Row [row].shared1.p)++] = col ;
		    Row [row].shared2.mark = col ;
		}
	    }
	}
    }
    else
    {
	/* if cols not jumbled, we don't need the mark (this is faster) */
	for (col = 0 ; col < n_col ; col++)
	{
	    cp = &A [p [col]] ;
	    cp_end = &A [p [col+1]] ;
	    while (cp < cp_end)
	    {
		A [(Row [*cp++].shared1.p)++] = col ;
	    }
	}
    }

    /* === Clear the row marks and set row degrees ========================== */

    for (row = 0 ; row < n_row ; row++)
    {
	Row [row].shared2.mark = 0 ;
	Row [row].shared1.degree = Row [row].length ;
    }

    /* === See if we need to re-create columns ============================== */

    if (stats [COLAMD_STATUS] == COLAMD_OK_BUT_JUMBLED)
    {
    	DEBUG0 (("colamd: reconstructing column form, matrix jumbled\n")) ;

#ifndef NDEBUG
	/* make sure column lengths are correct */
	for (col = 0 ; col < n_col ; col++)
	{
	    p [col] = Col [col].length ;
	}
	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		p [*rp++]-- ;
	    }
	}
	for (col = 0 ; col < n_col ; col++)
	{
	    ASSERT (p [col] == 0) ;
	}
	/* now p is all zero (different than when debugging is turned off) */
#endif /* NDEBUG */

	/* === Compute col pointers ========================================= */

	/* col form of the matrix starts at A [0]. */
	/* Note, we may have a gap between the col form and the row */
	/* form if there were duplicate entries, if so, it will be */
	/* removed upon the first garbage collection */
	Col [0].start = 0 ;
	p [0] = Col [0].start ;
	for (col = 1 ; col < n_col ; col++)
	{
	    /* note that the lengths here are for pruned columns, i.e. */
	    /* no duplicate row indices will exist for these columns */
	    Col [col].start = Col [col-1].start + Col [col-1].length ;
	    p [col] = Col [col].start ;
	}

	/* === Re-create col form =========================================== */

	for (row = 0 ; row < n_row ; row++)
	{
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		A [(p [*rp++])++] = row ;
	    }
	}
    }

    /* === Done.  Matrix is not (or no longer) jumbled ====================== */

    return (TRUE) ;
}


/* ========================================================================== */
/* === init_scoring ========================================================= */
/* ========================================================================== */

/*
    Kills dense or empty columns and rows, calculates an initial score for
    each column, and places all columns in the degree lists.  Not user-callable.
*/

PRIVATE void init_scoring
(
    /* === Parameters ======================================================= */

    integer n_row,			/* number of rows of A */
    integer n_col,			/* number of columns of A */
    mbdyn_Colamd_Row Row [],		/* of size n_row+1 */
    mbdyn_Colamd_Col Col [],		/* of size n_col+1 */
    integer A [],			/* column form and row form of A */
    integer head [],		/* of size n_col+1 */
    double knobs [COLAMD_KNOBS],/* parameters */
    integer *p_n_row2,		/* number of non-dense, non-empty rows */
    integer *p_n_col2,		/* number of non-dense, non-empty columns */
    integer *p_max_deg		/* maximum row degree */
)
{
    /* === Local variables ================================================== */

    integer c ;			/* a column index */
    integer r, row ;		/* a row index */
    integer *cp ;			/* a column pointer */
    integer deg ;			/* degree of a row or column */
    integer *cp_end ;		/* a pointer to the end of a column */
    integer *new_cp ;		/* new column pointer */
    integer col_length ;		/* length of pruned column */
    integer score ;			/* current column score */
    integer n_col2 ;		/* number of non-dense, non-empty columns */
    integer n_row2 ;		/* number of non-dense, non-empty rows */
    integer dense_row_count ;	/* remove rows with more entries than this */
    integer dense_col_count ;	/* remove cols with more entries than this */
    integer min_score ;		/* smallest column score */
    integer max_deg ;		/* maximum row degree */
    integer next_col ;		/* Used to add to degree list.*/

#ifndef NDEBUG
    integer debug_count ;		/* debug only. */
#endif /* NDEBUG */

    /* === Extract knobs ==================================================== */

    dense_row_count = MAX (0, MIN (knobs [COLAMD_DENSE_ROW] * n_col, n_col)) ;
    dense_col_count = MAX (0, MIN (knobs [COLAMD_DENSE_COL] * n_row, n_row)) ;
    DEBUG1 (("colamd: densecount: %d %d\n", dense_row_count, dense_col_count)) ;
    max_deg = 0 ;
    n_col2 = n_col ;
    n_row2 = n_row ;

    /* === Kill empty columns =============================================== */

    /* Put the empty columns at the end in their natural order, so that LU */
    /* factorization can proceed as far as possible. */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	deg = Col [c].length ;
	if (deg == 0)
	{
	    /* this is a empty column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: null columns killed: %d\n", n_col - n_col2)) ;

    /* === Kill dense columns =============================================== */

    /* Put the dense columns at the end, in their natural order */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip any dead columns */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	deg = Col [c].length ;
	if (deg > dense_col_count)
	{
	    /* this is a dense column, kill and order it last */
	    Col [c].shared2.order = --n_col2 ;
	    /* decrement the row degrees */
	    cp = &A [Col [c].start] ;
	    cp_end = cp + Col [c].length ;
	    while (cp < cp_end)
	    {
		Row [*cp++].shared1.degree-- ;
	    }
	    KILL_PRINCIPAL_COL (c) ;
	}
    }
    DEBUG1 (("colamd: Dense and null columns killed: %d\n", n_col - n_col2)) ;

    /* === Kill dense and empty rows ======================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	deg = Row [r].shared1.degree ;
	ASSERT (deg >= 0 && deg <= n_col) ;
	if (deg > dense_row_count || deg == 0)
	{
	    /* kill a dense or empty row */
	    KILL_ROW (r) ;
	    --n_row2 ;
	}
	else
	{
	    /* keep track of max degree of remaining rows */
	    max_deg = MAX (max_deg, deg) ;
	}
    }
    DEBUG1 (("colamd: Dense and null rows killed: %d\n", n_row - n_row2)) ;

    /* === Compute initial column scores ==================================== */

    /* At this point the row degrees are accurate.  They reflect the number */
    /* of "live" (non-dense) columns in each row.  No empty rows exist. */
    /* Some "live" columns may contain only dead rows, however.  These are */
    /* pruned in the code below. */

    /* now find the initial matlab score for each column */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* skip dead column */
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	score = 0 ;
	cp = &A [Col [c].start] ;
	new_cp = cp ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    /* skip if dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    /* compact the column */
	    *new_cp++ = row ;
	    /* add row's external degree */
	    score += Row [row].shared1.degree - 1 ;
	    /* guard against integer overflow */
	    score = MIN (score, n_col) ;
	}
	/* determine pruned column length */
	col_length = (int) (new_cp - &A [Col [c].start]) ;
	if (col_length == 0)
	{
	    /* a newly-made null column (all rows in this col are "dense" */
	    /* and have already been killed) */
	    DEBUG2 (("Newly null killed: %d\n", c)) ;
	    Col [c].shared2.order = --n_col2 ;
	    KILL_PRINCIPAL_COL (c) ;
	}
	else
	{
	    /* set column length and set score */
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    Col [c].length = col_length ;
	    Col [c].shared2.score = score ;
	}
    }
    DEBUG1 (("colamd: Dense, null, and newly-null columns killed: %d\n",
    	n_col-n_col2)) ;

    /* At this point, all empty rows and columns are dead.  All live columns */
    /* are "clean" (containing no dead rows) and simplicial (no supercolumns */
    /* yet).  Rows may contain dead columns, but all live rows contain at */
    /* least one live column. */

#ifndef NDEBUG
    debug_structures (n_row, n_col, Row, Col, A, n_col2) ;
#endif /* NDEBUG */

    /* === Initialize degree lists ========================================== */

#ifndef NDEBUG
    debug_count = 0 ;
#endif /* NDEBUG */

    /* clear the hash buckets */
    for (c = 0 ; c <= n_col ; c++)
    {
	head [c] = EMPTY ;
    }
    min_score = n_col ;
    /* place in reverse order, so low column indices are at the front */
    /* of the lists.  This is to encourage natural tie-breaking */
    for (c = n_col-1 ; c >= 0 ; c--)
    {
	/* only add principal columns to degree lists */
	if (COL_IS_ALIVE (c))
	{
	    DEBUG4 (("place %d score %d minscore %d ncol %d\n",
		c, Col [c].shared2.score, min_score, n_col)) ;

	    /* === Add columns score to DList =============================== */

	    score = Col [c].shared2.score ;

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (score >= 0) ;
	    ASSERT (score <= n_col) ;
	    ASSERT (head [score] >= EMPTY) ;

	    /* now add this column to dList at proper score location */
	    next_col = head [score] ;
	    Col [c].shared3.prev = EMPTY ;
	    Col [c].shared4.degree_next = next_col ;

	    /* if there already was a column with the same score, set its */
	    /* previous pointer to this new column */
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = c ;
	    }
	    head [score] = c ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, score) ;

#ifndef NDEBUG
	    debug_count++ ;
#endif /* NDEBUG */

	}
    }

#ifndef NDEBUG
    DEBUG1 (("colamd: Live cols %d out of %d, non-princ: %d\n",
	debug_count, n_col, n_col-debug_count)) ;
    ASSERT (debug_count == n_col2) ;
    debug_deg_lists (n_row, n_col, Row, Col, head, min_score, n_col2, max_deg) ;
#endif /* NDEBUG */

    /* === Return number of remaining columns, and max row degree =========== */

    *p_n_col2 = n_col2 ;
    *p_n_row2 = n_row2 ;
    *p_max_deg = max_deg ;
}


/* ========================================================================== */
/* === find_ordering ======================================================== */
/* ========================================================================== */

/*
    Order the principal columns of the supercolumn form of the matrix
    (no supercolumns on input).  Uses a minimum approximate column minimum
    degree ordering method.  Not user-callable.
*/

PRIVATE integer find_ordering	/* return the number of garbage collections */
(
    /* === Parameters ======================================================= */

    integer n_row,			/* number of rows of A */
    integer n_col,			/* number of columns of A */
    integer Alen,			/* size of A, 2*nnz + n_col or larger */
    mbdyn_Colamd_Row Row [],		/* of size n_row+1 */
    mbdyn_Colamd_Col Col [],		/* of size n_col+1 */
    integer A [],			/* column form and row form of A */
    integer head [],		/* of size n_col+1 */
    integer n_col2,			/* Remaining columns to order */
    integer max_deg,		/* Maximum row degree */
    integer pfree			/* index of first free slot (2*nnz on entry) */
)
{
    /* === Local variables ================================================== */

    integer k ;			/* current pivot ordering step */
    integer pivot_col ;		/* current pivot column */
    integer *cp ;			/* a column pointer */
    integer *rp ;			/* a row pointer */
    integer pivot_row ;		/* current pivot row */
    integer *new_cp ;		/* modified column pointer */
    integer *new_rp ;		/* modified row pointer */
    integer pivot_row_start ;	/* pointer to start of pivot row */
    integer pivot_row_degree ;	/* number of columns in pivot row */
    integer pivot_row_length ;	/* number of supercolumns in pivot row */
    integer pivot_col_score ;	/* score of pivot column */
    integer needed_memory ;		/* free space needed for pivot row */
    integer *cp_end ;		/* pointer to the end of a column */
    integer *rp_end ;		/* pointer to the end of a row */
    integer row ;			/* a row index */
    integer col ;			/* a column index */
    integer max_score ;		/* maximum possible score */
    integer cur_score ;		/* score of current column */
    unsigned int hash ;		/* hash value for supernode detection */
    integer head_column ;		/* head of hash bucket */
    integer first_col ;		/* first column in hash bucket */
    integer tag_mark ;		/* marker value for mark array */
    integer row_mark ;		/* Row [row].shared2.mark */
    integer set_difference ;	/* set difference size of row with pivot row */
    integer min_score ;		/* smallest column score */
    integer col_thickness ;		/* "thickness" (no. of columns in a supercol) */
    integer max_mark ;		/* maximum value of tag_mark */
    integer pivot_col_thickness ;	/* number of columns represented by pivot col */
    integer prev_col ;		/* Used by Dlist operations. */
    integer next_col ;		/* Used by Dlist operations. */
    integer ngarbage ;		/* number of garbage collections performed */

#ifndef NDEBUG
    integer debug_d ;		/* debug loop counter */
    integer debug_step = 0 ;	/* debug loop counter */
#endif /* NDEBUG */

    /* === Initialization and clear mark ==================================== */

    max_mark = INT_MAX - n_col ;	/* INT_MAX defined in <limits.h> */
    tag_mark = clear_mark (n_row, Row) ;
    min_score = 0 ;
    ngarbage = 0 ;
    DEBUG1 (("colamd: Ordering, n_col2=%d\n", n_col2)) ;

    /* === Order the columns ================================================ */

    for (k = 0 ; k < n_col2 ; /* 'k' is incremented below */)
    {

#ifndef NDEBUG
	if (debug_step % 100 == 0)
	{
	    DEBUG2 (("\n...       Step k: %d out of n_col2: %d\n", k, n_col2)) ;
	}
	else
	{
	    DEBUG3 (("\n----------Step k: %d out of n_col2: %d\n", k, n_col2)) ;
	}
	debug_step++ ;
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k, max_deg) ;
	debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */

	/* === Select pivot column, and order it ============================ */

	/* make sure degree list isn't empty */
	ASSERT (min_score >= 0) ;
	ASSERT (min_score <= n_col) ;
	ASSERT (head [min_score] >= EMPTY) ;

#ifndef NDEBUG
	for (debug_d = 0 ; debug_d < min_score ; debug_d++)
	{
	    ASSERT (head [debug_d] == EMPTY) ;
	}
#endif /* NDEBUG */

	/* get pivot column from head of minimum degree list */
	while (head [min_score] == EMPTY && min_score < n_col)
	{
	    min_score++ ;
	}
	pivot_col = head [min_score] ;
	ASSERT (pivot_col >= 0 && pivot_col <= n_col) ;
	next_col = Col [pivot_col].shared4.degree_next ;
	head [min_score] = next_col ;
	if (next_col != EMPTY)
	{
	    Col [next_col].shared3.prev = EMPTY ;
	}

	ASSERT (COL_IS_ALIVE (pivot_col)) ;
	DEBUG3 (("Pivot col: %d\n", pivot_col)) ;

	/* remember score for defrag check */
	pivot_col_score = Col [pivot_col].shared2.score ;

	/* the pivot column is the kth column in the pivot order */
	Col [pivot_col].shared2.order = k ;

	/* increment order count by column thickness */
	pivot_col_thickness = Col [pivot_col].shared1.thickness ;
	k += pivot_col_thickness ;
	ASSERT (pivot_col_thickness > 0) ;

	/* === Garbage_collection, if necessary ============================= */

	needed_memory = MIN (pivot_col_score, n_col - k) ;
	if (pfree + needed_memory >= Alen)
	{
	    pfree = garbage_collection (n_row, n_col, Row, Col, A, &A [pfree]) ;
	    ngarbage++ ;
	    /* after garbage collection we will have enough */
	    ASSERT (pfree + needed_memory < Alen) ;
	    /* garbage collection has wiped out the Row[].shared2.mark array */
	    tag_mark = clear_mark (n_row, Row) ;

#ifndef NDEBUG
	    debug_matrix (n_row, n_col, Row, Col, A) ;
#endif /* NDEBUG */
	}

	/* === Compute pivot row pattern ==================================== */

	/* get starting location for this new merged row */
	pivot_row_start = pfree ;

	/* initialize new row counts to zero */
	pivot_row_degree = 0 ;

	/* tag pivot column as having been visited so it isn't included */
	/* in merged pivot row */
	Col [pivot_col].shared1.thickness = -pivot_col_thickness ;

	/* pivot row is the union of all rows in the pivot column pattern */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* get a row */
	    row = *cp++ ;
	    DEBUG4 (("Pivot col pattern %d %d\n", ROW_IS_ALIVE (row), row)) ;
	    /* skip if row is dead */
	    if (ROW_IS_DEAD (row))
	    {
		continue ;
	    }
	    rp = &A [Row [row].start] ;
	    rp_end = rp + Row [row].length ;
	    while (rp < rp_end)
	    {
		/* get a column */
		col = *rp++ ;
		/* add the column, if alive and untagged */
		col_thickness = Col [col].shared1.thickness ;
		if (col_thickness > 0 && COL_IS_ALIVE (col))
		{
		    /* tag column in pivot row */
		    Col [col].shared1.thickness = -col_thickness ;
		    ASSERT (pfree < Alen) ;
		    /* place column in pivot row */
		    A [pfree++] = col ;
		    pivot_row_degree += col_thickness ;
		}
	    }
	}

	/* clear tag on pivot column */
	Col [pivot_col].shared1.thickness = pivot_col_thickness ;
	max_deg = MAX (max_deg, pivot_row_degree) ;

#ifndef NDEBUG
	DEBUG3 (("check2\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Kill all rows used to construct pivot row ==================== */

	/* also kill pivot row, temporarily */
	cp = &A [Col [pivot_col].start] ;
	cp_end = cp + Col [pivot_col].length ;
	while (cp < cp_end)
	{
	    /* may be killing an already dead row */
	    row = *cp++ ;
	    DEBUG3 (("Kill row in pivot col: %d\n", row)) ;
	    KILL_ROW (row) ;
	}

	/* === Select a row index to use as the new pivot row =============== */

	pivot_row_length = pfree - pivot_row_start ;
	if (pivot_row_length > 0)
	{
	    /* pick the "pivot" row arbitrarily (first row in col) */
	    pivot_row = A [Col [pivot_col].start] ;
	    DEBUG3 (("Pivotal row is %d\n", pivot_row)) ;
	}
	else
	{
	    /* there is no pivot row, since it is of zero length */
	    pivot_row = EMPTY ;
	    ASSERT (pivot_row_length == 0) ;
	}
	ASSERT (Col [pivot_col].length > 0 || pivot_row_length == 0) ;

	/* === Approximate degree computation =============================== */

	/* Here begins the computation of the approximate degree.  The column */
	/* score is the sum of the pivot row "length", plus the size of the */
	/* set differences of each row in the column minus the pattern of the */
	/* pivot row itself.  The column ("thickness") itself is also */
	/* excluded from the column score (we thus use an approximate */
	/* external degree). */

	/* The time taken by the following code (compute set differences, and */
	/* add them up) is proportional to the size of the data structure */
	/* being scanned - that is, the sum of the sizes of each column in */
	/* the pivot row.  Thus, the amortized time to compute a column score */
	/* is proportional to the size of that column (where size, in this */
	/* context, is the column "length", or the number of row indices */
	/* in that column).  The number of row indices in a column is */
	/* monotonically non-decreasing, from the length of the original */
	/* column on input to colamd. */

	/* === Compute set differences ====================================== */

	DEBUG3 (("** Computing set differences phase. **\n")) ;

	/* pivot row is currently dead - it will be revived later. */

	DEBUG3 (("Pivot row: ")) ;
	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    DEBUG3 (("Col: %d\n", col)) ;

	    /* clear tags used to construct pivot row pattern */
	    col_thickness = -Col [col].shared1.thickness ;
	    ASSERT (col_thickness > 0) ;
	    Col [col].shared1.thickness = col_thickness ;

	    /* === Remove column from degree list =========================== */

	    cur_score = Col [col].shared2.score ;
	    prev_col = Col [col].shared3.prev ;
	    next_col = Col [col].shared4.degree_next ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (cur_score >= EMPTY) ;
	    if (prev_col == EMPTY)
	    {
		head [cur_score] = next_col ;
	    }
	    else
	    {
		Col [prev_col].shared4.degree_next = next_col ;
	    }
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = prev_col ;
	    }

	    /* === Scan the column ========================================== */

	    cp = &A [Col [col].start] ;
	    cp_end = cp + Col [col].length ;
	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row != pivot_row) ;
		set_difference = row_mark - tag_mark ;
		/* check if the row has been seen yet */
		if (set_difference < 0)
		{
		    ASSERT (Row [row].shared1.degree <= max_deg) ;
		    set_difference = Row [row].shared1.degree ;
		}
		/* subtract column thickness from this row's set difference */
		set_difference -= col_thickness ;
		ASSERT (set_difference >= 0) ;
		/* absorb this row if the set difference becomes zero */
		if (set_difference == 0)
		{
		    DEBUG3 (("aggressive absorption. Row: %d\n", row)) ;
		    KILL_ROW (row) ;
		}
		else
		{
		    /* save the new mark */
		    Row [row].shared2.mark = set_difference + tag_mark ;
		}
	    }
	}

#ifndef NDEBUG
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k-pivot_row_degree, max_deg) ;
#endif /* NDEBUG */

	/* === Add up set differences for each column ======================= */

	DEBUG3 (("** Adding set differences phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    /* get a column */
	    col = *rp++ ;
	    ASSERT (COL_IS_ALIVE (col) && col != pivot_col) ;
	    hash = 0 ;
	    cur_score = 0 ;
	    cp = &A [Col [col].start] ;
	    /* compact the column */
	    new_cp = cp ;
	    cp_end = cp + Col [col].length ;

	    DEBUG4 (("Adding set diffs for Col: %d.\n", col)) ;

	    while (cp < cp_end)
	    {
		/* get a row */
		row = *cp++ ;
		ASSERT(row >= 0 && row < n_row) ;
		row_mark = Row [row].shared2.mark ;
		/* skip if dead */
		if (ROW_IS_MARKED_DEAD (row_mark))
		{
		    continue ;
		}
		ASSERT (row_mark > tag_mark) ;
		/* compact the column */
		*new_cp++ = row ;
		/* compute hash function */
		hash += row ;
		/* add set difference */
		cur_score += row_mark - tag_mark ;
		/* integer overflow... */
		cur_score = MIN (cur_score, n_col) ;
	    }

	    /* recompute the column's length */
	    Col [col].length = (int) (new_cp - &A [Col [col].start]) ;

	    /* === Further mass elimination ================================= */

	    if (Col [col].length == 0)
	    {
		DEBUG4 (("further mass elimination. Col: %d\n", col)) ;
		/* nothing left but the pivot row in this column */
		KILL_PRINCIPAL_COL (col) ;
		pivot_row_degree -= Col [col].shared1.thickness ;
		ASSERT (pivot_row_degree >= 0) ;
		/* order it */
		Col [col].shared2.order = k ;
		/* increment order count by column thickness */
		k += Col [col].shared1.thickness ;
	    }
	    else
	    {
		/* === Prepare for supercolumn detection ==================== */

		DEBUG4 (("Preparing supercol detection for Col: %d.\n", col)) ;

		/* save score so far */
		Col [col].shared2.score = cur_score ;

		/* add column to hash table, for supercolumn detection */
		hash %= n_col + 1 ;

		DEBUG4 ((" Hash = %d, n_col = %d.\n", hash, n_col)) ;
		ASSERT (hash <= n_col) ;

		head_column = head [hash] ;
		if (head_column > EMPTY)
		{
		    /* degree list "hash" is non-empty, use prev (shared3) of */
		    /* first column in degree list as head of hash bucket */
		    first_col = Col [head_column].shared3.headhash ;
		    Col [head_column].shared3.headhash = col ;
		}
		else
		{
		    /* degree list "hash" is empty, use head as hash bucket */
		    first_col = - (head_column + 2) ;
		    head [hash] = - (col + 2) ;
		}
		Col [col].shared4.hash_next = first_col ;

		/* save hash function in Col [col].shared3.hash */
		Col [col].shared3.hash = (int) hash ;
		ASSERT (COL_IS_ALIVE (col)) ;
	    }
	}

	/* The approximate external column degree is now computed.  */

	/* === Supercolumn detection ======================================== */

	DEBUG3 (("** Supercolumn detection phase. **\n")) ;

	detect_super_cols (

#ifndef NDEBUG
		n_col, Row,
#endif /* NDEBUG */

		Col, A, head, pivot_row_start, pivot_row_length) ;

	/* === Kill the pivotal column ====================================== */

	KILL_PRINCIPAL_COL (pivot_col) ;

	/* === Clear mark =================================================== */

	tag_mark += (max_deg + 1) ;
	if (tag_mark >= max_mark)
	{
	    DEBUG2 (("clearing tag_mark\n")) ;
	    tag_mark = clear_mark (n_row, Row) ;
	}

#ifndef NDEBUG
	DEBUG3 (("check3\n")) ;
	debug_mark (n_row, Row, tag_mark, max_mark) ;
#endif /* NDEBUG */

	/* === Finalize the new pivot row, and column scores ================ */

	DEBUG3 (("** Finalize scores phase. **\n")) ;

	/* for each column in pivot row */
	rp = &A [pivot_row_start] ;
	/* compact the pivot row */
	new_rp = rp ;
	rp_end = rp + pivot_row_length ;
	while (rp < rp_end)
	{
	    col = *rp++ ;
	    /* skip dead columns */
	    if (COL_IS_DEAD (col))
	    {
		continue ;
	    }
	    *new_rp++ = col ;
	    /* add new pivot row to column */
	    A [Col [col].start + (Col [col].length++)] = pivot_row ;

	    /* retrieve score so far and add on pivot row's degree. */
	    /* (we wait until here for this in case the pivot */
	    /* row's degree was reduced due to mass elimination). */
	    cur_score = Col [col].shared2.score + pivot_row_degree ;

	    /* calculate the max possible score as the number of */
	    /* external columns minus the 'k' value minus the */
	    /* columns thickness */
	    max_score = n_col - k - Col [col].shared1.thickness ;

	    /* make the score the external degree of the union-of-rows */
	    cur_score -= Col [col].shared1.thickness ;

	    /* make sure score is less or equal than the max score */
	    cur_score = MIN (cur_score, max_score) ;
	    ASSERT (cur_score >= 0) ;

	    /* store updated score */
	    Col [col].shared2.score = cur_score ;

	    /* === Place column back in degree list ========================= */

	    ASSERT (min_score >= 0) ;
	    ASSERT (min_score <= n_col) ;
	    ASSERT (cur_score >= 0) ;
	    ASSERT (cur_score <= n_col) ;
	    ASSERT (head [cur_score] >= EMPTY) ;
	    next_col = head [cur_score] ;
	    Col [col].shared4.degree_next = next_col ;
	    Col [col].shared3.prev = EMPTY ;
	    if (next_col != EMPTY)
	    {
		Col [next_col].shared3.prev = col ;
	    }
	    head [cur_score] = col ;

	    /* see if this score is less than current min */
	    min_score = MIN (min_score, cur_score) ;

	}

#ifndef NDEBUG
	debug_deg_lists (n_row, n_col, Row, Col, head,
		min_score, n_col2-k, max_deg) ;
#endif /* NDEBUG */

	/* === Resurrect the new pivot row ================================== */

	if (pivot_row_degree > 0)
	{
	    /* update pivot row length to reflect any cols that were killed */
	    /* during super-col detection and mass elimination */
	    Row [pivot_row].start  = pivot_row_start ;
	    Row [pivot_row].length = (int) (new_rp - &A[pivot_row_start]) ;
	    Row [pivot_row].shared1.degree = pivot_row_degree ;
	    Row [pivot_row].shared2.mark = 0 ;
	    /* pivot row is no longer dead */
	}
    }

    /* === All principal columns have now been ordered ====================== */

    return (ngarbage) ;
}


/* ========================================================================== */
/* === order_children ======================================================= */
/* ========================================================================== */

/*
    The find_ordering routine has ordered all of the principal columns (the
    representatives of the supercolumns).  The non-principal columns have not
    yet been ordered.  This routine orders those columns by walking up the
    parent tree (a column is a child of the column which absorbed it).  The
    final permutation vector is then placed in p [0 ... n_col-1], with p [0]
    being the first column, and p [n_col-1] being the last.  It doesn't look
    like it at first glance, but be assured that this routine takes time linear
    in the number of columns.  Although not immediately obvious, the time
    taken by this routine is O (n_col), that is, linear in the number of
    columns.  Not user-callable.
*/

PRIVATE void order_children
(
    /* === Parameters ======================================================= */

    integer n_col,			/* number of columns of A */
    mbdyn_Colamd_Col Col [],		/* of size n_col+1 */
    integer p []			/* p [0 ... n_col-1] is the column permutation*/
)
{
    /* === Local variables ================================================== */

    integer i ;			/* loop counter for all columns */
    integer c ;			/* column index */
    integer parent ;		/* index of column's parent */
    integer order ;			/* column's order */

    /* === Order each non-principal column ================================== */

    for (i = 0 ; i < n_col ; i++)
    {
	/* find an un-ordered non-principal column */
	ASSERT (COL_IS_DEAD (i)) ;
	if (!COL_IS_DEAD_PRINCIPAL (i) && Col [i].shared2.order == EMPTY)
	{
	    parent = i ;
	    /* once found, find its principal parent */
	    do
	    {
		parent = Col [parent].shared1.parent ;
	    } while (!COL_IS_DEAD_PRINCIPAL (parent)) ;

	    /* now, order all un-ordered non-principal columns along path */
	    /* to this parent.  collapse tree at the same time */
	    c = i ;
	    /* get order of parent */
	    order = Col [parent].shared2.order ;

	    do
	    {
		ASSERT (Col [c].shared2.order == EMPTY) ;

		/* order this column */
		Col [c].shared2.order = order++ ;
		/* collaps tree */
		Col [c].shared1.parent = parent ;

		/* get immediate parent of this column */
		c = Col [c].shared1.parent ;

		/* continue until we hit an ordered column.  There are */
		/* guarranteed not to be anymore unordered columns */
		/* above an ordered column */
	    } while (Col [c].shared2.order == EMPTY) ;

	    /* re-order the super_col parent to largest order for this group */
	    Col [parent].shared2.order = order ;
	}
    }

    /* === Generate the permutation ========================================= */

    for (c = 0 ; c < n_col ; c++)
    {
	p [Col [c].shared2.order] = c ;
    }
}


/* ========================================================================== */
/* === detect_super_cols ==================================================== */
/* ========================================================================== */

/*
    Detects supercolumns by finding matches between columns in the hash buckets.
    Check amongst columns in the set A [row_start ... row_start + row_length-1].
    The columns under consideration are currently *not* in the degree lists,
    and have already been placed in the hash buckets.

    The hash bucket for columns whose hash function is equal to h is stored
    as follows:

	if head [h] is >= 0, then head [h] contains a degree list, so:

		head [h] is the first column in degree bucket h.
		Col [head [h]].headhash gives the first column in hash bucket h.

	otherwise, the degree list is empty, and:

		-(head [h] + 2) is the first column in hash bucket h.

    For a column c in a hash bucket, Col [c].shared3.prev is NOT a "previous
    column" pointer.  Col [c].shared3.hash is used instead as the hash number
    for that column.  The value of Col [c].shared4.hash_next is the next column
    in the same hash bucket.

    Assuming no, or "few" hash collisions, the time taken by this routine is
    linear in the sum of the sizes (lengths) of each column whose score has
    just been computed in the approximate degree computation.
    Not user-callable.
*/

PRIVATE void detect_super_cols
(
    /* === Parameters ======================================================= */

#ifndef NDEBUG
    /* these two parameters are only needed when debugging is enabled: */
    integer n_col,			/* number of columns of A */
    mbdyn_Colamd_Row Row [],		/* of size n_row+1 */
#endif /* NDEBUG */

    mbdyn_Colamd_Col Col [],		/* of size n_col+1 */
    integer A [],			/* row indices of A */
    integer head [],		/* head of degree lists and hash buckets */
    integer row_start,		/* pointer to set of columns to check */
    integer row_length		/* number of columns to check */
)
{
    /* === Local variables ================================================== */

    integer hash ;			/* hash value for a column */
    integer *rp ;			/* pointer to a row */
    integer c ;			/* a column index */
    integer super_c ;		/* column index of the column to absorb into */
    integer *cp1 ;			/* column pointer for column super_c */
    integer *cp2 ;			/* column pointer for column c */
    integer length ;		/* length of column super_c */
    integer prev_c ;		/* column preceding c in hash bucket */
    integer i ;			/* loop counter */
    integer *rp_end ;		/* pointer to the end of the row */
    integer col ;			/* a column index in the row to check */
    integer head_column ;		/* first column in hash bucket or degree list */
    integer first_col ;		/* first column in hash bucket */

    /* === Consider each column in the row ================================== */

    rp = &A [row_start] ;
    rp_end = rp + row_length ;
    while (rp < rp_end)
    {
	col = *rp++ ;
	if (COL_IS_DEAD (col))
	{
	    continue ;
	}

	/* get hash number for this column */
	hash = Col [col].shared3.hash ;
	ASSERT (hash <= n_col) ;

	/* === Get the first column in this hash bucket ===================== */

	head_column = head [hash] ;
	if (head_column > EMPTY)
	{
	    first_col = Col [head_column].shared3.headhash ;
	}
	else
	{
	    first_col = - (head_column + 2) ;
	}

	/* === Consider each column in the hash bucket ====================== */

	for (super_c = first_col ; super_c != EMPTY ;
	    super_c = Col [super_c].shared4.hash_next)
	{
	    ASSERT (COL_IS_ALIVE (super_c)) ;
	    ASSERT (Col [super_c].shared3.hash == hash) ;
	    length = Col [super_c].length ;

	    /* prev_c is the column preceding column c in the hash bucket */
	    prev_c = super_c ;

	    /* === Compare super_c with all columns after it ================ */

	    for (c = Col [super_c].shared4.hash_next ;
		 c != EMPTY ; c = Col [c].shared4.hash_next)
	    {
		ASSERT (c != super_c) ;
		ASSERT (COL_IS_ALIVE (c)) ;
		ASSERT (Col [c].shared3.hash == hash) ;

		/* not identical if lengths or scores are different */
		if (Col [c].length != length ||
		    Col [c].shared2.score != Col [super_c].shared2.score)
		{
		    prev_c = c ;
		    continue ;
		}

		/* compare the two columns */
		cp1 = &A [Col [super_c].start] ;
		cp2 = &A [Col [c].start] ;

		for (i = 0 ; i < length ; i++)
		{
		    /* the columns are "clean" (no dead rows) */
		    ASSERT (ROW_IS_ALIVE (*cp1))  ;
		    ASSERT (ROW_IS_ALIVE (*cp2))  ;
		    /* row indices will same order for both supercols, */
		    /* no gather scatter nessasary */
		    if (*cp1++ != *cp2++)
		    {
			break ;
		    }
		}

		/* the two columns are different if the for-loop "broke" */
		if (i != length)
		{
		    prev_c = c ;
		    continue ;
		}

		/* === Got it!  two columns are identical =================== */

		ASSERT (Col [c].shared2.score == Col [super_c].shared2.score) ;

		Col [super_c].shared1.thickness += Col [c].shared1.thickness ;
		Col [c].shared1.parent = super_c ;
		KILL_NON_PRINCIPAL_COL (c) ;
		/* order c later, in order_children() */
		Col [c].shared2.order = EMPTY ;
		/* remove c from hash bucket */
		Col [prev_c].shared4.hash_next = Col [c].shared4.hash_next ;
	    }
	}

	/* === Empty this hash bucket ======================================= */

	if (head_column > EMPTY)
	{
	    /* corresponding degree list "hash" is not empty */
	    Col [head_column].shared3.headhash = EMPTY ;
	}
	else
	{
	    /* corresponding degree list "hash" is empty */
	    head [hash] = EMPTY ;
	}
    }
}


/* ========================================================================== */
/* === garbage_collection =================================================== */
/* ========================================================================== */

/*
    Defragments and compacts columns and rows in the workspace A.  Used when
    all avaliable memory has been used while performing row merging.  Returns
    the index of the first free position in A, after garbage collection.  The
    time taken by this routine is linear is the size of the array A, which is
    itself linear in the number of nonzeros in the input matrix.
    Not user-callable.
*/

PRIVATE integer garbage_collection  /* returns the new value of pfree */
(
    /* === Parameters ======================================================= */

    integer n_row,			/* number of rows */
    integer n_col,			/* number of columns */
    mbdyn_Colamd_Row Row [],		/* row info */
    mbdyn_Colamd_Col Col [],		/* column info */
    integer A [],			/* A [0 ... Alen-1] holds the matrix */
    integer *pfree			/* &A [0] ... pfree is in use */
)
{
    /* === Local variables ================================================== */

    integer *psrc ;			/* source pointer */
    integer *pdest ;		/* destination pointer */
    integer j ;			/* counter */
    integer r ;			/* a row index */
    integer c ;			/* a column index */
    integer length ;		/* length of a row or column */

#ifndef NDEBUG
    integer debug_rows ;
    DEBUG2 (("Defrag..\n")) ;
    for (psrc = &A[0] ; psrc < pfree ; psrc++) ASSERT (*psrc >= 0) ;
    debug_rows = 0 ;
#endif /* NDEBUG */

    /* === Defragment the columns =========================================== */

    pdest = &A[0] ;
    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    psrc = &A [Col [c].start] ;

	    /* move and compact the column */
	    ASSERT (pdest <= psrc) ;
	    Col [c].start = (int) (pdest - &A [0]) ;
	    length = Col [c].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		r = *psrc++ ;
		if (ROW_IS_ALIVE (r))
		{
		    *pdest++ = r ;
		}
	    }
	    Col [c].length = (int) (pdest - &A [Col [c].start]) ;
	}
    }

    /* === Prepare to defragment the rows =================================== */

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    if (Row [r].length == 0)
	    {
		/* this row is of zero length.  cannot compact it, so kill it */
		DEBUG3 (("Defrag row kill\n")) ;
		KILL_ROW (r) ;
	    }
	    else
	    {
		/* save first column index in Row [r].shared2.first_column */
		psrc = &A [Row [r].start] ;
		Row [r].shared2.first_column = *psrc ;
		ASSERT (ROW_IS_ALIVE (r)) ;
		/* flag the start of the row with the one's complement of row */
		*psrc = ONES_COMPLEMENT (r) ;

#ifndef NDEBUG
		debug_rows++ ;
#endif /* NDEBUG */

	    }
	}
    }

    /* === Defragment the rows ============================================== */

    psrc = pdest ;
    while (psrc < pfree)
    {
	/* find a negative number ... the start of a row */
	if (*psrc++ < 0)
	{
	    psrc-- ;
	    /* get the row index */
	    r = ONES_COMPLEMENT (*psrc) ;
	    ASSERT (r >= 0 && r < n_row) ;
	    /* restore first column index */
	    *psrc = Row [r].shared2.first_column ;
	    ASSERT (ROW_IS_ALIVE (r)) ;

	    /* move and compact the row */
	    ASSERT (pdest <= psrc) ;
	    Row [r].start = (int) (pdest - &A [0]) ;
	    length = Row [r].length ;
	    for (j = 0 ; j < length ; j++)
	    {
		c = *psrc++ ;
		if (COL_IS_ALIVE (c))
		{
		    *pdest++ = c ;
		}
	    }
	    Row [r].length = (int) (pdest - &A [Row [r].start]) ;

#ifndef NDEBUG
	    debug_rows-- ;
#endif /* NDEBUG */

	}
    }
    /* ensure we found all the rows */
    ASSERT (debug_rows == 0) ;

    /* === Return the new value of pfree ==================================== */

    return ((int) (pdest - &A [0])) ;
}


/* ========================================================================== */
/* === clear_mark =========================================================== */
/* ========================================================================== */

/*
    Clears the Row [].shared2.mark array, and returns the new tag_mark.
    Return value is the new tag_mark.  Not user-callable.
*/

PRIVATE integer clear_mark	/* return the new value for tag_mark */
(
    /* === Parameters ======================================================= */

    integer n_row,		/* number of rows in A */
    mbdyn_Colamd_Row Row []	/* Row [0 ... n_row-1].shared2.mark is set to zero */
)
{
    /* === Local variables ================================================== */

    integer r ;

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    Row [r].shared2.mark = 0 ;
	}
    }
    return (1) ;
}


/* ========================================================================== */
/* === print_report ========================================================= */
/* ========================================================================== */

PRIVATE void print_report
(
    char *method,
    integer stats [COLAMD_STATS]
)
{

    integer i1, i2, i3 ;

    if (!stats)
    {
    	PRINTF ("%s: No statistics available.\n", method) ;
	return ;
    }

    i1 = stats [COLAMD_INFO1] ;
    i2 = stats [COLAMD_INFO2] ;
    i3 = stats [COLAMD_INFO3] ;

    if (stats [COLAMD_STATUS] >= 0)
    {
    	PRINTF ("%s: OK.  ", method) ;
    }
    else
    {
    	PRINTF ("%s: ERROR.  ", method) ;
    }

    switch (stats [COLAMD_STATUS])
    {

	case COLAMD_OK_BUT_JUMBLED:

	    PRINTF ("Matrix has unsorted or duplicate row indices.\n") ;

	    PRINTF ("%s: number of duplicate or out-of-order row indices: %d\n",
	    method, i3) ;

	    PRINTF ("%s: last seen duplicate or out-of-order row index:   %d\n",
	    method, INDEX (i2)) ;

	    PRINTF ("%s: last seen in column:                             %d",
	    method, INDEX (i1)) ;

	    /* no break - fall through to next case instead */

	case COLAMD_OK:

	    PRINTF ("\n") ;

 	    PRINTF ("%s: number of dense or empty rows ignored:           %d\n",
	    method, stats [COLAMD_DENSE_ROW]) ;

	    PRINTF ("%s: number of dense or empty columns ignored:        %d\n",
	    method, stats [COLAMD_DENSE_COL]) ;

	    PRINTF ("%s: number of garbage collections performed:         %d\n",
	    method, stats [COLAMD_DEFRAG_COUNT]) ;
	    break ;

	case COLAMD_ERROR_A_not_present:

	    PRINTF ("Array A (row indices of matrix) not present.\n") ;
	    break ;

	case COLAMD_ERROR_p_not_present:

	    PRINTF ("Array p (column pointers for matrix) not present.\n") ;
	    break ;

	case COLAMD_ERROR_nrow_negative:

	    PRINTF ("Invalid number of rows (%d).\n", i1) ;
	    break ;

	case COLAMD_ERROR_ncol_negative:

	    PRINTF ("Invalid number of columns (%d).\n", i1) ;
	    break ;

	case COLAMD_ERROR_nnz_negative:

	    PRINTF ("Invalid number of nonzero entries (%d).\n", i1) ;
	    break ;

	case COLAMD_ERROR_p0_nonzero:

	    PRINTF ("Invalid column pointer, p [0] = %d, must be zero.\n", i1) ;
	    break ;

	case COLAMD_ERROR_A_too_small:

	    PRINTF ("Array A too small.\n") ;
	    PRINTF ("        Need Alen >= %d, but given only Alen = %d.\n",
	    i1, i2) ;
	    break ;

	case COLAMD_ERROR_col_length_negative:

	    PRINTF
	    ("Column %d has a negative number of nonzero entries (%d).\n",
	    INDEX (i1), i2) ;
	    break ;

	case COLAMD_ERROR_row_index_out_of_bounds:

	    PRINTF
	    ("Row index (row %d) out of bounds (%d to %d) in column %d.\n",
	    INDEX (i2), INDEX (0), INDEX (i3-1), INDEX (i1)) ;
	    break ;

	case COLAMD_ERROR_out_of_memory:

	    PRINTF ("Out of memory.\n") ;
	    break ;

	case COLAMD_ERROR_internal_error:

	    /* if this happens, there is a bug in the code */
	    PRINTF
	    ("Internal error! Please contact authors (davis@cise.ufl.edu).\n") ;
	    break ;
    }
}




/* ========================================================================== */
/* === colamd debugging routines ============================================ */
/* ========================================================================== */

/* When debugging is disabled, the remainder of this file is ignored. */

#ifndef NDEBUG


/* ========================================================================== */
/* === debug_structures ===================================================== */
/* ========================================================================== */

/*
    At this point, all empty rows and columns are dead.  All live columns
    are "clean" (containing no dead rows) and simplicial (no supercolumns
    yet).  Rows may contain dead columns, but all live rows contain at
    least one live column.
*/

PRIVATE void debug_structures
(
    /* === Parameters ======================================================= */

    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A [],
    integer n_col2
)
{
    /* === Local variables ================================================== */

    integer i ;
    integer c ;
    integer *cp ;
    integer *cp_end ;
    integer len ;
    integer score ;
    integer r ;
    integer *rp ;
    integer *rp_end ;
    integer deg ;

    /* === Check A, Row, and Col ============================================ */

    for (c = 0 ; c < n_col ; c++)
    {
	if (COL_IS_ALIVE (c))
	{
	    len = Col [c].length ;
	    score = Col [c].shared2.score ;
	    DEBUG4 (("initial live col %5d %5d %5d\n", c, len, score)) ;
	    ASSERT (len > 0) ;
	    ASSERT (score >= 0) ;
	    ASSERT (Col [c].shared1.thickness == 1) ;
	    cp = &A [Col [c].start] ;
	    cp_end = cp + len ;
	    while (cp < cp_end)
	    {
		r = *cp++ ;
		ASSERT (ROW_IS_ALIVE (r)) ;
	    }
	}
	else
	{
	    i = Col [c].shared2.order ;
	    ASSERT (i >= n_col2 && i < n_col) ;
	}
    }

    for (r = 0 ; r < n_row ; r++)
    {
	if (ROW_IS_ALIVE (r))
	{
	    i = 0 ;
	    len = Row [r].length ;
	    deg = Row [r].shared1.degree ;
	    ASSERT (len > 0) ;
	    ASSERT (deg > 0) ;
	    rp = &A [Row [r].start] ;
	    rp_end = rp + len ;
	    while (rp < rp_end)
	    {
		c = *rp++ ;
		if (COL_IS_ALIVE (c))
		{
		    i++ ;
		}
	    }
	    ASSERT (i > 0) ;
	}
    }
}


/* ========================================================================== */
/* === debug_deg_lists ====================================================== */
/* ========================================================================== */

/*
    Prints the contents of the degree lists.  Counts the number of columns
    in the degree list and compares it to the total it should have.  Also
    checks the row degrees.
*/

PRIVATE void debug_deg_lists
(
    /* === Parameters ======================================================= */

    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer head [],
    integer min_score,
    integer should,
    integer max_deg
)
{
    /* === Local variables ================================================== */

    integer deg ;
    integer col ;
    integer have ;
    integer row ;

    /* === Check the degree lists =========================================== */

    if (n_col > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    have = 0 ;
    DEBUG4 (("Degree lists: %d\n", min_score)) ;
    for (deg = 0 ; deg <= n_col ; deg++)
    {
	col = head [deg] ;
	if (col == EMPTY)
	{
	    continue ;
	}
	DEBUG4 (("%d:", deg)) ;
	while (col != EMPTY)
	{
	    DEBUG4 ((" %d", col)) ;
	    have += Col [col].shared1.thickness ;
	    ASSERT (COL_IS_ALIVE (col)) ;
	    col = Col [col].shared4.degree_next ;
	}
	DEBUG4 (("\n")) ;
    }
    DEBUG4 (("should %d have %d\n", should, have)) ;
    ASSERT (should == have) ;

    /* === Check the row degrees ============================================ */

    if (n_row > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    for (row = 0 ; row < n_row ; row++)
    {
	if (ROW_IS_ALIVE (row))
	{
	    ASSERT (Row [row].shared1.degree <= max_deg) ;
	}
    }
}


/* ========================================================================== */
/* === debug_mark =========================================================== */
/* ========================================================================== */

/*
    Ensures that the tag_mark is less that the maximum and also ensures that
    each entry in the mark array is less than the tag mark.
*/

PRIVATE void debug_mark
(
    /* === Parameters ======================================================= */

    integer n_row,
    mbdyn_Colamd_Row Row [],
    integer tag_mark,
    integer max_mark
)
{
    /* === Local variables ================================================== */

    integer r ;

    /* === Check the Row marks ============================================== */

    ASSERT (tag_mark > 0 && tag_mark <= max_mark) ;
    if (n_row > 10000 && colamd_debug <= 0)
    {
	return ;
    }
    for (r = 0 ; r < n_row ; r++)
    {
	ASSERT (Row [r].shared2.mark < tag_mark) ;
    }
}


/* ========================================================================== */
/* === debug_matrix ========================================================= */
/* ========================================================================== */

/*
    Prints out the contents of the columns and the rows.
*/

PRIVATE void debug_matrix
(
    /* === Parameters ======================================================= */

    integer n_row,
    integer n_col,
    mbdyn_Colamd_Row Row [],
    mbdyn_Colamd_Col Col [],
    integer A []
)
{
    /* === Local variables ================================================== */

    integer r ;
    integer c ;
    integer *rp ;
    integer *rp_end ;
    integer *cp ;
    integer *cp_end ;

    /* === Dump the rows and columns of the matrix ========================== */

    if (colamd_debug < 3)
    {
	return ;
    }
    DEBUG3 (("DUMP MATRIX:\n")) ;
    for (r = 0 ; r < n_row ; r++)
    {
	DEBUG3 (("Row %d alive? %d\n", r, ROW_IS_ALIVE (r))) ;
	if (ROW_IS_DEAD (r))
	{
	    continue ;
	}
	DEBUG3 (("start %d length %d degree %d\n",
		Row [r].start, Row [r].length, Row [r].shared1.degree)) ;
	rp = &A [Row [r].start] ;
	rp_end = rp + Row [r].length ;
	while (rp < rp_end)
	{
	    c = *rp++ ;
	    DEBUG4 (("	%d col %d\n", COL_IS_ALIVE (c), c)) ;
	}
    }

    for (c = 0 ; c < n_col ; c++)
    {
	DEBUG3 (("Col %d alive? %d\n", c, COL_IS_ALIVE (c))) ;
	if (COL_IS_DEAD (c))
	{
	    continue ;
	}
	DEBUG3 (("start %d length %d shared1 %d shared2 %d\n",
		Col [c].start, Col [c].length,
		Col [c].shared1.thickness, Col [c].shared2.score)) ;
	cp = &A [Col [c].start] ;
	cp_end = cp + Col [c].length ;
	while (cp < cp_end)
	{
	    r = *cp++ ;
	    DEBUG4 (("	%d row %d\n", ROW_IS_ALIVE (r), r)) ;
	}
    }
}

PRIVATE void colamd_get_debug
(
    char *method
)
{
    colamd_debug = 0 ;		/* no debug printing */

    /* get "D" environment variable, which gives the debug printing level */
    if (getenv ("D"))
    {
    	colamd_debug = atoi (getenv ("D")) ;
    }

    DEBUG0 (("%s: debug version, D = %d (THIS WILL BE SLOW!)\n",
    	method, colamd_debug)) ;
}

#endif /* NDEBUG */

