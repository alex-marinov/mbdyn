/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2015
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
neq:		is the matrix size;
a:		is used to store matrix elements; a[row][col] is the element in row
  		row and column col;
nzr and nzc:	are vectors of size neq, used to count the nonzero elements of a particular row or
		column, respectively;
ri and ci:	are neq x neq matrices used to store nonzero element indices; ri[col][i]
		(resp. ci[row][i]), with i < nzr[col] (i < nzc[row]), is the row (column)
		index of one of the nzr[col] (nzc[row]) nonzero elements in
		column col (row row); note that indices in ri[col] and ci[row] are
		not ordered;
nz:		nz[row][col] is true if the element in row row and column col is
		nonzero, false otherwise;
piv:		is a vector of size neq. 


The subroutine naivfct perform the LU factorization, naivslv the back-solve.

*/


#ifndef mthrdslv_h
#define mthrdslv_h

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

enum {
	NOPIV,
	SPRSPIV,
	FULLPIV
};

/* #define  PIVMETH  NOPIV */
#define  PIVMETH  SPRSPIV
/* #define  PIVMETH  FULLPIV */

#define	NAIVE_MASK	(0xF0000000U)
#define NAIVE_MAX	(~NAIVE_MASK)

#define NAIVE_ENULCOL	(0x10000000U)
#define NAIVE_ENOPIV	(0x20000000U)
#define NAIVE_ERANGE	(0x40000000U)


typedef integer** IMAT;
typedef doublereal** RMAT;
typedef char** NZMAT;

extern int naivfct(RMAT a, integer neq, integer *nzr, IMAT ri,
		integer *nzc, IMAT ci, NZMAT nz, 
		integer *piv, doublereal minpiv);

extern int naivslv(RMAT a, integer neq, integer *nzc, IMAT ci,
		doublereal *rhs, doublereal *sol, integer *piv);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* mthrdslv_h */

