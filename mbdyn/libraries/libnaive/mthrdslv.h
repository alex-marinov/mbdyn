/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2004-2004
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

#ifndef mthrdslv_h
#define mthrdslv_h

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define  NOPIV    1
#define  SPRSPIV  2
#define  FULLPIV  3

/* #define  PIVMETH  NOPIV */
#define  PIVMETH  SPRSPIV
/* #define  PIVMETH  FULLPIV */

#define ENULCOL  0x10000000
#define ENOPIV   0x20000000

#define HIGH     0x80000000
#define LOW      0x7FFFFFFF

typedef integer** IMAT;
typedef doublereal** RMAT;

extern int naivfct(RMAT a, integer neq, integer *nzr, IMAT ri,
		integer *nzc, IMAT ci, integer *piv, doublereal minpiv);

extern void naivslv(RMAT a, integer neq, integer *nzc, IMAT ci,
		doublereal *rhs, integer *piv);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* mthrdslv_h */

