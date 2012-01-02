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

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */


#ifndef SHAPEFNC_H
#define SHAPEFNC_H

#include "myassert.h"
#include "except.h"

enum Order { ORD_ALG = 0, ORD_D1 = 1, ORD_D2 = 2 };

/* Funzioni di interpolazione */
extern const doublereal dN3[2][3];
extern const doublereal dN3P[2][3];
extern const doublereal dN3PP[2][3];

extern doublereal 
ShapeFunc2N(doublereal d, integer iNode, enum Order Ord = ORD_ALG);
extern doublereal 
DxDcsi2N(doublereal d, const Vec3& X1, const Vec3& X2);

extern const doublereal dN2[2];
extern const doublereal dN2P[2];
extern const doublereal dN2PP[2];

extern doublereal
ShapeFunc3N(doublereal d, integer iNode, enum Order Ord = ORD_ALG);
extern doublereal
DxDcsi3N(doublereal d, const Vec3& X1, const Vec3& X2, const Vec3& X3);

#endif /* SHAPEFNC_H */

