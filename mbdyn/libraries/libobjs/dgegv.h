/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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
 * copyright notice.
 * 
 * The original package can be found at Netlib:
 *
 *		http://www.netlib.org/lapack
 */

#ifndef DGEGV_H
#define DGEGV_H

#include <myf2c.h>

/* Subroutine */ extern int
__FC_DECL__(dgegv)(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar, 
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl, 
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork, 
	integer *info);

#endif /* DGEGV_H */

