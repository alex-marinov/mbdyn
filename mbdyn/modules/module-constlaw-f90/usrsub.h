/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef USRSUB_H
#define USRSUB_H

/* needs f2c.h or equivalent for Fortran types declaration */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern integer us1init_(integer *size, doublereal *vec, integer *err);
extern integer us1dstr_(integer *size, doublereal *vec, integer *err);
extern integer us1updt_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	const doublereal *const f,
	const doublereal *const fde,
	const doublereal *const fdep,
	integer *err);
extern integer us1aftc_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	integer *err);

extern integer us3init_(integer *size, doublereal *vec, integer *err);
extern integer us3dstr_(integer *size, doublereal *vec, integer *err);
extern integer us3updt_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	const doublereal *const f,
	const doublereal *const fde,
	const doublereal *const fdep,
	integer *err);
extern integer us3aftc_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	integer *err);

extern integer us6init_(integer *size, doublereal *vec, integer *err);
extern integer us6dstr_(integer *size, doublereal *vec, integer *err);
extern integer us6updt_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	const doublereal *const f,
	const doublereal *const fde,
	const doublereal *const fdep,
	integer *err);
extern integer us6aftc_(integer *size, doublereal *vec,
	const doublereal *const eps,
	const doublereal *const epsp,
	integer *err);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* USRSUB_H */
