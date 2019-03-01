/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#ifndef Y12LIB_H
#define Y12LIB_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef MBDYN_SINGLE_PRECISION
#define y12prefactor __FC_DECL__(y12mbe)
#define y12factor __FC_DECL__(y12mce)
#define y12solve __FC_DECL__(y12mde)
#else /* ! MBDYN_SINGLE_PRECISION */
#define y12prefactor __FC_DECL__(y12mbf)
#define y12factor __FC_DECL__(y12mcf)
#define y12solve __FC_DECL__(y12mdf)
#endif /* ! MBDYN_SINGLE_PRECISION */

extern int 
y12prefactor (integer *n, integer *z__, doublereal *const a, 
		   integer *const snr, integer *nn, 
		   integer *const rnr, integer *nn1, integer *ha,
		   integer *iha, doublereal *aflag, integer *iflag,
		   integer *ifail);
		   
extern int 
y12factor (integer *n, integer *z__, doublereal *const a, 
		   integer *const snr, integer *nn, 
		   integer *const rnr, integer *nn1, doublereal *pivot,
		   doublereal *const b, integer *ha, 
		   integer *iha, doublereal *aflag,
		   integer *iflag, integer *ifail);
		   
extern int 
y12solve (integer *n, doublereal *const a, 
		   integer *nn, doublereal *const b,
		   doublereal *pivot, integer *const snr, 
		   integer *ha, integer *iha,
		   integer *iflag, integer *ifail);

#ifdef __cplusplus 
}
#endif /* __cplusplus */

#endif /* Y12LIB_H */

