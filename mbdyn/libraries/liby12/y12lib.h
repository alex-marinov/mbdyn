/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

extern int 
__FC_DECL__(y12mbf) (integer *n, integer *z__, doublereal *a, integer *snr, 
		   integer *nn, integer *rnr, integer *nn1, integer *ha,
		   integer *iha, doublereal *aflag, integer *iflag,
		   integer *ifail);
		   
extern int 
__FC_DECL__(y12mcf) (integer *n, integer *z__, doublereal *a, integer *snr,
		   integer *nn, integer *rnr, integer *nn1, doublereal *pivot,
		   doublereal *b, integer *ha, integer *iha, doublereal *aflag,
		   integer *iflag, integer *ifail);
		   
extern int 
__FC_DECL__(y12mdf) (integer *n, doublereal *a, integer *nn, doublereal *b,
		   doublereal *pivot, integer *snr, integer *ha, integer *iha,
		   integer *iflag, integer *ifail);

#ifdef __cplusplus 
}
#endif /* __cplusplus */

#endif /* Y12LIB_H */

