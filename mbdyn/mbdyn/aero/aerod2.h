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

/* Funzione di calcolo delle forze aerodinamiche */


#ifndef AEROD_H
#define AEROD_H


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "ac/f2c.h"
   
extern int __FC_DECL__(aerod2)(doublereal* w,       /* velocita' nel sistema locale, 6 */
			       doublereal* vam,     /* dati, 6 */
			       doublereal* tng,     /* forze, 6 (Output) */
			       doublereal* outa,    /* vettore di lavoro, 20 */
			       integer* inst,       /* flag di instazionarieta', 0/1/2 */
			       doublereal* rspeed,  /* Omega */
			       integer* ipr);       /* profilo */
   
extern int __FC_DECL__(coeprd)(doublereal* da,      /* passo, delta t? */
			       doublereal* outa);   /* vettore di lavoro, 20 */

/* Da Max Lanz 2002/01/18 */
extern int __FC_DECL__(polcoe)(doublereal *x, 
			       doublereal *y, 
			       integer *n, 
			       doublereal *cof);

#ifdef __cplusplus
}
#endif /* __cplusplus */
   
#endif /* AEROD_H */
