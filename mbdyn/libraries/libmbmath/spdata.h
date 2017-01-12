/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef SPDATA_H
#define SPDATA_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <ac/f2c.h>


extern int __FC_DECL__(kd01a)(integer *iv, integer *maxkey, 
		integer *itable, integer *key);
extern int __FC_DECL__(kd01b)(integer *iv, integer *itable, 
		integer *key, integer *field, integer *ifree);
extern int __FC_DECL__(kd01c)(integer *iv, integer *itable, 
		integer *key, integer *field, integer *ifree);
extern logical __FC_DECL__(kd01h)(integer *iprime);

/* comlen kd01cm_ 12 */
extern struct ext_kd01cm_ {
    integer length, iprime, iempty;
} kd01cm_;

#define kd01cm_1 kd01cm_

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* SPDATA_H */

