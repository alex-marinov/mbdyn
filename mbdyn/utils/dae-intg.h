/*
 * MBDyn (C) is a multibody analysis code.
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

#ifndef DAE_INTG_H
#define DAE_INTG_H

#include <ac/iostream>

typedef int (*pread_t)(void**, const char*);
typedef int (*pinit_t)(void*, VectorHandler&);
typedef int (*psize_t)(void*);
typedef int (*pgrad_t)(void*, MatrixHandler&, MatrixHandler&, 
		     const VectorHandler&, const doublereal&);
typedef int (*pfunc_t)(void*, VectorHandler&,
		     const VectorHandler&, const doublereal&);
typedef std::ostream& (*pout_t)(void*, std::ostream&,
			 const VectorHandler&, const VectorHandler&);
typedef int (*pdestroy_t)(void**);

typedef struct _funcs {
	pread_t read;
	pinit_t init;
	psize_t size;
	pgrad_t grad;
	pfunc_t func;
	pout_t out;
	pdestroy_t destroy;
} funcs;

#endif /* DAE_INTG_H */

