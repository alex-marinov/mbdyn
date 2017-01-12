/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
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

#ifndef INTG_H
#define INTG_H

#include <iostream>

typedef int (pread_f)(void**, const char*);
typedef std::ostream& (phelp_f)(void *, std::ostream&);
typedef int (pinit_f)(void*, VectorHandler&);
typedef int (psize_f)(void*);
typedef int (pgrad_f)(void*, MatrixHandler&, const VectorHandler&, const doublereal&);
typedef int (pfunc_f)(void*, VectorHandler&, const VectorHandler&, const doublereal&);
typedef std::ostream& (pout_f)(void*, std::ostream&, const VectorHandler&, const VectorHandler&);
typedef int (pdestroy_f)(void**);

struct funcs {
	pread_f *read;
	phelp_f *help;
	pinit_f *init;
	psize_f *size;
	pgrad_f *grad;
	pfunc_f *func;
	pout_f *out;
	pdestroy_f *destroy;
};

#endif // INTG_H

