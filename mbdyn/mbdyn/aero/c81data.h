/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#ifndef C81DATA_H
#define C81DATA_H

#if defined(HAVE_IOSTREAM)
#include <iostream>
#elif defined(HAVE_IOSTREAM_H)
#include <iostream.h>
#endif

extern "C" {
extern c81_data* get_c81_data(long int jpro);
extern int set_c81_data(long int jpro, c81_data* data);
extern int read_c81_data(std::istream& in, c81_data* data);
extern int write_c81_data(std::ostream& out, c81_data* data);
extern double get_c81_coef(int nm, double* m, int na, double* a, double alpha, double mach);
}

#endif /* C81DATA_H */

