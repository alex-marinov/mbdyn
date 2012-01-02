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

#ifndef C81DATA_H
#define C81DATA_H

#ifdef __cplusplus
#include <iostream>

extern "C" {
#endif /* __cplusplus */
extern void c81_data_destroy(c81_data* data);
extern c81_data* c81_data_get(long int jpro);
extern int c81_data_set(long int jpro, c81_data* data);
extern int c81_data_merge(unsigned ndata, const c81_data **data, const doublereal *upper_bounds,
	doublereal dCsi, doublereal dcltol, c81_data *i_data);
extern int c81_data_do_stall(c81_data *data, const doublereal dcltol);
extern doublereal c81_data_get_coef(int nm, doublereal* m, int na, doublereal* a, doublereal alpha, doublereal mach);
extern int c81_data_flip(c81_data *data);

#ifdef __cplusplus
extern int c81_data_read(std::istream& in, c81_data* data, const doublereal dcptol, int *ff);
extern int c81_data_fc511_read(std::istream& in, c81_data* data, const doublereal dcptol);
extern int c81_data_nrel_read(std::istream& in, c81_data* data, const doublereal dcptol);
extern int c81_data_read_free_format(std::istream& in, c81_data* data, const doublereal dcptol);
extern int c81_data_write(std::ostream& out, c81_data* data);
extern int c81_data_write_free_format(std::ostream& out, c81_data* data);
}
#endif /* __cplusplus */

#endif /* C81DATA_H */

