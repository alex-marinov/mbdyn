/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dgeequ.h"
#ifdef USE_LAPACK
#include "ac/lapack.h"
#endif // USE_LAPACK

void
dgeequ_prepare(const MatrixHandler& mh,
	std::vector<doublereal>& r, std::vector<doublereal>& c,
	integer& nrows, integer& ncols)
{
	nrows = mh.iGetNumRows();
	ncols = mh.iGetNumCols();

	if (nrows <= 0 || ncols <= 0) {
		// error
	}

	if (r.empty()) {
		r.resize(nrows, 0.);
	} else {
		if (r.size() != unsigned(nrows)) {
			// error
		}
		r.assign(nrows, 0.);
	}

	if (c.empty()) {
		c.resize(ncols, 0.);
	} else {
		if (c.size() != unsigned(ncols)) {
			// error
		}
		c.assign(ncols, 0.);
	}

}

void
dgeequ(const FullMatrixHandler& mh, std::vector<doublereal>& r, std::vector<doublereal>& c,
	doublereal& rowcnd, doublereal& colcnd, doublereal& amax)
{
	integer nrows, ncols;
	dgeequ_prepare(mh, r, c, nrows, ncols);

	integer INFO;
	__FC_DECL__(dgeequ)(
		&nrows,
		&ncols,
		mh.pdGetMat(),
		&nrows,
		&r[0],
		&c[0],
		&rowcnd,
		&colcnd,
		&amax,
		&INFO);
	if (INFO != 0) {
		// error
	}
}

