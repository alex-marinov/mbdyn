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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <sstream>

#include "dgeequ.h"
#ifdef USE_LAPACK
#include "ac/lapack.h"
#endif // USE_LAPACK

// helper for dgeequ
void
dgeequ_prepare(const MatrixHandler& mh,
	std::vector<doublereal>& r, std::vector<doublereal>& c,
	integer& nrows, integer& ncols)
{
	nrows = mh.iGetNumRows();
	ncols = mh.iGetNumCols();

	if (nrows <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative row number in dgeequ_prepare()");
	}

	if (ncols <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative column number in dgeequ_prepare()");
	}

	if (r.empty()) {
		r.resize(nrows, 0.);
	} else {
		if (r.size() != unsigned(nrows)) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS,
				"row number mismatch in dgeequ_prepare()");
		}
		r.assign(nrows, 0.);
	}

	if (c.empty()) {
		c.resize(ncols, 0.);
	} else {
		if (c.size() != unsigned(ncols)) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS,
				"column number mismatch in dgeequ_prepare()");
		}
		c.assign(ncols, 0.);
	}

}

#ifdef USE_LAPACK
// computes scaling factors for a full matrix handler
// uses lapack's dgeequ
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
		std::ostringstream os;
		os << "LAPACK's dgeequ() failed: INFO=" << INFO;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, os.str());
	}
}
#else // !USE_LAPACK
void
dgeequ(const FullMatrixHandler& mh, std::vector<doublereal>& r, std::vector<doublereal>& c,
	doublereal& rowcnd, doublereal& colcnd, doublereal& amax)
{
	dgeequ<FullMatrixHandler, FullMatrixHandler::const_iterator>(mh, r, c, rowcnd, colcnd, amax);
}
#endif // !USE_LAPACK

// scales matrix for full matrix handler, in place
FullMatrixHandler&
dgeequ_scale(FullMatrixHandler& mh, std::vector<doublereal>& r, std::vector<doublereal>& c)
{
	integer nrows = mh.iGetNumRows();
	integer ncols = mh.iGetNumCols();

	doublereal *pd = mh.pdGetMat();

	for (integer ic = 0; ic < ncols; ic++) {
		for (integer ir = 0; ir < nrows; ir++) {
			pd[ir] *= r[ir]*c[ic];
		}
		pd += ncols;
	}

	return mh;
}

// scales vector, in place
// caller needs to guarantee the length of s is at least N
void
dgeequ_scale(integer N, doublereal *v_out, doublereal *v_in, doublereal *s)
{
	for (; N-- > 0; ) {
		v_out[N] = v_in[N]*s[N];
	}
}

// scales vector handler, in place
// caller needs to guarantee the length of s is at least N
VectorHandler&
dgeequ_scale(VectorHandler& v, doublereal *s)
{
	dgeequ_scale(v.iGetSize(), v.pdGetVec(), v.pdGetVec(), s);

	return v;
}

