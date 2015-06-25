/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

MatrixScaleBase::MatrixScaleBase(const SolutionManager::ScaleOpt& scale)
:dCondBefore(-1.),
 dCondAfter(-1.),
 uFlags(scale.uFlags),
 bOK(true)
{

}

MatrixScaleBase::~MatrixScaleBase()
{

}

std::ostream& MatrixScaleBase::Report(std::ostream& os) const
{
	if ((uFlags & SolutionManager::SCALEF_COND_NUM)
			&& dCondBefore > 0 && dCondAfter > 0) {
		os << "cond: (" << dCondBefore << ") " << dCondAfter << std::endl;
	}

	if ((uFlags & SolutionManager::SCALEF_VERBOSE)
		|| (!bOK && (uFlags & SolutionManager::SCALEF_WARN))) {
		vReport(os);
	}

	return os;
}

VectorHandler&
MatrixScaleBase::ScaleVector(VectorHandler& v, const std::vector<doublereal>& s)
{
	if (static_cast<size_t>(v.iGetSize()) != s.size()) {
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i < v.iGetSize(); ++i) {
		v(i + 1) *= s[i];
	}

	return v;
}

void MatrixScaleBase::PrepareRows(const MatrixHandler& mh, integer& nrows)
{
	nrows = mh.iGetNumRows();

	if (nrows <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative row number");
	}

	if (rowScale.empty()) {
		rowScale.resize(nrows, 0.);
	} else {
		if (rowScale.size() != unsigned(nrows)) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS,
				"row number mismatch");
		}
		rowScale.assign(nrows, 0.);
	}
}

void MatrixScaleBase::PrepareCols(const MatrixHandler& mh, integer& ncols)
{
	ncols = mh.iGetNumCols();

	if (ncols <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative column number");
	}

	if (colScale.empty()) {
		colScale.resize(ncols, 0.);
	} else {
		if (colScale.size() != unsigned(ncols)) {
			// error
			throw ErrGeneric(MBDYN_EXCEPT_ARGS,
				"column number mismatch");
		}
		colScale.assign(ncols, 0.);
	}
}
