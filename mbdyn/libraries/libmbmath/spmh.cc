/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2008
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "spmh.h"

SparseMatrixHandler::SparseMatrixHandler(const integer &n, const integer &nn)
:  NRows(n), NCols(nn == 0 ? n : nn), NZ(0)
{
	NO_OP;
}

SparseMatrixHandler::~SparseMatrixHandler(void)
{
	NO_OP;
}

CompactSparseMatrixHandler::CompactSparseMatrixHandler(const integer &n,
		const integer &nn,
		std::vector<doublereal>&x,
		const std::vector<integer>& i,
		const std::vector<integer>& p)
: SparseMatrixHandler(n, nn),
bMatDuplicate(false),
Ax(x),
Ai(i),
Ap(p)
{
	NZ = Ax.size();
}

CompactSparseMatrixHandler::~CompactSparseMatrixHandler()
{
	if (bMatDuplicate) {
		delete &Ax;
	}
}

/* used to sum CC matrices with identical indices */
void
CompactSparseMatrixHandler::AddUnchecked(const CompactSparseMatrixHandler& m)
{
	/* FIXME: put in stl-ish form;
	 * see if we can use something from optimized blas,
	 * e.g. ATLAS, goto or so... */

	/* checks - uncomment to enable */
#ifdef DEBUG
	ASSERT(Ax.size() == m.Ax.size());
	ASSERT(Ai.size() == m.Ai.size());
	ASSERT(Ap.size() == m.Ap.size());
	for (std::vector<doublereal>::size_type i = 0; i < Ai.size(); i++) {
		if (Ai[i] != m.Ai[i]) {
			silent_cerr("AddUnchecked: Ai[" << i << "] differs" << std::endl);
			throw ErrGeneric();
		}
	}

	for (std::vector<doublereal>::size_type i = 0; i < Ap.size(); i++) {
		if (Ap[i] != m.Ap[i]) {
			silent_cerr("AddUnchecked: Ap[" << i << "] differs" << std::endl);
			throw ErrGeneric();
		}
	}
#endif /* DEBUG */
		
	doublereal *d = &Ax[0], *s = &m.Ax[0];
	std::vector<doublereal>::size_type n = Ax.size();
	for (std::vector<doublereal>::size_type i = 0; i < n; i++) {
		d[i] += s[i];
	}
}

void
CompactSparseMatrixHandler::Reset(void)
{
	std::fill(Ax.begin(), Ax.end(), 0.);
}

integer
CompactSparseMatrixHandler::MakeCompressedColumnForm(doublereal *const Ax,
		integer *const Ai, integer *const Ap,
		int offset) const
{
	throw ErrGeneric();
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeCompressedColumnForm(std::vector<doublereal>& Ax,
		std::vector<integer>& Ai, std::vector<integer>& Ap,
		int offset) const
{
	throw ErrGeneric();
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeIndexForm(doublereal *const rAx,
		integer *const Arow, integer *const Acol,
		integer *const AcolSt, int offset) const
{
	throw ErrGeneric();
	return Nz();
}

integer
CompactSparseMatrixHandler::MakeIndexForm(std::vector<doublereal>& rAx,
                std::vector<integer>& Arow, std::vector<integer>& Acol,
		std::vector<integer>& AcolSt, int offset) const
{
	throw ErrGeneric();
	return Nz();
}
