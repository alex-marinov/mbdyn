/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2004
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

doublereal SparseMatrixHandler::zero = 0.;

SparseMatrixHandler::SparseMatrixHandler(const int &n, const int &nn)
:  NRows(n), NCols(nn == 0 ? n : nn), NZ(0)
{
	NO_OP;
}

SparseMatrixHandler::~SparseMatrixHandler(void)
{
	NO_OP;
}

void
SparseMatrixHandler::Init(const doublereal& c)
{
	Reset(c);
}

CompactSparseMatrixHandler::CompactSparseMatrixHandler(const int &n,
		const int &nn,
		std::vector<doublereal>&x,
		const std::vector<int>& i,
		const std::vector<int>& p)
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
	for (unsigned i = 0; i < Ai.size(); i++) {
		if (Ai[i] != m.Ai[i]) {
			std::cerr << "Ai[" << i << "]" << std::endl;
			THROW(ErrGeneric());
		}
	}

	for (unsigned i = 0; i < Ap.size(); i++) {
		if (Ap[i] != m.Ap[i]) {
			std::cerr << "Ap[" << i << "]" << std::endl;
			THROW(ErrGeneric());
		}
	}
#endif /* DEBUG */
		
	doublereal *d = &Ax[0], *s = &m.Ax[0];
	for (unsigned long i = 0; i < Ax.size(); i++) {
		d[i] += s[i];
	}
}

void
CompactSparseMatrixHandler::Reset(const doublereal &r)
{
	std::fill(Ax.begin(), Ax.end(), r);
}

