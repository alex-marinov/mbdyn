/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

#include <string.h>	/* for memset() */
#include <iostream>
#include <iomanip>

#include <fullmh.h>
#include <submat.h>

/* FullMatrixHandler - begin */

void
FullMatrixHandler::Reset(void)
{
#ifdef DEBUG
	IsValid();
#endif /* DEBUG */

	for (integer c = iNumCols; c > 0; c--) {
		for (integer r = iNumRows; r > 0; r--) {
			ppdColsm1[c][r] = 0.;
		}
	}
}

void
FullMatrixHandler::CreateColRow(integer iNR, integer iNC)
{
	ASSERT(iNC > 0);
	ASSERT(iNC > 0);
	ASSERT(pdRawm1);
	ASSERT(ppdColsm1);
	
	doublereal *pd = pdRawm1 + iNC*iNR;

	for (integer c = iNC; c > 0; c--) {
		pd -= iNR;
		ppdColsm1[c] = pd;
	}

	ASSERT(pd == pdRawm1);
}

FullMatrixHandler::FullMatrixHandler(doublereal* pd, doublereal** ppd,
		integer iSize, integer iNR, integer iNC,
		integer iMaxC)
: bOwnsMemory(false),
iNumRows(iNR), iNumCols(iNC), iRawSize(iSize), iMaxCols(iMaxC),
pdRaw(pd), pdRawm1(0),
ppdCols(ppd), ppdColsm1(0), m_end(*this, true)
{
	if (iMaxCols == 0) {
		iMaxCols = iNumCols;
	}

	ASSERT(pd);
	pdRawm1 = pd - 1;
	ASSERT(ppd);
	ppdColsm1 = ppd - 1;

	CreateColRow(iNumRows, iNumCols);

#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
}

/* costruttore che si alloca la memoria */
FullMatrixHandler::FullMatrixHandler(integer iNR, integer iNC)
: bOwnsMemory(true),
iNumRows(iNR), iNumCols(iNC ? iNC : iNumRows),
iRawSize(iNumRows*iNumCols), iMaxCols(iNumCols),
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL), m_end(*this, true)
{
	ASSERT(iNumRows > 0);
	ASSERT(iNumCols > 0);
	Resize(iNumRows, iNumCols);

#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
}


/* costruttore che non fa nulla */
FullMatrixHandler::FullMatrixHandler(void)
: bOwnsMemory(true),
iNumRows(0), iNumCols(0), iRawSize(0), iMaxCols(0),
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL), m_end(*this, true)
{
	NO_OP;
}


FullMatrixHandler::~FullMatrixHandler(void)
{
	Detach();
}

/* ridimensiona la matrice (se possiede la memoria) */
void
FullMatrixHandler::Resize(integer iNewRows, integer iNewCols)
{
	if (iNewRows < 0 || iNewCols < 0) {
		silent_cerr("Negative size!" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iSize = iNewRows*iNewCols;
	if (bOwnsMemory) {
		if (pdRaw != NULL) {
			if (iSize > iRawSize) {
				/* crea */
				doublereal* pd = NULL;
				SAFENEWARR(pd, doublereal, iSize);
				for (integer i = iRawSize; i-- > 0; ) {
					pd[i] = pdRaw[i];
				}
				iRawSize = iSize;
				SAFEDELETEARR(pdRaw);
				pdRaw = pd;
				pdRawm1 = pd - 1;
			}

			if (iNewCols > iMaxCols) {
				/* crea */
				SAFEDELETEARR(ppdCols);
				SAFENEWARR(ppdCols, doublereal*, iNewCols);
				iMaxCols = iNewCols;
				ppdColsm1 = ppdCols-1;
				CreateColRow(iNewRows, iNewCols);
			}

			/* aggiorna iNumCols, iNumRows */
			if (iNewRows != iNumRows || iNewCols != iNumCols) {
				CreateColRow(iNewRows, iNewCols);
			}

			iNumRows = iNewRows;
			iNumCols = iNewCols;

		} else {
			/* crea tutto */
			iRawSize = iSize;
			iMaxCols = iNewCols;
			iNumRows = iNewRows;
			iNumCols = iNewCols;
			SAFENEWARR(pdRaw, doublereal, iRawSize);
			pdRawm1 = pdRaw-1;
			SAFENEWARR(ppdCols, doublereal*, iNumCols);
			ppdColsm1 = ppdCols-1;
			CreateColRow(iNumRows, iNumCols);
		}

	} else {
		if (pdRaw != NULL) {
			if (iSize > iRawSize || iNewCols > iMaxCols) {
				if (iSize > iRawSize) {
					/* errore */
					silent_cerr("Can't resize to "
						<< iSize
						<< ": larger than max size "
						<< iRawSize << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);

				} else if (iNewCols > iMaxCols) {
					silent_cerr("Can't resize to "
						<< iNewCols
						<< " cols: larger than max "
						"num cols "
						<< iMaxCols << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else {
				/* aggiorna */
				if (iNewRows != iNumRows ||
						iNewCols != iNumCols) {
					CreateColRow(iNewRows, iNewCols);
				}
				iNumRows = iNewRows;
				iNumCols = iNewCols;
			}

		} else {
			/* errore */
			silent_cerr("internal error!" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}
}


/* si stacca dalla memoria a cui e' associato */
void
FullMatrixHandler::Detach(void)
{
	if (bOwnsMemory) {
		if (pdRaw != NULL) {
			SAFEDELETEARR(pdRaw);
		}
		if (ppdCols != NULL) {
			SAFEDELETEARR(ppdCols);
		}
	}
	iRawSize = iMaxCols = iNumRows = iNumCols = 0;
	pdRaw = pdRawm1 = NULL;
	ppdCols = ppdColsm1 = NULL;

	/* so Resize() can be safely invoked */
	bOwnsMemory = true;
}


/* Attacca un nuovo array, con n. righe, n. colonne e dim. massima;
 * se assente, assunta = nrighe*ncolonne */
void
FullMatrixHandler::Attach(integer iNewRows, integer iNewCols,
		doublereal* pd, doublereal** ppd,
		integer iMSize, integer iMaxC)
{
	if (bOwnsMemory || pdRaw != NULL) {
		Detach();
		bOwnsMemory = false;
	}

	iNumRows = iNewRows;
	iNumCols = iNewCols;
	if (iMSize > 0) {
		if (iMSize < iNumRows*iNumCols) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		iRawSize = iMSize;

	} else {
		iRawSize = iNumRows*iNumCols;
	}

	if (iMaxC > 0) {
		if (iMaxC < iNumCols) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		iMaxCols = iMaxC;

	} else {
		iMaxCols = iNumCols;
	}
	pdRaw = pd;
	pdRawm1 = pdRaw - 1;
	ppdCols = ppd;
	ppdColsm1 = ppdCols - 1;

#ifdef DEBUG
	IsValid();
#endif /* DEBUG */
}

#ifdef DEBUG
/* Usata per il debug */
void
FullMatrixHandler::IsValid(void) const
{
	ASSERT(pdRaw != NULL);
	ASSERT(pdRawm1 != NULL);
	ASSERT(ppdCols != NULL);
	ASSERT(ppdColsm1 != NULL);
	ASSERT(iNumRows > 0);
	ASSERT(iNumCols > 0);
	ASSERT(iRawSize >= iNumRows*iNumCols);

#ifdef DEBUG_MEMMANAGER
	ASSERT(defaultMemoryManager.fIsValid(pdRaw, iRawSize*sizeof(doublereal)));
	ASSERT(defaultMemoryManager.fIsValid(ppdCols, iNumCols*sizeof(doublereal*)));
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */


std::ostream&
operator << (std::ostream& out, const FullMatrixHandler& m)
{
#ifdef HAVE_FMTFLAGS_IN_IOS
	std::ios::fmtflags oldbits = out.setf(std::ios::scientific);
#else /* !HAVE_FMTFLAGS_IN_IOS */
	long oldbits = out.setf(ios::scientific);
#endif /* !HAVE_FMTFLAGS_IN_IOS */

	out << "<FullMatrixHandler> n. rows: " << m.iNumRows
		<< "; n. cols: " << m.iNumCols << std::endl;
	for (int i = 1; i <= m.iNumRows; i++) {
		out << "Row " << std::setw(8) << i;
		for (int j = 1; j <= m.iNumCols; j++) {
			out << std::setw(10) << std::setprecision(2)
				<< m.ppdColsm1[j][i];
		}
		out << std::endl;
	}

	out.flags(oldbits);
	return out;
}


/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
FullMatrixHandler::operator += (const SubMatrixHandler& SubMH)
{
	return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
FullMatrixHandler::operator -= (const SubMatrixHandler& SubMH)
{
	return SubMH.SubFrom(*this);
}

/* Esegue il prodotto tra due matrici e se lo memorizza */
void
FullMatrixHandler::MatMul(const FullMatrixHandler& m1,
		const FullMatrixHandler& m2)
{
	ASSERT(m1.iGetNumRows() == iGetNumRows());
	ASSERT(m2.iGetNumRows() == m1.iGetNumCols());
	ASSERT(m2.iGetNumCols() == iGetNumCols());

	integer iN = m1.iGetNumCols();

	for (integer iRow = 1; iRow <= iNumRows; iRow++) {
		for (integer iCol = 1; iCol <=iNumCols; iCol++) {
			ppdColsm1[iCol][iRow] = 0.;
			for (integer iK = 1; iK <= iN; iK++) {
				ppdColsm1[iCol][iRow] += m1(iRow,iK)*m2(iK,iCol);
			}
		}
	}
}

void
FullMatrixHandler::const_iterator::reset(bool is_end)
{
	if (is_end) {
		elem.iRow = m.iNumRows;
		elem.iCol = m.iNumCols;

	} else {
		i_idx = 0;
		elem.iRow = 0;
		elem.iCol = 0;
		elem.dCoef = m.pdRaw[i_idx];
	}
}

FullMatrixHandler::const_iterator::const_iterator(const FullMatrixHandler& m, bool is_end)
: m(m)
{
	reset(is_end);
}

FullMatrixHandler::const_iterator::~const_iterator(void)
{
	NO_OP;
}

const FullMatrixHandler::const_iterator&
FullMatrixHandler::const_iterator::operator ++ (void) const
{
#if 0
	// NOTE: this version only iterates on non-zero entries
	do {
		++i_idx;
		if (++elem.iRow == m.iNumRows) {
			if (++elem.iCol == m.iNumCols) {
				return *this;
			}
			elem.iRow = 0;
		}
	} while (m.pdRaw[i_idx] == 0.);
#else
	// NOTE: this version iterates on all coefficients
	++i_idx;
	if (++elem.iRow == m.iNumRows) {
		if (++elem.iCol == m.iNumCols) {
			return *this;
		}
		elem.iRow = 0;
	}
#endif

	elem.dCoef = m.pdRaw[i_idx];

	return *this;
}

const SparseMatrixHandler::SparseMatrixElement *
FullMatrixHandler::const_iterator::operator -> (void)
{
	return &elem;
}

const SparseMatrixHandler::SparseMatrixElement&
FullMatrixHandler::const_iterator::operator * (void)
{
	return elem;
}

bool
FullMatrixHandler::const_iterator::operator == (const FullMatrixHandler::const_iterator& op) const
{
	return elem == op.elem;
}

bool
FullMatrixHandler::const_iterator::operator != (const FullMatrixHandler::const_iterator& op) const
{
	return elem != op.elem;
}

/* FullMatrixHandler - end */

