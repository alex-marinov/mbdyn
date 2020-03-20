/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
#include <cstring>	/* for memset() */
#include <iostream>
#include <iomanip>
#include <ac/lapack.h>
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
	Reset();

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

/* costruttore di copia */
FullMatrixHandler::FullMatrixHandler(const FullMatrixHandler& MH)
: bOwnsMemory(true),
iNumRows(0), iNumCols(0), iRawSize(0), iMaxCols(0),
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL), m_end(*this, true)
{
	if (MH.pdRaw != 0) {
		Resize(MH.iNumRows, MH.iNumCols);

		for (integer i = 0; i < iNumRows*iNumCols; i++) {
			pdRaw[i] = MH.pdRaw[i];
		}
	}
}

/* Overload di = usato per l'assemblaggio delle matrici */
FullMatrixHandler&
FullMatrixHandler::operator = (const FullMatrixHandler& MH)
{
	if (pdRaw == 0) {
		Resize(MH.iNumRows, MH.iNumCols);

	} else if (iNumRows != MH.iNumRows || iNumCols != MH.iNumCols) {
		silent_cerr("FullMatrixHandler::operator = (const FullMatrixHandler&): "
			"incompatible size (" << iNumRows << ", " << iNumCols << ") "
			" <> (" << MH.iNumRows << ", " << MH.iNumCols << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (integer i = 0; i < iNumRows*iNumCols; i++) {
		pdRaw[i] = MH.pdRaw[i];
	}

	return *this;
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

extern std::ostream&
Write(std::ostream& out,
	const FullMatrixHandler& m,
	const char* s, 
	const char* s2)
{
#ifdef HAVE_FMTFLAGS_IN_IOS
	std::ios::fmtflags oldbits = out.setf(std::ios::scientific);
#else /* !HAVE_FMTFLAGS_IN_IOS */
	long oldbits = out.setf(ios::scientific);
#endif /* !HAVE_FMTFLAGS_IN_IOS */

	if (s2 == NULL) {
		s2 = s;
	}

	for (int i = 1; i <= m.iNumRows; i++) {
		for (int j = 1; j <= m.iNumCols; j++) {
			out << std::setw(20) << std::setprecision(12)
				<< m.ppdColsm1[j][i] << s;
		}
		out << s2;
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

/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler&
FullMatrixHandler::operator += (const VariableSubMatrixHandler& SubMH)
{
	return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler&
FullMatrixHandler::operator -= (const VariableSubMatrixHandler& SubMH)
{
	return SubMH.SubFrom(*this);
}

/* Esegue il prodotto tra due matrici e se lo memorizza */
void
FullMatrixHandler::MatMul(const FullMatrixHandler& m1,
		const FullMatrixHandler& m2)
{
	if (m1.iGetNumRows() != iGetNumRows()
		|| m2.iGetNumRows() != m1.iGetNumCols()
		|| m2.iGetNumCols() != iGetNumCols())
	{
		silent_cerr("FullMatrixHandler::MatMul: size mismatch "
			"this(" << iGetNumRows() << ", " << iGetNumCols() << ") "
			"= m1(" << m1.iGetNumRows() << ", " << m1.iGetNumCols() << ") "
			"* m2(" << m2.iGetNumRows() << ", " << m2.iGetNumCols() << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer nc = m1.iGetNumCols();

	for (integer iRow = 1; iRow <= iNumRows; iRow++) {
		for (integer iCol = 1; iCol <= iNumCols; iCol++) {
			ppdColsm1[iCol][iRow] = 0.;
			for (integer iK = 1; iK <= nc; iK++) {
				// ppdColsm1[iCol][iRow] += m1(iRow, iK)*m2(iK, iCol);
				ppdColsm1[iCol][iRow] += m1.ppdColsm1[iK][iRow]*m2.ppdColsm1[iCol][iK];
			}
		}
	}
}

#if defined HAVE_BLAS
#if defined HAVE_CBLAS
void FullMatrixHandler::LapackMatrixOp(const char* TRANSA,
                                       const char* TRANSB,
                                       doublereal* pout,
                                       integer out_nr,
                                       integer out_nc,
                                       const doublereal* pin,
                                       integer in_nr,
                                       doublereal ALPHA,
                                       doublereal BETA) const
{
        ASSERT((strcmp("N", TRANSA) == 0 && out_nr == iGetNumRows() && in_nr == iGetNumCols()) ||
               (strcmp("T", TRANSA) == 0 && out_nr == iGetNumCols() && in_nr == iGetNumRows()));
	CBLAS_TRANSPOSE TA, TB;
	if (strcmp("N", TRANSA) == 0) TA = CblasNoTrans;
	else if (strcmp("T", TRANSA) == 0) TA = CblasTrans;
	if (strcmp("N", TRANSB) == 0) TB = CblasNoTrans;
	else if (strcmp("T", TRANSB) == 0) TB = CblasTrans;
	cblas_dgemm(CblasColMajor, TA, TB, 
                           out_nr,
                           out_nc,
                           in_nr,
                           ALPHA,
                           pdGetMat(),
                           out_nr,
                           pin,
                           in_nr,
                           BETA,
                           pout,
                           out_nr);        
        
}
#else /* HAVE_BLAS && !HAVE_CBLAS */
/*
   Note: dgemm performs

       C := alpha*op( A )*op( B ) + beta*C,

   in MatMatMul_base and MatTMatMul_base context,

       A := *this
       B := in
       C := out

       op( B ) := B

   in MatMatMul_base context,

       op( A ) := A

   in MatTMatMul_base context,

       op( A ) := A^T

   out = op( A ) * B 	requires alpha = 1, beta = 0
   out += op( A ) * B 	requires alpha = 1, beta = 1
   out -= op( A ) * B 	requires alpha = -1, beta = 1
 */

void FullMatrixHandler::LapackMatrixOp(const char* TRANSA,
                                       const char* TRANSB,
                                       doublereal* pout,
                                       integer out_nr,
                                       integer out_nc,
                                       const doublereal* pin,
                                       integer in_nr,
                                       doublereal ALPHA,
                                       doublereal BETA) const
{
        ASSERT((strcmp("N", TRANSA) == 0 && out_nr == iGetNumRows() && in_nr == iGetNumCols()) ||
               (strcmp("T", TRANSA) == 0 && out_nr == iGetNumCols() && in_nr == iGetNumRows()));
        
        __FC_DECL__(dgemm)(TRANSA,
                           TRANSB,
                           &out_nr,
                           &out_nc,
                           &in_nr,
                           &ALPHA,
                           pdGetMat(),
                           &out_nr,
                           pin,
                           &in_nr,
                           &BETA,
                           pout,
                           &out_nr);        
}
#endif /* HAVE_CBLAS */
  
void FullMatrixHandler::LapackMatrixOp(const char* TRANSA,
                                       const char* TRANSB,
                                       doublereal* pout,
                                       integer out_nr,
                                       integer out_nc,
                                       const doublereal* pin,
                                       integer in_nr,
                                       VectorOperation op) const
{        
        doublereal ALPHA, BETA;
        
        if (op == &VectorHandler::IncCoef) {
                BETA = 1.;
                ALPHA = 1.;
        } else if (op == &VectorHandler::DecCoef) {
                BETA = 1.;
                ALPHA = -1.;

        } else if (op == &VectorHandler::PutCoef) {
                BETA = 0.;
                ALPHA = 1.;
        } else {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        LapackMatrixOp(TRANSA, TRANSB, pout, out_nr, out_nc, pin, in_nr, ALPHA, BETA);
}

void FullMatrixHandler::LapackMatrixOp(const char* TRANSA,
                                       const char* TRANSB,
                                       doublereal* pout,
                                       integer out_nr,
                                       integer out_nc,
                                       const doublereal* pin,
                                       integer in_nr,
                                       MatrixOperation op) const
{        
        doublereal ALPHA, BETA;
        
        if (op == &MatrixHandler::IncCoef) {
                BETA = 1.;
                ALPHA = 1.;
        } else if (op == &MatrixHandler::DecCoef) {
                BETA = 1.;
                ALPHA = -1.;
        } else if (op == &MatrixHandler::PutCoef) {
                BETA = 0.;
                ALPHA = 1.;
        } else {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        LapackMatrixOp(TRANSA, TRANSB, pout, out_nr, out_nc, pin, in_nr, ALPHA, BETA);
}
#endif // HAVE_BLAS

MatrixHandler&
FullMatrixHandler::MatMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol,
		const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	const FullMatrixHandler *pin = dynamic_cast<const FullMatrixHandler *>(&in);
	if (pin == 0) {
		// if input matrix is not FullMatrixHandler, use generic function
		return MatrixHandler::MatMatMul_base(op, out, in);
	}

	integer out_nc = out.iGetNumCols();
	integer out_nr = out.iGetNumRows();
	integer in_nr = in.iGetNumRows();

	if (out_nr != iGetNumRows()
		|| out_nc != in.iGetNumCols()
		|| in_nr != iGetNumCols())
	{
		const char *strop;

		if (op == &MatrixHandler::IncCoef) {
			strop = "+=";
		} else if (op == &MatrixHandler::DecCoef) {
			strop = "-=";
		} else if (op == &MatrixHandler::PutCoef) {
			strop = "=";
		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		silent_cerr("FullMatrixHandler::MatMatMul_base: size mismatch "
			"out(" << out_nr << ", " << out_nc << ") "
			<< strop << " this(" << iGetNumRows() << ", " << iGetNumCols() << ") "
			"* in(" << in_nr << ", " << in.iGetNumCols() << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	FullMatrixHandler *pout = dynamic_cast<FullMatrixHandler *>(&out);
	if (pout == 0) {
		// if output matrix is not FullMatrixHandler,
		// optimize coefficient computation
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(r, k) * in(k, c)
					d += ppdColsm1[k][r]*pin->ppdColsm1[c][k];
				}

				(out.*op)(r, c, d);
			}
		}

		return out;
	}

	// if all matrices are FullMatrixHandler,
	// either use LAPACK or optimize coefficient computation and store
#ifdef HAVE_BLAS
        LapackMatrixOp("N", "N", pout->pdGetMat(), out_nr, out_nc, pin->pdGetMat(), in_nr, op);
#else // ! HAVE_BLAS
	if (op == &MatrixHandler::IncCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(r, k) * in(k, c)
					d += ppdColsm1[k][r]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] += d;
			}
		}

	} else if (op == &MatrixHandler::DecCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(r, k) * in(k, c)
					d += ppdColsm1[k][r]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] -= d;
			}
		}

	} else if (op == &MatrixHandler::PutCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(r, k) * in(k, c)
					d += ppdColsm1[k][r]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] = d;
			}
		}

	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // ! HAVE_BLAS

	return out;
}

MatrixHandler&
FullMatrixHandler::MatTMatMul_base(
	void (MatrixHandler::*op)(integer iRow, integer iCol,
		const doublereal& dCoef),
	MatrixHandler& out, const MatrixHandler& in) const
{
	const FullMatrixHandler *pin = dynamic_cast<const FullMatrixHandler *>(&in);
	if (pin == 0) {
		// if input matrix is not FullMatrixHandler, use generic function
		return MatrixHandler::MatTMatMul_base(op, out, in);
	}

	integer out_nc = out.iGetNumCols();
	integer out_nr = out.iGetNumRows();
	integer in_nr = in.iGetNumRows();

	if (out_nr != iGetNumCols()
		|| out_nc != in.iGetNumCols()
		|| in_nr != iGetNumRows())
	{
		const char *strop;

		if (op == &MatrixHandler::IncCoef) {
			strop = "+=";
		} else if (op == &MatrixHandler::DecCoef) {
			strop = "-=";
		} else if (op == &MatrixHandler::PutCoef) {
			strop = "=";
		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		silent_cerr("FullMatrixHandler::MatTMatMul_base: size mismatch "
			"out(" << out_nr << ", " << out_nc << ") "
			<< strop << " this(" << iGetNumRows() << ", " << iGetNumCols() << ")^T "
			"* in(" << in_nr << ", " << in.iGetNumCols() << ")"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	FullMatrixHandler *pout = dynamic_cast<FullMatrixHandler *>(&out);
	if (pout == 0) {
		// if output matrix is not FullMatrixHandler,
		// optimize coefficient computation
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(k, r) * in(k, c)
					d += ppdColsm1[r][k]*pin->ppdColsm1[c][k];
				}

				(out.*op)(r, c, d);
			}
		}

		return out;
	}

	// if all matrices are FullMatrixHandler,
	// either use LAPACK or optimize coefficient computation and store
#ifdef HAVE_BLAS
        LapackMatrixOp("T", "N", pout->pdGetMat(), out_nr, out_nc, pin->pdGetMat(), in_nr, op);
#else // ! HAVE_BLAS
	if (op == &MatrixHandler::IncCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(k, r) * in(k, c)
					d += ppdColsm1[r][k]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] += d;
			}
		}

	} else if (op == &MatrixHandler::DecCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(k, r) * in(k, c)
					d += ppdColsm1[r][k]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] -= d;
			}
		}

	} else if (op == &MatrixHandler::PutCoef) {
		for (integer c = 1; c <= out_nc; c++) {
			for (integer r = 1; r <= out_nr; r++) {
				doublereal d = 0.;

				for (integer k = 1; k <= in_nr; k++) {
					// this(k, r) * in(k, c)
					d += ppdColsm1[r][k]*pin->ppdColsm1[c][k];
				}

				// (out.*op)(r, c, d);
				pout->ppdColsm1[c][r] = d;
			}
		}

	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif // ! HAVE_BLAS

	return out;
}

VectorHandler&
FullMatrixHandler::MatVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
        integer out_nr = out.iGetSize();
        integer in_nr = in.iGetSize();

        if (out_nr != iGetNumRows() || in_nr != iGetNumCols()) {
                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        
        // if all matrices are FullMatrixHandler,
        // either use LAPACK or optimize coefficient computation and store
#ifdef HAVE_BLAS               
        LapackMatrixOp("N", "N", out.pdGetVec(), out_nr, 1, in.pdGetVec(), in_nr, op);
#else
	integer nr = iGetNumRows();
	integer nc = iGetNumCols();

	ASSERT(nc == in.iGetSize());
	ASSERT(nr == out.iGetSize());

	const MyVectorHandler *pin = dynamic_cast<const MyVectorHandler *>(&in);
	if (pin == 0) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ir, ic) * in(ic)
				d += ppdColsm1[ic][ir]*in(ic);
			}
			(out.*op)(ir, d);
		}

		return out;
	}

	MyVectorHandler *pout = dynamic_cast<MyVectorHandler *>(&out);
	if (pout == 0) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ir, ic) * in(ic)
				d += ppdColsm1[ic][ir]*pin->pdVecm1[ic];
			}
			(out.*op)(ir, d);
		}

		return out;
	}

	if (op == &VectorHandler::IncCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ir, ic) * in(ic)
				d += ppdColsm1[ic][ir]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] += d;
		}

	} else if (op == &VectorHandler::DecCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ir, ic) * in(ic)
				d += ppdColsm1[ic][ir]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] -= d;
		}

	} else if (op == &VectorHandler::PutCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ir, ic) * in(ic)
				d += ppdColsm1[ic][ir]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] = d;
		}

	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif
	return out;
}

VectorHandler&
FullMatrixHandler::MatTVecMul_base(
	void (VectorHandler::*op)(integer iRow, const doublereal& dCoef),
	VectorHandler& out, const VectorHandler& in) const
{
        integer out_nr = out.iGetSize();
        integer in_nr = in.iGetSize();

        if (out_nr != iGetNumCols() || in_nr != iGetNumRows()) {
                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        
#ifdef HAVE_BLAS        
        LapackMatrixOp("T", "N", out.pdGetVec(), out_nr, 1, in.pdGetVec(), in_nr, op);
#else
	integer nr = iGetNumCols();
	integer nc = iGetNumRows();

	ASSERT(nc == in.iGetSize());
	ASSERT(nr == out.iGetSize());

	const MyVectorHandler *pin = dynamic_cast<const MyVectorHandler *>(&in);
	if (pin == 0) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ic, ir) * in(ic)
				d += ppdColsm1[ir][ic]*in(ic);
			}
			(out.*op)(ir, d);
		}

		return out;
	}

	MyVectorHandler *pout = dynamic_cast<MyVectorHandler *>(&out);
	if (pout == 0) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ic, ir) * in(ic)
				d += ppdColsm1[ir][ic]*pin->pdVecm1[ic];
			}
			(out.*op)(ir, d);
		}

		return out;
	}

	if (op == &VectorHandler::IncCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ic, ir) * in(ic)
				d += ppdColsm1[ir][ic]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] += d;
		}

	} else if (op == &VectorHandler::DecCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ic, ir) * in(ic)
				d += ppdColsm1[ir][ic]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] -= d;
		}

	} else if (op == &VectorHandler::PutCoef) {
		for (integer ir = 1; ir <= nr; ir++) {
			doublereal d = 0.;
			for (integer ic = 1; ic <= nc; ic++) {
				// out(ir) ? this(ic, ir) * in(ic)
				d += ppdColsm1[ir][ic]*pin->pdVecm1[ic];
			}
			pout->pdVecm1[ir] = d;
		}

	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif
	return out;
}

void
FullMatrixHandler::Add(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] += source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::Sub(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] -= source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::Put(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] = source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::Add(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] += dCoef*source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::Sub(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] -= dCoef*source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::Put(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumRows();
	integer nc = source.iGetNumCols();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ir, ic) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] = dCoef*source.ppdColsm1[ic][ir];
		}
	}

	return;
}

void
FullMatrixHandler::AddT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] += source.ppdColsm1[ir][ic];
		}
	}

	return;
}

void
FullMatrixHandler::SubT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] -= source.ppdColsm1[ir][ic];
		}
	}

	return;
}

void
FullMatrixHandler::PutT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] = source.ppdColsm1[ir][ic];
		}
	}

	return;
}

void
FullMatrixHandler::AddT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] += dCoef*source.ppdColsm1[ir][ic];
		}
	}

	return;
}

void
FullMatrixHandler::SubT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] -= dCoef*source.ppdColsm1[ir][ic];
		}
	}

	return;
}

void
FullMatrixHandler::PutT(integer iRow, 
	integer iCol, 
	const FullMatrixHandler & source,
	const doublereal dCoef)
{
	integer nr = source.iGetNumCols();
	integer nc = source.iGetNumRows();
	iRow--;
	iCol--;

#ifdef DEBUG
	IsValid();

	ASSERT(iRow >= 0);
	ASSERT(iCol >= 0);
	ASSERT(iRow + nr <= iGetNumRows());
	ASSERT(iCol + nc <= iGetNumCols());
#endif // DEBUG
#if 0
	if (iRow < 0 || iCol < 0
		|| iRow + nr > iGetNumRows()
		|| iCol + nc > iGetNumCols())
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif

	for (integer ir = 1; ir <= nr; ir++) {
		for (integer ic = 1; ic <= nc; ic++) {
			// IncCoef(ir + iRow, ic + iCol, source(ic, ir) * dCoef);
			ppdColsm1[ic + iCol][ir + iRow] = dCoef*source.ppdColsm1[ir][ic];
		}
	}

	return;
}

/* somma un vettore di tipo Vec3 in una data posizione */
void
FullMatrixHandler::Add(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] += pdFrom[V1];

	ppdColsm1[iCol][iRow + 1] += pdFrom[V2];

	ppdColsm1[iCol][iRow + 2] += pdFrom[V3];
}

/* sottrae un vettore di tipo Vec3 in una data posizione */
void
FullMatrixHandler::Sub(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] -= pdFrom[V1];

	ppdColsm1[iCol][iRow + 1] -= pdFrom[V2];

	ppdColsm1[iCol][iRow + 2] -= pdFrom[V3];
}

/* scrive un vettore di tipo Vec3 in una data posizione */
void
FullMatrixHandler::Put(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] = pdFrom[V1];

	ppdColsm1[iCol][iRow + 1] = pdFrom[V2];

	ppdColsm1[iCol][iRow + 2] = pdFrom[V3];
}

/* somma un vettore di tipo Vec3 trasposto in una data posizione */
void
FullMatrixHandler::AddT(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] += pdFrom[V1];

	ppdColsm1[iCol + 1][iRow] += pdFrom[V2];

	ppdColsm1[iCol + 2][iRow] += pdFrom[V3];
}

/* sottrae un vettore di tipo Vec3 trasposto in una data posizione */
void
FullMatrixHandler::SubT(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] -= pdFrom[V1];

	ppdColsm1[iCol + 1][iRow] -= pdFrom[V2];

	ppdColsm1[iCol + 2][iRow] -= pdFrom[V3];
}

/* scrive un vettore di tipo Vec3 trasposto in una data posizione */
void
FullMatrixHandler::PutT(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] = pdFrom[V1];

	ppdColsm1[iCol + 1][iRow] = pdFrom[V2];

	ppdColsm1[iCol + 2][iRow] = pdFrom[V3];
}

#if 0 /* FIXME: replace original? */
/* somma un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullMatrixHandler::AddDiag(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] += pdFrom[V1];

	ppdColsm1[iCol + 1][iRow + 1] += pdFrom[V2];

	ppdColsm1[iCol + 2][iRow + 2] += pdFrom[V3];
}

/* sottrae un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullMatrixHandler::SubDiag(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] -= pdFrom[V1];

	ppdColsm1[iCol + 1][iRow + 1] -= pdFrom[V2];

	ppdColsm1[iCol + 2][iRow + 2] -= pdFrom[V3];
}

/* scrive un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullMatrixHandler::PutDiag(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol][iRow] = pdFrom[V1];

	ppdColsm1[iCol + 1][iRow + 1] = pdFrom[V2];

	ppdColsm1[iCol + 2][iRow + 2] = pdFrom[V3];
}

/* somma un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullMatrixHandler::AddCross(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol + 1][iRow] -= pdFrom[V3];
	ppdColsm1[iCol + 2][iRow] += pdFrom[V2];

	ppdColsm1[iCol][iRow + 1] += pdFrom[V3];
	ppdColsm1[iCol + 2][iRow + 1] -= pdFrom[V1];

	ppdColsm1[iCol][iRow + 2] -= pdFrom[V2];
	ppdColsm1[iCol + 1][iRow + 2] += pdFrom[V1];
}

/* sottrae un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullMatrixHandler::SubCross(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol + 1][iRow] += pdFrom[V3];
	ppdColsm1[iCol + 2][iRow] -= pdFrom[V2];

	ppdColsm1[iCol][iRow + 1] -= pdFrom[V3];
	ppdColsm1[iCol + 2][iRow + 1] += pdFrom[V1];

	ppdColsm1[iCol][iRow + 2] += pdFrom[V2];
	ppdColsm1[iCol + 1][iRow + 2] -= pdFrom[V1];
}

/* scrive un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullMatrixHandler::PutCross(integer iRow, integer iCol, const Vec3& v)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN)
	 */
	const doublereal* pdFrom = v.pGetVec();

	ppdColsm1[iCol + 1][iRow] = -pdFrom[V3];
	ppdColsm1[iCol + 2][iRow] = pdFrom[V2];

	ppdColsm1[iCol][iRow + 1] = pdFrom[V3];
	ppdColsm1[iCol + 2][iRow + 1] = -pdFrom[V1];

	ppdColsm1[iCol][iRow + 2] = -pdFrom[V2];
	ppdColsm1[iCol + 1][iRow + 2] = pdFrom[V1];
}
#endif

/* somma una matrice di tipo Mat3x3 in una data posizione */
void
FullMatrixHandler::Add(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN) e assume
	 * che la matrice Mat3x3 sia organizzata anch'essa per colonne */
	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] += pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] += pdFrom[M21];
	ppdColsm1[iCol][iRow + 2] += pdFrom[M31];

	ppdColsm1[iCol + 1][iRow] += pdFrom[M12];
	ppdColsm1[iCol + 1][iRow + 1] += pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] += pdFrom[M32];

	ppdColsm1[iCol + 2][iRow] += pdFrom[M13];
	ppdColsm1[iCol + 2][iRow + 1] += pdFrom[M23];
	ppdColsm1[iCol + 2][iRow + 2] += pdFrom[M33];

#if 0
	for (unsigned c = 0; c < 2; c++) {
		for (unsigned r = 0; r < 2; r++) {
			ppdColsm1[iCol + c][iRow + r] += pdFrom[r];
		}
		pdFrom += 3;
	}
#endif
}

void
FullMatrixHandler::AddT(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	/* Nota: assume che la matrice sia organizzata
	 * "per colonne" (stile FORTRAN) e assume
	 * che la matrice Mat3x3 sia organizzata anch'essa per colonne */
	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] += pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] += pdFrom[M12];
	ppdColsm1[iCol][iRow + 2] += pdFrom[M13];

	ppdColsm1[iCol + 1][iRow] += pdFrom[M21];
	ppdColsm1[iCol + 1][iRow + 1] += pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] += pdFrom[M23];

	ppdColsm1[iCol + 2][iRow] += pdFrom[M31];
	ppdColsm1[iCol + 2][iRow + 1] += pdFrom[M32];
	ppdColsm1[iCol + 2][iRow + 2] += pdFrom[M33];
}


/* sottrae una matrice di tipo Mat3x3 in una data posizione
 * analoga ala precedente, con il meno (per evitare temporanei ecc) */
void
FullMatrixHandler::Sub(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1
	 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] -= pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] -= pdFrom[M21];
	ppdColsm1[iCol][iRow + 2] -= pdFrom[M31];

	ppdColsm1[iCol + 1][iRow] -= pdFrom[M12];
	ppdColsm1[iCol + 1][iRow + 1] -= pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] -= pdFrom[M32];

	ppdColsm1[iCol + 2][iRow] -= pdFrom[M13];
	ppdColsm1[iCol + 2][iRow + 1] -= pdFrom[M23];
	ppdColsm1[iCol + 2][iRow + 2] -= pdFrom[M33];

#if 0
	for (unsigned c = 0; c < 2; c++) {
		for (unsigned r = 0; r < 2; r++) {
			ppdColsm1[iCol + c][iRow + r] -= pdFrom[r];
		}
		pdFrom += 3;
	}
#endif
}

void
FullMatrixHandler::SubT(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1
	 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] -= pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] -= pdFrom[M12];
	ppdColsm1[iCol][iRow + 2] -= pdFrom[M13];

	ppdColsm1[iCol + 1][iRow] -= pdFrom[M21];
	ppdColsm1[iCol + 1][iRow + 1] -= pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] -= pdFrom[M23];

	ppdColsm1[iCol + 2][iRow] -= pdFrom[M31];
	ppdColsm1[iCol + 2][iRow + 1] -= pdFrom[M32];
	ppdColsm1[iCol + 2][iRow + 2] -= pdFrom[M33];
}


/* mette una matrice di tipo Mat3x3 in una data posizione;
 * analoga alle precedenti */
void
FullMatrixHandler::Put(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1
	 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] = pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] = pdFrom[M21];
	ppdColsm1[iCol][iRow + 2] = pdFrom[M31];

	ppdColsm1[iCol + 1][iRow] = pdFrom[M12];
	ppdColsm1[iCol + 1][iRow + 1] = pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] = pdFrom[M32];

	ppdColsm1[iCol + 2][iRow] = pdFrom[M13];
	ppdColsm1[iCol + 2][iRow + 1] = pdFrom[M23];
	ppdColsm1[iCol + 2][iRow + 2] = pdFrom[M33];

#if 0
	for (unsigned c = 0; c < 2; c++) {
		for (unsigned r = 0; r < 2; r++) {
			ppdColsm1[iCol + c][iRow + r] = pdFrom[r];
		}
		pdFrom += 3;
	}
#endif
}

void
FullMatrixHandler::PutT(integer iRow, integer iCol, const Mat3x3& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1
	 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	const doublereal* pdFrom = m.pGetMat();

	ppdColsm1[iCol][iRow] = pdFrom[M11];
	ppdColsm1[iCol][iRow + 1] = pdFrom[M12];
	ppdColsm1[iCol][iRow + 2] = pdFrom[M13];

	ppdColsm1[iCol + 1][iRow] = pdFrom[M21];
	ppdColsm1[iCol + 1][iRow + 1] = pdFrom[M22];
	ppdColsm1[iCol + 1][iRow + 2] = pdFrom[M23];

	ppdColsm1[iCol + 2][iRow] = pdFrom[M31];
	ppdColsm1[iCol + 2][iRow + 1] = pdFrom[M32];
	ppdColsm1[iCol + 2][iRow + 2] = pdFrom[M33];
}


/* somma una matrice di tipo Mat3xN in una data posizione */
void
FullMatrixHandler::Add(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - m.iGetNumCols() + 1);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = 3; r > 0; r--) {
		for (integer c = m.iGetNumCols(); c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] += m(r, c);
		}
	}
}

/* sottrae una matrice di tipo Mat3xN in una data posizione */
void
FullMatrixHandler::Sub(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - m.iGetNumCols() + 1);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = 3; r > 0; r--) {
		for (integer c = m.iGetNumCols(); c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] -= m(r, c);
		}
	}
}

/* setta una matrice di tipo Mat3xN in una data posizione */
void
FullMatrixHandler::Put(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - 2);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - m.iGetNumCols() + 1);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = 3; r > 0; r--) {
		for (integer c = m.iGetNumCols(); c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] = m(r, c);
		}
	}
}

/* somma una matrice di tipo Mat3xN in una data posizione, trasposta */
void
FullMatrixHandler::AddT(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumCols() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumCols(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] += m(c, r);
		}
	}
}

/* sottrae una matrice di tipo Mat3xN in una data posizione, trasposta */
void
FullMatrixHandler::SubT(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumCols() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumCols(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] -= m(c, r);
		}
	}
}

/* setta una matrice di tipo Mat3xN in una data posizione, trasposta */
void
FullMatrixHandler::PutT(integer iRow, integer iCol, const Mat3xN& m)
{
	/* iRow e iCol sono gli indici effettivi di riga e colonna
	 * es. per il primo coefficiente:
	 *     iRow = 1, iCol = 1 */

#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumCols() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumCols(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] = m(c, r);
		}
	}
}

/* somma una matrice di tipo MatNx3 in una data posizione */
void
FullMatrixHandler::Add(integer iRow, integer iCol, const MatNx3& m)
{
#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumRows() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumRows(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] += m(r, c);
		}
	}
}

/* sottrae una matrice di tipo MatNx3 in una data posizione */
void
FullMatrixHandler::Sub(integer iRow, integer iCol, const MatNx3& m)
{
#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumRows() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumRows(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] -= m(r, c);
		}
	}
}

/* setta una matrice di tipo MatNx3 in una data posizione */
void
FullMatrixHandler::Put(integer iRow, integer iCol, const MatNx3& m)
{
#ifdef DEBUG
	IsValid();

	ASSERT(iRow > 0);
	ASSERT(iRow <= iNumRows - m.iGetNumRows() + 1);
	ASSERT(iCol > 0);
	ASSERT(iCol <= iNumCols - 2);
#endif /* DEBUG */

	--iRow;
	--iCol;
	for (int r = m.iGetNumRows(); r > 0; r--) {
		for (integer c = 3; c > 0; c--) {
			ppdColsm1[iCol + c][iRow + r] = m(r, c);
		}
	}
}

void
FullMatrixHandler::CopyMatrixRow(integer dest_row,
	const FullMatrixHandler & source, integer source_row)
{
	integer nc = iGetNumCols();

	ASSERT(dest_row >= 1);
	ASSERT(dest_row <= iGetNumRows());
	ASSERT(source_row >= 1);
	ASSERT(source_row <= source.iGetNumRows());
	ASSERT(nc == source.iGetNumCols());

	for (integer ic = 1; ic <= nc; ic++) {
		// dest(dest_row, ic) = source(source_row, ic);
		ppdColsm1[ic][dest_row] = source.ppdColsm1[ic][source_row];
	}
}

void
FullMatrixHandler::CopyMatrixBlock(integer dest_row, integer dest_col,
	const FullMatrixHandler & source, 
	integer source_start_row, integer source_end_row,
	integer source_start_col, integer source_end_col)
{
	ASSERT(source_start_row >= 1);
	ASSERT(source_end_row <= source.iGetNumRows());
	ASSERT(source_start_row <= source_end_row);
	ASSERT(dest_row >= 1);
	ASSERT(dest_row + (source_end_row - source_start_row) <= iGetNumRows());

	ASSERT(source_start_col >= 1);
	ASSERT(source_end_col <= source.iGetNumCols());
	ASSERT(source_start_col <= source_end_col);
	ASSERT(dest_col >= 1);
	ASSERT(dest_col + (source_end_col - source_start_col) <= iGetNumCols());

	for (integer ir = source_start_row; ir <= source_end_row; ir++) {
		integer row = dest_row + (ir - source_start_row);
		for (integer ic = source_start_col; ic <= source_end_col; ic++) {
			integer col = dest_col + (ic - source_start_col);
			// dest(row, col) = source(ir, ic);
			ppdColsm1[col][row] = source.ppdColsm1[ic][ir];
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
FullMatrixHandler::const_iterator::operator -> (void) const
{
	return &elem;
}

const SparseMatrixHandler::SparseMatrixElement&
FullMatrixHandler::const_iterator::operator * (void) const
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

