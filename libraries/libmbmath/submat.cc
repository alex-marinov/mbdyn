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

/* sottomatrici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <string.h>	/* for memset() */
#include <iomanip>

#include <submat.h>

/* SubMatrixHandler - begin */

SubMatrixHandler::~SubMatrixHandler(void)
{
        NO_OP;
}

SubMatrixHandler* SubMatrixHandler::Copy() const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SubMatrixHandler::EnumerateNz(const std::function<EnumerateNzCallback>& func) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

/* SubMatrixHandler - end */


/* FullSubMatrixHandler - begin */

FullSubMatrixHandler::FullSubMatrixHandler(integer iIntSize,
                integer* piTmpVec,
                integer iDoubleSize,
                doublereal* pdTmpMat,
                integer iMaxCols,
                doublereal **ppdCols)
: FullMatrixHandler(pdTmpMat, ppdCols, iDoubleSize, 1, 1, iMaxCols),
iVecSize(iIntSize), piRow(0), piRowm1(0), piColm1(0)
{
#ifdef DEBUG
        IsValid();
#endif

        /*
         * Non posso fare controlli piu' precisi sulla consistenza
         * delle dimensioni delle aree di lavoro perche' e' legale
         * che la submatrix sia in seguito ridimensionata nei modi
         * piu' vari. */
        ASSERT(piTmpVec);
        piRow = piTmpVec;
        piRowm1 = piTmpVec - 1;
        piColm1 = piRowm1 + iNumRows;
}

FullSubMatrixHandler::FullSubMatrixHandler(integer iNR, integer iNC)
: FullMatrixHandler(iNR, iNC), iVecSize(iNR + iNC), piRow(0), piRowm1(0), piColm1(0)
{
        ASSERT(bOwnsMemory);
        SAFENEWARR(piRow, integer, iVecSize);

        piRowm1 = piRow - 1;
        piColm1 = piRowm1 + iNR;
}


FullSubMatrixHandler::~FullSubMatrixHandler(void)
{
        if (bOwnsMemory) {
                SAFEDELETEARR(piRow);
        }
}

#ifdef DEBUG
void
FullSubMatrixHandler::IsValid(void) const
{
        FullMatrixHandler::IsValid();

        ASSERT(iVecSize > 0);
        ASSERT(iNumRows + iNumCols <= iVecSize);
        ASSERT(piRowm1 != NULL);
        ASSERT(piColm1 != NULL);

#ifdef DEBUG_MEMMANAGER
        ASSERT(defaultMemoryManager.fIsValid(piRow,
                                iVecSize*sizeof(integer)));
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */


void
FullSubMatrixHandler::Reset(void)
{
#ifdef DEBUG
        IsValid();
#endif /* DEBUG */

        FullMatrixHandler::Reset();

#if 0
        /*
         * this is not strictly required, because all the indices should
         * be explicitly set before the matrix is used
         */
        for (integer i = iNumRows + iNumCols; i > 0; i--) {
                piRowm1[i] = 0;
        }
#endif
}

/*
 * Modifica le dimensioni correnti
 */
void
FullSubMatrixHandler::Resize(integer iNewRow, integer iNewCol) {
#ifdef DEBUG
        IsValid();
#endif /* DEBUG */

        ASSERT(iNewRow > 0);
        ASSERT(iNewCol > 0);
        ASSERT(iNewRow + iNewCol <= iVecSize);
        ASSERT(iNewRow*iNewCol <= iRawSize);

        if (iNewRow <= 0
                        || iNewCol <= 0
                        || iNewRow + iNewCol > iVecSize
                        || iNewRow*iNewCol > iRawSize
                        || iNewCol > iMaxCols) {
                silent_cerr("FullSubMatrixHandler::Resize() - error"
                                << std::endl);

                throw SubMatrixHandler::ErrResize(MBDYN_EXCEPT_ARGS);
        }

        iNumRows = iNewRow;
        iNumCols = iNewCol;

        /*
         * FIXME: don't call FullMatrixHandler::Resize, because
         * we know there's room enough and we don't want to
         * waste time regenerating support stuff; we simply
         * use a portion of the matrix as is.
         */
#if 0
        FullMatrixHandler::Resize(iNewRow, iNewCol);

        ASSERT(piRowm1 != NULL);
        piColm1 = piRowm1 + iNewRow;
#endif

#ifdef DEBUG
        IsValid();
#endif /* DEBUG */
}

/* Ridimensiona ed inizializza. */
void
FullSubMatrixHandler::ResizeReset(integer ir, integer ic)
{
        FullSubMatrixHandler::Resize(ir, ic);
        FullSubMatrixHandler::Reset();
}

/*
 * Collega la matrice Full alla memoria che gli viene passata
 * in ingresso
 */
void
FullSubMatrixHandler::Attach(int iRows, int iCols, integer* piTmpIndx)
{
        iNumRows = iRows;
        iNumCols = iCols;

        if (bOwnsMemory) {
                if (piRow != 0) {
                        SAFEDELETEARR(piRow);
                }
                piRow = piTmpIndx;
                piRowm1 = piTmpIndx - 1;
                piColm1 = piRowm1 + iNumRows;
        }
}

void
FullSubMatrixHandler::Add(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::Add(iRow, iCol, v);
}

void
FullSubMatrixHandler::Sub(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::Sub(iRow, iCol, v);
}

void
FullSubMatrixHandler::Put(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::Put(iRow, iCol, v);
}

void
FullSubMatrixHandler::AddT(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::AddT(iRow, iCol, v);
}

void
FullSubMatrixHandler::SubT(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::SubT(iRow, iCol, v);
}

void
FullSubMatrixHandler::PutT(integer iRow, integer iCol, const Vec3& v)
{
        FullMatrixHandler::PutT(iRow, iCol, v);
}

#if 0 /* FIXME: replace original? */
/* somma un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullSubMatrixHandler::AddDiag(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol][iRow] += pdFrom[V1];

        ppdColsm1[iCol + 1][iRow + 1] += pdFrom[V2];

        ppdColsm1[iCol + 2][iRow + 2] += pdFrom[V3];
}

/* sottrae un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullSubMatrixHandler::SubDiag(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol][iRow] -= pdFrom[V1];

        ppdColsm1[iCol + 1][iRow + 1] -= pdFrom[V2];

        ppdColsm1[iCol + 2][iRow + 2] -= pdFrom[V3];
}

/* scrive un vettore di tipo Vec3 in una data posizione in diagonale */
void
FullSubMatrixHandler::PutDiag(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol][iRow] = pdFrom[V1];

        ppdColsm1[iCol + 1][iRow + 1] = pdFrom[V2];

        ppdColsm1[iCol + 2][iRow + 2] = pdFrom[V3];
}

/* somma un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullSubMatrixHandler::AddCross(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol + 1][iRow] -= pdFrom[V3];
        ppdColsm1[iCol + 2][iRow] += pdFrom[V2];

        ppdColsm1[iCol][iRow + 1] += pdFrom[V3];
        ppdColsm1[iCol + 2][iRow + 1] -= pdFrom[V1];

        ppdColsm1[iCol][iRow + 2] -= pdFrom[V2];
        ppdColsm1[iCol + 1][iRow + 2] += pdFrom[V1];
}

/* sottrae un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullSubMatrixHandler::SubCross(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol + 1][iRow] += pdFrom[V3];
        ppdColsm1[iCol + 2][iRow] -= pdFrom[V2];

        ppdColsm1[iCol][iRow + 1] -= pdFrom[V3];
        ppdColsm1[iCol + 2][iRow + 1] += pdFrom[V1];

        ppdColsm1[iCol][iRow + 2] += pdFrom[V2];
        ppdColsm1[iCol + 1][iRow + 2] -= pdFrom[V1];
}

/* scrive un vettore di tipo Vec3 in una data posizione [ v x ] */
void
FullSubMatrixHandler::PutCross(integer iRow, integer iCol, const Vec3& v)
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
        doublereal* pdFrom = v.pGetVec();

        ppdColsm1[iCol + 1][iRow] = -pdFrom[V3];
        ppdColsm1[iCol + 2][iRow] = pdFrom[V2];

        ppdColsm1[iCol][iRow + 1] = pdFrom[V3];
        ppdColsm1[iCol + 2][iRow + 1] = -pdFrom[V1];

        ppdColsm1[iCol][iRow + 2] = -pdFrom[V2];
        ppdColsm1[iCol + 1][iRow + 2] = pdFrom[V1];
}
#endif

void
FullSubMatrixHandler::Add(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::Add(iRow, iCol, m);
}

void
FullSubMatrixHandler::AddT(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::AddT(iRow, iCol, m);
}

void
FullSubMatrixHandler::Sub(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::Sub(iRow, iCol, m);
}

void
FullSubMatrixHandler::SubT(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::SubT(iRow, iCol, m);
}

void
FullSubMatrixHandler::Put(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::Put(iRow, iCol, m);
}

void
FullSubMatrixHandler::PutT(integer iRow, integer iCol, const Mat3x3& m)
{
        FullMatrixHandler::PutT(iRow, iCol, m);
}

void
FullSubMatrixHandler::Add(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::Add(iRow, iCol, m);
}

void
FullSubMatrixHandler::Sub(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::Sub(iRow, iCol, m);
}

void
FullSubMatrixHandler::Put(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::Put(iRow, iCol, m);
}

void
FullSubMatrixHandler::AddT(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::AddT(iRow, iCol, m);
}

void
FullSubMatrixHandler::SubT(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::SubT(iRow, iCol, m);
}

void
FullSubMatrixHandler::PutT(integer iRow, integer iCol, const Mat3xN& m)
{
        FullMatrixHandler::PutT(iRow, iCol, m);
}

void
FullSubMatrixHandler::Add(integer iRow, integer iCol, const MatNx3& m)
{
        FullMatrixHandler::Add(iRow, iCol, m);
}

/* sottrae una matrice di tipo MatNx3 in una data posizione */
void
FullSubMatrixHandler::Sub(integer iRow, integer iCol, const MatNx3& m)
{
        FullMatrixHandler::Sub(iRow, iCol, m);
}

/* setta una matrice di tipo MatNx3 in una data posizione */
void
FullSubMatrixHandler::Put(integer iRow, integer iCol, const MatNx3& m)
{
        FullMatrixHandler::Put(iRow, iCol, m);
}

/* setta una matrice di tipo Mat3xN in una data posizione */
void
FullSubMatrixHandler::PutDiag(integer iFirstRow, integer iFirstCol,
                const Vec3& v)
{
        /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
         * es. per il primo coefficiente:
         *     iRow = 0, iCol = 0 */

#ifdef DEBUG
        IsValid();

        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstRow <= iNumRows - 3);
        ASSERT(iFirstCol >= 0);
        ASSERT(iFirstCol <= iNumCols - 3);
#endif /* DEBUG */

        const doublereal *pdv = v.pGetVec();

        ppdColsm1[iFirstCol + 1][iFirstRow + 1] = pdv[V1];
        ppdColsm1[iFirstCol + 2][iFirstRow + 2] = pdv[V2];
        ppdColsm1[iFirstCol + 3][iFirstRow + 3] = pdv[V3];
}


/* setta una matrice di tipo Mat3xN in una data posizione */
void
FullSubMatrixHandler::PutDiag(integer iFirstRow, integer iFirstCol,
                const doublereal& d)
{
        /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
         * es. per il primo coefficiente:
         *     iRow = 0, iCol = 0 */

#ifdef DEBUG
        IsValid();

        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstRow <= iNumRows - 3);
        ASSERT(iFirstCol >= 0);
        ASSERT(iFirstCol <= iNumCols - 3);
#endif /* DEBUG */

        ppdColsm1[iFirstCol + 1][iFirstRow + 1] = d;
        ppdColsm1[iFirstCol + 2][iFirstRow + 2] = d;
        ppdColsm1[iFirstCol + 3][iFirstRow + 3] = d;
}


/* setta una matrice di tipo Mat3xN in una data posizione */
void
FullSubMatrixHandler::PutCross(integer iFirstRow, integer iFirstCol,
                const Vec3& v)
{
        /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
         * es. per il primo coefficiente:
         *     iRow = 0, iCol = 0 */

#ifdef DEBUG
        IsValid();

        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstRow <= iNumRows - 3);
        ASSERT(iFirstCol >= 0);
        ASSERT(iFirstCol <= iNumCols - 3);
#endif /* DEBUG */

        const doublereal *pdv = v.pGetVec();

        ppdColsm1[iFirstCol + 1][ iFirstRow + 2] = pdv[V3];
        ppdColsm1[iFirstCol + 1][ iFirstRow + 3] = -pdv[V2];

        ppdColsm1[iFirstCol + 2][ iFirstRow + 1] = -pdv[V3];
        ppdColsm1[iFirstCol + 2][ iFirstRow + 3] = pdv[V1];

        ppdColsm1[iFirstCol + 3][ iFirstRow + 1] = pdv[V2];
        ppdColsm1[iFirstCol + 3][ iFirstRow + 2] = -pdv[V1];
}


/* somma la matrice ad un matrix handler usando i metodi generici */
MatrixHandler&
FullSubMatrixHandler::AddTo(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumRows);
        ASSERT(MH.iGetNumCols() >= iNumCols);

#ifdef USE_SPARSE_AUTODIFF
        sp_grad::SpGradient oRow, oItem;

#ifdef USE_MULTITHREAD
        std::vector<bool> rgStatus(iNumRows, false);
        bool bRepeatLoop;

        do {
             bRepeatLoop = false;
#endif
             for (integer r = 1; r <= iNumRows; ++r) {
#ifdef USE_MULTITHREAD
                  if (rgStatus[r - 1]) {
                       continue;
                  }
#endif
                  oRow.ResizeReset(0., iNumCols);

                  for (integer c = 1; c <= iNumCols; ++c) {
                       oItem.Reset(0., piColm1[c], ppdColsm1[c][r]);
                       oRow += oItem;
                  }

#ifdef USE_MULTITHREAD
                  bool bStatus =
#endif
                       MH.AddItem(piRowm1[r], oRow);

#ifdef USE_MULTITHREAD
                  if (bStatus) {
                       rgStatus[r - 1] = bStatus;
                  } else {
                       bRepeatLoop = true;
                  }
#endif
             }
#ifdef USE_MULTITHREAD
        } while (bRepeatLoop);
#endif
#else
        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumCols());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumRows());

                        MH(piRowm1[r], piColm1[c]) += ppdColsm1[c][r];
                }
        }
#endif
        return MH;
}


/* somma la matrice, trasposta, ad un matrix handler usando i metodi generici */
MatrixHandler&
FullSubMatrixHandler::AddToT(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumCols);
        ASSERT(MH.iGetNumCols() >= iNumRows);

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumRows());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumCols());

                        MH(piColm1[c], piRowm1[r]) += ppdColsm1[c][r];
                }
        }

        return MH;
}


/* somma la matrice ad un FullMatrixHandler */
MatrixHandler&
FullSubMatrixHandler::AddTo(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumRows);
        ASSERT(MH.iGetNumCols() >= iNumCols);

        doublereal **ppd = MH.ppdColsm1;

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumCols());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumRows());

                        ppd[piColm1[c]][piRowm1[r]] += ppdColsm1[c][r];
                }
        }

        return MH;
}


/* somma la matrice, trasposta, ad un FullMatrixHandler */
MatrixHandler&
FullSubMatrixHandler::AddToT(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumCols);
        ASSERT(MH.iGetNumCols() >= iNumRows);

        doublereal **ppd = MH.ppdColsm1;

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumRows());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumCols());

                        ppd[piRowm1[r]][piColm1[c]] += ppdColsm1[c][r];
                }
        }

        return MH;
}

#ifdef USE_SPARSE_AUTODIFF
VectorHandler& FullSubMatrixHandler::AddTo(VectorHandler& A, const VectorHandler& Y) const
{
#ifdef DEBUG
        IsValid();
        A.IsValid();
        Y.IsValid();
#endif
        ASSERT(A.iGetSize() == Y.iGetSize());
        
        for (integer c = 1; c <= iNumCols; ++c) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= Y.iGetSize());

                for (integer r = 1; r <= iNumRows; ++r) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= A.iGetSize());

                        A(piRowm1[r]) += ppdColsm1[c][r] * Y(piColm1[c]);
                }
        }
        
        return A;
}
#endif

/* sottrae la matrice da un matrix handler usando i metodi generici */
MatrixHandler&
FullSubMatrixHandler::SubFrom(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumRows);
        ASSERT(MH.iGetNumCols() >= iNumCols);

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumCols());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumRows());

                        MH(piRowm1[r], piColm1[c]) -= ppdColsm1[c][r];
                }
        }

        return MH;
}


/* sottrae la matrice, trasposta, da un matrix handler usando i metodi generici */
MatrixHandler&
FullSubMatrixHandler::SubFromT(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumCols);
        ASSERT(MH.iGetNumCols() >= iNumRows);

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumRows());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumCols());

                        MH(piColm1[c], piRowm1[r]) -= ppdColsm1[c][r];
                }
        }

        return MH;
}


/* sottrae la matrice da un FullMatrixHandler */
MatrixHandler&
FullSubMatrixHandler::SubFrom(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumRows);
        ASSERT(MH.iGetNumCols() >= iNumCols);

        doublereal **ppd = MH.ppdColsm1;

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumCols());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumRows());

                        ppd[piColm1[c]][piRowm1[r]] -= ppdColsm1[c][r];
                }
        }

        return MH;
};


/* sottrae la matrice, trasposta, da un FullMatrixHandler */
MatrixHandler&
FullSubMatrixHandler::SubFromT(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        ASSERT(MH.iGetNumRows() >= iNumCols);
        ASSERT(MH.iGetNumCols() >= iNumRows);

        doublereal **ppd = MH.ppdColsm1;

        for (integer c = iNumCols; c > 0; c--) {
                ASSERT(piColm1[c] > 0);
                ASSERT(piColm1[c] <= MH.iGetNumRows());

                for (integer r = iNumRows; r > 0; r--) {
                        ASSERT(piRowm1[r] > 0);
                        ASSERT(piRowm1[r] <= MH.iGetNumCols());

                        ppd[piRowm1[r]][piColm1[c]] -= ppdColsm1[c][r];
                }
        }

        return MH;
};


/* output, usato principalmente per debug */
std::ostream&
operator << (std::ostream& out, const FullSubMatrixHandler& m)
{
#ifdef DEBUG
        m.IsValid();
#endif /* DEBUG */

        ASSERT(m.iNumRows > 0);
        ASSERT(m.iNumCols > 0);
        ASSERT(m.piRowm1 != NULL);
        ASSERT(m.piColm1 != NULL);
        ASSERT(m.ppdColsm1 != NULL);

        out << std::setw(12) << "";
        for (integer c = 1; c <= m.iNumCols; c++) {
                out << std::setw(12) << m.piColm1[c];
        }
        out << std::endl << std::endl;

        for (integer r = 1; r <= m.iNumRows; r++) {
                out << std::setw(12) << m.piRowm1[r];
                for (integer c = 1; c <= m.iNumCols; c++) {
                        out << std::setw(12) << m.ppdColsm1[c][r];
                }
                out << std::endl;
        }

        return out << std::endl;
}

void
FullSubMatrixHandler::Add(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::Add(iRow, iCol, source);
}

void
FullSubMatrixHandler::Sub(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::Sub(iRow, iCol, source);
}

void
FullSubMatrixHandler::Put(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::Put(iRow, iCol, source);
}

void
FullSubMatrixHandler::Add(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::Add(iRow, iCol, source, dCoef);
}

void
FullSubMatrixHandler::Sub(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::Sub(iRow, iCol, source, dCoef);
}

void
FullSubMatrixHandler::Put(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::Put(iRow, iCol, source, dCoef);
}

void
FullSubMatrixHandler::AddT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::AddT(iRow, iCol, source);
}

void
FullSubMatrixHandler::SubT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::SubT(iRow, iCol, source);
}

void
FullSubMatrixHandler::PutT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source)
{
        FullMatrixHandler::PutT(iRow, iCol, source);
}

void
FullSubMatrixHandler::AddT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::AddT(iRow, iCol, source, dCoef);
}

void
FullSubMatrixHandler::SubT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::SubT(iRow, iCol, source, dCoef);
}

void
FullSubMatrixHandler::PutT(integer iRow,
        integer iCol,
        const FullMatrixHandler & source,
        const doublereal dCoef)
{
        FullMatrixHandler::PutT(iRow, iCol, source, dCoef);
}

/* FullSubMatrixHandler - end */


/* SparseSubMatrixHandler - begin */

SparseSubMatrixHandler::SparseSubMatrixHandler(integer iTmpInt,
                integer* piTmpIndex, integer iTmpDouble, doublereal* pdTmpMat)
: bOwnsMemory(false),
iIntSize(iTmpInt), iDoubleSize(iTmpDouble),
iNumItems(iTmpInt/2),
piRow(0), piRowm1(0), piColm1(0),
pdMat(0), pdMatm1(0)
{
        ASSERT(piTmpIndex);
        ASSERT(pdTmpMat);

        piRow = piTmpIndex;
        piRowm1 = piTmpIndex - 1;
        piColm1 = piRowm1 + iNumItems;

        pdMat = pdTmpMat;
        pdMatm1 = pdTmpMat - 1;

#ifdef DEBUG
        IsValid();
#endif /* DEBUG */
}

SparseSubMatrixHandler::SparseSubMatrixHandler(integer iTmpInt)
: iIntSize(2*iTmpInt), iDoubleSize(iTmpInt),
iNumItems(iTmpInt), piRow(0), piRowm1(0), piColm1(0), pdMat(0), pdMatm1(0)
{
        SAFENEWARR(pdMat, doublereal, iNumItems);
        pdMatm1 = pdMat - 1;

        SAFENEWARR(piRow, integer, iIntSize);
        piRowm1 = piRow - 1;
        piColm1 = piRowm1 + iNumItems;

        bOwnsMemory = true;
}

/* Distruttore banale.
 * Nota: dato che la classe non possiede la memoria,
 * non ne deve deallocare
 */
SparseSubMatrixHandler::~SparseSubMatrixHandler(void)
{
        if (bOwnsMemory) {
                SAFEDELETEARR(pdMat);

                SAFEDELETEARR(piRow);
        }
}

#ifdef DEBUG
void
SparseSubMatrixHandler::IsValid(void) const
{
        ASSERT(iIntSize > 0);
        ASSERT(iIntSize%2 == 0);
        ASSERT(iDoubleSize > 0);
        ASSERT(iIntSize >= 2*iDoubleSize);
        ASSERT(iNumItems >= 0);
        ASSERT(piRowm1 != NULL);
        ASSERT(piColm1 != NULL);
        ASSERT(pdMatm1 != NULL);

#ifdef DEBUG_MEMMANAGER
        ASSERT(defaultMemoryManager.fIsValid(piRow,
                                iIntSize*sizeof(integer)));
        ASSERT(defaultMemoryManager.fIsValid(pdMat,
                                iDoubleSize*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}
#endif /* DEBUG */


/*
 * Ridimensiona la matrice.
 * Nota: solo il primo argomento viene considerato,
 * e rappresenta il numero totale di entries.
 * Questo metodo deve essere chiamato prima di qualsiasi
 * operazione sulla matrice.
 */
void
SparseSubMatrixHandler::Resize(integer iNewRow, integer iNewCol)
{
#ifdef DEBUG
        IsValid();
#endif /* DEBUG */

        ASSERT(iNewRow > 0);
        ASSERT(2*iNewRow <= iIntSize);
        ASSERT(iNewRow <= iDoubleSize);
        ASSERT(piRowm1 != NULL);

        if (iNewRow <= 0
                        || 2*iNewRow > iIntSize
                        || iNewRow > iDoubleSize) {
                silent_cerr("SparseSubMatrixHandler::Resize() - error"
                                << std::endl);
                throw SparseSubMatrixHandler::ErrResize(MBDYN_EXCEPT_ARGS);
        }

        iNumItems = iNewRow;

#ifdef DEBUG
        IsValid();
#endif /* DEBUG */
}

/*
 * Ridimensiona ed inizializza.
 * Unione dei due metodi precedenti
 */
void
SparseSubMatrixHandler::ResizeReset(integer iNewRow, integer iNewCol)
{
        Resize(iNewRow, iNewCol);
        Reset();
}

/*
 * Collega la matrice sparsa alla memoria che gli viene passata
 * in ingresso
 */
void
SparseSubMatrixHandler::Attach(int iNumEntr, doublereal* pdTmpMat,
                integer* piTmpIndx)
{
        if (bOwnsMemory) {
                SAFEDELETEARR(piRow);

                SAFEDELETEARR(pdMat);

                bOwnsMemory = false;
        }

        ASSERT(iNumEntr > 0);
        ASSERT(piTmpIndx);
        ASSERT(pdTmpMat);

        iIntSize = iNumEntr*2;
        iDoubleSize = iNumEntr;
        iNumItems = iNumEntr;
        piRow = piTmpIndx;
        piRowm1 = piTmpIndx - 1;
        piColm1 = piRowm1 + iNumEntr;
        pdMat = pdTmpMat;
        pdMatm1 = pdTmpMat - 1;

#ifdef DEBUG
        IsValid();
#endif /* DEBUG */
}

void
SparseSubMatrixHandler::PutDiag(integer iSubIt, integer iFirstRow,
                integer iFirstCol, const Vec3& v)
{
#ifdef DEBUG
        IsValid();

        ASSERT(iNumItems >= 3);
        ASSERT(iSubIt > 0);
        ASSERT(iSubIt <= iNumItems - 2);
        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstCol >= 0);
#endif /* DEBUG */

        /*
         * Attenzione agli argomenti:
         * iSubIt e' il primo indice della matrice da utilizzare,
         * con 1 <= iSubit <= iCurSize;
         * iFirstRow e' il primo indice di riga -1, ovvero il primo indice
         * di riga della sottomatrice diag(v) piena e' iFirstRow + 1
         * iFirstCol e' il primo indice di colonna -1, ovvero il primo indice
         * di colonna della sottomatrice diag(v) piena e' iFirstCol + 1
         * v e' il vettore che genera diag(v)
         */

        /* Matrice diag(v) :
         *
         *         1   2   3
         *
         * 1    |  v1  0   0  |
         * 2    |  0   v2  0  |
         * 3    |  0   0   v3 |
         */

        /* assume che il Vec3 sia un'array di 3 reali */
        const doublereal* pdFrom = v.pGetVec();

        doublereal* pdm = pdMatm1 + iSubIt;
        integer* pir = piRowm1 + iSubIt;
        integer* pic = piColm1 + iSubIt;

        /* Coefficiente 1,1 */
        pdm[0] = pdFrom[V1];
        pir[0] = iFirstRow + 1;
        pic[0] = iFirstCol + 1;

        /* Coefficiente 2,2 */
        pdm[1] = pdFrom[V2];
        pir[1] = iFirstRow + 2;
        pic[1] = iFirstCol + 2;

        /* Coefficiente 3,3 */
        pdm[2] = pdFrom[V3];
        pir[2] = iFirstRow + 3;
        pic[2] = iFirstCol + 3;
}


void
SparseSubMatrixHandler::PutDiag(integer iSubIt, integer iFirstRow,
               integer iFirstCol, const doublereal& d)
{
#ifdef DEBUG
        IsValid();

        ASSERT(iNumItems >= 3);
        ASSERT(iSubIt > 0);
        ASSERT(iSubIt <= iNumItems - 2);
        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstCol >= 0);
#endif /* DEBUG */

        /* Attenzione agli argomenti:
         * iSubIt e' il primo indice della matrice da utilizzare,
         * con 1 <= iSubit <= iCurSize;
         * iFirstRow e' il primo indice di riga -1, ovvero il
         * primo indice di riga della sottomatrice I*d piena e' iFirstRow+1
         * iFirstCol e' il primo indice di colonna -1, ovvero il
         * primo indice di colonna della sottomatrice I*d piena e' iFirstCol+1
         * v e' il vettore che genera I*d */

        /* Matrice I*d :
         *
         *         1   2   3
         *
         * 1    |  d   0   0 |
         * 2    |  0   d   0 |
         * 3    |  0   0   d |
         */

        doublereal* pdm = pdMatm1 + iSubIt;
        integer* pir = piRowm1 + iSubIt;
        integer* pic = piColm1 + iSubIt;

        /* Coefficiente 1,1 */
        pdm[0] = d;
        pir[0] = iFirstRow + 1;
        pic[0] = iFirstCol + 1;

        /* Coefficiente 2,2 */
        pdm[1] = d;
        pir[1] = iFirstRow + 2;
        pic[1] = iFirstCol + 2;

        /* Coefficiente 3,3 */
        pdm[2] = d;
        pir[2] = iFirstRow + 3;
        pic[2] = iFirstCol + 3;
}


void
SparseSubMatrixHandler::PutCross(integer iSubIt, integer iFirstRow,
               integer iFirstCol, const Vec3& v)
{
#ifdef DEBUG
        IsValid();

        ASSERT(iNumItems >= 6);
        ASSERT(iSubIt > 0);
        ASSERT(iSubIt <= iNumItems - 5);
        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstCol >= 0);
#endif /* DEBUG */

        /* Attenzione agli argomenti:
         * iSubIt e' il primo indice della matrice da utilizzare,
         * con 1 <= iSubit <= iCurSize;
         * iFirstRow e' il primo indice di riga -1, ovvero il
         * primo indice di riga della sottomatrice v/\ piena e' iFirstRow+1
         * iFirstCol e' il primo indice di colonna -1, ovvero il
         * primo indice di colonna della sottomatrice v/\ piena e' iFirstCol+1
         * v e' il vettore che genera v/\ */

        /* Matrice v/\ :
         *
         *         1   2   3
         *
         * 1    |  0  -v3  v2 |
         * 2    |  v3  0  -v1 |
         * 3    | -v2  v1  0  |
         */

        /* assume che il Vec3 sia un'array di 3 reali */
        const doublereal* pdFrom = v.pGetVec();

        /* Coefficiente 1,2 */
        doublereal* pdm = pdMatm1 + iSubIt;
        integer* pir = piRowm1 + iSubIt;
        integer* pic = piColm1 + iSubIt;

        pdm[0] = -pdFrom[V3];               // -v.dGet(3);
        pir[0] = iFirstRow + 1;
        pic[0] = iFirstCol + 2;

        /* Coefficiente 1,3 */
        pdm[1] = pdFrom[V2];                // v.dGet(2);
        pir[1] = iFirstRow + 1;
        pic[1] = iFirstCol + 3;

        /* Coefficiente 2,1 */
        pdm[2] = pdFrom[V3];                // v.dGet(3);
        pir[2] = iFirstRow + 2;
        pic[2] = iFirstCol + 1;

        /* Coefficiente 2,3 */
        pdm[3] = -pdFrom[V1];               // -v.dGet(1);
        pir[3] = iFirstRow + 2;
        pic[3] = iFirstCol + 3;

        /* Coefficiente 3,1 */
        pdm[4] = -pdFrom[V2];                // -v.dGet(2);
        pir[4] = iFirstRow + 3;
        pic[4] = iFirstCol + 1;

        /* Coefficiente 3,2 */
        pdm[5] = pdFrom[V1];                 // v.dGet(1);
        pir[5] = iFirstRow + 3;
        pic[5] = iFirstCol + 2;
}


void
SparseSubMatrixHandler::Reset(void)
{
#ifdef DEBUG
        IsValid();
#endif /* DEBUG */

        ASSERT(iNumItems > 0);

#if defined HAVE_MEMSET
        memset(pdMatm1 + 1, 0, iNumItems*sizeof(doublereal));
#else /* !HAVE_MEMSET */
        for (integer i = iNumItems; i > 0; i--) {
                pdMatm1[i] = 0.;
        }
#endif /* HAVE_MEMSET */
}


/* Inserisce una matrice 3x3;
 * si noti che non ci sono Add, Sub, ecc. perche' la filosofia
 * della matrice sparsa prevede che ad ogni item (riga, colonna, valore)
 * corrisponda un termine che poi verra' sommato ad una matrice vera,
 * senza controlli su eventuali duplicazioni. */
void
SparseSubMatrixHandler::PutMat3x3(integer iSubIt, integer iFirstRow,
                integer iFirstCol, const Mat3x3& m)
{
#ifdef DEBUG
        IsValid();

        ASSERT(iNumItems >= 9);
        ASSERT(iSubIt > 0);
        ASSERT(iSubIt <= iNumItems - 8);
        ASSERT(iFirstRow >= 0);
        ASSERT(iFirstCol >= 0);
#endif /* DEBUG */

        /* Attenzione agli argomenti:
         * iSubIt e' il primo indice della matrice da utilizzare,
         * con 1 <= iSubit <= iCurSize;
         * iFirstRow e' il primo indice di riga -1, ovvero il
         * primo indice di riga della sottomatrice v/\ piena e' iFirstRow+1
         * iFirstCol e' il primo indice di colonna -1, ovvero il
         * primo indice di colonna della sottomatrice m piena e' iFirstCol+1
         */

        /* Per efficienza, vengono scritte esplicitamente tutte
         * le assegnazioni;
         * la funzione quindi non e' messa in linea intenzionalmente */

        /* Coefficienti 1,1-3,1 */
        doublereal* pdFrom = (doublereal*)m.pGetMat();
        doublereal* pdTmpMat = pdMatm1 + iSubIt;
        integer* piTmpRow = piRowm1 + iSubIt;
        integer* piTmpCol = piColm1 + iSubIt;

        iFirstRow++;
        iFirstCol++;

        /* Prima riga */
        pdTmpMat[0] = pdFrom[M11];
        piTmpRow[0] = iFirstRow++;
        piTmpCol[0] = iFirstCol;

        pdTmpMat[1] = pdFrom[M21];
        piTmpRow[1] = iFirstRow++;
        piTmpCol[1] = iFirstCol;

        pdTmpMat[2] = pdFrom[M31];
        piTmpRow[2] = iFirstRow;
        piTmpCol[2] = iFirstCol;

        /* Seconda riga */
        iFirstRow -= 2;
        iFirstCol++;

        pdTmpMat[3] = pdFrom[M12];
        piTmpRow[3] = iFirstRow++;
        piTmpCol[3] = iFirstCol;

        pdTmpMat[4] = pdFrom[M22];
        piTmpRow[4] = iFirstRow++;
        piTmpCol[4] = iFirstCol;

        pdTmpMat[5] = pdFrom[M32];
        piTmpRow[5] = iFirstRow;
        piTmpCol[5] = iFirstCol;

        /* Terza riga */
        iFirstRow -= 2;
        iFirstCol++;

        pdTmpMat[6] = pdFrom[M13];
        piTmpRow[6] = iFirstRow++;
        piTmpCol[6] = iFirstCol;

        pdTmpMat[7] = pdFrom[M23];
        piTmpRow[7] = iFirstRow++;
        piTmpCol[7] = iFirstCol;

        pdTmpMat[8] = pdFrom[M33];
        piTmpRow[8] = iFirstRow;
        piTmpCol[8] = iFirstCol;
}


/* somma la matrice ad un matrix handler usando i metodi generici */
MatrixHandler&
SparseSubMatrixHandler::AddTo(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

#ifdef USE_SPARSE_AUTODIFF
        sp_grad::SpGradient oItem;

#ifdef USE_MULTITHREAD
        std::vector<bool> rgStatus(iNumItems, false);
        bool bRepeatLoop;

        do {
             bRepeatLoop = false;
#endif
             for (integer i = 1; i <= iNumItems; ++i) {
#ifdef USE_MULTITHREAD
                  if (rgStatus[i - 1]) {
                       continue;
                  }
#endif
                  oItem.Reset(0., piColm1[i], pdMatm1[i]);

#ifdef USE_MULTITHREAD
                  bool bStatus =
#endif
                       MH.AddItem(piRowm1[i], oItem);

#ifdef USE_MULTITHREAD
                  if (bStatus) {
                       rgStatus[i - 1] = bStatus;
                  } else {
                       bRepeatLoop = true;
                  }
#endif
             }
#ifdef USE_MULTITHREAD
        } while (bRepeatLoop);
#endif
#else
        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumRows());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumCols());

                MH(piRowm1[i], piColm1[i]) += pdMatm1[i];
        }
#endif
        return MH;
}


/* somma la matrice, trasposta, ad un matrix handler usando i metodi generici */
MatrixHandler&
SparseSubMatrixHandler::AddToT(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumCols());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumRows());

                MH(piColm1[i], piRowm1[i]) += pdMatm1[i];
        }

        return MH;
}

/* somma la matrice ad un FullMatrixHandler */
MatrixHandler&
SparseSubMatrixHandler::AddTo(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        doublereal **ppd = MH.ppdColsm1;

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumRows());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumCols());

                ppd[piColm1[i]][piRowm1[i]] += pdMatm1[i];
        }

        return MH;
}


/* somma la matrice, trasposta, ad un FullMatrixHandler */
MatrixHandler&
SparseSubMatrixHandler::AddToT(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        doublereal **ppd = MH.ppdColsm1;

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumCols());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumRows());

                ppd[piRowm1[i]][piColm1[i]] += pdMatm1[i];
        }

        return MH;
}

#ifdef USE_SPARSE_AUTODIFF
VectorHandler& SparseSubMatrixHandler::AddTo(VectorHandler& A, const VectorHandler& Y) const
{
#ifdef DEBUG
        IsValid();
        A.IsValid();
        Y.IsValid();
#endif
        
        for (integer i = 1; i <= iNumItems; ++i) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= A.iGetSize());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= Y.iGetSize());

                A(piRowm1[i]) += pdMatm1[i] * Y(piColm1[i]);
        }
        
        return A;
}
#endif

/* sottrae la matrice da un matrix handler usando i metodi generici */
MatrixHandler&
SparseSubMatrixHandler::SubFrom(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumRows());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumCols());

                MH(piRowm1[i], piColm1[i]) -= pdMatm1[i];
        }

        return MH;
}


/* sottrae la matrice, trasposta, da un matrix handler usando i metodi generici */
MatrixHandler&
SparseSubMatrixHandler::SubFromT(MatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumCols());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumRows());

                MH(piColm1[i], piRowm1[i]) -= pdMatm1[i];
        }

        return MH;
}


/* sottrae la matrice da un FullMatrixHandler */
MatrixHandler&
SparseSubMatrixHandler::SubFrom(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        doublereal **ppd = MH.ppdColsm1;

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumRows());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumCols());

                ppd[piColm1[i]][piRowm1[i]] -= pdMatm1[i];
        }

        return MH;
}

/* sottrae la matrice, trasposta, da un FullMatrixHandler */
MatrixHandler&
SparseSubMatrixHandler::SubFromT(FullMatrixHandler& MH) const
{
#ifdef DEBUG
        IsValid();
        MH.IsValid();
#endif /* DEBUG */

        doublereal **ppd = MH.ppdColsm1;

        for (integer i = iNumItems; i > 0; i--) {
                ASSERT(piRowm1[i] > 0);
                ASSERT(piRowm1[i] <= MH.iGetNumCols());
                ASSERT(piColm1[i] > 0);
                ASSERT(piColm1[i] <= MH.iGetNumRows());

                ppd[piRowm1[i]][piColm1[i]] -= pdMatm1[i];
        }

        return MH;
}

/* SparseSubMatrixHandler - end */

#ifdef USE_SPARSE_AUTODIFF
SpGradientSubMatrixHandler::SpGradientSubMatrixHandler(integer iNumItemsMax) {
     oVec.reserve(iNumItemsMax);
}

SpGradientSubMatrixHandler::~SpGradientSubMatrixHandler()
{
}

void SpGradientSubMatrixHandler::Resize(integer iNumRows, integer)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSubMatrixHandler::Reset()
{
     oVec.clear();
}

void SpGradientSubMatrixHandler::ResizeReset(integer iNumRows, integer iNumCols)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

doublereal& SpGradientSubMatrixHandler::operator()(integer iRow, integer iCol)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

const doublereal& SpGradientSubMatrixHandler::operator()(integer iRow, integer iCol) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

integer SpGradientSubMatrixHandler::iGetNumRows() const
{
     return oVec.size();
}

integer SpGradientSubMatrixHandler::iGetNumCols() const
{
     return 0;
}

void SpGradientSubMatrixHandler::PutRowIndex(integer iSubIt, integer iRow)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void SpGradientSubMatrixHandler::PutColIndex(integer iSubIt, integer iRow)
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

integer SpGradientSubMatrixHandler::iGetRowIndex(integer) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

integer SpGradientSubMatrixHandler::iGetColIndex(integer) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

MatrixHandler& SpGradientSubMatrixHandler::AddTo(MatrixHandler& MH) const
{
#ifdef USE_MULTITHREAD
     for (const auto& oItem: oVec) {
          oItem.bInserted = false;
     }

     bool bRepeatLoop;

     do {
          bRepeatLoop = false;
#endif
          for (const auto& oItem: oVec) {
#ifdef USE_MULTITHREAD
               if (oItem.bInserted) {
                    continue;
               }
#endif

#ifdef USE_MULTITHREAD
               bool bStatus =
#endif
                    MH.AddItem(oItem.iEquationIdx, oItem.oResidual);

#ifdef USE_MULTITHREAD
               if (bStatus) {
                    oItem.bInserted = bStatus;
               } else {
                    bRepeatLoop = true;
               }
#endif
          }

#ifdef USE_MULTITHREAD
     } while (bRepeatLoop);
#endif
     return MH;
}

VectorHandler& SpGradientSubMatrixHandler::AddTo(VectorHandler& A, const VectorHandler& Y) const
{
#ifdef DEBUG
        IsValid();
        A.IsValid();
        Y.IsValid();
#endif

        ASSERT(A.iGetSize() == Y.iGetSize());
        
        for (const auto& oItem: oVec) {
                ASSERT(oItem.iEquationIdx >= 1);
                ASSERT(oItem.iEquationIdx <= A.iGetSize());
                
                for (const auto& oRec: oItem.oResidual) {
                        ASSERT(oRec.iDof >= 1);
                        ASSERT(oRec.iDof <= Y.iGetSize());
                        
                        A(oItem.iEquationIdx) += oRec.dDer * Y(oRec.iDof);
                }
        }

        return A;
}

MatrixHandler& SpGradientSubMatrixHandler::SubFrom(MatrixHandler& MH) const
{
     for (const auto& oItem: oVec) {
          for (const auto& oRec: oItem.oResidual) {
               if (oRec.dDer) {
                    MH(oItem.iEquationIdx, oRec.iDof) -= oRec.dDer;
               }
          }
     }

     return MH;
}

MatrixHandler& SpGradientSubMatrixHandler::AddToT(MatrixHandler& MH) const
{
     for (const auto& oItem: oVec) {
          for (const auto& oRec: oItem.oResidual) {
               if (oRec.dDer) {
                    MH(oRec.iDof, oItem.iEquationIdx) += oRec.dDer;
               }
          }
     }

     return MH;
}

MatrixHandler& SpGradientSubMatrixHandler::SubFromT(MatrixHandler& MH) const
{
     for (const auto& oItem: oVec) {
          for (const auto& oRec: oItem.oResidual) {
               if (oRec.dDer) {
                    MH(oRec.iDof, oItem.iEquationIdx) -= oRec.dDer;
               }
          }
     }

     return MH;
}

bool SpGradientSubMatrixHandler::AddItem(integer iEquationIdx, const sp_grad::SpGradient& oResidual)
{
     SP_GRAD_ASSERT(oVec.size() < oVec.capacity());

     oVec.emplace_back(iEquationIdx, oResidual);

     return true;
}

#ifdef DEBUG
void SpGradientSubMatrixHandler::IsValid(void) const
{
     for (const auto& oItem: oVec) {
          SP_GRAD_ASSERT(oItem.oResidual.bValid());
     }
}
#endif
#endif

/* MySubVectorHandler - begin */

MySubVectorHandler::MySubVectorHandler(integer iSize)
: MyVectorHandler(), piRow(NULL), piRowm1(NULL) {
        Resize(iSize);

#ifdef DEBUG
        for (integer i = 1; i <= iMaxSize; ++i) {
             piRowm1[i] = std::numeric_limits<integer>::max();
        }

        IsValid();
#endif
}

MySubVectorHandler::MySubVectorHandler(integer iSize, integer* piTmpRow,
               doublereal* pdTmpVec)
: MyVectorHandler(iSize, pdTmpVec), piRow(0), piRowm1(0)
{
        if (piTmpRow == 0) {
                /* ... because set by MyVectorHandler() */
                ASSERT(bOwnsMemory == true);

                if (pdTmpVec == 0) {
                        silent_cerr("MySubVectorHandler(): illegal args" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                SAFENEWARR(piRow, integer, iSize);
                piRowm1 = piRow;
        }

        piRowm1--;

#ifdef DEBUG
        for (integer i = 1; i <= iMaxSize; ++i) {
             piRowm1[i] = std::numeric_limits<integer>::max();
        }

        IsValid();
#endif /* DEBUG */
}

void
MySubVectorHandler::Resize(integer iSize)
{
#ifdef DEBUG
     IsValid();
#endif
        if (iSize < 0) {
                silent_cerr("Negative size!" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        ASSERT((piRowm1 == NULL && pdVecm1 == NULL)
                        || (piRowm1 != NULL && pdVecm1 != NULL));

        if (!bOwnsMemory && piRowm1 != NULL) {
                if (iSize > iMaxSize) {
                        silent_cerr("Can't resize to " << iSize
                                << ": larger than max size " << iMaxSize
                                << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
                iCurSize = iSize;

        } else {
                if (piRowm1 != NULL) {
                        if (iSize <= iMaxSize) {
                                iCurSize = iSize;

                        } else {
                                doublereal* pd = NULL;
                                SAFENEWARR(pd, doublereal, iSize);
                                //pd--;
                                for (integer i = iCurSize-1; i >= 0; i--) {
                                        pd[i] = pdVec[i];
                                }

                                integer* pi = NULL;
                                SAFENEWARR(pi, integer, iSize);
                                //pi--;
                                for (integer i = iCurSize-1; i >= 0; i--) {
                                        pi[i] = piRow[i];
                                }

#ifdef DEBUG
                                for (integer i = iCurSize; i < iSize; ++i) {
                                     pi[i] = std::numeric_limits<integer>::max();
                                }
#endif
                                //pdVecm1++;
                                SAFEDELETEARR(pdVec);

                                SAFEDELETEARR(piRow);

                                pdVec = pd;
                                pdVecm1 = pd - 1;
                                piRow = pi;
                                piRowm1 = pi - 1;
                                iMaxSize = iCurSize = iSize;
                        }

                } else {
                        if (iSize > 0) {
                                bOwnsMemory = true;

                                SAFENEWARR(pdVec, doublereal, iSize);
                                pdVecm1 = pdVec - 1;

                                SAFENEWARR(piRow, integer, iSize);
                                piRowm1 = piRow - 1;

                                iMaxSize = iCurSize = iSize;

#ifdef DEBUG
                                for (integer i = 1; i <= iMaxSize; ++i) {
                                     piRowm1[i] = std::numeric_limits<integer>::max();
                                }
#endif
                        }
                }
        }

#ifdef DEBUG
        IsValid();
#endif
}

void
MySubVectorHandler::Detach(void)
{
#ifdef DEBUG
     IsValid();
#endif
        if (bOwnsMemory) {
                if (pdVecm1 != NULL) {
                        //pdVecm1++;
                        SAFEDELETEARR(pdVec);

                        SAFEDELETEARR(piRow);
                }

                bOwnsMemory = false;
        }

        iMaxSize = iCurSize = 0;
        pdVec = NULL;
        pdVecm1 = NULL;
        piRow = NULL;
        piRowm1 = NULL;
}

void
MySubVectorHandler::Attach(integer iSize, doublereal* pd,
                integer* pi, integer iMSize)
{
#ifdef DEBUG
     IsValid();
#endif

        if (bOwnsMemory && pdVecm1 != NULL) {
                Detach();
                bOwnsMemory = false;
        }

        iMaxSize = iCurSize = iSize;
        if (iMSize >= iSize) {
                iMaxSize = iMSize;
        }

        pdVec = pd;
        pdVecm1 = pd - 1;
        piRow = pi;
        piRowm1 = pi - 1;

#ifdef DEBUG
        IsValid();
#endif
}

#ifdef DEBUG
void
MySubVectorHandler::IsValid(void) const
{
        MyVectorHandler::IsValid();

        ASSERT((piRowm1 != NULL && iMaxSize > 0) || (piRowm1 == NULL && iMaxSize == 0));

#ifdef DEBUG_MEMMANAGER
        ASSERT(defaultMemoryManager.fIsValid(piRow,
                                MyVectorHandler::iMaxSize*sizeof(integer)));
#endif /* DEBUG_MEMMANAGER */

        for (integer i = 1; i <= iMaxSize; ++i) {
             ASSERT(piRowm1[i] > 0);
        }
}
#endif /* DEBUG */

VectorHandler&
MySubVectorHandler::AddTo(VectorHandler& VH) const
{
#ifdef DEBUG
        IsValid();
        VH.IsValid();
#endif /* DEBUG */

        for (integer i = iGetSize(); i > 0; i--) {
#if 0
                /* FIXME: workaround for SchurVectorHandler... */
                VH(piRowm1[i]) += pdVecm1[i];
#endif
                VH.IncCoef(piRowm1[i], pdVecm1[i]);
        }

        return VH;
}

VectorHandler&
MySubVectorHandler::AddTo(MyVectorHandler& VH) const
{
#ifdef DEBUG
        IsValid();
        VH.IsValid();
#endif /* DEBUG */

        doublereal* pdm1 = VH.pdGetVec() - 1;
        for (integer i = iGetSize(); i > 0; i--) {
                pdm1[piRowm1[i]] += pdVecm1[i];
        }
        return VH;
}

VectorHandler&
MySubVectorHandler::AddAbsValuesTo(VectorHandler& VH) const
{
#ifdef DEBUG
        IsValid();
        VH.IsValid();
#endif /* DEBUG */

        for (integer i = iGetSize(); i > 0; i--) {
#if 0
                /* FIXME: workaround for SchurVectorHandler... */
                VH(piRowm1[i]) += pdVecm1[i];
#endif
                VH.IncCoef(piRowm1[i], std::abs(pdVecm1[i]));
        }

        return VH;
}

VectorHandler&
MySubVectorHandler::AddAbsValuesTo(MyVectorHandler& VH) const
{
#ifdef DEBUG
        IsValid();
        VH.IsValid();
#endif /* DEBUG */

        doublereal* pdm1 = VH.pdGetVec() - 1;
        for (integer i = iGetSize(); i > 0; i--) {
                pdm1[piRowm1[i]] += std::abs(pdVecm1[i]);
        }
        return VH;
}

std::ostream&
operator << (std::ostream& out, const SubVectorHandler& v)
{
#ifdef DEBUG
        v.IsValid();
#endif /* DEBUG */

        integer iRow = v.iGetSize();

        ASSERT(iRow > 0);

        for (integer i = 1; i <= iRow; i++) {
                out << std::setw(12) << v.iGetRowIndex(i)
                        << " " << std::setw(12) << v(i) << std::endl;
        }

        return out << std::endl;
}

/* MySubVectorHandler - end */
