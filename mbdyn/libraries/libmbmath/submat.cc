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

/* sottomatrici */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <string.h>	/* for memset() */
#include <ac/iomanip>

#include <submat.h>

/* SubMatrixHandler - begin */

SubMatrixHandler::~SubMatrixHandler(void)
{
   NO_OP;
}

/* SubMatrixHandler - end */


/* FullSubMatrixHandler - begin */

FullSubMatrixHandler::FullSubMatrixHandler(integer iIntSize,
					   integer iDoubleSize,
					   integer* piTmpVec, 
					   doublereal* pdTmpMat)
: iVecSize(iIntSize), iMatSize(iDoubleSize), 
iNumRows(1), iNumCols(1),
piRow(piTmpVec), piCol(piRow+iNumRows), pdMat(pdTmpMat) 
{
#ifdef DEBUG	
   IsValid();
#endif	
   
   /* Non posso fare controlli piu' precisi sulla consistenza 
    * delle dimensioni delle aree di lavoro perche' e' legale
    * che la submatrix sia in seguito ridimensionata nei modi
    * piu' vari. */
}


FullSubMatrixHandler::~FullSubMatrixHandler(void)
{
   NO_OP;
}
					   
					   
void FullSubMatrixHandler::IsValid(void) const
{
   ASSERT(iVecSize > 0);
   ASSERT(iMatSize > 0);
   ASSERT(iNumRows >= 0);
   ASSERT(iNumCols >= 0);
   ASSERT(iNumRows+iNumCols <= iVecSize);
   ASSERT(iNumRows*iNumCols <= iMatSize);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
   ASSERT(pdMat != NULL);
   
#ifdef DEBUG_MEMMANAGER       
   ASSERT(defaultMemoryManager.fIsValid(piRow, iVecSize*sizeof(integer)));
   ASSERT(defaultMemoryManager.fIsValid(pdMat, iMatSize*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}


void 
FullSubMatrixHandler::Init(const doublereal& dResetVal)
{
#ifdef DEBUG	
      IsValid();
#endif
      
#ifdef HAVE_MEMSET
      if (dResetVal == 0.) {
	      memset(pdMat, 0, iNumRows*iNumCols*sizeof(doublereal));
      } else {
#endif /* HAVE_MEMSET */
      for (integer i = iNumRows*iNumCols; i-- > 0; ) {
	 pdMat[i] = dResetVal;
      }
#ifdef HAVE_MEMSET
      }
#endif /* HAVE_MEMSET */
     
#if 0 
      /* 
       * this is not strictly required, because all the indices should
       * be explicitly set before the matrix is used
       */
      for (integer i = iNumRows+iNumCols; i-- > 0; ) {
	 piRow[i] = 0;
      }	
#endif
}   
   

/* somma una matrice di tipo Mat3x3 in una data posizione */
void FullSubMatrixHandler::Add(integer iRow, integer iCol, const Mat3x3& m)
{
   /* iRow e iCol sono gli indici effettivi di riga e colonna
    * es. per il primo coefficiente:
    *     iRow = 1, iCol = 1 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-2));
#endif	

   /* Nota: assume che la matrice sia organizzata "per colonne" (stile FORTRAN)
    * e assume che la matrice Mat3x3 si aorganizzata anch'essa per colonne */
   doublereal* pdTo = pdMat+(--iRow)+(--iCol)*iNumRows;
   doublereal* pdFrom = m.pGetMat();

   /* Prima colonna */
   pdTo[0] += pdFrom[M11];
   pdTo[1] += pdFrom[M21];
   pdTo[2] += pdFrom[M31];
   
   /* Seconda colonna */
   pdTo += iNumRows;
   pdTo[0] += pdFrom[M12];
   pdTo[1] += pdFrom[M22];
   pdTo[2] += pdFrom[M32];
   
   /* Terza colonna */
   pdTo += iNumRows;
   pdTo[0] += pdFrom[M13];
   pdTo[1] += pdFrom[M23];
   pdTo[2] += pdFrom[M33];
}


/* sottrae una matrice di tipo Mat3x3 in una data posizione
 * analoga ala precedente, con il meno (per evitare temporanei ecc) */
void FullSubMatrixHandler::Sub(integer iRow, integer iCol, const Mat3x3& m)
{
   /* iRow e iCol sono gli indici effettivi di riga e colonna
    * es. per il primo coefficiente:
    *     iRow = 1, iCol = 1
    */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-2));
#endif	
      
   doublereal* pdTo = pdMat+(--iRow)+(--iCol)*iNumRows;	
   doublereal* pdFrom = m.pGetMat();
   
   /* Prima colonna */
   pdTo[0] -= pdFrom[M11];
   pdTo[1] -= pdFrom[M21];
   pdTo[2] -= pdFrom[M31];
   
   /* Seconda colonna */
   pdTo += iNumRows;
   pdTo[0] -= pdFrom[M12];
   pdTo[1] -= pdFrom[M22];
   pdTo[2] -= pdFrom[M32];
   
   /* Terza colonna */
   pdTo += iNumRows;
   pdTo[0] -= pdFrom[M13];
   pdTo[1] -= pdFrom[M23];
   pdTo[2] -= pdFrom[M33];
}


/* mette una matrice di tipo Mat3x3 in una data posizione;
 * analoga alle precedenti */
void FullSubMatrixHandler::Put(integer iRow, integer iCol, const Mat3x3& m)
{
#ifdef DEBUG
   IsValid();

   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-2));
#endif	
   
   doublereal* pdTo = pdMat+(--iRow)+(--iCol)*iNumRows;	
   doublereal* pdFrom = m.pGetMat();
   
   /* Prima colonna */
   pdTo[0] = pdFrom[M11];
   pdTo[1] = pdFrom[M21];
   pdTo[2] = pdFrom[M31];
   
   /* Seconda colonna */
   pdTo += iNumRows;
   pdTo[0] = pdFrom[M12];
   pdTo[1] = pdFrom[M22];
   pdTo[2] = pdFrom[M32];
   
   /* Terza colonna */
   pdTo += iNumRows;
   pdTo[0] = pdFrom[M13];
   pdTo[1] = pdFrom[M23];
   pdTo[2] = pdFrom[M33];
}


/* somma una matrice di tipo Mat3xN in una data posizione */
void FullSubMatrixHandler::Add(integer iRow, integer iCol, const Mat3xN& m)
{
   /* iRow e iCol sono gli indici effettivi di riga e colonna
    * es. per il primo coefficiente:
    *     iRow = 1, iCol = 1 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-m.iGetNumCols()+1));
#endif	

   --iRow;
   --iCol;
   for (int i = 3; i > 0; i--) {
      for (integer j = m.iGetNumCols(); j > 0; j--) {
	 fIncCoef(i+iRow, j+iCol, m.dGet(i, j));
      }
   }
}

/* sottrae una matrice di tipo Mat3xN in una data posizione */
void FullSubMatrixHandler::Sub(integer iRow, integer iCol, const Mat3xN& m)
{
   /* iRow e iCol sono gli indici effettivi di riga e colonna
    * es. per il primo coefficiente:
    *     iRow = 1, iCol = 1 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-m.iGetNumCols()+1));
#endif	

   --iRow;
   --iCol;
   for (int i = 3; i > 0; i--) {
      for (integer j = m.iGetNumCols(); j > 0; j--) {
	 fDecCoef(i+iRow, j+iCol, m.dGet(i, j));
      }
   }
}

/* somma una matrice di tipo MatNx3 in una data posizione */
void FullSubMatrixHandler::Add(integer iRow, integer iCol, const MatNx3& m)
{
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-m.iGetNumRows()+1));
   ASSERT((iCol > 0) && (iCol <= iNumCols-2));
#endif	 

   --iRow;
   --iCol;
   for (int i = m.iGetNumRows(); i > 0; i--) {
      for (integer j = 3; j > 0; j--) {
	 fIncCoef(i+iRow, j+iCol, m.dGet(i, j));
      }
   }
}

/* sottrae una matrice di tipo MatNx3 in una data posizione */
void FullSubMatrixHandler::Sub(integer iRow, integer iCol, const MatNx3& m)
{
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-m.iGetNumRows()+1));
   ASSERT((iCol > 0) && (iCol <= iNumCols-2));
#endif	 

   --iRow;
   --iCol;
   for (int i = m.iGetNumRows(); i > 0; i--) {
      for (integer j = 3; j > 0; j--) {
	 fDecCoef(i+iRow, j+iCol, m.dGet(i, j));
      }
   }
}


/* setta una matrice di tipo Mat3xN in una data posizione */
void FullSubMatrixHandler::Put(integer iRow, integer iCol, const Mat3xN& m)
{
   /* iRow e iCol sono gli indici effettivi di riga e colonna
    * es. per il primo coefficiente:
    *     iRow = 1, iCol = 1 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iRow > 0) && (iRow <= iNumRows-2));
   ASSERT((iCol > 0) && (iCol <= iNumCols-m.iGetNumCols()+1));
#endif	

   --iRow;
   --iCol;
   for (int i = 3; i > 0; i--) {
      for (integer j = m.iGetNumCols(); j > 0; j--) {
	 fPutCoef(i+iRow, j+iCol, m.dGet(i, j));
      }
   }
}


/* setta una matrice di tipo Mat3xN in una data posizione */
flag
FullSubMatrixHandler::fPutDiag(integer iFirstRow, integer iFirstCol, 
		const Vec3& v)
{
   /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
    * es. per il primo coefficiente:
    *     iRow = 0, iCol = 0 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iFirstRow >= 0) && (iFirstRow <= iNumRows-3));
   ASSERT((iFirstCol >= 0) && (iFirstCol <= iNumCols-3));
#endif	

   const doublereal *pdv = v.pGetVec();

   fPutCoef(iFirstRow+1, iFirstCol+1, pdv[V1]);
   fPutCoef(iFirstRow+2, iFirstCol+2, pdv[V2]);
   fPutCoef(iFirstRow+3, iFirstCol+3, pdv[V3]);

   return flag(0);
}


/* setta una matrice di tipo Mat3xN in una data posizione */
flag
FullSubMatrixHandler::fPutDiag(integer iFirstRow, integer iFirstCol, 
		const doublereal& d)
{
   /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
    * es. per il primo coefficiente:
    *     iRow = 0, iCol = 0 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iFirstRow >= 0) && (iFirstRow <= iNumRows-3));
   ASSERT((iFirstCol >= 0) && (iFirstCol <= iNumCols-3));
#endif	

   fPutCoef(iFirstRow+1, iFirstCol+1, d);
   fPutCoef(iFirstRow+2, iFirstCol+2, d);
   fPutCoef(iFirstRow+3, iFirstCol+3, d);

   return flag(0);
}


/* setta una matrice di tipo Mat3xN in una data posizione */
flag
FullSubMatrixHandler::fPutCross(integer iFirstRow, integer iFirstCol, 
		const Vec3& v)
{
   /* iFirstRow e iFirstCol sono gli indici effettivi di riga e colonna -1
    * es. per il primo coefficiente:
    *     iRow = 0, iCol = 0 */
   
#ifdef DEBUG
   IsValid();
   
   ASSERT((iFirstRow >= 0) && (iFirstRow <= iNumRows-3));
   ASSERT((iFirstCol >= 0) && (iFirstCol <= iNumCols-3));
#endif	

   const doublereal *pdv = v.pGetVec();

   fPutCoef(iFirstRow+1, iFirstCol+2, -pdv[V3]);
   fPutCoef(iFirstRow+1, iFirstCol+3, pdv[V2]);

   fPutCoef(iFirstRow+2, iFirstCol+1, pdv[V3]);
   fPutCoef(iFirstRow+2, iFirstCol+3, -pdv[V1]);

   fPutCoef(iFirstRow+3, iFirstCol+1, -pdv[V2]);
   fPutCoef(iFirstRow+3, iFirstCol+2, pdv[V1]);

   return flag(0);
}


/* somma la matrice ad un matrix handler usando i metodi generici */
MatrixHandler& FullSubMatrixHandler::AddTo(MatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   ASSERT(MH.iGetNumRows() >= iNumRows);
   ASSERT(MH.iGetNumCols() >= iNumCols);
   
   doublereal* pd = pdMat+iNumRows*iNumCols;
   
   for (integer j = iNumCols; j-- > 0; ) {
      ASSERT(piCol[j] > 0 && piCol[j] <= MH.iGetNumCols());
      pd -= iNumRows;
      for (integer i = iNumRows; i-- > 0; ) {
	 ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
	 MH.fIncCoef(piRow[i], piCol[j], pd[i]);
      }
   }
   
   return MH;
}


/* somma la matrice ad un FullMatrixHandler */
MatrixHandler& FullSubMatrixHandler::AddTo(FullMatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   ASSERT(MH.iGetNumRows() >= iNumRows);
   ASSERT(MH.iGetNumCols() >= iNumCols);
   
   doublereal* pd = pdMat+iNumRows*iNumCols;
   
   for (integer j = iNumCols; j-- > 0; ) {
      ASSERT(piCol[j] > 0 && piCol[j] <= MH.iGetNumCols());
      
      pd -= iNumRows;
      doublereal* pdTo = MH.pdGetMat()+MH.iGetNumRows()*(piCol[j]-1)-1;
      
      for (integer i = iNumRows; i-- > 0; ) {
	 ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
	 pdTo[piRow[i]] += pd[i];
      }
   }
   
   return MH;
}


/* sottrae la matrice da un matrix handler usando i metodi generici */
MatrixHandler& FullSubMatrixHandler::SubFrom(MatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   ASSERT(MH.iGetNumRows() >= iNumRows);
   ASSERT(MH.iGetNumCols() >= iNumCols);
   
   doublereal* pd = pdMat+iNumRows*iNumCols;
   
   for (integer j = iNumCols; j-- > 0; ) {
      ASSERT(piCol[j] > 0 && piCol[j] <= MH.iGetNumCols());
      pd -= iNumRows;
      for (integer i = iNumRows; i-- > 0; ) {
	 ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
	 MH.fDecCoef(piRow[i], piCol[j], pd[i]);
      }
   }
   
   return MH;
}


/* sottrae la matrice da un FullMatrixHandler */
MatrixHandler& FullSubMatrixHandler::SubFrom(FullMatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   ASSERT(MH.iGetNumRows() >= iNumRows);
   ASSERT(MH.iGetNumCols() >= iNumCols);
   
   doublereal* pd = pdMat+iNumRows*iNumCols;
   
   for (integer j = iNumCols; j-- > 0; ) {
      ASSERT(piCol[j] > 0 && piCol[j] <= MH.iGetNumCols());
      
      pd -= iNumRows;
      doublereal* pdTo = MH.pdGetMat()+(MH.iGetNumRows()*piCol[j]-1)-1;
      
      for (integer i = iNumRows; i-- > 0; ) {
	 ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
	 pdTo[piRow[i]] -= pd[i];
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
#endif	
   
   integer iRow = (integer)m.iNumRows;
   integer iCol = (integer)m.iNumCols;
   integer* piRow = (integer*)m.piRow;
   integer* piCol = (integer*)m.piCol;
   doublereal* pd = (doublereal*)m.pdMat;
   
   ASSERT(iRow > 0);
   ASSERT(iCol > 0);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
   ASSERT(pd != NULL);
   
   integer* piCnt = piCol;
   out << std::setw(12) << "";
   while (piCnt < piCol+iCol) {      
      out << std::setw(12) << *piCnt++;
   }
   out << std::endl << std::endl;
     
   piCnt = piRow;
   while (piCnt < piRow+iRow) {	
      out << std::setw(12) << *piCnt++;
      for (integer iCnt = 0; iCnt < iRow; iCnt++) {	 
	 out << std::setw(12) << *(pd+iCnt*iRow);
      }
      out << std::endl;
      pd++;
   }
   
   return out << std::endl;   
}

/* FullSubMatrixHandler - end */


/* SparseSubMatrixHandler - begin */

void SparseSubMatrixHandler::IsValid(void) const
{
   ASSERT(iIntSize > 0);
   ASSERT(iIntSize%2 == 0);
   ASSERT(iDoubleSize > 0);
   ASSERT(iIntSize >= 2*iDoubleSize);
   ASSERT(iNumItems >= 0);
   ASSERT(piRow != NULL);
   ASSERT(piCol != NULL);
   ASSERT(pdMat != NULL);
   
#ifdef DEBUG_MEMMANAGER       
   ASSERT(defaultMemoryManager.fIsValid(piRow, iIntSize*sizeof(integer)));
   ASSERT(defaultMemoryManager.fIsValid(pdMat, iDoubleSize*sizeof(doublereal)));
#endif /* DEBUG_MEMMANAGER */
}


flag SparseSubMatrixHandler::fPutDiag(integer iSubIt, integer iFirstRow, 
				       integer iFirstCol, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   
   ASSERT(iNumItems >= 3);
   ASSERT(iSubIt > 0 && iSubIt <= iNumItems-2);
   ASSERT(iFirstRow >= 0);
   ASSERT(iFirstCol >= 0);
#endif	
   
   /* Attenzione agli argomenti:
    * iSubIt e' il primo indice della matrice da utilizzare,
    * con 1 <= iSubit <= iCurSize;
    * iFirstRow e' il primo indice di riga -1, ovvero il
    * primo indice di riga della sottomatrice diag(v) piena e' iFirstRow+1
    * iFirstCol e' il primo indice di colonna -1, ovvero il
    * primo indice di colonna della sottomatrice diag(v) piena e' iFirstCol+1
    * v e' il vettore che genera diag(v) */
   
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
   
   doublereal* pdm = pdMat+(--iSubIt);
   integer* pir = piRow+iSubIt;
   integer* pic = piCol+iSubIt;
   
   /* Coefficiente 1,1 */
   pdm[0] = pdFrom[V1];
   pir[0] = iFirstRow+1;
   pic[0] = iFirstCol+1;
   
   /* Coefficiente 2,2 */
   pdm[1] = pdFrom[V2];
   pir[1] = iFirstRow+2;
   pic[1] = iFirstCol+2;
   
   /* Coefficiente 3,3 */
   pdm[2] = pdFrom[V3];
   pir[2] = iFirstRow+3;
   pic[2] = iFirstCol+3;
   
   return flag(0);
}


flag SparseSubMatrixHandler::fPutDiag(integer iSubIt, integer iFirstRow, 
				       integer iFirstCol, const doublereal& d)
{
#ifdef DEBUG
   IsValid();
   
   ASSERT(iNumItems >= 3);
   ASSERT(iSubIt > 0 && iSubIt <= iNumItems-2);
   ASSERT(iFirstRow >= 0);
   ASSERT(iFirstCol >= 0);
#endif	
   
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

   doublereal* pdm = pdMat+(--iSubIt);
   integer* pir = piRow+iSubIt;
   integer* pic = piCol+iSubIt;
   
   /* Coefficiente 1,1 */
   pdm[0] = d;
   pir[0] = iFirstRow+1;
   pic[0] = iFirstCol+1;
   
   /* Coefficiente 2,2 */
   pdm[1] = d;
   pir[1] = iFirstRow+2;
   pic[1] = iFirstCol+2;
   
   /* Coefficiente 3,3 */
   pdm[2] = d;
   pir[2] = iFirstRow+3;
   pic[2] = iFirstCol+3;
   
   return flag(0);
}


flag SparseSubMatrixHandler::fPutCross(integer iSubIt, integer iFirstRow, 
				       integer iFirstCol, const Vec3& v)
{
#ifdef DEBUG
   IsValid();
   
   ASSERT(iNumItems >= 6);
   ASSERT(iSubIt > 0 && iSubIt <= iNumItems-5);
   ASSERT(iFirstRow >= 0);
   ASSERT(iFirstCol >= 0);
#endif	
   
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
   doublereal* pdm = pdMat+(--iSubIt);
   integer* pir = piRow+iSubIt;
   integer* pic = piCol+iSubIt;
   
   pdm[0] = -pdFrom[V3];               // -v.dGet(3);
   pir[0] = iFirstRow+1;
   pic[0] = iFirstCol+2;
   
   /* Coefficiente 1,3 */
   pdm[1] = pdFrom[V2];                // v.dGet(2);
   pir[1] = iFirstRow+1;
   pic[1] = iFirstCol+3;
   
   /* Coefficiente 2,1 */
   pdm[2] = pdFrom[V3];                // v.dGet(3);
   pir[2] = iFirstRow+2;
   pic[2] = iFirstCol+1;
   
   /* Coefficiente 2,3 */
   pdm[3] = -pdFrom[V1];               // -v.dGet(1);
   pir[3] = iFirstRow+2;
   pic[3] = iFirstCol+3;
   
   /* Coefficiente 3,1 */
   pdm[4] = -pdFrom[V2];                // -v.dGet(2);
   pir[4] = iFirstRow+3;
   pic[4] = iFirstCol+1;
   
   /* Coefficiente 3,2 */
   pdm[5] = pdFrom[V1];                 // v.dGet(1);
   pir[5] = iFirstRow+3;
   pic[5] = iFirstCol+2;
   
   return flag(0);
}


void 
SparseSubMatrixHandler::Init(const doublereal& dCoef)
{
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT(iNumItems > 0);	
      
#ifdef HAVE_MEMSET
      if (dCoef == 0.) {
	      memset(pdMat, 0, iNumItems*sizeof(doublereal));
      } else {
#endif /* HAVE_MEMSET */
      for (integer i = 0; i < iNumItems; i++) {
	      pdMat[i] = dCoef;
      }
#ifdef HAVE_MEMSET
      }
#endif /* HAVE_MEMSET */
}
   

/* Inserisce una matrice 3x3; 
 * si noti che non ci sono Add, Sub, ecc. perche' la filosofia 
 * della matrice sparsa prevede che ad ogni item (riga, colonna, valore)
 * corrisponda un termine che poi verra' sommato ad una matrice vera,
 * senza controlli su eventuali duplicazioni. */
flag SparseSubMatrixHandler::fPutMat3x3(integer iSubIt, integer iFirstRow, 
					integer iFirstCol, const Mat3x3& m)
{
#ifdef DEBUG
   IsValid();

   ASSERT(iNumItems >= 9);
   ASSERT(iSubIt > 0 && iSubIt <= iNumItems-8);
   ASSERT(iFirstRow >= 0);
   ASSERT(iFirstCol >= 0);
#endif	
   
   /* Attenzione agli argomenti:
    * iSubIt e' il primo indice della matrice da utilizzare,
    * con 1 <= iSubit <= iCurSize;
    * iFirstRow e' il primo indice di riga -1, ovvero il
    * primo indice di riga della sottomatrice v/\ piena e' iFirstRow+1
    * iFirstCol e' il primo indice di colonna -1, ovvero il
    * primo indice di colonna della sottomatrice m piena e' iFirstCol+1 
    */

/* Per efficienza, vengono scritte esplicitamente tutte le assegnazioni;
 * la funzione quindi non e' messa in linea intenzionalmente */
   
   /* Coefficienti 1,1-3,1 */
   doublereal* pdFrom = (doublereal*)m.pGetMat();
   doublereal* pdTmpMat = pdMat+(--iSubIt);
   integer* piTmpRow = piRow+iSubIt;
   integer* piTmpCol = piCol+iSubIt;
   
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

   return flag(0);
}


/* somma la matrice ad un matrix handler usando i metodi generici */
MatrixHandler& SparseSubMatrixHandler::AddTo(MatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   for (integer i = iNumItems; i-- > 0; ) {
      ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
      ASSERT(piCol[i] > 0 && piCol[i] <= MH.iGetNumCols());
      MH.fIncCoef(piRow[i], piCol[i], pdMat[i]);
   }
   
   return MH;
}


/* somma la matrice ad un FullMatrixHandler */
MatrixHandler& SparseSubMatrixHandler::AddTo(FullMatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   doublereal* pdTo = MH.pdGetMat();
   integer iNumR = MH.iGetNumRows();

   for (integer i = iNumItems; i-- > 0; ) {
      ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
      ASSERT(piCol[i] > 0 && piCol[i] <= MH.iGetNumCols());
      pdTo[(piRow[i]-1)+iNumR*(piCol[i]-1)] += pdMat[i];
   }
   
   return MH;
}


/* sottrae la matrice da un matrix handler usando i metodi generici */
MatrixHandler& SparseSubMatrixHandler::SubFrom(MatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   for (integer i = iNumItems; i-- > 0; ) {
      ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
      ASSERT(piCol[i] > 0 && piCol[i] <= MH.iGetNumCols());
      MH.fDecCoef(piRow[i], piCol[i], pdMat[i]);
   }
   
   return MH;
}


/* sottrae la matrice da un FullMatrixHandler */
MatrixHandler& SparseSubMatrixHandler::SubFrom(FullMatrixHandler& MH) const {
#ifdef DEBUG
   IsValid();
   MH.IsValid();
#endif // DEBUG
   
   doublereal* pdTo = MH.pdGetMat();
   integer iNumR = MH.iGetNumRows();
   
   for (integer i = iNumItems; i-- > 0; ) {
      ASSERT(piRow[i] > 0 && piRow[i] <= MH.iGetNumRows());
      ASSERT(piCol[i] > 0 && piCol[i] <= MH.iGetNumCols());
      pdTo[(piRow[i]-1)+iNumR*(piCol[i]-1)] -= pdMat[i];
   }
   
   return MH;
}

/* SparseSubMatrixHandler - end */


/* MySubVectorHandler - begin */

MySubVectorHandler::MySubVectorHandler(integer iSize) 
: MyVectorHandler(), piRow(NULL), piRowm1(NULL) {
   Resize(iSize);
}

MySubVectorHandler::MySubVectorHandler(integer iSize, integer* piTmpRow, 
				       doublereal* pdTmpVec)
: MyVectorHandler(iSize, pdTmpVec), piRow(piTmpRow), piRowm1(piTmpRow-1) {
#ifdef DEBUG
   IsValid();
#else
   NO_OP;
#endif
}
   
void MySubVectorHandler::Resize(integer iSize) 
{
   if (iSize < 0) {
      std::cerr << "Negative size!" << std::endl;
      THROW(ErrGeneric());
   }
   
   ASSERT((piRow == NULL && pdVec == NULL) 
	  || (piRow != NULL && pdVec != NULL));
   
   if (!fOwnsMemory && piRow != NULL) {
      if (iSize > iMaxSize) {
         std::cerr << "Can't resize to " << iSize 
	   << ": larger than max size " << iMaxSize << std::endl;
	 THROW(ErrGeneric());
      }
      iCurSize = iSize;
   } else {
      if (piRow != NULL) {
	 if (iSize < iMaxSize) {
	    iCurSize = iSize;
	 } else {
	    doublereal* pd = NULL;
	    SAFENEWARR(pd, doublereal, iSize);
	    integer* pi = NULL;
	    SAFENEWARR(pi, integer, iSize);
	    for (integer i = iCurSize; i-- > 0; ) {
	       pd[i] = pdVec[i];
	       pi[i] = piRow[i];
	    }
	    SAFEDELETEARR(pdVec);
	    SAFEDELETEARR(piRow);
	    pdVec = pd;
	    pdVecm1 = pdVec-1;
	    piRow = pi;
	    piRowm1 = piRow-1;
	    iMaxSize = iCurSize = iSize;
	 }
      } else {
	 if (iSize > 0) {
	    SAFENEWARR(pdVec, doublereal, iSize);
	    SAFENEWARR(piRow, integer, iSize);
	    pdVecm1 = pdVec-1;
	    piRowm1 = piRow-1;
	    iMaxSize = iCurSize = iSize;
	 }	 
      }
   }
}

void MySubVectorHandler::Detach(void) 
{
   if (fOwnsMemory) {
      if (pdVec != NULL) {
	 SAFEDELETEARR(pdVec);
	 SAFEDELETEARR(piRow);
      }
      fOwnsMemory = 0;
   }
   iMaxSize = iCurSize = 0;
   pdVec = pdVecm1 = NULL;
   piRow = piRowm1 = NULL;
}

void MySubVectorHandler::Attach(integer iSize, doublereal* pd, 
				integer* pi, integer iMSize) 
{
   if (fOwnsMemory && pdVec != NULL) {
      Detach();
      fOwnsMemory = 0;
   }
   iMaxSize = iCurSize = iSize;
   if (iMSize >= iSize) {
      iMaxSize = iMSize;
   }
   pdVec = pd;
   pdVecm1 = pd-1;
   piRow = pi;
   piRowm1 = piRow-1;
}

void MySubVectorHandler::IsValid(void) const
{
   MyVectorHandler::IsValid();
   
   ASSERT(piRow != NULL);
   ASSERT(piRowm1 < piRow);
   ASSERT(piRowm1+1 == piRow);	
   
#ifdef DEBUG_MEMMANAGER  
   ASSERT(defaultMemoryManager.fIsValid(piRow, 
			   MyVectorHandler::iMaxSize*sizeof(integer)));
#endif /* DEBUG_MEMMANAGER */
}

VectorHandler& MySubVectorHandler::AddTo(VectorHandler& VH) const {
#ifdef DEBUG 
   IsValid();
   VH.IsValid();
#endif
   
   for (integer i = iGetSize(); i > 0; i--) {
      VH.fIncCoef(piRowm1[i], pdVecm1[i]);
   }
   return VH;
}

VectorHandler& MySubVectorHandler::AddTo(MyVectorHandler& VH) const {
#ifdef DEBUG 
   IsValid();
   VH.IsValid();
#endif
   
   doublereal* pdm1 = VH.pdGetVec()-1;
   for (integer i = iGetSize(); i > 0; i--) {
      pdm1[piRowm1[i]] += pdVecm1[i];
   }
   return VH;
}


std::ostream& 
operator << (std::ostream& out, const SubVectorHandler& v)
{
#ifdef DEBUG
   v.IsValid();
#endif
   
   integer iRow = v.iGetSize();
   
   ASSERT(iRow > 0);
   
   for (integer i = 1; i <= iRow; i++) {
      out << std::setw(12) << v.iGetRowIndex(i) 
	      << std::setw(12) << v.dGetCoef(i) << std::endl;
   }   
      
   return out << std::endl;
}

/* MySubVectorHandler - end */
