/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
#include <ac/iostream>
#include <ac/iomanip>

#include <fullmh.h>
#include <submat.h>

/* FullMatrixHandler - begin */

void FullMatrixHandler::Init(const doublereal& dResetVal)
{
   IsValid();
#ifdef HAVE_MEMSET
   if (dResetVal == 0.) {
	   memset(pdRaw, 0, iNumCols*iNumRows*sizeof(doublereal));
   } else {
#endif /* HAVE_MEMSET */
   for (integer i = iNumCols*iNumRows; i-- > 0; ) {
      pdRaw[i] = dResetVal;
   }
#ifdef HAVE_MEMSET
   }
#endif /* HAVE_MEMSET */
}

void FullMatrixHandler::CreateColRow(integer iNR, integer iNC)
{
   for (integer i = iNC; i-- > 0; ) {
      ppdCols[i] = pdRawm1+i*iNR;
   }
}


FullMatrixHandler::FullMatrixHandler(doublereal* pd, doublereal** ppd, 
				     integer iSize, integer iNR, integer iNC,
				     integer iMaxC)
: fOwnsMemory(0), 
iNumRows(iNR), iNumCols(iNC), iRawSize(iSize), iMaxCols(iMaxC), 
pdRaw(pd), pdRawm1(pd-1),
ppdCols(ppd), ppdColsm1(ppd-1)
{   
   if (iMaxCols <= 0) {
      iMaxCols = iNumCols;
   }
   
   CreateColRow(iNumRows, iNumCols);
   IsValid();
}

/* costruttore che si alloca la memoria */
FullMatrixHandler::FullMatrixHandler(integer iNR, integer iNC)
: fOwnsMemory(1), 
iNumRows(iNR), iNumCols(iNC ? iNC : iNumRows), 
iRawSize(iNumRows*iNumCols), iMaxCols(iNumCols),
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL)
{
   ASSERT(iNumRows > 0);
   ASSERT(iNumCols > 0);
   Resize(iNumRows, iNumCols);
   IsValid();
}


/* costruttore che si alloca la memoria */
FullMatrixHandler::FullMatrixHandler(integer iNR, integer iNC, const doublereal& dVal)
: fOwnsMemory(1), 
iNumRows(iNR), iNumCols(iNC),
iRawSize(iNumRows*iNumCols), iMaxCols(iNumCols),
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL)
{
   ASSERT(iNumRows > 0);
   ASSERT(iNumCols > 0);
   Resize(iNumRows, iNumCols);
   Reset(dVal);
   IsValid();
}


/* costruttore che non fa nulla */
FullMatrixHandler::FullMatrixHandler(void) 
: fOwnsMemory(1), 
iNumRows(0), iNumCols(0), iRawSize(0), iMaxCols(0), 
pdRaw(NULL), pdRawm1(NULL),
ppdCols(NULL), ppdColsm1(NULL)
{
   NO_OP;
}


FullMatrixHandler::~FullMatrixHandler(void) 
{
   Detach();
}

/* ridimensiona la matrice (se possiede la memoria) */
void FullMatrixHandler::Resize(integer iNewRows, integer iNewCols)
{
   if (iNewRows < 0 || iNewCols < 0) {
      std::cerr << "Negative size!" << std::endl;
      THROW(ErrGeneric());
   }
   
   integer iSize = iNewRows*iNewCols;
   if (fOwnsMemory) {
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
	    pdRawm1 = pd-1;	    
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
	       std::cerr << "Can't resize to " << iSize 
		 << ": larger than max size " << iRawSize << std::endl;
	       THROW(ErrGeneric());
	    } else if (iNewCols > iMaxCols) {
	       std::cerr << "Can't resize to " << iNewCols
		 << " cols: larger than max num cols " 
		 << iMaxCols << std::endl;
	       THROW(ErrGeneric());
	    }
	 } else {
	    /* aggiorna */
	    if (iNewRows != iNumRows || iNewCols != iNumCols) {
	       CreateColRow(iNewRows, iNewCols);
	    }
	    iNumRows = iNewRows;
	    iNumCols = iNewCols;
	 }
      } else {
	 /* errore */
	 std::cerr << "internal error!" << std::endl;
	 THROW(ErrGeneric());
      }      
   }
}


/* si stacca dalla memoria a cui e' associato */
void FullMatrixHandler::Detach(void)
{
   if (fOwnsMemory) {
      if (pdRaw != NULL) {
	 SAFEDELETEARR(pdRaw);
      }
      if (ppdCols != NULL) {
	 SAFEDELETEARR(ppdCols);
      }
      fOwnsMemory = 0;
   }
   iRawSize = iMaxCols = iNumRows = iNumCols = 0;
   pdRaw = pdRawm1 = NULL;
   ppdCols = ppdColsm1 = NULL;
   
}


/* Attacca un nuovo array, con n. righe, n. colonne e dim. massima;
 * se assente, assunta = nrighe*ncolonne */
void FullMatrixHandler::Attach(integer iNewRows, integer iNewCols,
			       doublereal* pd, doublereal** ppd, 
			       integer iMSize, integer iMaxC)
{   
   if (fOwnsMemory || pdRaw != NULL) {
      Detach();
      fOwnsMemory = 0;
   }
   iNumRows = iNewRows;
   iNumCols = iNewCols;
   if (iMSize > 0) {
      if (iMSize < iNumRows*iNumCols) {
	 THROW(ErrGeneric());
      }
      iRawSize = iMSize;
   } else {
      iRawSize = iNumRows*iNumCols;
   }
   if (iMaxC > 0) {
      if (iMaxC < iNumCols) {
	 THROW(ErrGeneric());
      }
      iMaxCols = iMaxC;
   } else {
      iMaxCols = iNumCols;
   }
   pdRaw = pd;
   pdRawm1 = pdRaw-1;
   ppdCols = ppd;
   ppdColsm1 = ppdCols-1;
   
   IsValid();
}








/* Usata per il debug */
void FullMatrixHandler::IsValid(void) const
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

   
std::ostream& operator << (std::ostream& out, const FullMatrixHandler& m)
{
#ifdef HAVE_FMTFLAGS_IN_IOS
   std::ios::fmtflags oldbits = out.setf(std::ios::scientific);
#else /* !HAVE_FMTFLAGS_IN_IOS */
   long oldbits = out.setf(ios::scientific);
#endif /* !HAVE_FMTFLAGS_IN_IOS */
   
   out << "<FullMatrixHandler> n. rows: " << m.iNumRows 
     << "; n. cols: " << m.iNumCols << std::endl;  
   for (int i = 1; i <= m.iNumRows; i++) {
      out << "Row " << std::setw(4) << i;
      for (int j = 1; j <= m.iNumCols; j++) {
	 out << std::setw(10) << std::setprecision(2) << m.ppdColsm1[j][i];
      }
      out << std::endl;
   }
   
   out.flags(oldbits);
   return out;
}


/* Overload di += usato per l'assemblaggio delle matrici */
MatrixHandler& FullMatrixHandler::operator +=(const SubMatrixHandler& SubMH) {
   return SubMH.AddTo(*this);
}

/* Overload di -= usato per l'assemblaggio delle matrici */
MatrixHandler& FullMatrixHandler::operator -=(const SubMatrixHandler& SubMH) {
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
   
#if 0
   for (integer iRow = iNumRows; iRow-- > 0; ) {
      for (integer iCol = iNumCols; iCol-- > 0; ) {
	 ppdCols[iCol][iRow] = 0.;
	 for (integer iK = iN; iK-- > 0; ) {
	    ppdCols[iCol][iRow] += m1.ppdCols[iK][iRow]*m2.ppdCols[iCol][iK];
	 }
      }
   }
#endif
   
   for (integer iRow = 1; iRow <= iNumRows; iRow++) {
      for (integer iCol = 1; iCol <=iNumCols; iCol++) {
	 ppdColsm1[iCol][iRow] = 0.;
	 for (integer iK = 1; iK <= iN; iK++) {
	    ppdColsm1[iCol][iRow] += m1.dGetCoef(iRow,iK)*m2.dGetCoef(iK,iCol);
	 }
      }
   }
}


/* FullMatrixHandler - end */

