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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/*****************************************************************************
 *                                                                           *
 *                            HARWLIB C++ WRAPPER                            *
 *                                                                           *
 *****************************************************************************/

/* Pierangelo Masarati */

/* Wrapper for Harwell sparse matrix managing and LU solution 
 * C/C++ porting by means of f2c */

/* Uso:
 * l'HSLUSolutionManager rende disponibile lo spazio per la matrice, 
 * la soluzione ed il termine noto. Inoltre rende disponibili
 * metodi per l'handling della matrice e dei vettori 
 * (inserimento e lettura dei coefficienti). Infine consente la soluzione del 
 * problema lineare mediante fattorizzazione LU e successiva sostituzone.
 * 
 * Uso consigliato:
 * // creare un oggetto HSLUSolutionManager:
 *         HSLUSolutionManager SM(size, workspacesize, pivotfactor);
 * 
 * // ottenere l'handling alla matrice ed al termine noto
 *         SparseMatrixHandler* pMH = SM.pMatHdl();
 *         VectorHandler* pRH = SM.pResHdl();
 * 
 * // Ogni volta che si desidera riassemblare la matrice:
 * // inizializzare lo SparseMatrixHandler
 *         SM.MatrInit();
 * 
 * // Operare sulla matrice e sui vettori
 *         pMH->fPutCoef(row, col, coef);
 *         coef = pMH->dGetCoef(row, col);
 *         pRH->fPutCoef(row, coef);
 *         coef = pRH->dGetCoef(row);
 * 
 * // Risolvere il problema; in questa fase si assume che: 
 * //  - la matrice sia stata completata;
 * //  - il termine noto sia stato completato. 
 * // In tale caso viene eseguita la fattorizzazione e quindi la soluzione.
 * // Se la fattorizzazione e' gia' stata compiuta, la matrice non e' piu'
 * // stata modificata e si e' modificato il termine noto, la chiamata a 
 * // HSLUSolutionManager::Solve()
 * // semplicemente esegue una nuova sostituzione all'indietro.
 *         SM.Solve();
 * 
 * // Se si desidera modificare la matrice per una nuova soluzione, occorre
 * // inizializzare di nuovo lo SparseMatrixHandler, con:
 *         SM.MatrInit();
 * 
 * // Per i parametri di inizializzazione e per eventuali modifiche
 * // fare riferimento al codice sorgente ed alla libreria originaria
 * // in FORTRAN od al suo porting in C, files <harwlib.f>, <harwlib.c>
 * // Il pacchetto e' costituito dai files:
 * 
 *         <harwlib.f>  libreria originaria in FORTRAN
 *         <harwlib.c>  conversione in C mediante f2c
 *         <harwlib.h>  header per <harwlib.c> con dichiarazione delle funzioni 
 *                      e dei common
 *         <harwrap.cc> wrapper in C++
 *         <harwrap.h>  header per <harwrap.cc>. 
 *                      E' sufficiente includere questo per usare il wrapper.
*/


#ifndef HARWRAP_H
#define HARWRAP_H


/* include di debug */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

/* include generali */

/* include del programma */
#include "harwlib.h"
#include "solman.h"
#include "submat.h"



/* classi dichiarate in questo file */

class SparseData;           /* gestore di sparsita' */
class SparseMatrixHandler;      /* gestore matrice sparsa (assemblaggio) */
class VectorHandler;        /* gestore vettore (pieno, assemblaggio) */
class HarwellLUSolver;      /* solutore */
class HSLUSolutionManager;  /* gestore della soluzione */

#ifdef USE_SCHUR
class SchurSolutionManager;
#endif /* USE_SCHUSR */

/* Memory managers per le classi che allocano memoria */

#ifdef DEBUG_MEMMANAGER
extern clMemMan SMmm;  /* memory manager per SolutionManager (vedi mynewmem.h) */
extern clMemMan LUmm;  /* memory manager per HarwellLUSolver */
extern clMemMan MHmm;  /* memory manager per SparseMatrixHandler */
#endif


/* Unione usata per impaccare il vettore di indici di riga e colonna
 * nel gestore di sparsita' (i vettori sono integer = long int, mentre 
 * in questo modo gli indici massimi sono 65535) */			     

union uPacVec { 
   integer iInt;
   struct { 
      unsigned short int ir; 
      unsigned short int ic; 
   } sRC;
};


/* SparseData - begin */

/* Gestore di entita' sparse con hash + linked list. 
 * usa spazio di lavoro messo a disposizione da altri. */

class SparseData {
   friend class HSLUSolutionManager;
   friend class SparseMatrixHandler;

 public: 
   class ErrNoRoom {};
 private:
 
 protected:
   integer iMaxSize;   /* Dimensione della memoria allocata */
   integer iCurSize;   /* Dimensione corrente */
   integer iNumItem;   /* Numero di entries correnti */
   integer** ppiTable; /* Vettore della linked list */
   integer** ppiKeys;  /* Vettore delle keys */
   
   /* Inizializza i vettori con il valore -(iSize+1) ecc. */
   void ResetVec(void) {
      ASSERT(iMaxSize > 0);
      ASSERT(iCurSize > 0);
      ASSERT(ppiTable != NULL);
      ASSERT(ppiKeys != NULL);
      ASSERT(*ppiTable != NULL);
      ASSERT(*ppiKeys != NULL);
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiTable));
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiKeys));
#endif
      
      integer iSize = iCurSize;
      __FC_DECL__(kd01a)(&iSize, *ppiTable, *ppiKeys);	
   };
   
 public:
   
   /* Costruttore (banale) */
   /* Viene chiamato da chi lo utilizza; gli vengono passati:
    * iSize -   dimensione degli spazi di lavoro
    * piTable - puntatore a punt. ad un array di interi di dimensioni iSize, 
    *           che viene usato come tabella di stato del campo,
    * piKeys  - puntatore a punt. ad un array delle stesse dimensioni del 
    *           precedente, che viene usato per contenere le keys
    */
   SparseData(integer iSize, integer** ppiTmpTable, integer** ppiTmpKeys) 
     : iMaxSize(iSize), iCurSize(iSize), iNumItem(0),
     ppiTable(ppiTmpTable), ppiKeys(ppiTmpKeys) {
	ASSERT(iCurSize > 0);
	ASSERT(ppiTable != NULL);
	ASSERT(ppiKeys != NULL);
	ASSERT(*ppiTable != NULL);
	ASSERT(*ppiKeys != NULL);
#ifdef DEBUG_MEMMANAGER	
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiTable));
	ASSERT(SMmm.fIsPointerToBlock((void*)*ppiKeys));
#endif	
     };
   
   /* Distruttore (banalissimo): non c'e' nulla da distruggere  */
   ~SparseData(void) { 
      NO_OP;
   };
   
   /* Trova la posizione del termine dato da iKey */
   /* Ritorna un termine che rappresenta la posizione della key di ingresso
    * nel vettore piKey; 
    * - se il termine e' positivo, la key esiste ed e' valida, 
    *   oppure l'inserzione e' riuscita correttamente;
    * - se il termine e' negativo, la key e' stata definita e poi cancellata 
    *   (attualmente non supportato)
    * - se il termine e' zero l'inserzione non e' riuscita 
    *   (ad esempio perche' si e' esaurito lo spazio)
    */
   integer iGetIndex(integer iKey) {
      ASSERT(ppiTable != NULL);
      ASSERT(ppiKeys != NULL);
      ASSERT(*ppiTable != NULL);
      ASSERT(*ppiKeys != NULL);
      ASSERT(iKey > 0);
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiTable));
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiKeys));
#endif
      
      integer iFree = 0;
      __FC_DECL__(kd01b)(*ppiTable, *ppiKeys, &iKey, &iFree);
      if (iFree == 0) {
	 cerr << "SparseData: there's no room left in matrix" << endl;	 
	 THROW(SparseData::ErrNoRoom());
      }
      
      return iFree;
   };
};

/* SparseData - end */


/* SparseMatrixHandler - begin */

/* Gestore di matrici sparse; usa spazio messo a disposizione da altri;
 * usa il gestore di sparsita' SparseData */
 
class SparseMatrixHandler : public MatrixHandler {
   friend class HSLUSolutionManager;
 
 private:
   integer iWorkSpaceSize; /* dimensione del workspace */
   integer iCurSize; /* dimensione corrente del workspace */
   integer iNumItem; /* numero di elementi effettivamente contenuti */
   
   integer iMatSize; /* Ordine della matrice (supposta quadrata) */

   SparseData* pHS;  /* Puntatore all'oggetto SparseData da utilizzare */
   
   const doublereal dZero;
   
 protected:
   integer** ppiRow; /* puntatori a punt. ad array di interi di dimensione iMaxSize */
   integer** ppiCol; /* '' '' */
   
   doublereal** ppdMat; /* puntatore a punt. ad array di reali di dimensione iMaxSize */

 public:
   
   /* Costruttore; gli vengono passati:
    * iMSize   - ordine della matrice;
    * ppiRow   - puntatore ad array di interi di dimensione iSize,
    *            contiene la tabella usata da SparseData;
    * ppiCol   - puntatore ad array di interi di dimensione iSize,
    *            contiene le keys usate da SparseData
    * ppdMat   - puntatore ad array di reali, contiene la matrice;
    * iWSSize  - dimensione del workspace
    */
   SparseMatrixHandler(integer iMSize, 
		       integer** ppiTmpRow, 
		       integer** ppiTmpCol, 
		       doublereal** ppdTmpMat, 
		       integer iWSSize);
   
   /* Distruttore banale - non c'e' nulla da distruggere */
   ~SparseMatrixHandler(void);
   
   /* Usato per il debug */
   virtual void IsValid(void) const;

   /* Inizializzazione della matrice: la riempe con il valore desiderato
    * ed inizializza il gestore di sparsita' */
   void Init(const doublereal& dResetVal = 0.);

   /* Restituisce un puntatore all'array di reali della matrice */
   inline doublereal* pdGetMat(void) const { return *ppdMat; };

   /* Restituisce un puntatore al vettore delle righe */
   inline integer* piGetRows(void) const { return *ppiRow; };
   
   /* Restituisce un puntatore al vettore delle colonne */
   inline integer* piGetCols(void) const { return *ppiCol; };

   /* Impacchetta la matrice; restituisce il numero di elementi diversi da zero */
   integer PacMat(void);

   /* Inserisce un coefficiente */
   inline flag fPutCoef(integer iRow, integer iCol, const doublereal& dCoef);
   
   /* Incrementa un coefficiente - se non esiste lo crea */
   inline flag fIncCoef(integer iRow, integer iCol, const doublereal& dCoef);
   
   /* Incrementa un coefficiente - se non esiste lo crea */
   inline flag fDecCoef(integer iRow, integer iCol, const doublereal& dCoef);
   
   /* Restituisce un coefficiente - zero se non e' definito */
   inline const doublereal& dGetCoef(integer iRow, integer iCol) const;

   /* dimensioni */
   virtual integer iGetNumRows(void) const {
      return iCurSize;
   };
   
   virtual integer iGetNumCols(void) const {
      return iCurSize;
   };   
};



inline flag SparseMatrixHandler::fPutCoef(integer iRow, 
				      	  integer iCol, 
					  const doublereal& dCoef)
{
#ifdef DEBUG
   IsValid();
#endif
    
  
   ASSERT((iRow > 0) && (iRow <= iMatSize));
   ASSERT((iCol > 0) && (iCol <= iMatSize));

   if (dCoef != doublereal(0.) ) {
      union uPacVec uPV;
      uPV.sRC.ir = (unsigned short int)(iRow);
      uPV.sRC.ic = (unsigned short int)(iCol);
      
      integer iField = uPV.iInt;
      integer iReturnFlag = pHS->iGetIndex(iField);
      iReturnFlag = abs(iReturnFlag)-1;
      (*ppdMat)[iReturnFlag] = dCoef;
      
      return flag(0);	
   }
   
   return flag(1);
}


inline flag SparseMatrixHandler::fIncCoef(integer iRow, 
				      	  integer iCol, 
					  const doublereal& dCoef)
{
#ifdef DEBUG
   IsValid();
#endif
  
   ASSERT((iRow > 0) && (iRow <= iMatSize));
   ASSERT((iCol > 0) && (iCol <= iMatSize));

   if (dCoef != doublereal(0.)) {
      union uPacVec uPV;
      uPV.sRC.ir = (unsigned short int)(iRow);
      uPV.sRC.ic = (unsigned short int)(iCol);
      
      integer iField = uPV.iInt;
      integer iReturnFlag = pHS->iGetIndex(iField);
      iReturnFlag = abs(iReturnFlag)-1;
      (*ppdMat)[iReturnFlag] += dCoef;
      
      return flag(0);	
   }
   
   return flag(1);
}


inline flag SparseMatrixHandler::fDecCoef(integer iRow,
				      	  integer iCol, 
					  const doublereal& dCoef)
{
#ifdef DEBUG
   IsValid();
#endif

   ASSERT((iRow > 0) && (iRow <= iMatSize));
   ASSERT((iCol > 0) && (iCol <= iMatSize));

   if (dCoef != doublereal(0.)) {
      union uPacVec uPV;
      uPV.sRC.ir = (unsigned short int)(iRow);
      uPV.sRC.ic = (unsigned short int)(iCol);
      
      integer iField = uPV.iInt;
      integer iReturnFlag = pHS->iGetIndex(iField);
      iReturnFlag = abs(iReturnFlag)-1;
      (*ppdMat)[iReturnFlag] -= dCoef;
      
      return flag(0);	
   }
   
   return flag(1);
}


inline const doublereal& 
SparseMatrixHandler::dGetCoef(integer iRow, integer iCol) const
{
#ifdef DEBUG
   IsValid();
#endif

   ASSERT((iRow > 0) && (iRow <= iMatSize));
   ASSERT((iCol > 0) && (iCol <= iMatSize));

   union uPacVec uPV;
   uPV.sRC.ir = short(iRow);
   uPV.sRC.ic = short(iCol);
	
   integer iField = uPV.iInt;
   integer iReturnFlag = pHS->iGetIndex(iField);
   if (iReturnFlag != 0) {
      iReturnFlag = abs(iReturnFlag)-1;
      return (*ppdMat)[iReturnFlag];
   }

   return dZero; // zero!
}

/* SparseMatrixHandler - end */


/* VectorHandler */

/* Gestore di vettori. Usa spazio messo a disposizione da altri. */
/* Definito in <solman.h> */


/* HarwellLUSolver - begin*/

/* Solutore LU per matrici sparse. usa spazio messo a disposizione da altri 
 * e si aspetta le matrici gia' bell'e preparate */

static char sLUClassName[] = "HarwellLUSolver";

class HarwellLUSolver {
   friend class HSLUSolutionManager;

#ifdef USE_SCHUR
   friend class SchurSolutionManager;
#endif /* USE_SCHUR */



 public:
   class ErrFactorisation {
    private: 
      int iErrCode;
    public:
      ErrFactorisation(int i) : iErrCode(i) {};
      int iGetErrCode(void) const { return iErrCode; };
   };
   
 private:
   integer iMatSize;
   integer** ppiRow;
   integer** ppiCol;
   doublereal** ppdMat;
   
   integer iN;         // ordine della matrice
   integer iNonZeroes; // coeff. non nulli
   doublereal* pdRhs;  // Soluzione e termine noto
   
   doublereal dU;      // parametro di pivoting
   integer* piKeep;    // vettore di lavoro
   integer* piW;       // vettore di lavoro
   doublereal* pdW;    // vettore di lavoro
   
 protected:
   
   /* Costruttore: si limita ad allocare la memoria */
   HarwellLUSolver(integer iMatOrd, integer iSize,
		   integer** ppiTmpRow, integer** ppiTmpCol, 
		   doublereal** ppdTmpMat,
		   doublereal* pdTmpRhs, doublereal dPivotFact) :
   iMatSize(iSize), ppiRow(ppiTmpRow), ppiCol(ppiTmpCol), ppdMat(ppdTmpMat),
     iN(iMatOrd), iNonZeroes(0), pdRhs(pdTmpRhs), 
     dU(dPivotFact), piKeep(NULL), piW(NULL), pdW(NULL) {		
	ASSERT(iMatSize > 0);
	ASSERT(ppiRow != NULL);
	ASSERT(ppiCol != NULL);
	ASSERT(ppdMat != NULL);   
	ASSERT(*ppiRow != NULL);
	ASSERT(*ppiCol != NULL);
	ASSERT(*ppdMat != NULL);
	ASSERT(pdRhs != NULL);
	ASSERT(iN > 0);
	
#ifdef DEBUG_MEMMANAGER	
	ASSERT(SMmm.fIsValid((void*)*ppiRow, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppiCol, iMatSize*sizeof(integer)));
	ASSERT(SMmm.fIsValid((void*)*ppdMat, iMatSize*sizeof(doublereal)));
	ASSERT(SMmm.fIsValid((void*)pdRhs, iN*sizeof(doublereal)));
#endif

	SAFENEWARR(piKeep, integer, 5*iN, LUmm);
	SAFENEWARR(piW, integer, 8*iN, LUmm);
	SAFENEWARR(pdW, doublereal, iN, LUmm);
	
#ifdef DEBUG	
	for (int iCnt = 0; iCnt < 5*iN; iCnt++) {
	   *(piKeep+iCnt) = 0;
	}
	for (int iCnt = 0; iCnt < 8*iN; iCnt++) {
	   *(piW+iCnt) = 0;
	}
	for (int iCnt = 0; iCnt < 1*iN; iCnt++) {
	   *(pdW+iCnt) = 0.;
	}
#endif	
     };
   
   /* Distruttore */
   ~HarwellLUSolver(void) {
      if (pdW != NULL) {	     
	 SAFEDELETEARR(pdW, LUmm);
      }	
      if (piW != NULL) {	     
	 SAFEDELETEARR(piW, LUmm);
      }	
      if (piKeep != NULL) {	     
	 SAFEDELETEARR(piKeep, LUmm);
      }	
   };
   
   void IsValid(void) const {
      ASSERT(iMatSize > 0);
      ASSERT(ppiRow != NULL);
      ASSERT(ppiCol != NULL);
      ASSERT(ppdMat != NULL);   
      ASSERT(*ppiRow != NULL);
      ASSERT(*ppiCol != NULL);
      ASSERT(*ppdMat != NULL);
      ASSERT(pdRhs != NULL);
      ASSERT(iN > 0);
      
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsValid((void*)*ppiRow, iMatSize*sizeof(integer)));
      ASSERT(SMmm.fIsValid((void*)*ppiCol, iMatSize*sizeof(integer)));
      ASSERT(SMmm.fIsValid((void*)*ppdMat, iMatSize*sizeof(doublereal)));
      ASSERT(SMmm.fIsValid((void*)pdRhs, iN*sizeof(doublereal)));
#endif
      
      ASSERT(piKeep != NULL);
      ASSERT(piW != NULL);
      ASSERT(pdW != NULL);
      
#ifdef DEBUG_MEMMANAGER
      ASSERT(LUmm.fIsBlock((void*)piKeep, 5*iN*sizeof(integer)));
      ASSERT(LUmm.fIsBlock((void*)piW, 8*iN*sizeof(integer)));
      ASSERT(LUmm.fIsBlock((void*)pdW, 1*iN*sizeof(doublereal)));	
#endif
   };
   
   /* Fattorizza la matrice */
   flag fLUFactor(void) {
#ifdef DEBUG
      IsValid();
#endif
      
      /*
      ASSERT(iMatSize > 0);
      ASSERT(ppiRow != NULL);
      ASSERT(ppiCol != NULL);
      ASSERT(ppdMat != NULL);
      ASSERT(*ppiRow != NULL);
      ASSERT(*ppiCol != NULL);
      ASSERT(*ppdMat != NULL);   
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiRow));
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiCol));
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppdMat));
#endif
      ASSERT(iN > 0);
       */
	
      ASSERT(iNonZeroes > 0);
      
      /*
      ASSERT(piKeep != NULL);
      ASSERT(piW != NULL);
      ASSERT(pdW != NULL);
#ifdef DEBUG_MEMMANAGER	
      ASSERT(LUmm.fIsPointerToBlock((void*)piKeep));
      ASSERT(LUmm.fIsPointerToBlock((void*)piW));
      ASSERT(LUmm.fIsPointerToBlock((void*)pdW));
#endif
      */
      
      integer iLicn = iMatSize;
      integer iLirn = iMatSize;
      integer iFlag = 0;
      
      DEBUGCOUT("Calling ma28ad_()," << endl
		<< "iN         = " << iN << endl
		<< "iNonZeroes = " << iNonZeroes << endl
		<< "pdMat      = " << *ppdMat << endl
		<< "iLicn      = " << iLicn << endl
		<< "piRow      = " << *ppiRow << endl
		<< "iLirn      = " << iLirn << endl
		<< "piCol      = " << *ppiCol << endl
		<< "dU         = " << dU << endl
		<< "piKeep     = " << piKeep << endl
		<< "piW        = " << piW << endl
		<< "iFlag      = " << iFlag << endl);
      
      __FC_DECL__(ma28ad)(&iN, &iNonZeroes, *ppdMat, &iLicn, *ppiRow, &iLirn, 
	      *ppiCol, &dU, piKeep, piW, pdW, &iFlag);
      
      if (iFlag < 0) { 
	 cerr << sLUClassName 
	   << ": error during factorisation, code " << iFlag << endl;	 
	 THROW(HarwellLUSolver::ErrFactorisation(iFlag));
      }  
      
      return iFlag;		
   };      
   
   /* Risolve */
   void Solve(void) {
#ifdef DEBUG
      IsValid();
#endif
      
      /*
      ASSERT(iMatSize > 0);
      ASSERT(ppiCol != NULL);
      ASSERT(ppdMat != NULL);
      ASSERT(*ppiCol != NULL);
      ASSERT(*ppdMat != NULL);
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppiCol));
      ASSERT(SMmm.fIsPointerToBlock((void*)*ppdMat));
#endif
      ASSERT(iN > 0);
      ASSERT(pdRhs != NULL);
      ASSERT(piKeep != NULL);
#ifdef DEBUG_MEMMANAGER	
      ASSERT(SMmm.fIsPointerToBlock((void*)pdRhs));
      ASSERT(LUmm.fIsPointerToBlock((void*)piKeep));
#endif
       */
   		
      integer iLicn = iMatSize;
      integer iMtype = 1;
      __FC_DECL__(ma28cd)(&iN, *ppdMat, &iLicn, *ppiCol, piKeep, pdRhs, pdW, &iMtype);
   };
   
};


/* HSLUSolutionManager - begin */

/* Gestisce la soluzione del problema; alloca le matrici occorrenti
 * e gli oggetti dedicati alla gestione delle matrici ed alla soluzione */

class HSLUSolutionManager : public SolutionManager {
 public: 
   class ErrGeneric {};
   
 private:
   
 protected:
   integer iMatMaxSize;  /* Dimensione massima della matrice (per resize) */
   integer iMatSize;     /* ordine della matrice */
   integer* piRow;       /* puntatore ad array di interi con: tabella di SparseData/indici di riga di HarwellLUSolver */
   integer* piCol;       /* puntatore ad array di interi con: keys di SparseData/indici di colonna di HarwellLUSolver */
   doublereal* pdMat;    /* puntatore ad array di reali con la matrice */
   doublereal* pdVec;    /* punattore ad array di reali con residuo/soluzione */
   
   SparseMatrixHandler* pMH; /* puntatore a SparseMatrixHandler */
   VectorHandler* pVH;   /* puntatore a VectorHandler */
   HarwellLUSolver* pLU; /* puntatore a HarwellLUSolver */
   
   flag fHasBeenReset;   /* flag di matrice resettata */
   
   /* Prepara i vettori e la matrice per il solutore */
   void PacVec(void); 
   
   /* Fattorizza la matrice (non viene mai chiamato direttamente, 
    * ma da Solve se la matrice ancora non e' stata fattorizzata) */
   void Factor(void);
   
 public:   
   /* Costruttore: usa e mette a disposizione matrici che gli sono date */
   HSLUSolutionManager(integer iSize, integer iWorkSpaceSize = 0,
		       const doublereal& dPivotFactor = 1.0);
   
   /* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
   ~HSLUSolutionManager(void);
   
   /* Usata per il debug */
   void IsValid(void) const;
   
   /* Ridimensiona le matrici */
   void MatrResize(integer iNewSize) {
   	THROW(ErrNotImplementedYet());
   };
   
   /* Inizializza il gestore delle matrici */
   void MatrInit(const doublereal& dResetVal = 0.);
   
   /* Risolve il sistema */
   void Solve(void);
   
   /* Rende disponibile l'handler per la matrice */
   MatrixHandler* pMatHdl(void) const {
      ASSERT(pMH != NULL);	
      return pMH;
   };
   
   /* Rende disponibile l'handler per il termine noto */
   VectorHandler* pResHdl(void) const {
      ASSERT(pVH != NULL);	
      return pVH;
   };

   /* Rende disponibile l'handler per la soluzione (e' lo stesso 
    * del termine noto, ma concettualmente sono separati) */
   VectorHandler* pSolHdl(void) const {
      ASSERT(pVH != NULL);	
      return pVH;
   };  
};

/* HSLUSolutionManager - end */

#endif
