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

/* Sottomatrici */


#ifndef SUBMAT_H
#define SUBMAT_H


#include "myassert.h"
#include "except.h"


/* include del programma */
/* si assume che solman.h includa piu' o meno direttamnte f2c.h, 
 * che contiene le dichiarazioni dei tipi derivati dal fortran. */

#include "solman.h"
#include "fullmh.h"
#include "matvec3.h"
#include "matvec3n.h"

/* SubMatrixHandler - begin */    

/**
 Classe virtuale delle sottomatrici. 
 Le SubMatrixHandler sono matrici dotate di vettori di incidenza
 per righe e colonne. 
 Sono usate per scrivere le sottomatrici di ogni elemento, che poi si
 sommano con apposite routines alla matrice jacobiana.
 */

class SubMatrixHandler : public MatrixHandler {   
 protected:
   
 public:
   /**@name Costruttori */
   //@{
   
   /** 
    Distruttore virtuale. 
    Necessario per la corretta distruzione degli oggetti derivati 
    */
   virtual ~SubMatrixHandler(void);
   //@}

   /**@name Metodi di servizio */
   //@{
   /**
    Routine di verifica della validita' dell'oggetto.
    Usata per il debug.
    */
   virtual void IsValid(void) const = 0;
   //@}
   
   /**@name Inizializzazione */
   //@{
   
   /**
    Inizializza la matrice con d 
    */
   virtual void Init(const doublereal&) = 0;
      
   /** 
    Ridimensiona la matrice.
    Nota: nell'implementazione corrente le dimensioni possono essere solamente
    inferiori alle dimensioni massime con cui e' stata dimensionata.
    */
   virtual void Resize(integer, integer) = 0;
      
   /**
    Ridimensiona ed inizializza.
    Combina le due funzioni precedenti in una chiamata.
    */
   virtual void ResizeInit(integer, integer, const doublereal&) = 0;
   //@}
   
   /**@name Gestione dei vettori di incidenza */
   //@{
      
   /** 
    Scrive l'indice di riga.
    */
   virtual inline flag fPutRowIndex(integer, integer) = 0;
      
   /** 
    Scrive l'indice di colonna.
    */
   virtual inline flag fPutColIndex(integer, integer) = 0;
   
   /** 
    Ottiene l'indice di riga.
    */
   virtual inline integer iGetRowIndex(integer) = 0;

   /** 
    Ottiene l'indice di colonna.
    */
   virtual inline integer iGetColIndex(integer) = 0;
   //@}

   /**@name Funzioni di interazione con le matrici */
   //@{
      
   /**
    Si somma ad una matrice.
    Nota: le dimensioni devono essere compatibili.
    */
   virtual MatrixHandler& AddTo(MatrixHandler& MH) const = 0;
      
   /** 
    Si sottrae da una matrice.
    Nota: le dimensioni devono essere compatibili.
    */
   virtual MatrixHandler& SubFrom(MatrixHandler& MH) const = 0;   
};
    
/* SubMatrixHandler - end */    
    
    
/* FullSubMatrixHandler */
    
/** 
 Sottomatrice piena. E' costituita da un vettore di interi, piRow, che 
 contiene i due vettori di incidenza, e da un vettore di reali, pdMat, 
 che contiene la matrice.
 Il vettore di incidenza delle righe e' piRow, mentre quello delle colonne
 e' piCol = piRow+iNumRows. Le dimensioni possono essere cambiate con
 Resize(), con i vincoli che la somma di righe e colonne non deve eccedere
 la lunghezza del vettore di interi, ed il loro prodotto non deve eccedere
 la lunghezza del vettore di reali. 
 */

class FullSubMatrixHandler : public SubMatrixHandler {
   friend ostream& operator << (ostream& out, const FullSubMatrixHandler& m);
   
 public:
   /**@name Errori */
   //@{
   
   /** 
    Errore di ridimensionamento illegale.
    */
   class ErrResize {};
   //@}
   
 protected:
   /** Dimensione totale del vettore di incidenza */
   integer iVecSize;  
   /** Dimensione totale del vettore dei coefficienti */
   integer iMatSize;
   /** Numero di righe correnti */
   integer iNumRows;
   /** Numero di colonne correnti */
   integer iNumCols;
   /** Puntatore al vettore di incidenza delle righe.
    Nota: coincide con il puntatore al vettore di incidenza */
   integer* piRow;
   /** Puntatore al vettore di incidenza delle colonne */
   integer* piCol;
   /** Puntatore al vettore dei coefficienti */
   doublereal* pdMat; 
   
 public:

   /**@name Costruttori */
   //@{
   
   /** 
    Riceve i vettori per gli indici ed i coefficienti, con relative 
    dimensioni
    @param iIntSize    dimensione del vettore degli indici
    @param iDoubleSize dimensione del vettore dei coefficienti
    @param piTmpVec    puntatore al vettore degli indici
    @param pdTmpMat    puntatore al vettore dei coefficienti
    */
   FullSubMatrixHandler(integer iIntSize, integer iDoubleSize,
			integer* piTmpVec, doublereal* pdTmpMat);
   
   /** 
    Distruttore banale. 
    Nota: l'elemento non possiede memoria e quindi non ne dealloca.
    */
   virtual ~FullSubMatrixHandler(void);
   //@}

   /**@name Metodi di servizio */
   //@{
      
   /**
    Routine di verifica della validita' dell'oggetto.
    Usata per il debug.
    */
   virtual void IsValid(void) const;

   /**
    Numero di righe della sottomatrice 
    */
   integer iGetNumRows(void) const {
      return iNumRows;
   };

   /**
    Numero di colonne della sottomatrice 
    */
   integer iGetNumCols(void) const {
      return iNumCols;
   };
   //@}
   
   
   /**@name Metodi di inizializzazione */
   //@{
      
   /**
    Inizializza la porzione utilizzata con il valore desiderato
    */
   void Init(const doublereal& dResetVal = 0.) {
#ifdef DEBUG	
      IsValid();
#endif
      integer i;
      
      for (i = iNumRows*iNumCols; i-- > 0; ) {
	 pdMat[i] = dResetVal;
      }
      
      for (i = iNumRows+iNumCols; i-- > 0; ) {
	 piRow[i] = 0;
      }	
   };   
   
   /**
    Modifica le dimensioni correnti
    */
   void Resize(integer iNewRow, integer iNewCol) {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT(iNewRow > 0);
      ASSERT(iNewCol > 0);
      ASSERT(iNewRow+iNewCol <= iVecSize);
      ASSERT(iNewRow*iNewCol <= iMatSize);		
      
      if (iNewRow <= 0 
	  || iNewCol <= 0 
	  || iNewRow+iNewCol > iVecSize 
	  || iNewRow*iNewCol > iMatSize) {
	 cerr << "FullSubMatrixHandler::Resize() - error" << endl;
	 
	 THROW(FullSubMatrixHandler::ErrResize());
      }
      
      iNumRows = iNewRow;
      iNumCols = iNewCol;
      
      ASSERT(piRow != NULL);
      piCol = piRow+iNewRow;
      
#ifdef DEBUG	
      IsValid();
#endif		
   };
   
   /**
    Ridimensiona ed inizializza.
    */
   void ResizeInit(integer iNewRow, integer iNewCol, const doublereal& dResetVal = 0.) {
      this->Resize(iNewRow, iNewCol);	
      this->Init(dResetVal);
   };

   /**
   Collega la matrice Full alla memoria che gli viene passata in ingresso
    */
   void Attach(int iRows, int iCols, doublereal* pdTmpMat, integer* piTmpIndx) {
     iNumRows = iRows;
     iNumCols = iCols;
     piRow = piTmpIndx;
     piCol =(piRow+iNumRows);
     pdMat = pdTmpMat;
   };
   //@}
   
   /**@name Gestione dei coefficienti */
   //@{
      
   /** 
    Scrive un coefficiente in base ai sottoindici.
    */
   inline flag fPutCoef(integer iSubRow, integer iSubCol, const doublereal& dCoef) {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));
      
      pdMat[(--iSubCol)*iNumRows+(--iSubRow)] = dCoef;
      
      return flag(0);
   };
   
   /** 
    Incrementa un coefficiente in base ai sottoindici.
    */
   inline flag fIncCoef(integer iSubRow, integer iSubCol, const doublereal& dCoef) {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));	
      
      pdMat[(--iSubCol)*iNumRows+(--iSubRow)] += dCoef;
      
      return flag(0);
   };
   
   /** 
    Decrementa un coefficiente in base ai sottoindici.
    */
   inline flag fDecCoef(integer iSubRow, integer iSubCol, const doublereal& dCoef) {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));	
      
      pdMat[(--iSubCol)*iNumRows+(--iSubRow)] -= dCoef;
      
      return flag(0);
   };
   
   /** 
    Ottiene un coefficiente in base ai sottoindici.
    */
   inline const doublereal& dGetCoef(integer iSubRow, integer iSubCol) const {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));	
      
      return pdMat[(--iSubCol)*iNumRows+(--iSubRow)];
   };      
  
   /**
    Scrive un indice di riga 
    */
   inline flag fPutRowIndex(integer iSubRow, integer iRow) {	
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      
      *(piRow+(--iSubRow)) = iRow;
      
      return flag(0);
   };   
   
   /** 
    Scrive un indice di colonna 
    */
   inline flag fPutColIndex(integer iSubCol, integer iCol) {
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));
      
      *(piCol+(--iSubCol)) = iCol;
      
      return flag(0);
   };   
   
   /** 
    Legge un indice di riga 
    */
   inline integer iGetRowIndex(integer iSubRow) {	
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubRow > 0) && (iSubRow <= iNumRows));
      
      return *(piRow+(--iSubRow));
   };   
   
   /** 
    Legge un indice di colonna 
    */
   inline integer iGetColIndex(integer iSubCol) {	
#ifdef DEBUG	
      IsValid();
#endif	
      
      ASSERT((iSubCol > 0) && (iSubCol <= iNumCols));
      
      return *(piCol+(--iSubCol));
   };
   
   /**
    Somma una matrice di tipo Mat3x3 in una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3x3 da sommare
    */
   void Add(integer iRow, integer iCol, const Mat3x3& m);
      
   /**
    Sottrae una matrice di tipo Mat3x3 da una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3x3 da sottrarre
    */
   void Sub(integer iRow, integer iCol, const Mat3x3& m);

   /**
    Scrive una matrice di tipo Mat3x3 in una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3x3 da scrivere
    */
   void Put(integer iRow, integer iCol, const Mat3x3& m);
 
   /**
    Somma una matrice di tipo Mat3xN in una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.   
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3xN
    */
   void Add(integer iRow, integer iCol, const Mat3xN& m);
   
   /**
    Sottrae una matrice di tipo Mat3xN da una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3xN
    */ 
   void Sub(integer iRow, integer iCol, const Mat3xN& m);
   
   /**
    Scrive una matrice di tipo Mat3xN in una data posizione.
    Nota: si assume che nella sottomatrice vi sia spazio per la matrice 3x3.
    Nota: gli indici sono a base 1, in stile FORTRAN.  
    @param iRow indice di riga della sottomatrice da cui iniziare
    @param iCol indice di colonna della sottomatrice da cui iniziare
    @param m    Mat3xN
    */
   void Put(integer iRow, integer iCol, const Mat3xN& m);
   //@}

   /** come sopra, ma per matrici Nx3 **/
    void Add(integer iRow, integer iCol, const MatNx3& m);
    void Sub(integer iRow, integer iCol, const MatNx3& m);

   /**@name Interazione con le matrici */
   //@{
   /** 
    Somma la matrice ad un matrix handler usando i metodi generici 
    */
   MatrixHandler& AddTo(MatrixHandler& MH) const;
   
   /**
    Somma la matrice ad un FullMatrixHandler 
    */
   MatrixHandler& AddTo(FullMatrixHandler& MH) const;

   /**
    Sottrae la matrice da un matrix handler usando i metodi generici 
    */
   MatrixHandler& SubFrom(MatrixHandler& MH) const;

   /**
    Sottrae la matrice da un FullMatrixHandler 
    */
   MatrixHandler& SubFrom(FullMatrixHandler& MH) const;
   //@}
};

/* FullSubMatrixHandler - end */


/* SparseSubMatrixHandler */
    
/** 
 Gestore di sottomatrici sparse, piuttosto rozzo, va usato con cautela.
 E' formato da due vettori di interi e da uno di reali, tutti della stessa 
 lunghezza. 
 Ad ogni indice del sotto-vettore corrisponde un coefficiente con i suoi 
 due indici nella matrice completa.
 La scrittura non e' strutturata, per cui l'utilizzatore deve badare a non
 lasciare vuoti e a non ripetere i coefficienti (in realta' non succede
 nulla se la si usa per assemblare, solo uno stesso coefficiente puo' essere 
 dato da piu' contributi). 
 */

class SparseSubMatrixHandler : public SubMatrixHandler {
   friend class SparseMatrixHandler;
   friend class FullMatrixHandler;
   
 public:
   /**@name Errori */
   //@{
   class ErrResize {};
   //@}
   
 private:
   /** Dimensioni dell'array degli indici */
   integer iIntSize;
   /** Dimensioni dell'array dei coefficienti */
   integer iDoubleSize;
   /** Numero di entries definite */
   integer iNumItems;
   /** Puntatore all'array degli indici di riga.
    Coincide con il puntatore all'array degli inidici, che e' unico */
   integer* piRow;
   /** Puntatore all'array degli indici di colonna */
   integer* piCol;
   /** Puntatore all'array dei coefficienti */
   doublereal* pdMat;
   
 protected:
   
 public:
   /**@name Costruttori */
   //@{
   
   /** Costruttore.
    @param iTmpInt    dimensione dell'array degli indici
    @param iTmpDouble dimensione dell'array dei coefficienti
    @param piTmpIndex puntatore all'array degli indici
    @param pdTmpMat   puntatore all'array dei coefficienti
    */
   SparseSubMatrixHandler(integer iTmpInt, integer iTmpDouble,
			integer* piTmpIndex, doublereal* pdTmpMat)
     : iIntSize(iTmpInt), iDoubleSize(iTmpDouble),
     iNumItems(iTmpInt/2), piRow(piTmpIndex), piCol(piTmpIndex+iNumItems), 
     pdMat(pdTmpMat) {
#ifdef DEBUG
	IsValid();
#endif	
     };
   
   /** Distruttore banale.
    Nota: dato che la classe non possiede la memoria, non ne deve deallocare
    */
   virtual ~SparseSubMatrixHandler(void) { 
      NO_OP;
   };
   //@}
   
   /**@name Metodi di servizio */
   //@{
      
   /**
    Routine di verifica della validita' dell'oggetto.
    Usata per il debug.
    */
   virtual void IsValid(void) const;

   /**
    Numero di righe della sottomatrice.
    Nota: rappresenta il numero totale di entries della sottomatrice.
    */
   integer iGetNumRows(void) const {
      return iNumItems;
   };

   /**
    Numero di colonne della sottomatrice.
    Nota: e' sempre 1, ovvero la matrice e' interpretata come un vettore.
    */
   integer iGetNumCols(void) const {
      return 1;
   };
   //@}
   
   /**@name Metodi di inizializzazione */
   //@{
      
   /**
    Inizializza la matrice con d.
    */
   void Init(const doublereal& dCoef = 0.) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT(iNumItems > 0);	
      
      integer* piTmpRow = piRow;
      integer* piTmpCol = piCol;
      doublereal* pdTmpMat = pdMat;
      
      while (piTmpRow < piRow+iNumItems) {
	 *piTmpRow++ = 0;
	 *piTmpCol++ = 0;
	 *pdTmpMat++ = dCoef;
      }	
   };
   
   /** 
    Ridimensiona la matrice.
    Nota: solo il primo argomento viene considerato, e rappresenta il numero
    totale di entries.
    Questo metodo deve essere chiamato prima di qualsiasi operazione sulla
    matrice.
    */
   void Resize(integer iNewRow, integer /* iNewCol */ ) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT(iNewRow > 0);
      ASSERT(2*iNewRow <= iIntSize);
      ASSERT(iNewRow <= iDoubleSize);
      ASSERT(piRow != NULL);
      
      if (iNewRow <= 0 || 2*iNewRow > iIntSize || iNewRow > iDoubleSize) {
	 cerr << "SparseSubMatrixHandler::Resize() - error" << endl;	 
	 THROW(SparseSubMatrixHandler::ErrResize());
      }		       
      
      iNumItems = iNewRow;
      piCol = piRow+iNumItems;
      
#ifdef DEBUG
      IsValid();
#endif		
   };
   
   /**
    Ridimensiona ed inizializza.
    Unione dei due metodi precedenti 
    */
   void ResizeInit(integer iNewRow, integer iNewCol, const doublereal& dCoef = 0.) {
      this->Resize(iNewRow, iNewCol);
      this->Init(dCoef);
   };

   /** 
    Collega la matrice sparsa alla memoria che gli viene passata in ingresso
    */
    void Attach(int iNumEntr, doublereal* pdTmpMat, integer* piTmpIndx) { 
    
     iIntSize = iNumEntr*2;
     iDoubleSize = iNumEntr;
     iNumItems = iNumEntr;
     piRow = piTmpIndx;
     piCol = piTmpIndx+iNumEntr;
     pdMat = pdTmpMat;
#ifdef DEBUG
     IsValid();
#endif
    };   

   //@}
   
   /**@name Gestione dei coefficienti */
   //@{
         
   /** 
    Scrive un coefficiente in base ai sottoindici.
    */
   inline flag fPutCoef(integer iSubIt, integer /* iDmy */ , const doublereal& dCoef) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));
      
      pdMat[--iSubIt] = dCoef;
      
      return flag(0);
   };
   
   /** 
    Incrementa un coefficiente in base ai sottoindici.
    */
   inline flag fIncCoef(integer iSubIt, integer /* iDmy */ , const doublereal& dCoef) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));   
      pdMat[--iSubIt] += dCoef;
      
      return flag(0);
   };
   
   /** 
    Decrementa un coefficiente in base ai sottoindici.
    */
   inline flag fDecCoef(integer iSubIt, integer /* iDmy */ , const doublereal& dCoef) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));   
      pdMat[--iSubIt] -= dCoef;
      
      return flag(0);
   };
   
   /** 
    Ottiene un coefficiente in base ai sottoindici.
    */
   inline const doublereal& dGetCoef(integer iSubIt, integer /* iDmy */ ) const {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));
      
      return pdMat[--iSubIt];
   };
   
   /**
    Scrive un indice di riga 
    */
   inline flag fPutRowIndex(integer iSubIt, integer iRow) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));      
      *(piRow+(--iSubIt)) = iRow;
      
      return flag(0);
   };   
   
   /**
    Scrive un indice di colonna
    */
   inline flag fPutColIndex(integer iSubIt, integer iCol) {
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));      
      *(piCol+(--iSubIt)) = iCol;
      
      return flag(0);
   };   
   
   /**
    Ottiene un indice di riga 
    */
   inline integer iGetRowIndex(integer iSubIt) {	
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));
      
      return *(piRow+(--iSubIt));
   };   
   
   /**
    Ottiene un indice di colonna
    */
   inline integer iGetColIndex(integer iSubIt) {	
#ifdef DEBUG
      IsValid();
#endif	
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));
      
      return *(piCol+(--iSubIt));
   };   
   
   /**
    Scrive un'entry completa.
    @param iSubIt   sottoindice (numero della entry)
    @param iRow     indice di riga
    @param iCol     indice di colonna
    @param dCoef    coefficiente
    */
   inline flag fPutItem(integer iSubIt, integer iRow, integer iCol, const doublereal& dCoef) {
#ifdef DEBUG
      IsValid();
#endif
      
      ASSERT((iSubIt > 0) && (iSubIt <= iNumItems));
      ASSERT(iRow > 0);
      ASSERT(iCol > 0);
      
      *(pdMat+(--iSubIt)) = dCoef;
      *(piRow+iSubIt) = iRow;
      *(piCol+iSubIt) = iCol;
      
      return flag(0);
   }
   
   /**
    Scrive una matrice prodotto vettore nella posizione assegnata.
    @param iSubIt     sottoindice iniziale (numero della prima entry)
    @param iFirstRow  indice della prima riga della matrice completa
    @param iFirstCol  indice della prima colonna della matrice completa
    @param v          vettore da cui viene calcolata la matrice prodotto 
                      vettore
    */
   flag fPutCross(integer iSubIt, integer iFirstRow, 
		  integer iFirstCol, const Vec3& v);
   
   /**
    Scrive una Mat3x3 nella posizione assegnata.
    @param iSubIt     sottoindice iniziale (numero della prima entry)
    @param iFirstRow  indice della prima riga della matrice completa
    @param iFirstCol  indice della prima colonna della matrice completa
    @param m          matrice da inserire
    */
   flag fPutMat3x3(integer iSubIt, integer iFirstRow, 
		   integer iFirstCol, const Mat3x3& m);
   //@}

   /**@name Interazione con le matrici */
   //@{
      
   /** 
    Somma la matrice ad un matrix handler usando i metodi generici 
    */ 
   MatrixHandler& AddTo(MatrixHandler& MH) const;
   
   /**
    Somma la matrice ad un FullMatrixHandler 
    */
   MatrixHandler& AddTo(FullMatrixHandler& MH) const;
   
   /** 
    Sottrae la matrice da un matrix handler usando i metodi generici 
    */
   MatrixHandler& SubFrom(MatrixHandler& MH) const;
   
   /** 
    Sottrae la matrice da un FullMatrixHandler 
    */
   MatrixHandler& SubFrom(FullMatrixHandler& MH) const;
   //@}
};
            
/* SparseSubMatrixHandler - end */


/* VariableSubMatrixHandler - begin */

/**
 Matrice che puo' diventare via via una sottomatrice piena o una sottomatrice
 sparsa, condividendo la memoria.
 Viene passata agli elementi che, a seconda della loro convenienza, 
 la configurano nel modo piu' opportuno.
 Quindi, con metodi opportuni, viene sommata alla matrice completa.
 */
   
class VariableSubMatrixHandler 
: public FullSubMatrixHandler, public SparseSubMatrixHandler {
 private:
   /**
    Stato della matrice.
    */
   enum { NULLMATRIX, FULL, SPARSE } eStatus;
   
 public:
   /**@name Costruttori */
   //@{
   
   /** 
    Costruttore: riceve gli spazi di lavoro con le loro dimensioni ed 
    inizializza le matrici piena e sparsa.
    @param iIntSize    dimensioni dell'array degli indici
    @param iDoubleSize dimensioni dell'array dei coefficienti
    @param piInt       array degli indici
    @param pdDouble    array dei coefficienti
    */
   VariableSubMatrixHandler(integer iIntSize, integer iDoubleSize, 
			    integer* piInt, doublereal* pdDouble)
     : FullSubMatrixHandler(iIntSize, iDoubleSize, piInt, pdDouble), 
     SparseSubMatrixHandler(iIntSize, iDoubleSize, piInt, pdDouble),
     eStatus(NULLMATRIX) { 
	NO_OP;
     };
   
   /** 
    Distruttore banale 
    */
#if 0
   virtual ~VariableSubMatrixHandler(void) { 
      NO_OP;
   };
#endif
   //@}
   
   /** Metodi di servizio */
   //@{
   
   /**
    Setta la matrice come vuota.
    Di conseguenza non viene assemblata.
    */
   void SetNullMatrix(void) { 
      eStatus = NULLMATRIX; 
   };
   
   /** 
    Setta la matrice come piena.
    Ritorna un riferimento a matrice piena, che puo' essere usato per le 
    normali operazioni di scrittura delle matrici piene.
    */
   FullSubMatrixHandler& SetFull(void) { 
      eStatus = FULL; 
      return (*(FullSubMatrixHandler*)this);
   };
   
   /** 
    Setta la matrice come sparsa.
    Ritorna un riferimento a matrice sparsa, che puo' essere usato per le 
    normali operazioni di scrittura delle matrici sparse.
    */
   SparseSubMatrixHandler& SetSparse(void) { 
      eStatus = SPARSE; 
      return (*(SparseSubMatrixHandler*)this);
   };
   
   /**
    Verifica se la matrice e' vuota.
    */
   flag fIsNullMatrix(void) const { 
      return flag(eStatus == NULLMATRIX);
   };
   
   /**
    Verifica se la matrice e' piena.
    */
   flag fIsFull(void) const { 
      return flag(eStatus == FULL);
   };
   
   /**
    Verifica se la matrice e' sparsa.
    */
   flag fIsSparse(void) const { 
      return flag(eStatus == SPARSE); 
   };

   /**
    Numero di righe della sottomatrice
   */ 
   integer iGetNumRows(void) const {
     if (fIsFull()) { 
       return FullSubMatrixHandler::iGetNumRows();
     }
     else {
       return SparseSubMatrixHandler::iGetNumRows();
     }
   };

   /**
    Numero di colonne della sottomatrice
   */ 
   integer iGetNumCols(void) const {
     if (fIsFull()) { 
       return FullSubMatrixHandler::iGetNumCols();
     }
     else {
       return SparseSubMatrixHandler::iGetNumCols();
     }
   };
      

   /** 
    Collega una matrice sparsa con della memoria già assegnata 
   */
   void Attach(int iNumEntr, doublereal* pdTmpMat, integer* piTmpIndx) { 
     SetSparse();
     SparseSubMatrixHandler::Attach(iNumEntr, pdTmpMat, piTmpIndx);
   };
 
   /** 
    Collega una matrice Full con della memoria già assegnata 
   */
   void Attach(int iNumRows,int iNumCols, doublereal* pdTmpMat, integer* piTmpIndx) { 
     SetFull();
     FullSubMatrixHandler::Attach(iNumRows, iNumCols, pdTmpMat, piTmpIndx);
   };

   //@}
   
   /**@name Interazione con le matrici */
   //@{
      
   /** 
    Si somma ad una matrice completa con metodi generici.
    */
   MatrixHandler& AddTo(MatrixHandler& MH) const {
      if (fIsSparse()) {
	 return SparseSubMatrixHandler::AddTo(MH);
      }
      if (fIsFull()) {
	 return FullSubMatrixHandler::AddTo(MH);
      }
      return MH;
   };

   /** 
    Si sottrae da una matrice completa con metodi generici.
    */
   MatrixHandler& SubFrom(MatrixHandler& MH) const {
      if (fIsSparse()) {
	 return SparseSubMatrixHandler::SubFrom(MH);
      }
      if (fIsFull()) {
	 return FullSubMatrixHandler::SubFrom(MH);
      }
      return MH;
   };
   //@}
};

/* VariableSubMatrixHandler - end */


/* SubVectorHandler - begin */

/**
 Classe virtuale dei sottovettori.
 */

class SubVectorHandler : public VectorHandler {
 public:
   
   /**@name Costruttori */
   //@{
   
   /**
    Distruttore virtuale.
    */
   virtual ~SubVectorHandler(void) { 
      NO_OP;
   };
   //@}
   
   /**@name Metodi di servizio */
   //@{
      
   /**
    Routine di verifica della validita' dell'oggetto.
    Usata per il debug.
    */
   virtual void IsValid(void) const = 0;
   //@}
   
   /**@name Operazioni su indici e coefficienti */
   //@{
      
   /**
    Scrive un indice di riga 
    */
   virtual inline flag fPutRowIndex(integer iSubRow, integer iRow) = 0;
   
   /**
    Ottiene un indice di riga
    */
   virtual inline integer iGetRowIndex(integer iSubRow) const = 0;
   
   /**
    Scrive una entry completa.
    @param iSubRow numero della entry (indice del sotto-vettore)
    @param iRow    indice della entry
    @param dCoef   coefficiente della entry
    */
   virtual inline flag fPutItem(integer iSubRow, integer iRow,
				const doublereal& dCoef) {
      fPutRowIndex(iSubRow, iRow);
      return fPutCoef(iSubRow, dCoef);
   };
   //@}
   
   /**@name Interazione con i vettori */
   //@{
  
   /** 
    Si somma ad un vettore con metodi generici 
    */
   virtual VectorHandler& AddTo(VectorHandler& VH) const = 0;
   //@}
};

   
/**
 Sottovettore standard, formato da un vettore di reali con associato 
 un vettore di interi che contiene gli indici di ogni coefficiente.
 Per il vettore di reali viene usato un VectorHandler, da cui la classe e' 
 derivata. 
 */
class MySubVectorHandler : public SubVectorHandler, public MyVectorHandler {
   friend ostream& operator << (ostream& out, const SubVectorHandler& v);
   
 protected:
   /** Puntatore all'array degli indici */   
   integer* piRow;
   /** Puntatore alla posizione precedente a piRow.
    Usato per rendere piu' efficiente l'accesso, dato che gli indici sono 
    a base 1, in stile FORTRAN
    */
   integer* piRowm1;
   
 protected:
 public:
   /**@name Costruttori */
   //@{
   
   /**
    Costruttore per memoria posseduta.
    Specifica solo la dimensione dell'array, che deve essere non-nulla.
    La memoria viene allocata e gestita dal VectorHandler.
    */
   MySubVectorHandler(integer iSize);
      
   /**
    Costruttore per memoria in prestito.
    Riceve la dimensione dell'array e i puntatori alla memoria.
    */      
   MySubVectorHandler(integer iSize, integer* piTmpRow, doublereal* pdTmpVec);
   
   /**
    Distruttore.
    */
   virtual ~MySubVectorHandler(void) { 
      Detach();
   };
   //@}


   /** 
    Tutti questi metodi sono richiesti perche' la classe MySubVectorHandler
    dipende due volte da VectorHandler e quindi vi e' un'ambiguita' 
    che va risolta (in realta' solo le funzioni di MyVectorHandler
    sono definite, tutte le altre sono virtuali pure!) 
    */
   /**@name Metodi di servizio */
   //@{
      
   /** 
    Puntatore alla base del vettore (deprecato)
    */
   virtual doublereal* pdGetVec(void) const {
      return MyVectorHandler::pdGetVec();
   };

   /**
    Dimensioni del vettore 
    */
   virtual integer iGetSize(void) const {
      return MyVectorHandler::iGetSize();
   };
      
   /**
    Ridimensiona il vettore.
    Nota: se il vettore possiede la memoria a cui punta,
    la nuova dimensione puo' eccedere la massima dimensione corrente.
    */
   virtual void Resize(integer iSize);
      
   /**
    Inizializza il vettore con d
    */
   virtual void Reset(doublereal d = 0.) {
      MyVectorHandler::Reset(d);
   };
   
   /**
    Scollega il vettore dalla memoria ad esso associata.
    Se il vettore possiede la memoria, viene delallocata.
    */
   void Detach(void);
   
   /**
    Collega il vettore alla memoria che gli viene passata.
    La memoria a cui il vettore era collegato viene deallocata se
    era posseduta dal vettore.
    */
   void Attach(integer iSize, doublereal* pd, integer* pi, integer iMSize = 0);   
   
   /**
    Verifica la validita' del vettore.
    Usata per debug
    */
   virtual void IsValid(void) const;
   //@}


   /**@name Operazioni sugli indici e sui coefficienti */
   //@{
      
   /** 
    Scrive un coefficiente in base al sottoindice.
    */
   virtual flag fPutCoef(integer i, const doublereal& d) {
      return MyVectorHandler::fPutCoef(i, d);
   };

   /** 
    Incrementa un coefficiente in base al sottoindice.
    */
   virtual flag fIncCoef(integer i, const doublereal& d) {
      return MyVectorHandler::fIncCoef(i, d);
   };

   /** 
    Decrementa un coefficiente in base al sottoindice.
    */
   virtual flag fDecCoef(integer i, const doublereal& d) {
      return MyVectorHandler::fDecCoef(i, d);
   };

   /** 
    Ottiene un coefficiente in base al sottoindice.
    */
   virtual const doublereal& dGetCoef(integer i) const {
      return MyVectorHandler::dGetCoef(i);
   };
   
   /** 
    Scrive un indice di riga in base al sottoindice.
    */
   virtual inline flag fPutRowIndex(integer iSubRow, integer iRow);
      
   /** 
    Ottiene un indice di riga in base al sottoindice.
    */
   virtual inline integer iGetRowIndex(integer iSubRow) const;
      
   /** 
    Scrive una entry completa.
    @param iSubRow numero della entry (indice del sotto-vettore)
    @param iRow    indice della entry
    @param dCoef   coefficiente della entry    
    */
   virtual inline flag fPutItem(integer iSubRow, integer iRow, const doublereal& dCoef);
   //@}
   
   /**@name Interazione con i vettori */
   //@{
  
   /** 
    Si somma ad un vettore con metodi generici 
    */
   virtual VectorHandler& AddTo(VectorHandler& VH) const;
      
   /** 
    Si somma ad un MyVectorHandler
    */
   virtual VectorHandler& AddTo(MyVectorHandler& VH) const;
   //@}
};

inline flag MySubVectorHandler::fPutRowIndex(integer iSubRow, integer iRow) {
#ifdef DEBUG
   IsValid();
   ASSERT((iSubRow > 0) && (iSubRow <= iCurSize));
   ASSERT(iRow > 0);
#endif
   
   piRowm1[iSubRow] = iRow;
   
   return flag(0);
}
   
inline integer MySubVectorHandler::iGetRowIndex(integer iSubRow) const {
#ifdef DEBUG
   IsValid();
   ASSERT((iSubRow > 0) && (iSubRow <= iCurSize));
#endif
   
   return piRowm1[iSubRow];
}

inline flag MySubVectorHandler::fPutItem(integer iSubRow, integer iRow, 
					 const doublereal& dCoef) {	
#ifdef DEBUG
   IsValid();
   ASSERT((iSubRow > 0) && (iSubRow <= iCurSize));
   ASSERT(iRow > 0);
#endif	
   
   piRowm1[iSubRow] = iRow;
   pdVecm1[iSubRow] = dCoef;
   
   return flag(0);
}


/**@name Operazioni esterne su SubMatrixHandler e su SubVectorHandler */
//@{
   
/**
 Operatore per scrittura di SubVectorHandler su ostream.
 Usato principalmente per debug
 */
extern ostream& operator << (ostream& out, const SubVectorHandler& v);
   
/**
 Operatore per scrittura di FullSubMatrixHandler su ostream.
 Usato principalmente per debug
 */
extern ostream& operator << (ostream& out, const FullSubMatrixHandler& m);
//@}

/* SubVectorHandler - end */

#endif
