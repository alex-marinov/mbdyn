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

/* vettori e matrici 3x3 - operazioni collegate */


#ifndef MATVEC3_H
#define MATVEC3_H

#include <ac/iostream>
#include <ac/f2c.h>

#include <myassert.h>
#include <except.h>
#include <solman.h>

enum {
   V1 = 0,
   V2 = 1,
   V3 = 2
};

enum {
   M11 = 0,
   M12 = 3,
   M13 = 6,
   M21 = 1,
   M22 = 4,
   M23 = 7,
   M31 = 2,
   M32 = 5,
   M33 = 8
};

/* classi principali definite in questo file */
class Vec3;      /* vettore 3x1 */
class Mat3x3;    /* matrice 3x3 */

class _Mat3x3_Manip;  /* manipolatori per la matrice 3x3 */


/* Vec3 - begin */

/// Vettori di dimensione 3
class Vec3 {
   friend class Mat3x3;   
   friend Vec3 operator - (const Vec3& v);
   friend class _Mat3x3_Manip;
   friend class VecN;
   friend class Mat3xN;
   friend class MatNx3;
   

 private:
   ///vettore di tre reali che contiene i coefficienti
   doublereal pdVec[3];

 protected:
   
 public:
   /**@name Costruttori */
   //@{

   /** 
    Costruttore banalissimo: non fa nulla. Attenzione che cosi' il
    valore di Vec3 e' unpredictable. Se si desidera un vettore nullo
    usare Vec3(0.);
    */   
   Vec3(void) { 
      NO_OP; 
   };
   
   /**
    Assegna i tre valori. Se gli ultimi due non sono presenti li setta a zero.
    Per azzerare il vettore, basta usare Vec3(0.)
    */
   Vec3(const doublereal& v1, 
	const doublereal& v2 = 0., 
	const doublereal& v3 = 0.) {
      pdVec[V1] = v1; 
      pdVec[V2] = v2; 
      pdVec[V3] = v3; 
   };
   
   /** 
    Costruttore di copia.
    */
   Vec3(const Vec3& v) {
      pdVec[V1] = v.pdVec[V1];
      pdVec[V2] = v.pdVec[V2];
      pdVec[V3] = v.pdVec[V3];
   };
   
   /** 
    Costruttore da array di reali. Copia i primi 3 valori dell'array.
    Si assume che l'array da pd sia lungo almeno 3. 
    */
   Vec3(const doublereal* pd) {
      ASSERT(pd != NULL);
      
      pdVec[V1] = pd[0];
      pdVec[V2] = pd[1];
      pdVec[V3] = pd[2];
   };
   
   /**
    Costruttore da VectorHandler. Prende i valori da iFirstIndex 
    a iFirstIndex+1. Nota: gli indici del VectorHandler partono da 1, 
    in stile FORTRAN.
    */
   Vec3(const VectorHandler& vh, integer iFirstIndex) {
      ASSERT(iFirstIndex > 0 && iFirstIndex <= vh.iGetSize()-2);
      pdVec[V1] = vh.dGetCoef(iFirstIndex);
      pdVec[V2] = vh.dGetCoef(++iFirstIndex);
      pdVec[V3] = vh.dGetCoef(++iFirstIndex);
   };
   
   /**
    Distruttore banale: non fa nulla.
    */
   ~Vec3(void) { 
      NO_OP; 
   };
      
   //@}

   /**@name Metodi di servizio */
   //@{
      
   /**
    Dirty job: restituisce il puntatore al vettore (deprecato).
    */
   doublereal* pGetVec(void) const { 
      return (doublereal*)pdVec; 
   };
   //@}
      
   /**@name Operatori su vettori e matrici */
   //@{
      
   /** 
    Prodotto "tensore". 
    Restituisce se stesso per v trasposto.
    */
   Mat3x3 Tens(const Vec3& v) const;
      
   /**
    Prodotto vettore. 
    Restituisce il prodotto vettore tra se stesso e v in un temporaneo.
    */
   Vec3 Cross(const Vec3& v) const {     
      return Vec3(pdVec[V2]*v.pdVec[V3]-pdVec[V3]*v.pdVec[V2],
		  pdVec[V3]*v.pdVec[V1]-pdVec[V1]*v.pdVec[V3],
		  pdVec[V1]*v.pdVec[V2]-pdVec[V2]*v.pdVec[V1]);
   };
   
   /** 
    Prodotto vettore per matrice. 
    Restituisce il prodotto vettore tra se stesso e m in un temporaneo.
    */
   Mat3x3 Cross(const Mat3x3& m) const;
   
   /** 
    Prodotto scalare. 
    Restituisce il prodotto scalare tra se e v
    */
   doublereal Dot(const Vec3& v) const {
      return 
	pdVec[V1]*v.pdVec[V1]+
	pdVec[V2]*v.pdVec[V2]+
	pdVec[V3]*v.pdVec[V3];
   };
   
   /** 
    Prodotto scalare per se stesso.
    */
   doublereal Dot(void) const { 
      return 
	pdVec[V1]*pdVec[V1]+
	pdVec[V2]*pdVec[V2]+
	pdVec[V3]*pdVec[V3];
   };
   
   /**
    Norma: sqrt(Dot())
    */
   doublereal Norm(void) const { 
      return 
	sqrt(pdVec[V1]*pdVec[V1]+
	     pdVec[V2]*pdVec[V2]+
	     pdVec[V3]*pdVec[V3]);
   };
   //@}

   /**@name Operazioni sui coefficienti */
   //@{
      
   /** 
    Assegnazione di un coefficiente. 
    Nota: l'indice ha base 1, in stile FORTRAN.
    */
   inline void Put(unsigned short int iRow, const doublereal& dCoef) {
      ASSERT(iRow >= 1 && iRow <= 3);   
      pdVec[--iRow] = dCoef;
   };
   
   /**
    Lettura di un coefficiente.
    Nota: l'indice ha base 1, in stile FORTRAN.
    */
   inline const doublereal& dGet(unsigned short int iRow) const {
      ASSERT(iRow >= 1 && iRow <= 3);
      return pdVec[--iRow];
   };

   inline doublereal& operator () (unsigned short int iRow) {
      ASSERT(iRow >= 1 && iRow <= 3);
      return pdVec[--iRow];
   };

   inline const doublereal& operator () (unsigned short int iRow) const {
      ASSERT(iRow >= 1 && iRow <= 3);
      return pdVec[--iRow];
   };

   inline doublereal& operator [] (unsigned short int iRow) {
      ASSERT(iRow <= 2);
      return pdVec[iRow];
   };

   inline const doublereal& operator [] (unsigned short int iRow) const {
      ASSERT(iRow <= 2);
      return pdVec[iRow];
   };
   //@}
   
   
   /**@name Operazioni con array di reali */
   //@{
      
   /**
    Somma se stesso all'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void AddTo(doublereal* pd) const {	
      ASSERT(pd != NULL);
      pd[0] += pdVec[V1];
      pd[1] += pdVec[V2];
      pd[2] += pdVec[V3];
   };
   
   /** 
    Sottrae se stesso dall'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void SubFrom(doublereal* pd) const {	
      ASSERT(pd != NULL);
      pd[0] -= pdVec[V1];
      pd[1] -= pdVec[V2];
      pd[2] -= pdVec[V3];
   };
   
   /**
    Scrive se stesso sull'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void PutTo(doublereal* pd) const {
      ASSERT(pd != NULL);
      pd[0] = pdVec[V1];
      pd[1] = pdVec[V2];
      pd[2] = pdVec[V3];
   };   
   
   /**
    Si legge dall'array pd.
    Si assume che l'array pd sia lungo almeno 3
    */
   void GetFrom(const doublereal* pd) {
      ASSERT(pd != NULL);
      pdVec[V1] = pd[0];
      pdVec[V2] = pd[1];
      pdVec[V3] = pd[2];
   };
   //@}

   /**@name Operatori */
   //@{
      
   /** 
    Operatore di assegnazione 
    */
   const Vec3& operator = (const Vec3& v) {
      pdVec[V1] = v.pdVec[V1];
      pdVec[V2] = v.pdVec[V2];
      pdVec[V3] = v.pdVec[V3];
      
      return *this;
   };
      
   /**
    Operatore somma. 
    Restituisce v sommato a se stesso in un temporaneo.
    */
   Vec3 operator + (const Vec3& v) const {
      return Vec3(pdVec[V1]+v.pdVec[V1],
		  pdVec[V2]+v.pdVec[V2],
		  pdVec[V3]+v.pdVec[V3]);
   };
   
   /** 
    Operatore somma e assegnazione. 
    Somma v a se stesso in loco.
    */
   const Vec3& operator += (const Vec3& v) {
      pdVec[V1] += v.pdVec[V1];
      pdVec[V2] += v.pdVec[V2];
      pdVec[V3] += v.pdVec[V3];
      return *this;
   };
   
   /** 
    Operatore sottrazione. 
    Sottrae v da se stesso in un temporaneo.
    */
   Vec3 operator - (const Vec3& v) const {
      return Vec3(pdVec[V1]-v.pdVec[V1],
		  pdVec[V2]-v.pdVec[V2],
		  pdVec[V3]-v.pdVec[V3]);
   };
   
   /**
    Operatore sottrazione e assegnazione. 
    Sottrae v da se stesso in loco.
    */
   const Vec3& operator -= (const Vec3& v) {
      pdVec[V1] -= v.pdVec[V1];
      pdVec[V2] -= v.pdVec[V2];
      pdVec[V3] -= v.pdVec[V3];
      return *this;
   };
   
   /** 
    Operatore prodotto per scalare. 
    Moltiplica se stesso per d in un temporaneo.
    */
   Vec3 operator * (const doublereal& d) const {
      return Vec3(pdVec[V1]*d,
		  pdVec[V2]*d,
		  pdVec[V3]*d);    
   };
   
   /**
    Operatore prodotto e assegnazione per scalare.
    Moltiplica se stesso per d in loco.
    */
   const Vec3& operator *= (const doublereal& d) {
      pdVec[V1] *= d;
      pdVec[V2] *= d;
      pdVec[V3] *= d;
      return *this;
   };   

   /**
    Prodotto scalare. 
    Moltiplica se stesso per v.
    */
   doublereal operator * (const Vec3& v) const {
      return pdVec[V1]*v.pdVec[V1]
	+pdVec[V2]*v.pdVec[V2]
	+pdVec[V3]*v.pdVec[V3];
   };   
   
   /**
    Operatore divisione per scalare. 
    Divide se stesso per d in un temporaneo.
    */
   Vec3 operator / (const doublereal& d) const {
      ASSERT(d != 0.);      
      return Vec3(pdVec[V1]/d,
		  pdVec[V2]/d,
		  pdVec[V3]/d);    
   };
   
   /**
    Operatore divisione e assegnazione per scalare. 
    Divide se stesso per d in loco.
    */
   const Vec3& operator /= (const doublereal& d) {
      ASSERT(d != 0.);	
      pdVec[V1] /= d;
      pdVec[V2] /= d;
      pdVec[V3] /= d;
      return *this;
   };   
   
   /**
    Operatore booleano di uguaglianza tra vettori.
    */
      bool operator == (const Vec3& v) const {
	 return (pdVec[V1] == v.pdVec[V1] 
		 && pdVec[V2] == v.pdVec[V2] 
		 && pdVec[V3] == v.pdVec[V3]);
      };

   /**
    Operatore booleano di disuguaglianza tra vettori.
    */
      bool operator != (const Vec3& v) const {
	 return (pdVec[V1] != v.pdVec[V1] 
		 || pdVec[V2] != v.pdVec[V2] 
		 || pdVec[V3] != v.pdVec[V3]);
      };
   //@}
     
   /**@name Input/Output */
   //@{
      
   /**
    Scrive se stesso sull'ostream out.
    I coefficienti sono separati dalla stringa sFill (spazio di default).
    */
   std::ostream& Write(std::ostream& out, const char* sFill = " ") const;
   //@}
};
   
/* Vec3 - end */


/* _Mat3x3_Manip - begin */

/* classe virtuale dei manipolatori */
/* 
 Manipolatori per la costruzione di matrici 3x3.
 Sono usati per ottenere le matrici di rotazione e lel loro derivate
 dai parametri di rotazione
 */
class _Mat3x3_Manip {
 public:
   /**
    Operatore deprecato (e' poco efficiente).
    L'uso e': m = _Mat3x3_Manip << v;
    dove m e' una Mat3x3, v e' un Vec3 e _Mat3x3_Manip e' un oggetto
    derivato da questa classe che implementa la trasformazione da
    vettore a matrice.
    */
   virtual inline Mat3x3 operator << (const Vec3& v) const = 0;
   
   /**
    Metodo che trasforma il vettore v nella matrice m.
    Viene usato da un costruttore di Mat3x3 che riceve come
    argomenti un manipolatore e un vettore di parametri di rotazione.
    */
   virtual inline void Make(Mat3x3& m, const Vec3& v) const = 0;
};

/* _Mat3x3_Manip - end */


/* Mat3x3 - begin */
/// Matrici 3x3
class Mat3x3 {
   friend class Vec3;
   friend class SparseSubMatrixHandler;
   friend class _Mat3x3_Manip;   
   friend class Mat3xN;
   friend class MatNx3;
   
 protected:
   ///Vettore di 9 reali che contiene i coefficienti
   doublereal pdMat[9];
   
 public:

   /**@name Costruttori */
   //@{
   /**
    Costruttore banale: non inizializza i coefficienti.
    Per azzerare la matrice usare Mat3x3(0.)
    */
   Mat3x3(void) {
      NO_OP;
   };
   
   /* Costruttore di matrice diagonale.
    Scrive il valore in ingresso sulla diagonale principale. 
    Di conseguenza, se e' nullo automaticamente azzera la matrice.
    */
   Mat3x3(const doublereal& d) {
      pdMat[M11] = d;
      pdMat[M21] = 0.;
      pdMat[M31] = 0.;
      pdMat[M12] = 0;
      pdMat[M22] = d;
      pdMat[M32] = 0.;
      pdMat[M13] = 0.;
      pdMat[M23] = 0.;
      pdMat[M33] = d;
   };
         
   /** 
    Costrutture completo.
    */
   Mat3x3(const doublereal& m11, const doublereal& m21, const doublereal& m31,
	  const doublereal& m12, const doublereal& m22, const doublereal& m32,
	  const doublereal& m13, const doublereal& m23, const doublereal& m33) {
      pdMat[M11] = m11;
      pdMat[M21] = m21;
      pdMat[M31] = m31;
      pdMat[M12] = m12;
      pdMat[M22] = m22;
      pdMat[M32] = m32;
      pdMat[M13] = m13;
      pdMat[M23] = m23;
      pdMat[M33] = m33;      
   };
   
   /**
    Costruttore di copia.
    */
   Mat3x3(const Mat3x3& m) {
      pdMat[M11] = m.pdMat[M11];
      pdMat[M21] = m.pdMat[M21];
      pdMat[M31] = m.pdMat[M31];
      pdMat[M12] = m.pdMat[M12];
      pdMat[M22] = m.pdMat[M22];
      pdMat[M32] = m.pdMat[M32];
      pdMat[M13] = m.pdMat[M13];
      pdMat[M23] = m.pdMat[M23];
      pdMat[M33] = m.pdMat[M33];
   };
   
   /**
    Costruttore prodotto vettore.
    La matrice viene inizializzata con il vettore v disposto a dare 
    la matrice prodotto vettore 
    (ovvero la matrice che, moltiplicata per w, da' v vettor w).
    */
   Mat3x3(const Vec3& v) {
      pdMat[M11] = 0.;
      pdMat[M21] = v.pdVec[V3];
      pdMat[M31] = -v.pdVec[V2];
      pdMat[M12] = -v.pdVec[V3];
      pdMat[M22] = 0.;      
      pdMat[M32] = v.pdVec[V1];
      pdMat[M13] = v.pdVec[V2];
      pdMat[M23] = -v.pdVec[V1];
      pdMat[M33] = 0.;
   };
   
   /**
    Costruttore doppio prodotto vettore.
    restituisce la matrice data dal prodotto di due matrici prodotto vettore.
    */
   Mat3x3(const Vec3& a, const Vec3& b) {
      pdMat[M11] = -a.pdVec[V2]*b.pdVec[V2]-a.pdVec[V3]*b.pdVec[V3];
      pdMat[M21] = a.pdVec[V1]*b.pdVec[V2];
      pdMat[M31] = a.pdVec[V1]*b.pdVec[V3];
      pdMat[M12] = a.pdVec[V2]*b.pdVec[V1];
      pdMat[M22] = -a.pdVec[V3]*b.pdVec[V3]-a.pdVec[V1]*b.pdVec[V1];
      pdMat[M32] = a.pdVec[V2]*b.pdVec[V3];
      pdMat[M13] = a.pdVec[V3]*b.pdVec[V1];
      pdMat[M23] = a.pdVec[V3]*b.pdVec[V2];
      pdMat[M33] = -a.pdVec[V1]*b.pdVec[V1]-a.pdVec[V2]*b.pdVec[V2];
   };
   
   
   /**
    Costruttore che piglia tre vettori e li affianca a dare la matrice
    */
   Mat3x3(const Vec3& v1, const Vec3& v2, const Vec3& v3) {
      pdMat[M11] = v1.pdVec[V1];
      pdMat[M21] = v1.pdVec[V2];
      pdMat[M31] = v1.pdVec[V3];

      pdMat[M12] = v2.pdVec[V1];
      pdMat[M22] = v2.pdVec[V2];
      pdMat[M32] = v2.pdVec[V3];

      pdMat[M13] = v3.pdVec[V1];
      pdMat[M23] = v3.pdVec[V2];
      pdMat[M33] = v3.pdVec[V3];
   };
   
   
   /**
    Costruttore di copia da array.
    si assume che l'array pd sia lungo almeno 9 
    */
   Mat3x3(const doublereal* pd, integer iSize) {
      ASSERT(pd != NULL);      
      GetFrom(pd, iSize);
   };   
   
   /**
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dal vettore v.
    */
   Mat3x3(const _Mat3x3_Manip& Manip, const Vec3& v) {
      Manip.Make(*this, v);
   };
   
   
   /**
    Costruttore che genera la matrice identita' moltiplicata per d e sommata
    alla matrice prodotto vettore ottenuta da v.
    */
   Mat3x3(const doublereal& d, const Vec3& v) {
      pdMat[M11] = d;
      pdMat[M21] = v.pdVec[V3];
      pdMat[M31] = -v.pdVec[V2];
      pdMat[M12] = -v.pdVec[V3];
      pdMat[M22] = d;
      pdMat[M32] = v.pdVec[V1];
      pdMat[M13] = v.pdVec[V2];
      pdMat[M23] = -v.pdVec[V1];
      pdMat[M33] = d;
   };
   
   /** 
    Distruttore banale.
    */
   ~Mat3x3(void) { 
      NO_OP;
   };
   //@}


   /**@name Metodi di servizio */
   //@{
      
   /**
    Dirty job: restituisce il puntatore alla matrice (deprecato).
    */
   doublereal* pGetMat(void) const { 
      return (doublereal*)pdMat;
   };
   //@}

      
   /**@name Operazioni sui coefficienti */
   //@{
      
   /**
    Assegnazione di un coefficiente.
    Nota: gli indici hanno base 1, in stile FORTRAN.
    */
   inline void Put(unsigned short int iRow,
		   unsigned short int iCol, 
		   const doublereal& dCoef) {
      ASSERT(iRow >= 1 && iRow <= 3);
      ASSERT(iCol >= 1 && iCol <= 3);      
      pdMat[--iRow+3*--iCol] = dCoef;
   };
   
   /**
    Lettura di un coefficiente.
    Nota: gli indici hanno base 1, in stile FORTRAN.
    */
   inline const doublereal& dGet(unsigned short int iRow, 
				 unsigned short int iCol) const {
      ASSERT(iRow >= 1 && iRow <= 3);
      ASSERT(iCol >= 1 && iCol <= 3);      
      return pdMat[--iRow+3*--iCol];
   };

   inline doublereal& operator () (unsigned short int iRow, 
		   unsigned short int iCol) {
       ASSERT(iRow >= 1 && iRow <= 3);
       ASSERT(iCol >= 1 && iCol <= 3);
       return pdMat[--iRow+3*--iCol];
   };

   inline const doublereal& operator () (unsigned short int iRow, 
		   unsigned short int iCol) const {
       ASSERT(iRow >= 1 && iRow <= 3);
       ASSERT(iCol >= 1 && iCol <= 3);
       return pdMat[--iRow+3*--iCol];
   };

   //@}
   
   /**@name Operazioni su matrici e vettori */
   //@{
      
   /**
    Prodotto di matrici prodotto vettore.
    Setta se stessa pari al prodotto delle matrici prodotto vettore 
    ottenute dai vettori a e b.
    */
   const Mat3x3& Tens(const Vec3& a, const Vec3& b) {
      pdMat[M11] = a.pdVec[V1]*b.pdVec[V1];   /* m(1,1) = a(1)*b(1) */
      pdMat[M21] = a.pdVec[V2]*b.pdVec[V1];   /* m(2,1) = a(2)*b(1) */
      pdMat[M31] = a.pdVec[V3]*b.pdVec[V1];   /* m(3,1) = a(3)*b(1) */
      pdMat[M12] = a.pdVec[V1]*b.pdVec[V2];   /* m(1,2) = a(1)*b(2) */
      pdMat[M22] = a.pdVec[V2]*b.pdVec[V2];   /* m(2,2) = a(2)*b(2) */
      pdMat[M32] = a.pdVec[V3]*b.pdVec[V2];   /* m(3,2) = a(3)*b(2) */
      pdMat[M13] = a.pdVec[V1]*b.pdVec[V3];   /* m(1,3) = a(1)*b(3) */
      pdMat[M23] = a.pdVec[V2]*b.pdVec[V3];   /* m(2,3) = a(2)*b(3) */
      pdMat[M33] = a.pdVec[V3]*b.pdVec[V3];   /* m(3,3) = a(3)*b(3) */
      return *this;
   };
   
   /**
    Traspone la matrice.
    Restituisce la trasposta di se stessa in un temporaneo (poco efficiente).
    */
   Mat3x3 Transpose(void) const {
      /* La funzione non e' ottimizzata. Genera una nuova matrice che 
       * costruisce in modo opportuno. Se si ritiene di dover usare 
       * ripetutamente la matrice trasposta, conviene memorizzarla in
       * una variabile di servizio */
      
      return Mat3x3(pdMat[M11], pdMat[M12], pdMat[M13],
		    pdMat[M21], pdMat[M22], pdMat[M23],
		    pdMat[M31], pdMat[M32], pdMat[M33]);
   };

   Vec3 Ax(void) const {
      /* assumendo che sia una matrice di rotazione ?!? */
      return Vec3(
		      .5*(pdMat[M32]-pdMat[M23]),
		      .5*(pdMat[M13]-pdMat[M31]),
		      .5*(pdMat[M21]-pdMat[M12])
		      );
   };

   doublereal Trace(void) const {
      return pdMat[M11]+pdMat[M22]+pdMat[M33];
   };

   Mat3x3 Symm(void) const {
      doublereal m12 = .5*(pdMat[M21]+pdMat[M12]);
      doublereal m13 = .5*(pdMat[M31]+pdMat[M13]);
      doublereal m23 = .5*(pdMat[M32]+pdMat[M23]);

      return Mat3x3(
		      pdMat[M11], m12, m13,
		      m12, pdMat[M22], m23,
		      m13, m23, pdMat[M33]
		      );
   };

   Mat3x3 Symm2(void) const {
      doublereal m12 = pdMat[M21]+pdMat[M12];
      doublereal m13 = pdMat[M31]+pdMat[M13];
      doublereal m23 = pdMat[M32]+pdMat[M23];

      return Mat3x3(
		      2.*pdMat[M11], m12, m13,
		      m12, 2.*pdMat[M22], m23,
		      m13, m23, 2.*pdMat[M33]
		      );
   };

   const Mat3x3& Symm(const Mat3x3& m) {
      pdMat[M11] = m.pdMat[M11];
      pdMat[M22] = m.pdMat[M22];
      pdMat[M33] = m.pdMat[M33];
      pdMat[M12] = pdMat[M21] = .5*(m.pdMat[M21]+m.pdMat[M12]);
      pdMat[M13] = pdMat[M31] = .5*(m.pdMat[M31]+m.pdMat[M13]);
      pdMat[M23] = pdMat[M32] = .5*(m.pdMat[M32]+m.pdMat[M23]);

      return *this;
   };

   Mat3x3 Skew(void) const {
      return Mat3x3(this->Ax());
   };

   const Mat3x3& Skew(const Mat3x3& m) {
      pdMat[M11] = pdMat[M22] = pdMat[M13] = 0.;
      pdMat[M12] = m.pdMat[M12] - m.pdMat[M21];
      pdMat[M21] = -pdMat[M12];
      pdMat[M13] = m.pdMat[M13] - m.pdMat[M31];
      pdMat[M31] = -pdMat[M13];
      pdMat[M23] = m.pdMat[M23] - m.pdMat[M32];
      pdMat[M32] = -pdMat[M23];

      return *this;
   };

   /** 
    Ottiene un sottovettore dalla matrice.
    Nota: l'indice e' a base 1, in stile FORTRAN.
    */
   Vec3 GetVec(unsigned short int i) const {
      ASSERT(i >= 1 && i <= 3);
      return Vec3(pdMat+3*--i);
   };   

   /** 
    Inversione. 
    Restituisce l'inversa di se stessa in un temporaneo.
    */
   Mat3x3 Inv(void) const;
   
   /**
    Soluzione.
    Restituisce l'inversa di se stessa per v in un temporaneo.
    */
   Vec3 Inv(const Vec3& v) const;      
   //@}
   
   /**@name Operazioni su arrays di reali */
   //@{
      
   /**
    Si legge da un array.
    Si assume che l'array pd sia lungo almeno 9.
    @param iSize e' il numero di righe di dimensionamento dell'array, 
    in stile FORTRAN
    */
   void GetFrom(const doublereal* pd, integer iSize) {
      ASSERT(pd != NULL);
      ASSERT(iSize >= 3);

      doublereal* pdFrom = (doublereal*)pd;
      pdMat[M11] = pdFrom[0];
      pdMat[M21] = pdFrom[1];
      pdMat[M31] = pdFrom[2];
      
      pdFrom += iSize;
      pdMat[M12] = pdFrom[0];
      pdMat[M22] = pdFrom[1];
      pdMat[M32] = pdFrom[2];
      
      pdFrom += iSize;
      pdMat[M13] = pdFrom[0];
      pdMat[M23] = pdFrom[1];
      pdMat[M33] = pdFrom[2];
   };
   
   /**
    Somma se stesso ad un array con iNRows righe.
    si assume che l'array sia lungo almeno 2*iNRows+3 a partire da pd.
    @param iNRows e' il numero di righe di dimensionamento dell'array, 
    in stile FORTRAN
    */
   void AddTo(doublereal* pd, integer iNRows) const {
      ASSERT(pd != NULL);
      ASSERT(iNRows >= 3);
      
      doublereal* pdTo = pd;
      
      pdTo[0] += pdMat[M11];
      pdTo[1] += pdMat[M21];
      pdTo[2] += pdMat[M31];
      
      pdTo += iNRows;
      pdTo[0] += pdMat[M12];
      pdTo[1] += pdMat[M22];
      pdTo[2] += pdMat[M32];
      
      pdTo += iNRows;
      pdTo[0] += pdMat[M13];
      pdTo[1] += pdMat[M23];
      pdTo[2] += pdMat[M33];      
   };
   
   /**
    Copia se stesso su un array con iNRows righe.
    Si assume che l'array sia lungo almeno 2*iNRows+3.
    @param iNRows e' il numero di righe di dimensionamento dell'array, 
    in stile FORTRAN
    */
   void PutTo(doublereal* pd, integer iNRows) const {
      ASSERT(pd != NULL);
      ASSERT(iNRows >= 3);
      
      doublereal* pdTo = pd;
      
      pdTo[0] = pdMat[M11];
      pdTo[1] = pdMat[M21];
      pdTo[2] = pdMat[M31];
      
      pdTo += iNRows;
      pdTo[0] = pdMat[M12];
      pdTo[1] = pdMat[M22];
      pdTo[2] = pdMat[M32];
      
      pdTo += iNRows;
      pdTo[0] = pdMat[M13];
      pdTo[1] = pdMat[M23];
      pdTo[2] = pdMat[M33];
   };
   //@}
   
   /**@name Operatori */
   //@{
   
   /**
    Operatore di assegnazione.
    */
   const Mat3x3& operator = (const Mat3x3& m) {
      
      pdMat[M11] = m.pdMat[M11];
      pdMat[M21] = m.pdMat[M21];
      pdMat[M31] = m.pdMat[M31];
      pdMat[M12] = m.pdMat[M12];
      pdMat[M22] = m.pdMat[M22];
      pdMat[M32] = m.pdMat[M32];
      pdMat[M13] = m.pdMat[M13];
      pdMat[M23] = m.pdMat[M23];
      pdMat[M33] = m.pdMat[M33];
      
      return *this;
   };
   
   /**
    Operatore somma. 
    Restituisce v sommato a se stesso in un temporaneo.
    */
   Mat3x3 operator + (const Mat3x3& m) const {
      return Mat3x3(pdMat[M11]+m.pdMat[M11],
		    pdMat[M21]+m.pdMat[M21],
		    pdMat[M31]+m.pdMat[M31],
		    pdMat[M12]+m.pdMat[M12],
		    pdMat[M22]+m.pdMat[M22],
		    pdMat[M32]+m.pdMat[M32],
		    pdMat[M13]+m.pdMat[M13],
		    pdMat[M23]+m.pdMat[M23],
		    pdMat[M33]+m.pdMat[M33]);
   };
   
   /**
    Operatore somma e assegnazione.
    Somma v a se stesso in loco.
    */
   const Mat3x3& operator += (const Mat3x3& m) {
      pdMat[M11] += m.pdMat[M11];
      pdMat[M21] += m.pdMat[M21];
      pdMat[M31] += m.pdMat[M31];
      pdMat[M12] += m.pdMat[M12];
      pdMat[M22] += m.pdMat[M22];
      pdMat[M32] += m.pdMat[M32];
      pdMat[M13] += m.pdMat[M13];
      pdMat[M23] += m.pdMat[M23];
      pdMat[M33] += m.pdMat[M33];
      
      return *this;
   };
   
   /**
    Operatore differenza.
    Restituisce v sottratto da se stesso in un temporaneo.
    */
   Mat3x3 operator - (const Mat3x3& m) const {
      return Mat3x3(pdMat[M11]-m.pdMat[M11],
		    pdMat[M21]-m.pdMat[M21],
		    pdMat[M31]-m.pdMat[M31],
		    pdMat[M12]-m.pdMat[M12],
		    pdMat[M22]-m.pdMat[M22],
		    pdMat[M32]-m.pdMat[M32],
		    pdMat[M13]-m.pdMat[M13],
		    pdMat[M23]-m.pdMat[M23],
		    pdMat[M33]-m.pdMat[M33]);
   };
   
   /**
    Operatore differenza e assegnazione.
    Sottrae v da se stesso in loco.
    */
   const Mat3x3& operator -= (const Mat3x3& m) {
      pdMat[M11] -= m.pdMat[M11];
      pdMat[M21] -= m.pdMat[M21];
      pdMat[M31] -= m.pdMat[M31];
      pdMat[M12] -= m.pdMat[M12];
      pdMat[M22] -= m.pdMat[M22];
      pdMat[M32] -= m.pdMat[M32];
      pdMat[M13] -= m.pdMat[M13];
      pdMat[M23] -= m.pdMat[M23];
      pdMat[M33] -= m.pdMat[M33];

      return *this;
   };
   
   /**
    Operatore moltiplicazione per scalare.
    Restituisce se stesso moltiplicato per d in un temporaneo.
    */
   Mat3x3 operator * (const doublereal& d) const {
      if (d != 1.) {
	 return Mat3x3(pdMat[M11]*d,
		       pdMat[M21]*d,
		       pdMat[M31]*d,
		       pdMat[M12]*d,
		       pdMat[M22]*d,
		       pdMat[M32]*d,
		       pdMat[M13]*d,
		       pdMat[M23]*d,
		       pdMat[M33]*d);
      }
      return *this;
   };
   
   /**
    Operatore moltiplicazione per scalare e assegnazione.
    Moltiplica se stesso per d in loco.
    */
   const Mat3x3& operator *= (const doublereal& d) {
      if (d != 1.) {
	 pdMat[M11] *= d;
	 pdMat[M21] *= d;
	 pdMat[M31] *= d;
	 pdMat[M12] *= d;
	 pdMat[M22] *= d;
	 pdMat[M32] *= d;
	 pdMat[M13] *= d;
	 pdMat[M23] *= d;
	 pdMat[M33] *= d;
      }
      
      return *this;
   };
   
   /**
    Operatore divisione per scalare.
    Restituisce se stesso diviso per d in un temporaneo.
    */
   Mat3x3 operator / (const doublereal& d) const {
      ASSERT(d != 0.);
      if (d != 1.) {
	 return Mat3x3(pdMat[M11]/d,
		       pdMat[M21]/d,
		       pdMat[M31]/d,
		       pdMat[M12]/d,
		       pdMat[M22]/d,
		       pdMat[M32]/d,
		       pdMat[M13]/d,
		       pdMat[M23]/d,
		       pdMat[M33]/d);
      }
      
      return *this;
   };
   
   /**
    Operatore divisione per scalare e assegnazione.
    Divide se stesso per d in loco
    */
   Mat3x3& operator /= (const doublereal& d) {
      ASSERT(d != 0.);
      if (d != 1.) {
	 pdMat[M11] /= d;
	 pdMat[M21] /= d;
	 pdMat[M31] /= d;
	 pdMat[M12] /= d;
	 pdMat[M22] /= d;
	 pdMat[M32] /= d;
	 pdMat[M13] /= d;
	 pdMat[M23] /= d;
	 pdMat[M33] /= d;
      }
      
      return *this;
   };
   
   
   /**
    Operatore prodotto matrice vettore.
    Restituisce se stesso moltiplicato per v in un temporaneo.
    */
   Vec3 operator * (const Vec3& v) const {

      return Vec3(pdMat[M11]*v.pdVec[V1]+pdMat[M12]*v.pdVec[V2]+pdMat[M13]*v.pdVec[V3],
		  pdMat[M21]*v.pdVec[V1]+pdMat[M22]*v.pdVec[V2]+pdMat[M23]*v.pdVec[V3],
		  pdMat[M31]*v.pdVec[V1]+pdMat[M32]*v.pdVec[V2]+pdMat[M33]*v.pdVec[V3]);
   };

   /**
    */
    bool operator == (const Mat3x3& m) const {
      return pdMat[M11] == m.pdMat[M11]
	      && pdMat[M12] == m.pdMat[M12]
	      && pdMat[M13] == m.pdMat[M13]
	      && pdMat[M21] == m.pdMat[M21]
	      && pdMat[M22] == m.pdMat[M22]
	      && pdMat[M23] == m.pdMat[M23]
	      && pdMat[M31] == m.pdMat[M31]
	      && pdMat[M32] == m.pdMat[M32]
	      && pdMat[M33] == m.pdMat[M33];
   };
   
   /**
    */
    bool operator != (const Mat3x3& m) const {
      return pdMat[M11] != m.pdMat[M11]
	      || pdMat[M12] != m.pdMat[M12]
	      || pdMat[M13] != m.pdMat[M13]
	      || pdMat[M21] != m.pdMat[M21]
	      || pdMat[M22] != m.pdMat[M22]
	      || pdMat[M23] != m.pdMat[M23]
	      || pdMat[M31] != m.pdMat[M31]
	      || pdMat[M32] != m.pdMat[M32]
	      || pdMat[M33] != m.pdMat[M33];
   };
   
   /** 
    Prodotto matrice per matrice.
    Restituisce il prodotto di se stessa per m in un temporaneo.
    */
   Mat3x3 operator * (const Mat3x3& m) const;      
   //@}

   doublereal Tr(void) const {
      return pdMat[M11]+pdMat[M22]+pdMat[M33];
   };
      
   /**@name Input/Output */
   //@{
   
   /**
    Scrittura su ostream della matrice.
    Scrive se stessa sull'ostream out usando come separatore tra colonne sFill
    e come separatore tra le righe sFill2.
    */
   std::ostream& Write(std::ostream& out, 
		  const char* sFill = " ", 
		  const char* sFill2 = NULL) const;
   //@}
};

extern const Mat3x3 Eye3;
extern const Mat3x3 Zero3x3;
extern const Vec3 Zero3;

/* Mat3x3 - end */

/**@name Operazioni esterne su Vec3 e Mat3x3 */
//@{

/**
 Operatore "meno" unario su Vec3.
 Restituisce l'opposto di se stesso in un temporaneo.
 */
extern Vec3 operator - (const Vec3& v);

/**
 Operatore "meno" unario su Mat3x3. 
 Restituisce l'opposto di se stesso in un temporaneo.
 */
extern Mat3x3 operator - (const Mat3x3& v);

/**
 Operatore di scrittura di Vec3 su ostream.
 Nota: i coefficienti sono separati da spazi. Non c'e' endl al termine
 */
extern std::ostream& operator << (std::ostream& out, const Vec3& v);
      
/**
 Operatore di scrittura di Mat3x3 su ostream.
 Nota: i coefficienti sono separati da spazi e sono scritti consecutivamente,
 per righe. Non c'e' endl al termine.
 */
extern std::ostream& operator << (std::ostream& out, const Mat3x3& m);


/** 
 Funzione di Output di reali su ostream.
 Necessarie per poter generare i templates dei legami costitutivi.
 Il terzo parametro e' definito solo per compatibilita'.
 @param out   ostream su cui avviene la scrittura.
 @param d     valore da scrivere.
 */
extern std::ostream& Write(std::ostream& out, const doublereal& d, const char*);
   
/**
 Funzione di Output di Vec3 su ostream.
 Necessarie per poter generare i templates dei legami costitutivi.
 @param out   ostream su cui avviene la scrittura.
 @param v     vettore da scrivere.
 @param s     separatore tra i valori.
 */
extern std::ostream& Write(std::ostream& out, const Vec3& v, const char* s = " ");
   
/** 
 Funzione di Output di Mat3x3 su ostream.
 Necessarie per poter generare i templates dei legami costitutivi.
 @param out   ostream su cui avviene la scrittura.
 @param d     valore da scrivere.
 @param s     separatore tra i valori (colonne).
 @param s2    separatore tra i valori (righe); se nullo, usa s.
 */
extern std::ostream& Write(std::ostream& out,
		      const Mat3x3& m,
		      const char* s = " ", 
		      const char* s2 = NULL);


/**
 Calcola i parametri di Rodriguez g a partire dalla matrice di rotazione R.
 Nota: i parametri devono essere definiti, ovvero R non deve rappresentare 
 una rotazione a cui corrispondono parametri singolari.
 */
extern Vec3 gparam(const Mat3x3& R);
   
/** 
 Calcola la matrice di rotazione a partire da due vettori sghembi. 
 @param ia indice del vettore va nel sistema locale.
 @param va vettore ia nel sistema locale.
 @param ib indice del vettore che risulta dall'elaborazione di vb.
 @param vb vettore che viene trasformato come segue:
 il vettore va viene assunto come asse ia nel sistema locale,
 mentre la componente di vb normale a va e giacente nel piano 
 normale al prodotto vettore tra va e vb viene assunta come
 componente ib nel sistema locale.
 */
extern Mat3x3 MatR2vec(unsigned short int ia, 
		       const Vec3& va, 
		       unsigned short int ib, 
		       const Vec3& vb);


/**
 Calcola gli angoli di Eulero a partire dalla matrice di rotazione R.
 Nota: gli angoli di Eulero sono ritornati in gradi.
 */
extern Vec3 EulerAngles(const Mat3x3& R);
extern void EulerParams(const Mat3x3& R, doublereal& e0, Vec3& e);

/**
 Calcola la matrice di rotazione corrispondente agli angoli di Eulero v.
 Nota: gli angoli di Eulero vengono letti in radianti.
 */
extern Mat3x3 RFromEulerAngles(const Vec3& v);
//@}


/* _MatR_Manip - begin */

/// Manipolatore per matrice R con parametri di Rodriguez.
class _MatR_Manip : public _Mat3x3_Manip {   
 public:
   /**
    Operatore deprecato (e' poco efficiente).
    */
   inline Mat3x3 operator << (const Vec3& g) const {
      doublereal d = (4./(4.+g.Dot()));
      return Eye3+Mat3x3(g*d)+Mat3x3(g, g*(d/2.));
   };
   
   /**
    Crea in m la matrice R corrispondente ai parametri g.
    */
   inline void Make(Mat3x3& m, const Vec3& g) const {
      doublereal d = (4./(4.+g.Dot()));
      
      /*
       m = Eye3;
       m += Mat3x3(g*d);
       */
      
      /* E' piu' efficiente se creo contemporaneamente I+d*g/\ */
      m = Mat3x3(1., g*d);
      
      /* Alla fine sommo il termine d/2*g/\g/\, che e' una matrice piena */
      m += Mat3x3(g, g*(d/2.));
   };
};

/* _MatR_Manip - end */


/* _MatG_Manip - begin */

/// Manipolatore per matrice G con parametri di Rodriguez.
class _MatG_Manip : public _Mat3x3_Manip {
 public:
   /**
    Operatore deprecato (e' poco efficiente).
    */
   inline Mat3x3 operator << (const Vec3& g) const {
      doublereal d = (4./(4.+g.Dot()));
      return Eye3*d+Mat3x3(g*(d/2.));
   };

   /**
    Crea in m la matrice G corrispondente ai parametri g.
    */
   inline void Make(Mat3x3& m, const Vec3& g) const {
      doublereal d = (4./(4.+g.Dot()));
      m = Mat3x3(d, g*(d/2.));
   };
};

/* _MatG_Manip - end */


/* _MatGm1_Manip - begin */

/// Manipolatore per inversa della matrice G con parametri di Rodriguez */
class _MatGm1_Manip : public _Mat3x3_Manip {
 public:
   /**
    Operatore deprecato (e' poco efficiente).
    */
   inline Mat3x3 operator << (const Vec3& g) const {
      return Eye3+g.Tens(g/4.)-Mat3x3(g/2.);
   };
   
   /**
    Crea in m l'inversa della matrice G corrispondente ai parametri g.
    */
   inline void Make(Mat3x3& m, const Vec3& g) const {
      m = Mat3x3(1., g/(-2.));
      m += g.Tens(g/4.);
   };
};

/* _MatGm1_Manip - end */


/**@name Costanti */
//@{
   
/**
 Matrice identita' 3x3 
 */
extern const Mat3x3 Eye3;
   
/**
 Matrice nulla 3x3
 */
extern const Mat3x3 Zero3x3;
   
/**
 Vettore nullo 3
 */
extern const Vec3 Zero3;


/** 
 Manipolatore per matrice R 
 */
extern _MatR_Manip MatR;

/** 
 Manipolatore per matrice G
 */
extern _MatG_Manip MatG;

/** 
 Manipolatore per inversa della matrice G 
 */
extern _MatGm1_Manip MatGm1;
//@}

#endif /* MATVEC3_H */

