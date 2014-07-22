/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#include <iostream>
#include <limits>
#include "ac/f2c.h"

#include "myassert.h"
#include "except.h"
#include "solman.h"

#include "tpls.h"

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

class Mat3x3_Manip;  /* manipolatori per la matrice 3x3 */

/* Vec3_Manip - begin */

class Vec3_Manip {
 public:
   /*
    Metodo che trasforma la matrice m nel vettore v.
    Viene usato da un costruttore di Vec3 che riceve come
    argomenti un manipolatore e una matrice di rotazione.
    
    NOTE: renamed from Make() to Manipulate() May 2009
    */
   virtual void Manipulate(Vec3& v, const Mat3x3& m) const = 0;

   virtual ~Vec3_Manip(void) { 
      NO_OP;
   };
};

/* Vec3_Manip - end */




/* Vec3 - begin */

// Vettori di dimensione 3
class Vec3 {
   friend class Mat3x3;   
   friend Vec3 operator - (const Vec3& v);
   // friend class Mat3x3_Manip;
   friend class VecN;
   friend class Mat3xN;
   friend class MatNx3;
   

 private:
   //vettore di tre reali che contiene i coefficienti
   doublereal pdVec[3];

 protected:
   
 public:
   /*Costruttori */

   /* 
    Costruttore banalissimo: non fa nulla. Attenzione che cosi' il
    valore di Vec3 e' unpredictable. Se si desidera un vettore nullo
    usare Zero3.
    */   
   Vec3(void) { 
      NO_OP; 
   };
   
   /*
    Assegna i tre valori. 
    Per azzerare il vettore, usare Zero3.
    */
   Vec3(const doublereal& v1, 
	const doublereal& v2, 
	const doublereal& v3) {
      pdVec[V1] = v1; 
      pdVec[V2] = v2; 
      pdVec[V3] = v3; 
   };
   
   /*
    Costruttore di copia.
    */
   Vec3(const Vec3& v) {
      pdVec[V1] = v.pdVec[V1];
      pdVec[V2] = v.pdVec[V2];
      pdVec[V3] = v.pdVec[V3];
   };
   
   /*
    Costruttore da array di reali. Copia i primi 3 valori dell'array.
    Si assume che l'array da pd sia lungo almeno 3. 
    */
   Vec3(const doublereal* pd) {
      ASSERT(pd != NULL);
      
      pdVec[V1] = pd[0];
      pdVec[V2] = pd[1];
      pdVec[V3] = pd[2];
   };
   
   /*
    Costruttore da VectorHandler. Prende i valori da iFirstIndex 
    a iFirstIndex+2. Nota: gli indici del VectorHandler partono da 1, 
    in stile FORTRAN.
    */
   Vec3(const VectorHandler& vh, integer iFirstIndex) {
      ASSERT(iFirstIndex > 0 && iFirstIndex <= vh.iGetSize()-2);
      pdVec[V1] = vh(iFirstIndex);
      pdVec[V2] = vh(++iFirstIndex);
      pdVec[V3] = vh(++iFirstIndex);
   };
   
   /*
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dal vettore v.
    */
   Vec3(const Vec3_Manip& Manip, const Mat3x3& m) {
      Manip.Manipulate(*this, m);
   };
   
   /*
    Distruttore banale: non fa nulla.
    */
   ~Vec3(void) { 
      NO_OP; 
   };
      

   /*Metodi di servizio */
      
   /*
    Dirty job: restituisce il puntatore al vettore (deprecato).
    */
   const doublereal* pGetVec(void) const { 
      return pdVec;
   };
      
   doublereal* pGetVec(void) { 
      return pdVec;
   };
      
   /*Operatori su vettori e matrici */
      
   /*
    Prodotto "tensore". 
    Restituisce se stesso per v trasposto.
    */
   Mat3x3 Tens(const Vec3& v) const;
      
   /*
    Prodotto "tensore". 
    Restituisce se stesso per se stesso
    */
   Mat3x3 Tens(void) const;
 
   /*
    Prodotto vettore. 
    Restituisce il prodotto vettore tra se stesso e v in un temporaneo.
    */
   Vec3 Cross(const Vec3& v) const {     
      return Vec3(pdVec[V2]*v.pdVec[V3]-pdVec[V3]*v.pdVec[V2],
		  pdVec[V3]*v.pdVec[V1]-pdVec[V1]*v.pdVec[V3],
		  pdVec[V1]*v.pdVec[V2]-pdVec[V2]*v.pdVec[V1]);
   };
   
   /*
    * element by element product of vectors
    */
   Vec3 EBEMult(const Vec3& v) const {
		return Vec3(pdVec[V1]*v.pdVec[V1],
		pdVec[V2]*v.pdVec[V2],
		pdVec[V3]*v.pdVec[V3]);
   };

   /*
    Prodotto vettore per matrice. 
    Restituisce il prodotto vettore tra se stesso e m in un temporaneo.
    */
   Mat3x3 Cross(const Mat3x3& m) const;
   
   /*
    Prodotto scalare. 
    Restituisce il prodotto scalare tra se e v
    */
   doublereal Dot(const Vec3& v) const {
      return 
	pdVec[V1]*v.pdVec[V1]+
	pdVec[V2]*v.pdVec[V2]+
	pdVec[V3]*v.pdVec[V3];
   };
   
   /*
    Prodotto scalare per se stesso.
    */
   doublereal Dot(void) const { 
      return 
	pdVec[V1]*pdVec[V1]+
	pdVec[V2]*pdVec[V2]+
	pdVec[V3]*pdVec[V3];
   };
   
   /*
    Norma: sqrt(Dot())
    */
   doublereal Norm(void) const { 
      return 
	sqrt(pdVec[V1]*pdVec[V1]+
	     pdVec[V2]*pdVec[V2]+
	     pdVec[V3]*pdVec[V3]);
   };

   /*Operazioni sui coefficienti */
      
   /*
    Assegnazione di un coefficiente. 
    Nota: l'indice ha base 1, in stile FORTRAN.
    */
   inline void Put(unsigned short int iRow, const doublereal& dCoef) {
      ASSERT(iRow >= 1 && iRow <= 3);   
      pdVec[--iRow] = dCoef;
   };
   
   /*
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
   
   
   /*Operazioni con array di reali */
      
   /*
    Somma se stesso all'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void AddTo(doublereal* pd) const {	
      ASSERT(pd != NULL);
      pd[0] += pdVec[V1];
      pd[1] += pdVec[V2];
      pd[2] += pdVec[V3];
   };
   
   /*
    Sottrae se stesso dall'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void SubFrom(doublereal* pd) const {	
      ASSERT(pd != NULL);
      pd[0] -= pdVec[V1];
      pd[1] -= pdVec[V2];
      pd[2] -= pdVec[V3];
   };
   
   /*
    Scrive se stesso sull'array pd.
    Si assume che l'array pd sia lungo almeno 3 
    */
   void PutTo(doublereal* pd) const {
      ASSERT(pd != NULL);
      pd[0] = pdVec[V1];
      pd[1] = pdVec[V2];
      pd[2] = pdVec[V3];
   };   
   
   /*
    Si legge dall'array pd.
    Si assume che l'array pd sia lungo almeno 3
    */
   void GetFrom(const doublereal* pd) {
      ASSERT(pd != NULL);
      pdVec[V1] = pd[0];
      pdVec[V2] = pd[1];
      pdVec[V3] = pd[2];
   };

   /*Operatori */
      
   /*
    Operatore di assegnazione 
    */
   const Vec3& operator = (const Vec3& v) {
      pdVec[V1] = v.pdVec[V1];
      pdVec[V2] = v.pdVec[V2];
      pdVec[V3] = v.pdVec[V3];
      
      return *this;
   };
      
   /*
    Operatore somma. 
    Restituisce v sommato a se stesso in un temporaneo.
    */
   Vec3 operator + (const Vec3& v) const {
      return Vec3(pdVec[V1]+v.pdVec[V1],
		  pdVec[V2]+v.pdVec[V2],
		  pdVec[V3]+v.pdVec[V3]);
   };
   
   /*
    Operatore somma e assegnazione. 
    Somma v a se stesso in loco.
    */
   const Vec3& operator += (const Vec3& v) {
      pdVec[V1] += v.pdVec[V1];
      pdVec[V2] += v.pdVec[V2];
      pdVec[V3] += v.pdVec[V3];
      return *this;
   };
   
   /*
    Operatore sottrazione. 
    Sottrae v da se stesso in un temporaneo.
    */
   Vec3 operator - (const Vec3& v) const {
      return Vec3(pdVec[V1]-v.pdVec[V1],
		  pdVec[V2]-v.pdVec[V2],
		  pdVec[V3]-v.pdVec[V3]);
   };
   
   /*
    Operatore sottrazione e assegnazione. 
    Sottrae v da se stesso in loco.
    */
   const Vec3& operator -= (const Vec3& v) {
      pdVec[V1] -= v.pdVec[V1];
      pdVec[V2] -= v.pdVec[V2];
      pdVec[V3] -= v.pdVec[V3];
      return *this;
   };
   
   /*
    Operatore prodotto per scalare. 
    Moltiplica se stesso per d in un temporaneo.
    */
   Vec3 operator * (const doublereal& d) const {
      return Vec3(pdVec[V1]*d,
		  pdVec[V2]*d,
		  pdVec[V3]*d);    
   };
   
   /*
    Operatore prodotto e assegnazione per scalare.
    Moltiplica se stesso per d in loco.
    */
   const Vec3& operator *= (const doublereal& d) {
      pdVec[V1] *= d;
      pdVec[V2] *= d;
      pdVec[V3] *= d;
      return *this;
   };   

   /*
    Prodotto scalare. 
    Moltiplica se stesso per v.
    */
   doublereal operator * (const Vec3& v) const {
      return pdVec[V1]*v.pdVec[V1]
	+pdVec[V2]*v.pdVec[V2]
	+pdVec[V3]*v.pdVec[V3];
   };   
   
	/**
	 * Scalar product
	 * multiplies self by matrix m; equivalent to m.Transpose() * this.
	 */
	Vec3 operator * (const Mat3x3& m) const;
 
   /*
    Operatore divisione per scalare. 
    Divide se stesso per d in un temporaneo.
    */
   Vec3 operator / (const doublereal& d) const {
      ASSERT(d != 0.);      
      return Vec3(pdVec[V1]/d,
		  pdVec[V2]/d,
		  pdVec[V3]/d);    
   };
   
   /*
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

   bool IsNull(void) const {
      return (pdVec[V1] == 0. && pdVec[V2] == 0. && pdVec[V3] == 0.);
   };
     
   bool IsExactlySame(const Vec3& v) const {
      return (pdVec[V1] == v.pdVec[V1]
           && pdVec[V2] == v.pdVec[V2]
           && pdVec[V3] == v.pdVec[V3]);
   };
 
   bool IsSame(const Vec3& v, const doublereal& dTol) const {
      doublereal d2 = 0.;

      for (int i = 0; i < 3; i++) {
         doublereal d = pdVec[i] - v.pdVec[i];
         d2 += d*d;
      }

      return sqrt(d2) <= dTol;
   };

   void Reset(void);
 
   /*Input/Output */
      
   /*
    Scrive se stesso sull'ostream out.
    I coefficienti sono separati dalla stringa sFill (spazio di default).
    */
   std::ostream& Write(std::ostream& out, const char* sFill = " ") const;
};
   
/* Vec3 - end */


/* Mat3x3_Manip - begin */

/* classe virtuale dei manipolatori */
/* 
 Manipolatori per la costruzione di matrici 3x3.
 Sono usati per ottenere le matrici di rotazione e lel loro derivate
 dai parametri di rotazione
 */
class Mat3x3_Manip {
 public:
   virtual void Manipulate(Mat3x3& m) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   };

   virtual void Manipulate(Mat3x3& m, const doublereal d) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   };

   /*
    Metodo che trasforma il vettore v nella matrice m.
    Viene usato da un costruttore di Mat3x3 che riceve come
    argomenti un manipolatore e un vettore di parametri di rotazione.
    
    NOTE: renamed from Make() to Manipulate() May 2009
    */
   virtual void Manipulate(Mat3x3& m, const Vec3& v) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   };

   virtual void Manipulate(Mat3x3& m, const Vec3& v1, const Vec3& v2) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   };

   virtual ~Mat3x3_Manip(void) { 
      NO_OP;
   };
};

/* Mat3x3_Manip - end */


/* Mat3x3 - begin */
// Matrici 3x3
class Mat3x3 {
   friend class Vec3;
   friend class SparseSubMatrixHandler;
   friend class Mat3x3_Manip;   
   friend class Mat3xN;
   friend class MatNx3;
   
 protected:
   //Vettore di 9 reali che contiene i coefficienti
   doublereal pdMat[9];
   
 public:

   /*Costruttori */
   /*
    Costruttore banale: non inizializza i coefficienti.
    Per azzerare la matrice usare Zero3x3.
    */
   Mat3x3(void) {
      NO_OP;
   };

   // replaced by Mat3x3Zero_Manip when passed 0.
   // replaced by Mat3x3DEye_Manip when passed a nonzero.
#if 0
   Mat3x3(const doublereal& d);
#endif

   /*
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
   
   /*
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

#if 0
   // will be replaced by MatCross_Manip 
   /*
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
#endif

#if 0 
   // will be replaced by MatCrossCross_Manip 
   /*
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
#endif 
   
   /*
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
   
   
   /*
    Costruttore di copia da array.
    si assume che l'array pd sia lungo almeno 9 
    */
   Mat3x3(const doublereal* pd, integer iSize) {
      ASSERT(pd != NULL);      
      GetFrom(pd, iSize);
   };   
   
   /*
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dal nulla.
    */
   Mat3x3(const Mat3x3_Manip& Manip) {
      Manip.Manipulate(*this);
   };

   /*
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dal reale d.
    */
   Mat3x3(const Mat3x3_Manip& Manip, const doublereal d) {
      Manip.Manipulate(*this, d);
   };

   /*
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dal vettore v.
    */
   Mat3x3(const Mat3x3_Manip& Manip, const Vec3& v) {
      Manip.Manipulate(*this, v);
   };

   /*
    Costruttore con manipolatore. 
    Invoca un metodo del manipolatore Manip che restituisce una matrice 
    a partire dai vettori v1 e v2.
    */
   Mat3x3(const Mat3x3_Manip& Manip, const Vec3& v1, const Vec3& v2) {
      Manip.Manipulate(*this, v1, v2);
   };

   // FIXME: is this "safe"?
   /*
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
   
   /*
    Distruttore banale.
    */
   ~Mat3x3(void) { 
      NO_OP;
   };


   /*Metodi di servizio */
      
   /*
    Dirty job: restituisce il puntatore alla matrice (deprecato).
    */
   const doublereal* pGetMat(void) const { 
      return pdMat;
   };

   doublereal* pGetMat(void) { 
      return pdMat;
   };

      
   /*Operazioni sui coefficienti */
      
   /*
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
   
   /*
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

   
   /*Operazioni su matrici e vettori */
      
   /*
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

   /*
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

   Mat3x3 Skew(void) const;

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

   /*
    Ottiene un sottovettore dalla matrice.
    Nota: l'indice e' a base 1, in stile FORTRAN.
    */
   Vec3 GetVec(unsigned short int i) const {
      ASSERT(i >= 1 && i <= 3);
      return Vec3(pdMat+3*--i);
   };

   /*
    Ottiene un sottovettore dalla matrice.
    Nota: l'indice e' a base 1, in stile FORTRAN.
    Alias di GetVec()
    */
   Vec3 GetCol(unsigned short int i) const {
      ASSERT(i >= 1 && i <= 3);
      return Vec3(pdMat + 3*--i);
   };

   /*
    Ottiene un sottovettore dalla matrice.
    Nota: l'indice e' a base 1, in stile FORTRAN.
    */
   Vec3 GetRow(unsigned short int i) const {
      ASSERT(i >= 1 && i <= 3);
      --i;
      return Vec3(pdMat[i], pdMat[3 + i], pdMat[6 + i]);
   };

   void PutVec(unsigned short int i, const Vec3& v) {
      ASSERT(i >= 1 && i <= 3);

      i--; i = 3*i;
      pdMat[i++] = v.pdVec[V1];
      pdMat[i++] = v.pdVec[V2];
      pdMat[i] = v.pdVec[V3];
   };

   void AddVec(unsigned short int i, const Vec3& v) {
      ASSERT(i >= 1 && i <= 3);

      i--; i = 3*i;
      pdMat[i++] += v.pdVec[V1];
      pdMat[i++] += v.pdVec[V2];
      pdMat[i] += v.pdVec[V3];
   };

   void SubVec(unsigned short int i, const Vec3& v) {
      ASSERT(i >= 1 && i <= 3);

      i--; i = 3*i;
      pdMat[i++] -= v.pdVec[V1];
      pdMat[i++] -= v.pdVec[V2];
      pdMat[i] -= v.pdVec[V3];
   };

   /*
    Inversione. 
    Restituisce l'inversa di se stessa in un temporaneo.
    */
   doublereal dDet(void) const;
   Mat3x3 Inv(const doublereal& ddet) const;
   Mat3x3 Inv(void) const;

   /*
    Soluzione.
    Restituisce l'inversa di se stessa per v in un temporaneo.
    */
   Vec3 Solve(const Vec3& v) const;      
   Vec3 Solve(const doublereal& d, const Vec3& v) const;      
   
   Vec3 LDLSolve(const Vec3& v) const;      

   /*
    Eigenvalues
    */
   bool EigSym(Vec3& EigenValues) const;
   bool EigSym(Vec3& EigenValues, Mat3x3& EigenVectors) const;
 
   /*Operazioni su arrays di reali */
      
   /*
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
   
   /*
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
   
   /*
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
   
   /*Operatori */
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   /*
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
   
   
   /*
    Operatore prodotto matrice vettore.
    Restituisce se stesso moltiplicato per v in un temporaneo.
    */
   Vec3 operator * (const Vec3& v) const {

      return Vec3(pdMat[M11]*v.pdVec[V1]+pdMat[M12]*v.pdVec[V2]+pdMat[M13]*v.pdVec[V3],
		  pdMat[M21]*v.pdVec[V1]+pdMat[M22]*v.pdVec[V2]+pdMat[M23]*v.pdVec[V3],
		  pdMat[M31]*v.pdVec[V1]+pdMat[M32]*v.pdVec[V2]+pdMat[M33]*v.pdVec[V3]);
   };

   bool IsNull(void) const {
      return (pdMat[M11] == 0. && pdMat[M12] == 0. && pdMat[M13] == 0.
           && pdMat[M21] == 0. && pdMat[M22] == 0. && pdMat[M23] == 0.
           && pdMat[M31] == 0. && pdMat[M32] == 0. && pdMat[M33] == 0.);
   };
     
   bool IsExactlySame(const Mat3x3& m) const {
      return (pdMat[M11] == m.pdMat[M11]
           && pdMat[M12] == m.pdMat[M12]
           && pdMat[M13] == m.pdMat[M13]
           && pdMat[M21] == m.pdMat[M21]
           && pdMat[M22] == m.pdMat[M22]
           && pdMat[M23] == m.pdMat[M23]
           && pdMat[M31] == m.pdMat[M31]
           && pdMat[M32] == m.pdMat[M32]
           && pdMat[M33] == m.pdMat[M33]);
   };
 
   bool IsSame(const Mat3x3& m, const doublereal& dTol) const {
      doublereal d2 = 0.;

      for (int i = 0; i < 9; i++) {
         doublereal d = pdMat[i] - m.pdMat[i];
         d2 += d*d;
      }

      return sqrt(d2) <= dTol;
   };

	bool IsSymmetric(void) const {
		if (pdMat[M12] != pdMat[M21]
			|| pdMat[M13] != pdMat[M31]
			|| pdMat[M23] != pdMat[M32])
		{
			return false;
		}

		return true;
	};
 
	bool IsSymmetric(const doublereal& dTol) const {
		ASSERT(dTol > 0.);

		if (fabs(pdMat[M12] - pdMat[M21]) > dTol
			|| fabs(pdMat[M13] - pdMat[M31]) > dTol
			|| fabs(pdMat[M23] - pdMat[M32]) > dTol)
		{
			return false;
		}

		return true;
	};

	bool IsDiag(void) const {
		if (pdMat[M12] != 0.
			|| pdMat[M21] != 0.
			|| pdMat[M13] != 0.
			|| pdMat[M31] != 0.
			|| pdMat[M23] != 0.
			|| pdMat[M32] != 0.)
		{
			return false;
		}

		return true;
	};
 
	bool IsDiag(const doublereal& dTol) const {
		ASSERT(dTol > 0.);

		if (fabs(pdMat[M12]) > dTol
			|| fabs(pdMat[M21]) > dTol
			|| fabs(pdMat[M13]) > dTol
			|| fabs(pdMat[M31]) > dTol
			|| fabs(pdMat[M23]) > dTol
			|| fabs(pdMat[M32]) > dTol)
		{
			return false;
		}

		return true;
	};
 
   /*
    Prodotto matrice per matrice.
    Restituisce il prodotto di se stessa per m in un temporaneo.
    */
   Mat3x3 operator * (const Mat3x3& m) const
   {
       return Mat3x3(pdMat[M11]*m.pdMat[M11]+pdMat[M12]*m.pdMat[M21]+pdMat[M13]*m.pdMat[M31],
		 pdMat[M21]*m.pdMat[M11]+pdMat[M22]*m.pdMat[M21]+pdMat[M23]*m.pdMat[M31],
		 pdMat[M31]*m.pdMat[M11]+pdMat[M32]*m.pdMat[M21]+pdMat[M33]*m.pdMat[M31],
		 pdMat[M11]*m.pdMat[M12]+pdMat[M12]*m.pdMat[M22]+pdMat[M13]*m.pdMat[M32],
		 pdMat[M21]*m.pdMat[M12]+pdMat[M22]*m.pdMat[M22]+pdMat[M23]*m.pdMat[M32],
		 pdMat[M31]*m.pdMat[M12]+pdMat[M32]*m.pdMat[M22]+pdMat[M33]*m.pdMat[M32],
		 pdMat[M11]*m.pdMat[M13]+pdMat[M12]*m.pdMat[M23]+pdMat[M13]*m.pdMat[M33],
		 pdMat[M21]*m.pdMat[M13]+pdMat[M22]*m.pdMat[M23]+pdMat[M23]*m.pdMat[M33],
		 pdMat[M31]*m.pdMat[M13]+pdMat[M32]*m.pdMat[M23]+pdMat[M33]*m.pdMat[M33]);
    };      

	/**
	 * multiply by another matrix, transposed: this * m^T
	 */
	Mat3x3 MulMT(const Mat3x3& m) const;

	/**
	 * multiply self transposed by a vector: this^T * v
	 */
	Vec3 MulTV(const Vec3& v) const;

	/**
	 * multiply self transposed by another matrix: this^T * m
	 */
	Mat3x3 MulTM(const Mat3x3& m) const;

	/**
	 * multiply self transposed by another matrix, transposed: this^T * m^T
	 */
	Mat3x3 MulTMT(const Mat3x3& m) const;

	/**
	 * multiply self times vCross
	 */
	Mat3x3 MulVCross(const Vec3& v) const;

	/**
	 * multiply self transposed times vCross
	 */
	Mat3x3 MulTVCross(const Vec3& v) const;

   doublereal Tr(void) const {
      return pdMat[M11]+pdMat[M22]+pdMat[M33];
   };
      
   void Reset(void);

   /*Input/Output */
   
   /*
    Scrittura su ostream della matrice.
    Scrive se stessa sull'ostream out usando come separatore tra colonne sFill
    e come separatore tra le righe sFill2.
    */
   std::ostream& Write(std::ostream& out, 
		  const char* sFill = " ", 
		  const char* sFill2 = NULL) const;
};

/* Mat3x3 - end */

/*Operazioni esterne su Vec3 e Mat3x3 */

/*
 Operatore "meno" unario su Vec3.
 Restituisce l'opposto di se stesso in un temporaneo.
 */
extern Vec3 operator - (const Vec3& v);

/*
 Operatore "meno" unario su Mat3x3. 
 Restituisce l'opposto di se stesso in un temporaneo.
 */
extern Mat3x3 operator - (const Mat3x3& v);

/*
 Operatore di scrittura di Vec3 su ostream.
 Nota: i coefficienti sono separati da spazi. Non c'e' endl al termine
 */
extern std::ostream& operator << (std::ostream& out, const Vec3& v);
      
/*
 Operatore di scrittura di Mat3x3 su ostream.
 Nota: i coefficienti sono separati da spazi e sono scritti consecutivamente,
 per righe. Non c'e' endl al termine.
 */
extern std::ostream& operator << (std::ostream& out, const Mat3x3& m);


/*
 Funzione di Output di reali su ostream.
 Necessarie per poter generare i templates dei legami costitutivi.
 Il terzo parametro e' definito solo per compatibilita'.
 @param out   ostream su cui avviene la scrittura.
 @param d     valore da scrivere.
 */
extern std::ostream& Write(std::ostream& out, const doublereal& d, const char*);
   
/*
 Funzione di Output di Vec3 su ostream.
 Necessarie per poter generare i templates dei legami costitutivi.
 @param out   ostream su cui avviene la scrittura.
 @param v     vettore da scrivere.
 @param s     separatore tra i valori.
 */
extern std::ostream& Write(std::ostream& out, const Vec3& v, const char* s = " ");
   
/*
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

// replaces Mat3x3(const doublereal&) used as Mat3x3(0.)
class Mat3x3Zero_Manip : public Mat3x3_Manip {
public:
	Mat3x3Zero_Manip() {};
	inline void Manipulate(Mat3x3& m) const {
		doublereal *pdm = m.pGetMat();

		pdm[M11] = 0.;
		pdm[M12] = 0.;
		pdm[M13] = 0.;
		pdm[M21] = 0.;
		pdm[M22] = 0.;
		pdm[M23] = 0.;
		pdm[M31] = 0.;
		pdm[M32] = 0.;
		pdm[M33] = 0.;
	};
};

// replaces Mat3x3(const doublereal&) used as Mat3x3(d) with d != 0.
class Mat3x3DEye_Manip : public Mat3x3_Manip {
public:
	Mat3x3DEye_Manip() {};
	inline void Manipulate(Mat3x3& m, const doublereal d) const {
		doublereal *pdm = m.pGetMat();

		pdm[M11] = d;
		pdm[M12] = 0.;
		pdm[M13] = 0.;
		pdm[M21] = 0.;
		pdm[M22] = d;
		pdm[M23] = 0.;
		pdm[M31] = 0.;
		pdm[M32] = 0.;
		pdm[M33] = d;
	};
};

class Mat3x3Diag_Manip : public Mat3x3_Manip {
public:
	Mat3x3Diag_Manip() {};
	inline void Manipulate(Mat3x3& m, const Vec3& v) const {
		doublereal *pdm = m.pGetMat();
		const doublereal *pdv = v.pGetVec();

		pdm[M11] = pdv[V1];
		pdm[M12] = 0.;
		pdm[M13] = 0.;
		pdm[M21] = 0.;
		pdm[M22] = pdv[V2];
		pdm[M23] = 0.;
		pdm[M31] = 0.;
		pdm[M32] = 0.;
		pdm[M33] = pdv[V3];
	};
};

// will replace Mat3x3(const Vec3&)
class MatCross_Manip : public Mat3x3_Manip {
public:
	MatCross_Manip() {};
	inline void Manipulate(Mat3x3& m, const Vec3& v) const {
		doublereal *pdm = m.pGetMat();
		const doublereal *pdv = v.pGetVec();

		pdm[M11] = 0.;
		pdm[M12] = -pdv[V3];
		pdm[M13] = pdv[V2];
		pdm[M21] = pdv[V3];
		pdm[M22] = 0.;
		pdm[M23] = -pdv[V1];
		pdm[M31] = -pdv[V2];
		pdm[M32] = pdv[V1];
		pdm[M33] = 0.;
	};
};

// will replace Mat3x3(const Vec3&, const Vec3&)
class MatCrossCross_Manip : public Mat3x3_Manip {
public:
	MatCrossCross_Manip() {};
	inline void Manipulate(Mat3x3& m, const Vec3& v1, const Vec3& v2) const {
		doublereal *pdm = m.pGetMat();
		const doublereal *pdv1 = v1.pGetVec();
		const doublereal *pdv2 = v2.pGetVec();

		double d11 = pdv1[V1]*pdv2[V1];
		double d22 = pdv1[V2]*pdv2[V2];
		double d33 = pdv1[V3]*pdv2[V3];

		pdm[M11] = -d22 - d33;
		pdm[M12] = pdv2[V1]*pdv1[V2];
		pdm[M13] = pdv2[V1]*pdv1[V3];
		pdm[M21] = pdv2[V2]*pdv1[V1];
		pdm[M22] = -d33 - d11;
		pdm[M23] = pdv2[V2]*pdv1[V3];
		pdm[M31] = pdv2[V3]*pdv1[V1];
		pdm[M32] = pdv2[V3]*pdv1[V2];
		pdm[M33] = -d11 - d22;
	};
};

extern const Mat3x3Zero_Manip Mat3x3Zero;
extern const Mat3x3DEye_Manip Mat3x3DEye;
extern const Mat3x3Diag_Manip Mat3x3Diag;
extern const MatCross_Manip MatCross;
extern const MatCrossCross_Manip MatCrossCross;

extern const doublereal Zero1;
extern const Vec3 Zero3;
extern const Mat3x3 Zero3x3;
extern const Mat3x3 Eye3;

template <>
inline const doublereal& mb_zero<doublereal>(void)
{
	return ::Zero1;
}

template <>
inline const Vec3& mb_zero<Vec3>(void)
{
	return ::Zero3;
}

template <>
inline const Mat3x3& mb_zero<Mat3x3>(void)
{
	return ::Zero3x3;
}

template <>
inline doublereal mb_deye<doublereal>(const doublereal d)
{
	return d;
}

template <>
inline doublereal& mb_deye<doublereal>(doublereal& out, const doublereal d)
{
	return out = d;
}

template <>
inline Mat3x3 mb_deye<Mat3x3>(const doublereal d)
{
	return Mat3x3(Mat3x3DEye, d);
}

template <>
inline Mat3x3& mb_deye<Mat3x3>(Mat3x3& out, const doublereal d)
{
	Mat3x3DEye.Manipulate(out, d);
	return out;
}

/* known orientation descriptions */
enum OrientationDescription {
	UNKNOWN_ORIENTATION_DESCRIPTION	= -1,
	EULER_123			= 0,
	EULER_313,
	EULER_321,
	ORIENTATION_VECTOR,
	ORIENTATION_MATRIX,

	LAST_ORIENTATION_DESCRIPTION
};

/*
 Calcola i parametri di Rodrigues g a partire dalla matrice di rotazione R.
 Nota: i parametri devono essere definiti, ovvero R non deve rappresentare 
 una rotazione a cui corrispondono parametri singolari.
 */
extern Vec3 MatR2gparam(const Mat3x3& R);

/*
 Computes so-called linear parametrization
 */
extern Vec3 MatR2LinParam(const Mat3x3& R);
   
/*
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


/*
 Calcola gli angoli di Eulero a partire dalla matrice di rotazione R.
 Nota: gli angoli di Eulero sono ritornati in gradi.
 */
extern const doublereal dRaDegr;
extern Vec3 MatR2EulerAngles(const Mat3x3& R);
extern Vec3 MatR2EulerAngles123(const Mat3x3& R);
extern Vec3 MatR2EulerAngles313(const Mat3x3& R);
extern Vec3 MatR2EulerAngles321(const Mat3x3& R);
extern void MatR2EulerParams(const Mat3x3& R, doublereal& e0, Vec3& e);

extern Vec3 Unwrap(const Vec3& vPrev, const Vec3& v);

/*
 Calcola la matrice di rotazione corrispondente agli angoli di Eulero v.
 Nota: gli angoli di Eulero vengono letti in radianti.
 */
extern Mat3x3 EulerAngles2MatR(const Vec3& v);
extern Mat3x3 EulerAngles123_2MatR(const Vec3& v);
extern Mat3x3 EulerAngles313_2MatR(const Vec3& v);
extern Mat3x3 EulerAngles321_2MatR(const Vec3& v);

/*
 * Cayley-Gibbs-Rodrigues rotation manipulation namespace
 */
namespace CGR_Rot {

class Param_Manip : public Vec3_Manip {
public:
	Param_Manip() {};
	inline void Manipulate(Vec3& v, const Mat3x3& m) const {
		// singularity test
		doublereal d = 1. + m.Trace();
   
		if (fabs(d) < std::numeric_limits<doublereal>::epsilon()) {
			silent_cerr("Param_Manip(): divide by zero, "
			"probably due to singularity in rotation parameters" << std::endl);
			throw ErrDivideByZero(MBDYN_EXCEPT_ARGS);
		}
   
		v = m.Ax()*(4./d);
	};
};

/* MatR_Manip - begin */

// Manipolatore per matrice R con parametri di Rodrigues.
class MatR_Manip : public Mat3x3_Manip {   
 public:
   MatR_Manip() {}
   /*
    Crea in m la matrice R corrispondente ai parametri g.
    */
   inline void Manipulate(Mat3x3& m, const Vec3& g) const {
      doublereal d = (4./(4. + g.Dot()));
      
      /*
       m = Eye3;
       m += Mat3x3(g*d);
       */
      
      /* E' piu' efficiente se creo contemporaneamente I+d*g/\ */
      m = Mat3x3(1., g*d);
      
      /* Alla fine sommo il termine d/2*g/\g/\, che e' una matrice piena */
      m += Mat3x3(MatCrossCross, g, g*(d/2.));
   };
};

/* MatR_Manip - end */


/* MatG_Manip - begin */

// Manipolatore per matrice G con parametri di Rodrigues.
class MatG_Manip : public Mat3x3_Manip {
 public:
   MatG_Manip() {};
   /*
    Crea in m la matrice G corrispondente ai parametri g.
    */
   inline void Manipulate(Mat3x3& m, const Vec3& g) const {
      doublereal d = (4./(4.+g.Dot()));
      m = Mat3x3(d, g*(d/2.));
   };
};

/* MatG_Manip - end */


/* MatGm1_Manip - begin */

// Manipolatore per inversa della matrice G con parametri di Rodrigues */
class MatGm1_Manip : public Mat3x3_Manip {
 public:
   MatGm1_Manip() {};
   /*
    Crea in m l'inversa della matrice G corrispondente ai parametri g.
    */
   inline void Manipulate(Mat3x3& m, const Vec3& g) const {
      m = Mat3x3(1., g/(-2.));
      m += g.Tens()/4.;
   };
};

/* MatGm1_Manip - end */


/*
 Manipolatore per parametri di Cayley-Gibbs-Rodrigues
 */
extern const Param_Manip Param;

/*
 Manipolatore per matrice R 
 */
extern const MatR_Manip MatR;

/*
 Manipolatore per matrice G
 */
extern const MatG_Manip MatG;

/*
 Manipolatore per inversa della matrice G 
 */
extern const MatGm1_Manip MatGm1;

} // end of namespace CGR_Rot

/* test */
template <class T>
bool
IsNull(const T& t)
{
	return t.IsNull();
}

template <class T>
bool
IsExactlySame(const T& t1, const T& t2)
{
	return t1.IsExactlySame(t2);
}

template <class T>
bool
IsSame(const T& t1, const T& t2, const doublereal& dTol)
{
	return t1.IsSame(t2, dTol);
}

template <>
bool
IsNull(const doublereal& d);
 
template <>
bool
IsExactlySame(const doublereal& d1, const doublereal& d2);

template <>
bool
IsSame(const doublereal& d1, const doublereal& d2, const doublereal& dTol);

extern Vec3 MultRV(const Vec3& v, const Mat3x3& R);

extern Mat3x3 MultRM(const Mat3x3& m, const Mat3x3& R);
extern Mat3x3 MultMRt(const Mat3x3& m, const Mat3x3& R);
extern Mat3x3 MultRMRt(const Mat3x3& m, const Mat3x3& R);

#endif /* MATVEC3_H */

