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

#ifndef TPLDRIVE_H
#define TPLDRIVE_H

#include "drive.h"
#include "mbpar.h"


/* TplDriveCaller - begin */

template <class T>
class TplDriveCaller {
 protected:
   T TReturnValue;
      
 public:
   virtual ~TplDriveCaller(void) { 
      NO_OP;
   };
   
   /* copia */
   virtual TplDriveCaller<T>* pCopy(void) const = 0;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual ostream& Restart(ostream& out) const = 0;
   
   /* Restituisce il valore del driver */
   virtual inline const T& Get(void) const = 0;
};

/* TplDriveCaller - end */


/* TplDriveOwner - begin */

template <class T>
class TplDriveOwner {
 protected:
   TplDriveCaller<T>* pTplDriveCaller;
   
 public:
   TplDriveOwner(const TplDriveCaller<T>* pDC = NULL)
     : pTplDriveCaller((TplDriveCaller<T>*)pDC) { 
	NO_OP;
     };
   
   virtual ~TplDriveOwner(void) { 
      NO_OP;
   };
   
   void Set(const TplDriveCaller<T>* pDC) {
      ASSERT(pDC != NULL);
      ASSERT(pTplDriveCaller == NULL);
      pTplDriveCaller = (TplDriveCaller<T>*)pDC;
   };
   
   TplDriveCaller<T>* pGetDriveCaller(void) const { 
      return pTplDriveCaller;
   };
   
   const T& Get(void) const { 
      return pTplDriveCaller->Get(); 
   };
};

/* TplDriveOwner - end */


/* SingleTplDriveCaller - begin */

template <class T>
class SingleTplDriveCaller : public TplDriveCaller<T>, public DriveOwner {
 protected:
   T t;
   
 public:
   SingleTplDriveCaller(const DriveCaller* pDC, const T& x)
     : DriveOwner(pDC), t((T&)x) {
      NO_OP;
   };
   
   ~SingleTplDriveCaller(void) {
      NO_OP;
   };
   
   /* copia */
   virtual TplDriveCaller<T>* pCopy(void) const {
      typedef SingleTplDriveCaller<T> dc;
      TplDriveCaller<T>* pDC = NULL;
      
      SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy(), t), DMmm);
      
      return pDC;
   };

   /* Scrive il contributo del DriveCaller al file di restart */
   virtual ostream& Restart(ostream& out) const {
      out << "single, ",
        Write(out, t, ", ") << ", ";
      return pGetDriveCaller()->Restart(out);
   }; 
   
   inline const T& Get(void) const {
      return ((T&)TplDriveCaller<T>::TReturnValue = t*dGet());
   };   
};

/* Nota: in caso di reale, viene semplificata la classe in modo da
 *       usare solo il drive senza pesatura che viene assunta unitaria */

class SingleTplDriveCaller<doublereal>
: public TplDriveCaller<doublereal>, public DriveOwner {
 public:
   SingleTplDriveCaller(const DriveCaller* pDC, const doublereal& = 0.)
     : DriveOwner(pDC) {
      NO_OP;
   };
   
   ~SingleTplDriveCaller(void) {
      NO_OP;
   };
   
   /* copia */
   virtual TplDriveCaller<doublereal>* pCopy(void) const {
      TplDriveCaller<doublereal>* pDC = NULL;
      
      typedef SingleTplDriveCaller<doublereal> dc;
      SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy()), DMmm);
      
      return pDC;
   };

   /* Scrive il contributo del DriveCaller al file di restart */
   virtual ostream& Restart(ostream& out) const {
      out << "single, ";
      return pGetDriveCaller()->Restart(out);
   };
   
   inline const doublereal& Get(void) const {
      return dGet();
   };   
};

/* SingleTplDriveCaller - end */


/* ArrayTplDriveCaller - begin */

template <class T>
struct DrivesArray {
   DriveCaller* pDriveCaller;
   T t;
};


template <class T>
class ArrayTplDriveCaller : public TplDriveCaller<T> {
 protected:   
   DrivesArray<T>* pDrivesArray;
   unsigned short int iNumDrives;
   
 public:
   ArrayTplDriveCaller(unsigned short int i, 
		       DrivesArray<T>* pDA)
     : pDrivesArray(pDA), iNumDrives(i) {
	ASSERT(i > 0);
	ASSERT(pDA != NULL);
     };
   
   ~ArrayTplDriveCaller(void) {
      for (int i = 0; i < iNumDrives; i++) {
	 SAFEDELETE(pDrivesArray[i].pDriveCaller, DMmm);
      }
      SAFEDELETEARR(pDrivesArray, DMmm);
   };
   
   /* copia */
   virtual TplDriveCaller<T>* pCopy(void) const {      
      typedef DrivesArray<T> da;
      da* pDA = NULL;
      
      SAFENEWARR(pDA, da, iNumDrives, DMmm);
      
      for (int i = 0; i < iNumDrives; i++) {
	 pDA[i].pDriveCaller = pDrivesArray[i].pDriveCaller->pCopy();
	 pDA[i].t = pDrivesArray[i].t;
      }

      typedef ArrayTplDriveCaller<T> dc;
      TplDriveCaller<T>* pDC = NULL;
      
      SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(iNumDrives, pDA), DMmm);
      
      return pDC;
   };

   /* Scrive il contributo del DriveCaller al file di restart */
   virtual ostream& Restart(ostream& out) const {
      out << "array, " << iNumDrives;
      for (int i = 0; i < iNumDrives; i++) {
         out << ", ",
           Write(out, pDrivesArray[i].t, ", ") << ", ",
           pDrivesArray[i].pDriveCaller->Restart(out);
      }
      return out;
   };
   
   inline const T& Get(void) const {
      (T&)TplDriveCaller<T>::TReturnValue = 0.;
      for (int i = 0; i < iNumDrives; i++) {
	 (T&)TplDriveCaller<T>::TReturnValue += (pDrivesArray[i].t)*(pDrivesArray[i].pDriveCaller->dGet());
      }      
      return TplDriveCaller<T>::TReturnValue;
   };   
};

/* Nota: di questa classe non viene scritta esplicitamente la versione
 *       per reali in quanto il coefficiente moltiplicativo
 *       puo' essere usato per un'ulteriore pesatura del drive */

/* ArrayTplDriveCaller - end */


extern doublereal GetT(MBDynParser& HP, const doublereal& t);

template <class T> T GetT(MBDynParser& HP, const T& t)
{
	return HP.Get(t);
}
  

/* enum delle parole chiave */
template <class T>
class ReadTplDriveKeyWords {
 public:
   enum KeyWords {
      UNKNOWNTPLDRIVE = -1,
	SINGLE = 0,
	ARRAY,
	LASTKEYWORD
   };
};


template <class T> 
TplDriveCaller<T>* ReadTplDrive(const DataManager* pDM,
				MBDynParser& HP,
				const DriveHandler* pDH, 
				const T& t)
{
   DEBUGCOUT("Entering ReadTplDrive" << endl);
   
   const char* sKeyWords[] = {
      "single",
	"array"
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)ReadTplDriveKeyWords<T>::LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);   
   
   TplDriveCaller<T>* pTplDC = NULL;

   /* Nota: siccome capita che non ci sia la deformazione imposta, 
    * faccio chiedere se e' desiderata */
   if (!HP.fIsArg()) {
      DriveCaller* pDC = NULL;
      SAFENEWWITHCONSTRUCTOR(pDC,
			     NullDriveCaller,
			     NullDriveCaller(pDH),
			     DMmm);

      T t(0.);
      SAFENEWWITHCONSTRUCTOR(pTplDC,
			     SingleTplDriveCaller<T>,
			     SingleTplDriveCaller<T>(pDC, t),
			     DMmm);
      
   } else { /* c'e' la deformazione */
      int CurrKeyWord = HP.IsKeyWord();

      /* valore di default (pericoloso!) */
      if (CurrKeyWord == ReadTplDriveKeyWords<T>::UNKNOWNTPLDRIVE) {
	 CurrKeyWord = ReadTplDriveKeyWords<T>::SINGLE;
      }
      
restart:   
   
      switch (CurrKeyWord) {
       case ReadTplDriveKeyWords<T>::SINGLE: {	  
	  T t(0.);
	  	  
	  t = GetT(HP, t);
	  
	  DriveCaller* pDC = ReadDriveData(pDM, HP, pDH);
	  HP.PutKeyTable(K);
	  
	  SAFENEWWITHCONSTRUCTOR(pTplDC, 
				 SingleTplDriveCaller<T>,
				 SingleTplDriveCaller<T>(pDC, t),
				 DMmm);
	  
	  break;
       }
	 
       case ReadTplDriveKeyWords<T>::ARRAY: {
	  unsigned short int iNumDr = HP.GetInt();
	  if (iNumDr == 0) {
	     cerr << "At least one drive is required in array template drive" << endl;
	     THROW(ErrGeneric());
	  } else if (iNumDr == 1) {
	     CurrKeyWord = ReadTplDriveKeyWords<T>::SINGLE;
	     goto restart;
	  } /* else */
	  
	  
	  DrivesArray<T>* pDA = NULL;
	  SAFENEWARR(pDA, DrivesArray<T>, iNumDr, DMmm);
	  
	  for (unsigned short int i = 0; i < iNumDr; i++) {
	     T t(0.);	     	   
	     pDA[i].t = GetT(HP, t);
	     pDA[i].pDriveCaller = ReadDriveData(pDM, HP, pDH);
	  }
	  HP.PutKeyTable(K);
	  
	  SAFENEWWITHCONSTRUCTOR(pTplDC, 
				 ArrayTplDriveCaller<T>,
				 ArrayTplDriveCaller<T>(iNumDr, pDA),
				 DMmm);
	  
	  break;
       }
	 
       default: {
	  cerr << "Unknown template drive type" << endl;
	  THROW(ErrGeneric());
       }
      }
   }
   
   if (pTplDC == NULL) {
      cerr << "Error in allocation of template drive" << endl;
      THROW(ErrMemory());
   }
   
   return pTplDC;
}

#endif /* TPLDRIVE_H */
