/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* drivers */

#ifndef DRIVE_H
#define DRIVE_H


/* include generali */
#include <time.h>
#include <ac/f2c.h>
#include <ac/pthread.h>

/* include per il debug */
#include <myassert.h>
#include <mynewmem.h>

/* include del programma */
#include <mathp.h>
#include <output.h>
#include <solman.h>
#include <withlab.h>
#include <llist.h>


/* Tipi di Drive */
class DriveType {
 public:
   enum Func {
      UNKNOWNFUNC = -1,
	
	CONST = 0,
	STEP,
	DOUBLESTEP,
	RAMP,	
	DOUBLERAMP,
	SINE,
	COSINE,
	
	LASTFUNCTYPE
   };
};      
   
extern const char* psDriveNames[];
extern const char* psReadControlDrivers[];

extern const char* psFuncNames[];


/* Drive - begin */

/* Classe dei drivers, ovvero oggetti che sono in grado di restituire
 * il valore di certi parametri utilizzati da elementi in funzione
 * di altri parametri indipendenti, principalmente del tempo.
 * Per ora la principale classe derivata da questi e' costituita
 * da funzioni ottenute per interpolazione di dati forniti dall'utente.
 * Esistono altri due tipi di drivers che, per varie ragioni, non
 * sono stati derivati da questa classe, ma incorporati
 * direttamente nel DriveHandler (vedi sotto). Questi due tipi restituiscono
 * il valore mediante valutazione di funzioni prestabilite in dipendenza da
 * parametri forniti dal valutatore, dato dalla classe DriveCaller.
 * Nel primo caso la funzione e' "built-in" nel DriveHandler, nel secondo 
 * si ottiene mediante valutazione nel MathParser di una stringa fornita 
 * dall'utente. La tabella dei simboli del MathParser contiene di default
 * la variabile Time, che viene mantenuta aggiornata e puo' essere usata
 * per valutare il valore delle espressioni fornite dall'utente. 
 * 
 * Quindi: 
 * - ogni elemento che usa drivers contiene gli opportuni DriveCaller.
 * - il DriveCaller contiene i dati necessari per chiamare il DriveHandler.
 * - il DriveHandler restituisce il valore richiesto:
 *   - valutando le funzioni "built-in"
 *   - interrogando i singoli Drive.
 * */

class DriveHandler;

class Drive : public WithLabel {
   /* Tipi di drive */
 public:
   enum Type {
      UNKNOWN = -1,
	
	FILEDRIVE = 0,
	
	LASTDRIVETYPE
   };
   
 protected:
   const DriveHandler* pDrvHdl;   
   static doublereal dReturnValue;
   
 public:
   Drive(unsigned int uL, const DriveHandler* pDH);
   virtual ~Drive(void);
   
   /* Tipo del drive (usato solo per debug ecc.) */
   virtual Drive::Type GetDriveType(void) const = 0;
   
   virtual void ServePending(const doublereal& t) = 0;
};

/* Drive - end */


/* DriveHandler - begin */

/* Classe che gestisce la valutazione del valore dei Drive. L'oggetto e' 
 * posseduto dal DataManager, che lo costruisce e lo mette a disposizione
 * dei DriveCaller. Contiene i Drive definiti dall'utente piu' i due 
 * "built-in" che valutano funzioni prestabilite o fornite dall'utente
 * sotto forma di stringa */

class DriveHandler {
   friend class DataManager;
   friend class RandDriveCaller;
   
 private:
   doublereal dTime;
#ifdef USE_MULTITHREAD
   mutable pthread_mutex_t parser_mutex;
#endif /* USE_MULTITHREAD */
   mutable MathParser Parser;
   
   /* variabili predefinite: tempo e variabile generica */
   Var* pTime;
   Var* pVar;
   
   static doublereal dDriveHandlerReturnValue; /* Usato per ritornare un reference */
   
   const VectorHandler* pXCurr;
   const VectorHandler* pXPrimeCurr;
   
   integer iCurrStep;
   
   /* For random drivers */
   class MyRand : public WithLabel {
    protected:
      integer iSteps;
      integer iRand;
      
    public:
      MyRand(unsigned int uLabel, integer iS = 1, integer iR = 0);
      virtual ~MyRand(void);
      
      inline integer iGetSteps(void) const;
      inline integer iGetRand(void) const;
      inline void SetRand(integer& iR);
   };
   
   HardDestructor<MyRand> MyRandD;
   MyLList<MyRand> MyRandLL;
   integer iRandDriveSize;
   MyRand** ppMyRand;
   
   
 protected:
   void SetTime(const doublereal& dt, flag fNewStep = 1);
   void LinkToSolution(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr);
   integer iRandInit(integer iSteps);
   
 public:
   DriveHandler(Table& SymbolTable);
   ~DriveHandler(void);
   
   void PutSymbolTable(Table& T);
   void SetVar(const doublereal& dVar);
   
   doublereal dGet(InputStream& InStr) const;
      
   inline const doublereal& dGetTime(void) const;
   inline long int iGetRand(integer iNumber) const;
};


inline integer 
DriveHandler::MyRand::iGetSteps(void) const 
{
   return iSteps;
}


inline integer 
DriveHandler::MyRand::iGetRand(void) const
{
   return iRand;
}


inline void 
DriveHandler::MyRand::SetRand(integer& iR)
{
   iRand = iR;
}


inline const doublereal& 
DriveHandler::dGetTime(void) const
{
   return dTime;
}


inline long int 
DriveHandler::iGetRand(integer iNumber) const
{
   return ppMyRand[iNumber]->iGetRand();
}

/* DriveHandler - end */



/* DriveCaller - begin */

/* Classe che chiama il drive handler per ottenere il valore del driver ad
 * un dato istante di tempo. Gli oggetti sono posseduti dai singoli elementi.
 * Ogni oggetto contiene i parametri che gli occorrono per la chiamata. */

class DriveCaller : public WithLabel {
 protected:
   const DriveHandler* pDrvHdl;   
   
 public:
   DriveCaller(const DriveHandler* pDH);
   virtual ~DriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const = 0;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const = 0;
   
   /* Restituisce il valore del driver */
   virtual doublereal dGet(const doublereal& dVal) const = 0;
   virtual inline doublereal dGet(void) const;
   
   virtual void SetDrvHdl(const DriveHandler* pDH);
};

inline doublereal 
DriveCaller:: dGet(void) const 
{
	return dGet(pDrvHdl->dGetTime());
}

/* DriveCaller - end */


/* NullDriveCaller - begin */

class NullDriveCaller : public DriveCaller {   
 public:
   NullDriveCaller(void);
   virtual ~NullDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Restituisce il valore del driver */
   virtual inline doublereal dGet(const doublereal& dVal) const;
   virtual inline doublereal dGet(void) const;
};


inline doublereal 
NullDriveCaller::dGet(const doublereal& /* dVal */ ) const
{
   return 0.;
}


inline doublereal
NullDriveCaller::dGet(void) const
{
   return 0.;
}

/* NullDriveCaller - end */


/* OneDriveCaller - begin */

class OneDriveCaller : public DriveCaller {   
 public:
   OneDriveCaller(void);
   virtual ~OneDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Restituisce il valore del driver */
   virtual inline doublereal dGet(const doublereal& dVal) const;
   virtual inline doublereal dGet(void) const;
};


inline doublereal 
OneDriveCaller::dGet(const doublereal& /* dVal */ ) const
{
   return 1.;
}


inline doublereal 
OneDriveCaller::dGet(void) const
{
   return 1.;
}

/* OneDriveCaller - end */


/* ConstDriveCaller - begin */

class ConstDriveCaller : public DriveCaller {
 private:
   doublereal dConst;
   
 public:
   ConstDriveCaller(doublereal d);
   virtual ~ConstDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& /* dVar */ ) const;
   inline doublereal dGet(void) const;
};


inline doublereal 
ConstDriveCaller::dGet(const doublereal& /* dVar */ ) const 
{
   return dConst; 
}
 

inline doublereal
ConstDriveCaller::dGet(void) const
{
   return dConst;
}

/* ConstDriveCaller - end */


/* DriveOwner - begin */

/* Possessore di DriveCaller, ne garantisce la corretta distruzione */

class DriveOwner {
 protected:
   DriveCaller* pDriveCaller;
   
 public:
   DriveOwner(const DriveCaller* pDC = NULL);
   virtual ~DriveOwner(void);
   
   void Set(const DriveCaller* pDC);
   DriveCaller* pGetDriveCaller(void) const;
   doublereal dGet(const doublereal& dVal) const;
   doublereal dGet(void) const;
};

/* DriveOwner - end */


class DataManager;
class MBDynParser;

extern DriveCaller* 
ReadDriveData(const DataManager* pDM, MBDynParser& HP);

#endif /* DRIVE_H */

