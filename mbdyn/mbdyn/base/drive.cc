/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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
 
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <drive.h>

doublereal Drive::dReturnValue = 0.;
doublereal DriveHandler::dDriveHandlerReturnValue = 0.;

/* Drive - begin */

Drive::Drive(unsigned int uL, const DriveHandler* pDH)
: WithLabel(uL), pDrvHdl(pDH) 
{
   NO_OP; 
}


Drive::~Drive(void) {
   NO_OP; 
}

/* Drive - end */


/* DriveHandler - begin */

DriveHandler::DriveHandler(MathParser& mp)
: dTime(0.), 
Parser(mp),
pTime(NULL), 
pVar(NULL),
pXCurr(NULL), 
pXPrimeCurr(NULL),
iCurrStep(0),
MyMeterD(),
MyMeterLL(MyMeterD),
iMeterDriveSize(0),
ppMyMeter(NULL),
MyRandD(),
MyRandLL(MyRandD),
iRandDriveSize(0),
ppMyRand(NULL)
{	
#ifdef USE_MULTITHREAD
   if (pthread_mutex_init(&parser_mutex, NULL)) {
      silent_cerr("DriveHandler::DriveHandler(): mutex init failed"
		      << std::endl);
      throw ErrGeneric();
   }
#endif /* USE_MULTITHREAD */

   /* Inserisce la variabile Time nella tabella dei simboli; sara'
    * mantenuta aggiornata dal DriveHandler */
   NamedValue *v = Parser.GetSymbolTable().Get("Time");
   if (v == NULL) {
      pTime = Parser.GetSymbolTable().Put("Time", Real(0));  
      if (pTime == NULL) {
	 silent_cerr("DriveHandler::DriveHandler(): "
			 "error in Time symbol insertion" << std::endl); 
	 throw ErrGeneric();
      }
   } else {
      if (!v->IsVar()) {
	 silent_cerr("Symbol 'Time' must be a variable" << std::endl);
         throw ErrGeneric();
      }
      pTime = (Var *)v;
   }
   
   /* Inserisce la variabile Var nella tabella dei simboli; sara'
    * mantenuta aggiornata dai DriveCaller attraverso il DriveHandler */
   v = Parser.GetSymbolTable().Get("Var");
   if (v == NULL) {
      pVar = Parser.GetSymbolTable().Put("Var", Real(0)); 
      if (pVar == NULL) {
	 silent_cerr("DriveHandler::DriveHandler(): "
			 "error in Var symbol insertion" << std::endl);
	 throw ErrGeneric();
      }

   } else {
      if (!v->IsVar()) {
	 silent_cerr("Symbol 'Var' must be a variable" << std::endl);
	 throw ErrGeneric();
      }
      pVar = (Var *)v;
   }
   
   /* Calcola il seed di riferimento per i numeri random */       
   srand(time(NULL));
}

DriveHandler::~DriveHandler(void) 
{
#ifdef USE_MULTITHREAD
   pthread_mutex_destroy(&parser_mutex);
#endif /* USE_MULTITHREAD */

   if (iMeterDriveSize > 0) { 
      if (ppMyMeter != NULL) {
	 SAFEDELETEARR(ppMyMeter);
      } else {
	 silent_cerr("Error, meter drive data array should exist" << std::endl);
      }
   }

   if (iRandDriveSize > 0) { 
      if (ppMyRand != NULL) {
	 SAFEDELETEARR(ppMyRand);
      } else {
	 silent_cerr("Error, random drive data array should exist" << std::endl);
      }
   }
}


void DriveHandler::SetTime(const doublereal& dt, flag fNewStep)
{      
   dTime = dt;
   
   /* in case of new step */
   if (fNewStep) {
      iCurrStep++;
      
      /* update the meter drivers */
      for (long int iCnt = 0; iCnt < iMeterDriveSize; iCnt++) {
	 MyMeter* pmm = ppMyMeter[iCnt];
	 pmm->SetMeter();
      }      
      
      /* update the random drivers */
      for (long int iCnt = 0; iCnt < iRandDriveSize; iCnt++) {
	 MyRand* pmr = ppMyRand[iCnt];
	 pmr->SetMeter();
      }      
   }
}


void DriveHandler::LinkToSolution(const VectorHandler& XCurr,
				  const VectorHandler& XPrimeCurr) 
{
   (VectorHandler*&)pXCurr = (VectorHandler*)&XCurr;
   (VectorHandler*&)pXPrimeCurr = (VectorHandler*)&XPrimeCurr;
}
   

/*
 * se iSteps == 0 inizializza la lista dei random drivers;
 * se iSteps != 0 allora alloca un nuovo gestore dei dati del
 * random driver, e ritorna il numero d'ordine
 */
integer DriveHandler::iRandInit(integer iSteps) 
{
   if (iSteps == 0) {
      /* initialises the structure */
      iRandDriveSize = MyRandLL.iGetSize();
      if (iRandDriveSize == 0) {
	 return 0;
      }
      
      SAFENEWARR(ppMyRand, MyRand*, iRandDriveSize);
      
      MyRand** ppmr = ppMyRand;
      MyRand* pmr = NULL;
      if (!MyRandLL.GetFirst(pmr)) {
	 silent_cerr("Error in getting first random drive data" << std::endl);
	 
	 throw ErrGeneric();
      }
      
#ifdef DEBUG
      long int iCnt = 0;
#endif      
      do {
	 ASSERT(++iCnt <= iRandDriveSize); 	  
	 *ppmr++ = pmr;	 
      } while (MyRandLL.GetNext(pmr));
	                   
      return 0;
   }
   
   /* else, adds a driver */
   MyRand* pmr = NULL;
   integer iNumber = MyRandLL.iGetSize();
   SAFENEWWITHCONSTRUCTOR(pmr, 
			  MyRand, 
			  MyRand((unsigned int)iNumber, iSteps, rand()));
   
   if (MyRandLL.Add(pmr)) {
      silent_cerr("Error in insertion of random driver data" << std::endl);
      throw ErrGeneric();
   }
   
   return iNumber;
}


/*
 * se iSteps == 0 inizializza la lista dei meter drivers;
 * se iSteps != 0 allora alloca un nuovo gestore dei dati del
 * meter driver, e ritorna il numero d'ordine
 */
integer DriveHandler::iMeterInit(integer iSteps) 
{
   if (iSteps == 0) {
      /* initialises the structure */
      iMeterDriveSize = MyMeterLL.iGetSize();
      if (iMeterDriveSize == 0) {
	 return 0;
      }
      
      SAFENEWARR(ppMyMeter, MyMeter*, iMeterDriveSize);
      
      MyMeter** ppmm = ppMyMeter;
      MyMeter* pmm = NULL;
      if (!MyMeterLL.GetFirst(pmm)) {
	 silent_cerr("Error in getting first meter drive data" << std::endl);
	 
	 throw ErrGeneric();
      }
      
#ifdef DEBUG
      long int iCnt = 0;
#endif      
      do {
	 ASSERT(++iCnt <= iMeterDriveSize);
	 *ppmm++ = pmm;	 
      } while (MyMeterLL.GetNext(pmm));
	                   
      return 0;
   }
   
   /* else, adds a driver */
   MyMeter* pmm = NULL;
   integer iNumber = MyMeterLL.iGetSize();
   SAFENEWWITHCONSTRUCTOR(pmm, 
			  MyMeter, 
			  MyMeter((unsigned int)iNumber, iSteps));
   
   if (MyMeterLL.Add(pmm)) {
      silent_cerr("Error in insertion of meter driver data" << std::endl);
      throw ErrGeneric();
   }
   
   return iNumber;
}


void DriveHandler::PutSymbolTable(Table& T) 
{
   Parser.PutSymbolTable(T);
}


void DriveHandler::SetVar(const doublereal& dVar)
{
   ASSERT(pVar != NULL);
   pVar->SetVal(dVar);
}
 

doublereal DriveHandler::dGet(InputStream& InStr) const 
{
   doublereal d;
#ifdef USE_MULTITHREAD
   pthread_mutex_lock(&parser_mutex);
#endif /* USE_MULTITHREAD */
   d = Parser.GetLastStmt(InStr);
#ifdef USE_MULTITHREAD
   pthread_mutex_unlock(&parser_mutex);
#endif /* USE_MULTITHREAD */
   return d;
}


DriveHandler::MyRand::MyRand(unsigned int uLabel, integer iS, integer iR)
: MyMeter(uLabel, iS), iRand(iR) 
{
   NO_OP;
}


DriveHandler::MyRand::~MyRand(void)
{
   NO_OP;
}

 
DriveHandler::MyMeter::MyMeter(unsigned int uLabel, integer iS)
: WithLabel(uLabel), iSteps(iS), iCurr(0)
{
   NO_OP;
}


DriveHandler::MyMeter::~MyMeter(void)
{
   NO_OP;
}
 
/* DriveHandler - end */


/* DriveCaller - begin */

DriveCaller::DriveCaller(const DriveHandler* pDH)
: pDrvHdl(pDH)
{
   NO_OP;
}


DriveCaller::~DriveCaller(void)
{
   NO_OP;
}
 

void DriveCaller::SetDrvHdl(const DriveHandler* pDH)
{
   (DriveHandler*&)pDrvHdl = (DriveHandler*)pDH;
}

/* DriveCaller - end */


/* NullDriveCaller - begin */

NullDriveCaller::NullDriveCaller(void)
: DriveCaller(0)
{
   NO_OP;
}


NullDriveCaller::~NullDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* NullDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEW(pDC, NullDriveCaller);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& NullDriveCaller::Restart(std::ostream& out) const
{      
   return out << "null";
}
 
/* NullDriveCaller - end */


/* OneDriveCaller - begin */

OneDriveCaller::OneDriveCaller(void)
: DriveCaller(0)
{
   NO_OP;
}


OneDriveCaller::~OneDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* OneDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEW(pDC, OneDriveCaller);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& OneDriveCaller::Restart(std::ostream& out) const
{      
   return out << "one";
}
 
/* OneDriveCaller - end */


/* DriveOwner - begin */

DriveOwner::DriveOwner(const DriveCaller* pDC)
: pDriveCaller((DriveCaller*)pDC) 
{
   NO_OP;
}
 

DriveOwner::~DriveOwner(void)
{ 
   ASSERT(pDriveCaller != NULL);
   
   if (pDriveCaller != NULL) {
      SAFEDELETE(pDriveCaller);
   }
}


void
DriveOwner::Set(const DriveCaller* pDC)
{
	ASSERT(pDC != NULL);
	if (pDriveCaller != NULL) {
		DEBUGCERR("warning: the original pointer to a drive caller is not null!" << std::endl);
		SAFEDELETE(pDriveCaller);
	}
	pDriveCaller = (DriveCaller *)pDC;
}


DriveCaller* DriveOwner::pGetDriveCaller(void) const
{
   return pDriveCaller;
}


doublereal DriveOwner::dGet(const doublereal& dVal) const
{
   return pDriveCaller->dGet(dVal);
}


doublereal DriveOwner::dGet(void) const
{
   return pDriveCaller->dGet();
}

/* DriveOwner - end */

