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

/* drivers */
 
#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <drive.h>

doublereal DriveCaller::dDriveReturnValue = 0.;
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

DriveHandler::DriveHandler(Table& SymbolTable)
: dTime(0.), 
Parser(SymbolTable),
pTime(NULL), 
pVar(NULL),
pXCurr(NULL), 
pXPrimeCurr(NULL),
iCurrStep(0),
#ifdef DEBUG_MEMMANAGER
MyRandD(DMmm),
#else
MyRandD(),
#endif
MyRandLL(MyRandD),
iRandDriveSize(0),
ppMyRand(NULL)
{	
   /* Inserisce la variabile Time nella tabella dei simboli; sara'
    * mantenuta aggiornata dal DriveHandler */
   NamedValue *v = SymbolTable.Get("Time");
   if (v == NULL) {
      pTime = SymbolTable.Put("Time", Real(0));  
      if (pTime == NULL) {
	 cerr << "DriveHandler::DriveHandler(): error in insertion Time symbol" << endl;	 
	 THROW(ErrGeneric());
      }
   } else {
      if (!v->IsVar()) {
         cerr << "Symbol 'Time' must be a variable" << endl;
         THROW(ErrGeneric());
      }
      pTime = (Var *)v;
   }
   
   /* Inserisce la variabile Var nella tabella dei simboli; sara'
    * mantenuta aggiornata dai DriveCaller attraverso il DriveHandler */
   v = SymbolTable.Get("Var");
   if (v == NULL) {
      pVar = SymbolTable.Put("Var", Real(0)); 
      if (pVar == NULL) {
	 cerr << "DriveHandler::DriveHandler(): error in insertion Var symbol" << endl;	 
	 THROW(ErrGeneric());      
      }
      // pVar->SetVal(Real(0));
   } else {
      if (!v->IsVar()) {
 	 cerr << "Symbol 'Var' must be a variable" << endl;
	 THROW(ErrGeneric());
      }
      pVar = (Var *)v;
   }
   
   /* Calcola il seed di riferimento per i numeri random */       
   srand(time(NULL));
}

DriveHandler::~DriveHandler(void) 
{
   if (iRandDriveSize > 0) { 
      if (ppMyRand != NULL) {
	 SAFEDELETEARR(ppMyRand, DMmm);
      } else {
	 cerr << "Error, random drive data array should exist" << endl;
      }
   }
}


void DriveHandler::SetTime(const doublereal& dt, flag fNewStep)
{      
   dTime = dt;
   
   /* in case of new step */
   if (fNewStep) {
      iCurrStep++;
      
      /* update the random drivers */
      for (long int iCnt = 0; iCnt < iRandDriveSize; iCnt++) {
	 MyRand* pmr = ppMyRand[iCnt];
	 if (iCurrStep%pmr->iGetSteps() == 0) {
	    integer iR = rand();
	    pmr->SetRand(iR);
	 }
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
      
      SAFENEWARR(ppMyRand, MyRand*, iRandDriveSize, DMmm);
      
      MyRand** ppmr = ppMyRand;
      MyRand* pmr = NULL;
      if (!MyRandLL.iGetFirst(pmr)) {
	 cerr << "Error in getting first random drive data" << endl;
	 
	 THROW(ErrGeneric());
      }
      
#ifdef DEBUG
      long int iCnt = 0;
#endif      
      do {
	 ASSERT(++iCnt <= iRandDriveSize); 	  
	 *ppmr++ = pmr;	 
      } while (MyRandLL.iGetNext(pmr));
	                   
      return 0;
   }
   
   /* else, adds a driver */
   MyRand* pmr = NULL;
   integer iNumber = MyRandLL.iGetSize();
   SAFENEWWITHCONSTRUCTOR(pmr, 
			  MyRand, 
			  MyRand((unsigned int)iNumber, iSteps, rand()),
			  DMmm);
   
   if (MyRandLL.iAdd(pmr)) {
      cerr << "Error in insertion of random driver data" << endl;
      THROW(ErrGeneric());
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
 

const doublereal& DriveHandler::dGet(InputStream& InStr) const 
{
   return (dDriveHandlerReturnValue = ((MathParser&)Parser).GetLastStmt(InStr));
}


DriveHandler::MyRand::MyRand(unsigned int uLabel, integer iS, integer iR)
: WithLabel(uLabel), iSteps(iS), iRand(iR) 
{
   NO_OP;
}


DriveHandler::MyRand::~MyRand(void)
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

NullDriveCaller::NullDriveCaller(const DriveHandler* pDH)
: DriveCaller(pDH)
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
   SAFENEWWITHCONSTRUCTOR(pDC, NullDriveCaller, NullDriveCaller(pDrvHdl), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& NullDriveCaller::Restart(ostream& out) const
{      
   return out << " null";
}
 
/* NullDriveCaller - end */


/* OneDriveCaller - begin */

OneDriveCaller::OneDriveCaller(const DriveHandler* pDH)
: DriveCaller(pDH)
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
   SAFENEWWITHCONSTRUCTOR(pDC, OneDriveCaller, OneDriveCaller(pDrvHdl), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& OneDriveCaller::Restart(ostream& out) const
{      
   return out << " one";
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
      SAFEDELETE(pDriveCaller, DMmm);
   }
}


void DriveOwner::Set(const DriveCaller* pDC)
{
   ASSERT(pDC != NULL);
#ifdef DEBUG
   if (pDriveCaller != NULL) {
      DEBUGCERR("warning: the original pointer to a drive caller is not null!" << endl);
   }
#endif
   pDriveCaller = (DriveCaller*)pDC;
}


DriveCaller* DriveOwner::pGetDriveCaller(void) const
{
   return pDriveCaller;
}


const doublereal& DriveOwner::dGet(void) const
{
   return pDriveCaller->dGet();
}

/* DriveOwner - end */
