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

/* file driver */

#include <mbconfig.h>

#include <dataman.h>
#include <filedrv.h>
#include <sockdrv.h>

/* FileDrive - begin */

FileDrive::FileDrive(unsigned int uL, const DriveHandler* pDH,
		     const char* const s, integer nd)
: Drive(uL, pDH), sFileName(NULL), iNumDrives(nd)
{   
   ASSERT(s != NULL);
   SAFENEWARR(sFileName, char, strlen(s)+1, DMmm);
   strcpy(sFileName, s);
}


FileDrive::~FileDrive(void)
{
   SAFEDELETEARR(sFileName, DMmm);
}


DriveType::Type FileDrive::GetDriveType(void) const
{
   return DriveType::FILEDRIVE;
}

/* FileDrive - end */


/* FileDriveCaller - begin */

FileDriveCaller::FileDriveCaller(const DriveHandler* pDH, 
				 const FileDrive* p, integer i)
: DriveCaller(pDH), pFileDrive((FileDrive*)p), iNumDrive(i)
{
   ASSERT(pFileDrive != NULL);
   ASSERT(iNumDrive > 0 && iNumDrive <= pFileDrive->iGetNumDrives());
}


FileDriveCaller::~FileDriveCaller(void)
{
   NO_OP;
}

/* Copia */
DriveCaller* FileDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  FileDriveCaller,
			  FileDriveCaller(pDrvHdl,  pFileDrive, iNumDrive),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& FileDriveCaller::Restart(ostream& out) const
{
   return out << " file, " 
     << pFileDrive->GetLabel() << ", " 
     << iNumDrive << endl;
}
   
/* FileDriveCaller - end */


/* FixedStepFileDrive - begin */

FixedStepFileDrive::FixedStepFileDrive(unsigned int uL, 
				       const DriveHandler* pDH,
				       const char* const sFileName,
				       integer is, integer nd,
				       doublereal t0, doublereal dt)
: FileDrive(uL, pDH, sFileName, nd),
dT0(t0), dDT(dt), iNumSteps(is), pd(NULL), pvd(NULL)
{
   ASSERT(iNumSteps > 0);
   ASSERT(iNumDrives > 0);
   ASSERT(sFileName != NULL);
   ASSERT(dDT > 0.);
   
   ifstream in(sFileName);
   if (!in) {
      cerr << "can't open file \""
	<< sFileName << "\"" << endl;
      THROW(ErrGeneric());
   }
   
   SAFENEWARR(pd, doublereal, iNumDrives*iNumSteps, DMmm);
   SAFENEWARR(pvd, doublereal*, iNumDrives+1, DMmm);
   
   /* Attenzione: il primo puntatore e' vuoto 
    * (ne e' stato allocato uno in piu'),
    * cosi' i drives possono essere numerati da 1 a n */
   for (integer i = iNumDrives; i-- > 0; ) {
      pvd[i+1] = pd+i*iNumSteps;
   }
   
   for (integer j = 0; j < iNumSteps; j++) {
      for (integer i = 1; i <= iNumDrives; i++) {
	 in >> pvd[i][j];
	 if (in.eof()) {
	    cerr << "unexpected end of file" << endl;
	    THROW(ErrGeneric());
	 }
      }
   }
}


FixedStepFileDrive::~FixedStepFileDrive(void)
{
   SAFEDELETEARR(pd, DMmm);
   SAFEDELETEARR(pvd, DMmm);
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& FixedStepFileDrive::Restart(ostream& out) const
{
   return out << "FixedStepFileDrive::Restart(): not implemented yet!" << endl;
}


const doublereal& FixedStepFileDrive::dGet(const doublereal& t, int i) const
{
   ASSERT(i > 0 && i <= iNumDrives);

   integer j1 = integer(floor((t-dT0)/dDT));
   integer j2 = j1+1;
   
   if (j2 < 0 || j1 > iNumSteps) {
      return (Drive::dReturnValue = 0.);
   }
   
   doublereal dt1 = dT0+j1*dDT;
   doublereal dt2 = dt1+dDT;
   
   return (Drive::dReturnValue = (pvd[i][j2]*(t-dt1)-pvd[i][j1]*(t-dt2))/dDT);
}


void FixedStepFileDrive::ServePending(void)
{
   NO_OP;
}

/* FixedStepFileDrive - end */


/* legge i drivers tipo file */

Drive* ReadFileDriver(DataManager* pDM, 
		      MBDynParser& HP,
		      unsigned int uLabel)
{
   Drive* pDr = NULL;
   
   const char* sKeyWords[] = {
      "fixedstep",
	"socket"
   };
   
   enum KeyWords {
      UNKNOWN = -1,
	
	FIXEDSTEP = 0,
	SOCKET,
	
	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);   
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord = KeyWords(HP.GetWord());
   
   switch (CurrKeyWord) {
      
    case FIXEDSTEP: {
     
       integer isteps = HP.GetInt();
       integer idrives = HP.GetInt();
       doublereal t0 = HP.GetReal();
       doublereal dt = HP.GetReal();
       const char* filename = HP.GetFileName();
       
       SAFENEWWITHCONSTRUCTOR(pDr,
			      FixedStepFileDrive,
			      FixedStepFileDrive(uLabel, pDM->pGetDrvHdl(),
						 filename, isteps, idrives,
						 t0, dt), 
			      DMmm);
       break;
    } 
      
    case SOCKET: {
#ifdef USE_SOCKET_DRIVES
       
       integer idrives = HP.GetInt();
      
       unsigned short int port = MBDynSocketDrivePort;
       if (HP.IsKeyWord("port")) {
	  port = HP.GetInt();      
       }
       
       AuthMethod* pAuth = ReadAuthMethod(pDM, HP);
       if (pAuth == NULL) {
	  cerr << "need an authentication method at line " << HP.GetLineData() << endl;
	  THROW(ErrGeneric());
       }
       
       SAFENEWWITHCONSTRUCTOR(pDr,
			      SocketDrive,
			      SocketDrive(uLabel, pDM->pGetDrvHdl(),
					  port, pAuth, idrives), 
			      DMmm);    
             
#else /* USE_SOCKET_DRIVES */
       cerr << "Sorry, socket drives not supported." << endl;
       THROW(ErrGeneric());
#endif /* USE_SOCKET_DRIVES */
       break;
    }
      
    default:     
      cerr << "unknown file drive at line " << HP.GetLineData() << endl;
      THROW(ErrGeneric());
   }
      
   /* scrittura dei dati specifici */

   return pDr;
} /* End of ReadFileDriver */
