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

/* file driver */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/fstream>

#include <dataman.h>
#include <filedrv.h>
#include <sockdrv.h>
#ifdef USE_RTAI
#include <rtai_in_drive.h>
#endif /* USE_RTAI */

/* FileDrive - begin */

FileDrive::FileDrive(unsigned int uL, const DriveHandler* pDH,
		const char* const s, integer nd)
: Drive(uL, pDH), sFileName(NULL), iNumDrives(nd), pdVal(NULL)
{   
	ASSERT(s != NULL);
	SAFESTRDUP(sFileName, s);
   	SAFENEWARR(pdVal, doublereal, nd + 1);
   	for (int iCnt = 0; iCnt <= nd; iCnt++) {
      		pdVal[iCnt] = 0.;
   	}
}


FileDrive::~FileDrive(void)
{
	if (sFileName != NULL) {
		SAFEDELETEARR(sFileName);
	}

   	if (pdVal != NULL) {
      		SAFEDELETEARR(pdVal);
   	}
}


Drive::Type
FileDrive::GetDriveType(void) const
{
   return Drive::FILEDRIVE;
}

doublereal
FileDrive::dGet(const doublereal& /* t */ , int i) const
{
   	ASSERT(i > 0 && i <= iNumDrives);
   	return pdVal[i];
}

/* FileDrive - end */


/* FileDriveCaller - begin */

FileDriveCaller::FileDriveCaller(const DriveHandler* pDH, 
				 const FileDrive* p, integer i,
				 const doublereal& da)
: DriveCaller(pDH), pFileDrive((FileDrive*)p), iNumDrive(i), dAmplitude(da)
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
			  FileDriveCaller(pDrvHdl,  pFileDrive, 
				  iNumDrive, dAmplitude));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
FileDriveCaller::Restart(std::ostream& out) const
{
	out << " file, " << pFileDrive->GetLabel()
		<< ", " << iNumDrive;
	if (dAmplitude != 1.) {
		out << ", amplitude, " << dAmplitude;
	}
	return out;
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
   
   std::ifstream in(sFileName);
   if (!in) {
      std::cerr << "can't open file \""
	<< sFileName << "\"" << std::endl;
      THROW(ErrGeneric());
   }
   
   SAFENEWARR(pd, doublereal, iNumDrives*iNumSteps);
   SAFENEWARR(pvd, doublereal*, iNumDrives+1);
   
   /* Attenzione: il primo puntatore e' vuoto 
    * (ne e' stato allocato uno in piu'),
    * cosi' i drives possono essere numerati da 1 a n */
   for (integer i = iNumDrives; i-- > 0; ) {
      pvd[i+1] = pd+i*iNumSteps;
   }

   /*
    * Mangia gli eventuali commenti iniziali
    */
   char c = '\0';
   while (in.get(c), c == '#') {
      char buf[1024];

      do {
	 in.getline(buf, sizeof(buf));
      } while (strlen(buf) == sizeof(buf) - 1 && buf[sizeof(buf) - 1] != '\n');
   }

   if (c != '#') {
      in.putback(c);
   }
   
   for (integer j = 0; j < iNumSteps; j++) {
      for (integer i = 1; i <= iNumDrives; i++) {
	 in >> pvd[i][j];
	 if (in.eof()) {
	    std::cerr << "unexpected end of file" << std::endl;
	    THROW(ErrGeneric());
	 }
      }
   }
}


FixedStepFileDrive::~FixedStepFileDrive(void)
{
   SAFEDELETEARR(pd);
   SAFEDELETEARR(pvd);
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& FixedStepFileDrive::Restart(std::ostream& out) const
{
   return out << "FixedStepFileDrive::Restart(): not implemented yet!" << std::endl;
}


void
FixedStepFileDrive::ServePending(const doublereal& t)
{
	integer j1 = integer(floor((t-dT0)/dDT));
	integer j2 = j1+1;
   
	if (j2 < 0 || j1 > iNumSteps) {
		for (int i = 1; i <= iNumDrives; i++) {
			pdVal[i] = 0.;
		}
	} else {
		doublereal dt1 = dT0+j1*dDT;
		doublereal dt2 = dt1+dDT;

		for (int i = 1; i <= iNumDrives; i++) {
   			pdVal[i] = (pvd[i][j2]*(t-dt1)-pvd[i][j1]*(t-dt2))/dDT;
		}
	}
}

/* FixedStepFileDrive - end */


/* legge i drivers tipo file */

Drive* ReadFileDriver(DataManager* pDM, 
		      MBDynParser& HP,
		      unsigned int uLabel)
{
   Drive* pDr = NULL;
   
   const char* sKeyWords[] = {
      "fixed" "step",
	"socket",
	"rtai" "input",

	NULL
   };
   
   enum KeyWords {
      UNKNOWN = -1,
	
	FIXEDSTEP = 0,
	SOCKET,
	RTAIINPUT,
	
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
						 t0, dt));
       break;
    } 
      
    case SOCKET: {
#ifdef USE_SOCKET_DRIVES
       
       integer idrives = HP.GetInt();
       unsigned short int port = MBDynSocketDrivePort;
       const char *path = NULL;
      
       if (HP.IsKeyWord("local")) {
	  path = HP.GetFileName();
	  ASSERT(path != NULL);
       } else if (HP.IsKeyWord("port")) {
          port = HP.GetInt();      
       }
      
       if (path == NULL) {
          AuthMethod* pAuth = ReadAuthMethod(pDM, HP);
          if (pAuth == NULL) {
	     std::cerr << "need an authentication method at line " << HP.GetLineData() << std::endl;
	     THROW(ErrGeneric());
	  }
          SAFENEWWITHCONSTRUCTOR(pDr,
			  SocketDrive,
			  SocketDrive(uLabel, pDM->pGetDrvHdl(),
				  port, pAuth, idrives));
       } else {
	  SAFENEWWITHCONSTRUCTOR(pDr,
			  SocketDrive,
			  SocketDrive(uLabel, pDM->pGetDrvHdl(), 
				  path, idrives));
       }
             
#else /* USE_SOCKET_DRIVES */
       std::cerr << "Sorry, socket drives not supported." << std::endl;
       THROW(ErrGeneric());
#endif /* USE_SOCKET_DRIVES */
       break;
    }

    case RTAIINPUT: {
#ifdef USE_RTAI
       pDr = ReadRTAIInDrive(pDM, HP, uLabel);
#else /* ! USE_RTAI */
       std::cerr << "Sorry, RTAI input requires configure --with-rtai"
	       << std::endl;
       THROW(ErrGeneric());
#endif /* ! USE_RTAI */
       break;
    }
     
    default:     
      std::cerr << "unknown file drive at line " << HP.GetLineData() << std::endl;
      THROW(ErrGeneric());
   }
      
   /* scrittura dei dati specifici */

   return pDr;
} /* End of ReadFileDriver */
