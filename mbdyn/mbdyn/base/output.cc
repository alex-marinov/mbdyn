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

/* gestore dell'output */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <output.h>

/* OutputHandler - begin */

#ifndef OUTPUT_PRECISION
#define OUTPUT_PRECISION 6
#endif /* OUTPUT_PRECISION */

const int iDefaultWidth = OUTPUT_PRECISION+6;
const int iDefaultPrecision = OUTPUT_PRECISION;

const char* psExt[] = {
   ".out",		/*  0 */
   ".mov",
   ".ele",
   ".abs",
   ".ine",
   ".jnt",		/*  5 */
   ".frc",
   ".act",
   ".rot",
   ".rst",
   ".rst.X",		/* 10 */
   ".aer",  
   ".hyd",
   ".prs",
   ".usr",
   ".gen",		/* 15 */
   ".par",
   ".res",
   ".ada",
   ".amd",
   ".rfm",		/* 20 */
   ".log",
   ".air",
   ".prm",
   ".ext",
   NULL
};


/* Costruttore senza inizializzazione */
OutputHandler::OutputHandler(void) 
: FileName(NULL),
ofOutput(),
ofStrNodes(),
ofElectric(),
ofAbstract(),
ofInertia(),
ofJoints(),
ofForces(),
ofBeams(),
ofRotors(),
ofRestart(),
ofRestartXSol(),
ofAerodynamic(),
ofHydraulic(),
ofPresNodes(),
ofLoadable(),
ofGenels(),
ofPartition(),
ofAdamsRes(),
ofAdamsCmd(),
ofAeroModals(),
ofReferenceFrames(),
ofLog(),
ofAirProps(),
ofParameters(),
ofExternals(),
iCurrWidth(iDefaultWidth), iCurrPrecision(iDefaultPrecision), nCurrRestartFile(0)
{
   OutData[OUTPUT].UseDefaultPrecision = false;
   OutData[OUTPUT].UseScientific = false;
   OutData[OUTPUT].pof = &ofOutput;
   
   OutData[STRNODES].UseDefaultPrecision = true;
   OutData[STRNODES].UseScientific = true;
   OutData[STRNODES].pof = &ofStrNodes;
   
   OutData[ELECTRIC].UseDefaultPrecision = true;
   OutData[ELECTRIC].UseScientific = true;
   OutData[ELECTRIC].pof= &ofElectric;
   
   OutData[ABSTRACT].UseDefaultPrecision = true;
   OutData[ABSTRACT].UseScientific = true;
   OutData[ABSTRACT].pof = &ofAbstract;
   
   OutData[INERTIA].UseDefaultPrecision = true;
   OutData[INERTIA].UseScientific = true;
   OutData[INERTIA].pof = &ofInertia;
   
   OutData[JOINTS].UseDefaultPrecision = true;
   OutData[JOINTS].UseScientific = true;
   OutData[JOINTS].pof = &ofJoints;
   
   OutData[FORCES].UseDefaultPrecision = true;
   OutData[FORCES].UseScientific = true;
   OutData[FORCES].pof = &ofForces;
   
   OutData[BEAMS].UseDefaultPrecision = true;
   OutData[BEAMS].UseScientific = true;
   OutData[BEAMS].pof = &ofBeams;
   
   OutData[ROTORS].UseDefaultPrecision = true;
   OutData[ROTORS].UseScientific = true;
   OutData[ROTORS].pof = &ofRotors;
   
   OutData[RESTART].UseDefaultPrecision = true;
   OutData[RESTART].UseScientific = true;
   OutData[RESTART].pof = &ofRestart;

   OutData[RESTARTXSOL].UseDefaultPrecision = true;
   OutData[RESTARTXSOL].UseScientific = true;
   OutData[RESTARTXSOL].pof = &ofRestartXSol;

   OutData[AERODYNAMIC].UseDefaultPrecision = true;
   OutData[AERODYNAMIC].UseScientific = true;
   OutData[AERODYNAMIC].pof = &ofAerodynamic;
   
   OutData[HYDRAULIC].UseDefaultPrecision = true;
   OutData[HYDRAULIC].UseScientific = true;
   OutData[HYDRAULIC].pof = &ofHydraulic;

   OutData[PRESNODES].UseDefaultPrecision = true;
   OutData[PRESNODES].UseScientific = true;
   OutData[PRESNODES].pof = &ofPresNodes;

   OutData[LOADABLE].UseDefaultPrecision = true;
   OutData[LOADABLE].UseScientific = true;
   OutData[LOADABLE].pof = &ofLoadable;
   
   OutData[GENELS].UseDefaultPrecision = true;
   OutData[GENELS].UseScientific = true;
   OutData[GENELS].pof = &ofGenels;

   OutData[PARTITION].UseDefaultPrecision = true;
   OutData[PARTITION].UseScientific = true;
   OutData[PARTITION].pof = &ofPartition;
   
   OutData[ADAMSRES].UseDefaultPrecision = true;
   OutData[ADAMSRES].UseScientific = true;
   OutData[ADAMSRES].pof = &ofAdamsRes;

   OutData[ADAMSCMD].UseDefaultPrecision = true;
   OutData[ADAMSCMD].UseScientific = true;
   OutData[ADAMSCMD].pof = &ofAdamsCmd;

   OutData[AEROMODALS].UseDefaultPrecision = true;
   OutData[AEROMODALS].UseScientific = true;
   OutData[AEROMODALS].pof = &ofAeroModals;
   
   OutData[REFERENCEFRAMES].UseDefaultPrecision = true;
   OutData[REFERENCEFRAMES].UseScientific = true;
   OutData[REFERENCEFRAMES].pof = &ofReferenceFrames;
   
   OutData[LOG].UseDefaultPrecision = false;
   OutData[LOG].UseScientific = false;
   OutData[LOG].pof = &ofLog;
   
   OutData[AIRPROPS].UseDefaultPrecision = true;
   OutData[AIRPROPS].UseScientific = true;
   OutData[AIRPROPS].pof = &ofAirProps;

   OutData[PARAMETERS].UseDefaultPrecision = true;
   OutData[PARAMETERS].UseScientific = true;
   OutData[PARAMETERS].pof = &ofParameters;

   OutData[EXTERNALS].UseDefaultPrecision = true;
   OutData[EXTERNALS].UseScientific = true;
   OutData[EXTERNALS].pof = &ofExternals;

   for (int iCnt = 0; iCnt < LASTFILE; iCnt++) {
      OutData[iCnt].IsOpen = false;
   }
}



/* Costruttore con inizializzazione */
OutputHandler::OutputHandler(const char* sFName, int iExtNum)
: FileName(sFName, iExtNum),
ofOutput(_sPutExt((char*)(psExt[OUTPUT]))),
ofStrNodes(),
ofElectric(),
ofAbstract(),
ofInertia(),
ofJoints(),
ofForces(),
ofBeams(),
ofRotors(),
ofRestart(),
ofRestartXSol(),
ofAerodynamic(),
ofHydraulic(),
ofPresNodes(),
ofLoadable(),
ofGenels(),
ofPartition(),
ofAdamsRes(),
ofAdamsCmd(),
ofAeroModals(),
ofReferenceFrames(),
ofLog(),
ofAirProps(),
ofParameters(),
ofExternals(),
iCurrWidth(iDefaultWidth), iCurrPrecision(iDefaultPrecision), nCurrRestartFile(0)
{
   OutData[OUTPUT].UseDefaultPrecision = false;
   OutData[OUTPUT].UseScientific = false;
   OutData[OUTPUT].pof = &ofOutput;
   ofOutput.width(iCurrWidth);
   ofOutput.precision(iCurrPrecision);
   
   OutData[STRNODES].UseDefaultPrecision = true;
   OutData[STRNODES].UseScientific = true;
   OutData[STRNODES].pof = &ofStrNodes;
   
   OutData[ELECTRIC].UseDefaultPrecision = true;
   OutData[ELECTRIC].UseScientific = true;
   OutData[ELECTRIC].pof = &ofElectric;
   
   OutData[ABSTRACT].UseDefaultPrecision = true;
   OutData[ABSTRACT].UseScientific = true;
   OutData[ABSTRACT].pof = &ofAbstract;
   
   OutData[INERTIA].UseDefaultPrecision = true;
   OutData[INERTIA].UseScientific = true;
   OutData[INERTIA].pof = &ofInertia;
   
   OutData[JOINTS].UseDefaultPrecision = true;
   OutData[JOINTS].UseScientific = true;
   OutData[JOINTS].pof = &ofJoints;
   
   OutData[FORCES].UseDefaultPrecision = true;
   OutData[FORCES].UseScientific = true;
   OutData[FORCES].pof = &ofForces;
   
   OutData[BEAMS].UseDefaultPrecision = true;
   OutData[BEAMS].UseScientific = true;
   OutData[BEAMS].pof = &ofBeams;
   
   OutData[ROTORS].UseDefaultPrecision = true;
   OutData[ROTORS].UseScientific = true;
   OutData[ROTORS].pof = &ofRotors;
   
   OutData[RESTART].UseDefaultPrecision = true;
   OutData[RESTART].UseScientific = true;
   OutData[RESTART].pof = &ofRestart;
   
   OutData[RESTARTXSOL].UseDefaultPrecision = true;
   OutData[RESTARTXSOL].UseScientific = true;
   OutData[RESTARTXSOL].pof = &ofRestartXSol;

   OutData[AERODYNAMIC].UseDefaultPrecision = true;
   OutData[AERODYNAMIC].UseScientific = true;
   OutData[AERODYNAMIC].pof = &ofAerodynamic;
   
   OutData[HYDRAULIC].UseDefaultPrecision = true;
   OutData[HYDRAULIC].UseScientific = true;
   OutData[HYDRAULIC].pof = &ofHydraulic;

   OutData[PRESNODES].UseDefaultPrecision = true;
   OutData[PRESNODES].UseScientific = true;
   OutData[PRESNODES].pof = &ofPresNodes;

   OutData[LOADABLE].UseDefaultPrecision = true;
   OutData[LOADABLE].UseScientific = true;
   OutData[LOADABLE].pof = &ofLoadable;
   
   OutData[GENELS].UseDefaultPrecision = true;
   OutData[GENELS].UseScientific = true;
   OutData[GENELS].pof = &ofGenels;

   OutData[PARTITION].UseDefaultPrecision = true;
   OutData[PARTITION].UseScientific = true;
   OutData[PARTITION].pof = &ofPartition;
   
   OutData[ADAMSRES].UseDefaultPrecision = true;
   OutData[ADAMSRES].UseScientific = true;
   OutData[ADAMSRES].pof = &ofAdamsRes;

   OutData[ADAMSCMD].UseDefaultPrecision = true;
   OutData[ADAMSCMD].UseScientific = true;
   OutData[ADAMSCMD].pof = &ofAdamsCmd;
   
   OutData[AEROMODALS].UseDefaultPrecision = true;
   OutData[AEROMODALS].UseScientific = true;
   OutData[AEROMODALS].pof = &ofAeroModals;
   
   OutData[REFERENCEFRAMES].UseDefaultPrecision = true;
   OutData[REFERENCEFRAMES].UseScientific = true;
   OutData[REFERENCEFRAMES].pof = &ofReferenceFrames;

   OutData[LOG].UseDefaultPrecision = false;
   OutData[LOG].UseScientific = false;
   OutData[LOG].pof = &ofLog;
   ofLog.width(iCurrWidth);
   ofLog.precision(iCurrPrecision);
   
   OutData[AIRPROPS].UseDefaultPrecision = true;
   OutData[AIRPROPS].UseScientific = true;
   OutData[AIRPROPS].pof = &ofAirProps;

   OutData[PARAMETERS].UseDefaultPrecision = true;
   OutData[PARAMETERS].UseScientific = true;
   OutData[PARAMETERS].pof = &ofParameters;

   OutData[EXTERNALS].UseDefaultPrecision = true;
   OutData[EXTERNALS].UseScientific = true;
   OutData[EXTERNALS].pof = &ofExternals;

   for (int iCnt = 0; iCnt < LASTFILE; iCnt++) {
      OutData[iCnt].IsOpen = false;
   }
}


/* Inizializzazione */
void OutputHandler::Init(const char* sFName, int iExtNum)
{
   FileName::iInit(sFName, iExtNum);
   OutputOpen();
   LogOpen();
}


/* Distruttore */
OutputHandler::~OutputHandler(void) 
{
   for(int iCnt = 0; iCnt < LASTFILE; iCnt++) {
      if (OutData[iCnt].IsOpen) {
#ifdef HAVE_ISOPEN
	 ASSERT(OutData[iCnt].pof->is_open());
#endif /* HAVE_ISOPEN */
	 OutData[iCnt].pof->close();
	 OutData[iCnt].IsOpen = false;
      }
   }
}
   
   
/* Aggiungere qui le funzioni che aprono i singoli stream */
bool
OutputHandler::Open(const OutputHandler::OutFiles out)
{
   if(!OutData[out].IsOpen) {
#ifdef HAVE_ISOPEN
      ASSERT(!OutData[out].pof->is_open());
#endif /* HAVE_ISOPEN */

      /* Apre lo stream */
      OutData[out].pof->open(_sPutExt((char*)(psExt[out])));

      if(!(*OutData[out].pof)) {
	 silent_cerr("Unable to open file \"" << _sPutExt((char*)(psExt[out]))
	   << "\"" << std::endl);
	 throw ErrFile();
      }

      OutData[out].IsOpen = true;
      
      /* Setta la formattazione dei campi */
      if(OutData[out].UseDefaultPrecision) {		        
	 OutData[out].pof->precision(iCurrPrecision);
      }
      
      /* Setta la notazione */
      if(OutData[out].UseScientific) {
	 OutData[out].pof->setf(std::ios::scientific);
      }

      return false;
   } /* else */

   return true;
}

bool
OutputHandler::IsOpen(const OutputHandler::OutFiles out) const
{
	return OutData[out].IsOpen;
}
   
bool
OutputHandler::Close(const OutputHandler::OutFiles out)
{
   if (!OutData[out].IsOpen) {
#ifdef HAVE_ISOPEN
      ASSERT(!OutData[out].pof->is_open());
#endif /* HAVE_ISOPEN */

      return false;
   }

   /* Chiude lo stream */
   OutData[out].pof->close();
   OutData[out].IsOpen = false;

   return true;
}

bool
OutputHandler::OutputOpen(void) 
{
   return Open(OUTPUT);
}

bool
OutputHandler::RestartOpen(bool openResXSol)
{
	if (!OutData[RESTART].IsOpen) {
	
#ifdef HAVE_ISOPEN
		ASSERT(!OutData[RESTART].pof->is_open());
#endif /* HAVE_ISOPEN */
		char *resExt = NULL;
		int n = nCurrRestartFile > 0 ?
			(int)log10(nCurrRestartFile) + 1 : 1;
		int lenExt = sizeof(".") - 1
			+ n
			+ sizeof(".rst") - 1
			+ sizeof("\0") - 1;
	
		SAFENEWARR(resExt, char, lenExt);
		snprintf(resExt, lenExt, ".%.*d.rst", n, nCurrRestartFile);
		/* Apre lo stream */
	      	OutData[RESTART].pof->open(_sPutExt(resExt));
	
	      	if(!(*OutData[RESTART].pof)) {
		 	std::cerr << "Unable to open file '" << _sPutExt(resExt)
		   		<< '\'' << std::endl;
			throw ErrFile();
		}
		SAFEDELETEARR(resExt);
		OutData[RESTART].IsOpen = true;
      
		/* Setta la formattazione dei campi */
		if(OutData[RESTART].UseDefaultPrecision) {		        
			OutData[RESTART].pof->precision(iCurrPrecision);
		}

		/* Setta la notazione */
		if(OutData[RESTART].UseScientific) {
			OutData[RESTART].pof->setf(std::ios::scientific);
		}
		
		if (openResXSol) {
			char *resXSolExt = NULL;
			int n = nCurrRestartFile > 0 ?
				(int)log10(nCurrRestartFile) + 1 : 1;
			int lenXSolExt = sizeof(".") - 1
				+ n
				+ sizeof(".rst.X") - 1
				+ sizeof("\0") - 1;
		
			SAFENEWARR(resXSolExt, char, lenXSolExt);
			snprintf(resXSolExt, lenXSolExt, ".%.*d.rst.X", n, nCurrRestartFile);
			/* Apre lo stream */
		      	OutData[RESTARTXSOL].pof->open(_sPutExt(resXSolExt));
		      	if(!(*OutData[RESTARTXSOL].pof)) {
			 	std::cerr << "Unable to open file '" << _sPutExt(resExt)
			   		<< '\'' << std::endl;
				throw ErrFile();
			}
			SAFEDELETEARR(resXSolExt);
			OutData[RESTARTXSOL].IsOpen = true;
			/*non occorre settare la precisone e il formato
			perchè il file è binario*/		
		}

		nCurrRestartFile++;
		
 	     	return false;
	}
	return true;
}


bool
OutputHandler::PartitionOpen(void)
{ 
#ifdef HAVE_ISOPEN
  ASSERT(!ofPartition.is_open());
#endif /* HAVE_ISOPEN */
  return Open(PARTITION);
}

bool
OutputHandler::AdamsResOpen(void)
{ 
#ifdef HAVE_ISOPEN
  ASSERT(!ofAdamsRes.is_open());
#endif /* HAVE_ISOPEN */
  return Open(ADAMSRES);
}

bool
OutputHandler::AdamsCmdOpen(void)
{
#ifdef HAVE_ISOPEN
   ASSERT(!ofAdamsCmd.is_open());
#endif /* HAVE_ISOPEN */
   return Open(ADAMSCMD);
}

bool
OutputHandler::LogOpen(void) 
{
#ifdef HAVE_ISOPEN
   ASSERT(!ofLog.is_open());
#endif /* HAVE_ISOPEN */
   return Open(LOG);
}


/* Setta precisione e dimensioni campo */
const int iWidth = 7; /* Caratteri richiesti dalla notazione esponenziale */

void OutputHandler::SetWidth(int iNewWidth)
{
   ASSERT(iNewWidth > iWidth);
   if (iNewWidth > iWidth) {
      iCurrWidth = iNewWidth;
      iCurrPrecision = iCurrWidth-iWidth;
      for (int iCnt = 0; iCnt < LASTFILE; iCnt++) {
	 if(OutData[iCnt].UseDefaultPrecision && *OutData[iCnt].pof) {
	    OutData[iCnt].pof->width(iCurrWidth);
	    OutData[iCnt].pof->precision(iCurrPrecision);
	 }		  
      }
   }
}

void OutputHandler::SetPrecision(int iNewPrecision)
{
   ASSERT(iNewPrecision > 0);
   if (iNewPrecision > 0) {
      iCurrPrecision = iNewPrecision;
      iCurrWidth = iNewPrecision+iWidth;
      for (int iCnt = 0; iCnt < LASTFILE; iCnt++) {
	 if (OutData[iCnt].UseDefaultPrecision && *OutData[iCnt].pof) {
	    OutData[iCnt].pof->width(iCurrWidth);
	    OutData[iCnt].pof->precision(iCurrPrecision);
	 }
      }
   }
}

/* OutputHandler - end */


/* ToBeOutput - begin */

ToBeOutput::ToBeOutput(flag fOut) 
: fOutput(fOut)
{
	NO_OP;
}

ToBeOutput::~ToBeOutput(void)
{
	NO_OP;
}

/* Regular output */
void
ToBeOutput::Output(OutputHandler& OH) const
{
	NO_OP;
}

/* Output of perturbed solution (modes ...) */
void
ToBeOutput::Output(OutputHandler& OH,
		const VectorHandler& X, const VectorHandler& XP) const
{
	NO_OP;
}

/* Output of modes in NASTRAN's pch/f06 format */
void
ToBeOutput::Output_pch(std::ostream &pch) const
{
	NO_OP;
}

void
ToBeOutput::Output_f06(std::ostream &f06, const VectorHandler& X) const
{
	NO_OP;
}

void
ToBeOutput::Output_f06(std::ostream &f06,
		const VectorHandler& Xr, const VectorHandler& Xi) const
{
	NO_OP;
}

flag
ToBeOutput::fToBeOutput(void) const
{
  	return fOutput;
}
   
void
ToBeOutput::SetOutputFlag(flag f)
{
  	fOutput = f;
}

/* ToBeOutput - end */
