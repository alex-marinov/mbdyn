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
   ".out",  /*  0 */
   ".mov",
   ".ele",
   ".abs",
   ".ine",
   ".jnt",  /*  5 */
   ".frc",
   ".act",
   ".rot",
   ".rst",
   ".aer",  /* 10 */
   ".hyd",
   ".prs",
   ".usr",
   ".gen",
   ".par",  /* 15 */
   ".res",
   ".ada",
   ".amd",
   ".rfm",
   ".log",  /* 20 */
   ".air",
   
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
iCurrWidth(iDefaultWidth), iCurrPrecision(iDefaultPrecision)
{
   OutData[OUTPUT].fToUseDefaultPrecision = flag(0);
   OutData[OUTPUT].fToUseScientific = flag(0);
   OutData[OUTPUT].pof = &ofOutput;
   
   OutData[STRNODES].fToUseDefaultPrecision = flag(1);
   OutData[STRNODES].fToUseScientific = flag(1);
   OutData[STRNODES].pof = &ofStrNodes;
   
   OutData[ELECTRIC].fToUseDefaultPrecision = flag(1);
   OutData[ELECTRIC].fToUseScientific = flag(1);
   OutData[ELECTRIC].pof= &ofElectric;
   
   OutData[ABSTRACT].fToUseDefaultPrecision = flag(1);
   OutData[ABSTRACT].fToUseScientific = flag(1);
   OutData[ABSTRACT].pof = &ofAbstract;
   
   OutData[INERTIA].fToUseDefaultPrecision = flag(1);
   OutData[INERTIA].fToUseScientific = flag(1);
   OutData[INERTIA].pof = &ofInertia;
   
   OutData[JOINTS].fToUseDefaultPrecision = flag(1);
   OutData[JOINTS].fToUseScientific = flag(1);
   OutData[JOINTS].pof = &ofJoints;
   
   OutData[FORCES].fToUseDefaultPrecision = flag(1);
   OutData[FORCES].fToUseScientific = flag(1);
   OutData[FORCES].pof = &ofForces;
   
   OutData[BEAMS].fToUseDefaultPrecision = flag(1);
   OutData[BEAMS].fToUseScientific = flag(1);
   OutData[BEAMS].pof = &ofBeams;
   
   OutData[ROTORS].fToUseDefaultPrecision = flag(1);
   OutData[ROTORS].fToUseScientific = flag(1);
   OutData[ROTORS].pof = &ofRotors;
   
   OutData[RESTART].fToUseDefaultPrecision = flag(1);
   OutData[RESTART].fToUseScientific = flag(1);
   OutData[RESTART].pof = &ofRestart;

   OutData[AERODYNAMIC].fToUseDefaultPrecision = flag(1);
   OutData[AERODYNAMIC].fToUseScientific = flag(1);
   OutData[AERODYNAMIC].pof = &ofAerodynamic;
   
   OutData[HYDRAULIC].fToUseDefaultPrecision = flag(1);
   OutData[HYDRAULIC].fToUseScientific = flag(1);
   OutData[HYDRAULIC].pof = &ofHydraulic;

   OutData[PRESNODES].fToUseDefaultPrecision = flag(1);
   OutData[PRESNODES].fToUseScientific = flag(1);
   OutData[PRESNODES].pof = &ofPresNodes;

   OutData[LOADABLE].fToUseDefaultPrecision = flag(1);
   OutData[LOADABLE].fToUseScientific = flag(1);
   OutData[LOADABLE].pof = &ofLoadable;
   
   OutData[GENELS].fToUseDefaultPrecision = flag(1);
   OutData[GENELS].fToUseScientific = flag(1);
   OutData[GENELS].pof = &ofGenels;

   OutData[PARTITION].fToUseDefaultPrecision = flag(1);
   OutData[PARTITION].fToUseScientific = flag(1);
   OutData[PARTITION].pof = &ofPartition;
   
   OutData[ADAMSRES].fToUseDefaultPrecision = flag(1);
   OutData[ADAMSRES].fToUseScientific = flag(1);
   OutData[ADAMSRES].pof = &ofAdamsRes;

   OutData[ADAMSCMD].fToUseDefaultPrecision = flag(1);
   OutData[ADAMSCMD].fToUseScientific = flag(1);
   OutData[ADAMSCMD].pof = &ofAdamsCmd;

   OutData[AEROMODALS].fToUseDefaultPrecision = flag(1);
   OutData[AEROMODALS].fToUseScientific = flag(1);
   OutData[AEROMODALS].pof = &ofAeroModals;
   
   OutData[REFERENCEFRAMES].fToUseDefaultPrecision = flag(1);
   OutData[REFERENCEFRAMES].fToUseScientific = flag(1);
   OutData[REFERENCEFRAMES].pof = &ofReferenceFrames;
   
   OutData[LOG].fToUseDefaultPrecision = flag(0);
   OutData[LOG].fToUseScientific = flag(0);
   OutData[LOG].pof = &ofLog;
   
   OutData[AIRPROPS].fToUseDefaultPrecision = flag(1);
   OutData[AIRPROPS].fToUseScientific = flag(1);
   OutData[AIRPROPS].pof = &ofAirProps;
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
iCurrWidth(iDefaultWidth), iCurrPrecision(iDefaultPrecision)
{
   OutData[OUTPUT].fToUseDefaultPrecision = flag(0);
   OutData[OUTPUT].fToUseScientific = flag(0);
   OutData[OUTPUT].pof = &ofOutput;
   ofOutput.width(iCurrWidth);
   ofOutput.precision(iCurrPrecision);
   
   OutData[STRNODES].fToUseDefaultPrecision = flag(1);
   OutData[STRNODES].fToUseScientific = flag(1);
   OutData[STRNODES].pof = &ofStrNodes;
   
   OutData[ELECTRIC].fToUseDefaultPrecision = flag(1);
   OutData[ELECTRIC].fToUseScientific = flag(1);
   OutData[ELECTRIC].pof = &ofElectric;
   
   OutData[ABSTRACT].fToUseDefaultPrecision = flag(1);
   OutData[ABSTRACT].fToUseScientific = flag(1);
   OutData[ABSTRACT].pof = &ofAbstract;
   
   OutData[INERTIA].fToUseDefaultPrecision = flag(1);
   OutData[INERTIA].fToUseScientific = flag(1);
   OutData[INERTIA].pof = &ofInertia;
   
   OutData[JOINTS].fToUseDefaultPrecision = flag(1);
   OutData[JOINTS].fToUseScientific = flag(1);
   OutData[JOINTS].pof = &ofJoints;
   
   OutData[FORCES].fToUseDefaultPrecision = flag(1);
   OutData[FORCES].fToUseScientific = flag(1);
   OutData[FORCES].pof = &ofForces;
   
   OutData[BEAMS].fToUseDefaultPrecision = flag(1);
   OutData[BEAMS].fToUseScientific = flag(1);
   OutData[BEAMS].pof = &ofBeams;
   
   OutData[ROTORS].fToUseDefaultPrecision = flag(1);
   OutData[ROTORS].fToUseScientific = flag(1);
   OutData[ROTORS].pof = &ofRotors;
   
   OutData[RESTART].fToUseDefaultPrecision = flag(1);
   OutData[RESTART].fToUseScientific = flag(1);
   OutData[RESTART].pof = &ofRestart;
   
   OutData[AERODYNAMIC].fToUseDefaultPrecision = flag(1);
   OutData[AERODYNAMIC].fToUseScientific = flag(1);
   OutData[AERODYNAMIC].pof = &ofAerodynamic;
   
   OutData[HYDRAULIC].fToUseDefaultPrecision = flag(1);
   OutData[HYDRAULIC].fToUseScientific = flag(1);
   OutData[HYDRAULIC].pof = &ofHydraulic;

   OutData[PRESNODES].fToUseDefaultPrecision = flag(1);
   OutData[PRESNODES].fToUseScientific = flag(1);
   OutData[PRESNODES].pof = &ofPresNodes;

   OutData[LOADABLE].fToUseDefaultPrecision = flag(1);
   OutData[LOADABLE].fToUseScientific = flag(1);
   OutData[LOADABLE].pof = &ofLoadable;
   
   OutData[GENELS].fToUseDefaultPrecision = flag(1);
   OutData[GENELS].fToUseScientific = flag(1);
   OutData[GENELS].pof = &ofGenels;

   OutData[PARTITION].fToUseDefaultPrecision = flag(1);
   OutData[PARTITION].fToUseScientific = flag(1);
   OutData[PARTITION].pof = &ofPartition;
   
   OutData[ADAMSRES].fToUseDefaultPrecision = flag(1);
   OutData[ADAMSRES].fToUseScientific = flag(1);
   OutData[ADAMSRES].pof = &ofAdamsRes;

   OutData[ADAMSCMD].fToUseDefaultPrecision = flag(1);
   OutData[ADAMSCMD].fToUseScientific = flag(1);
   OutData[ADAMSCMD].pof = &ofAdamsCmd;
   
   OutData[AEROMODALS].fToUseDefaultPrecision = flag(1);
   OutData[AEROMODALS].fToUseScientific = flag(1);
   OutData[AEROMODALS].pof = &ofAeroModals;
   
   OutData[REFERENCEFRAMES].fToUseDefaultPrecision = flag(1);
   OutData[REFERENCEFRAMES].fToUseScientific = flag(1);
   OutData[REFERENCEFRAMES].pof = &ofReferenceFrames;

   OutData[LOG].fToUseDefaultPrecision = flag(0);
   OutData[LOG].fToUseScientific = flag(0);
   OutData[LOG].pof = &ofLog;
   ofLog.width(iCurrWidth);
   ofLog.precision(iCurrPrecision);
   
   OutData[AIRPROPS].fToUseDefaultPrecision = flag(1);
   OutData[AIRPROPS].fToUseScientific = flag(1);
   OutData[AIRPROPS].pof = &ofAirProps;
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
#if HAVE_ISOPEN
      if (OutData[iCnt].pof->is_open()) {
#endif /* HAVE_ISOPEN */
	 OutData[iCnt].pof->close();
#if HAVE_ISOPEN
      }
#endif /* HAVE_ISOPEN */
   }
}
   
   
/* Aggiungere qui le funzioni che aprono i singoli stream */
int OutputHandler::Open(OutFiles out)
{
#if HAVE_ISOPEN
   if(!OutData[out].pof->is_open()) {
#endif /* HAVE_ISOPEN */
      /* Apre lo stream */
      OutData[out].pof->open(_sPutExt((char*)(psExt[out])));

      if(!(*OutData[out].pof)) {
	 std::cerr << "Unable to open file <" << _sPutExt((char*)(psExt[out]))
	   << '>' << std::endl;
	 THROW(ErrFile());
      }
      
      /* Setta la formattazione dei campi */
      if(OutData[out].fToUseDefaultPrecision) {		        
	 OutData[out].pof->precision(iCurrPrecision);
      }
      
      /* Setta la notazione */
      if(OutData[out].fToUseScientific) {
	 OutData[out].pof->setf(std::ios::scientific);
      }

      return 0;
#if HAVE_ISOPEN
   }
   return 1;
#endif /* HAVE_ISOPEN */
}
   
   
int OutputHandler::OutputOpen(void) 
{
#if HAVE_ISOPEN
   ASSERT(!ofOutput.is_open());
#endif /* HAVE_ISOPEN */
   return Open(OUTPUT);
}

#if 0
int OutputHandler::StrNodesOpen(void) 
{ 
#if HAVE_ISOPEN
   ASSERT(!ofStrNodes.is_open());
#endif /* HAVE_ISOPEN */
   return Open(STRNODES);
}


int OutputHandler::ElectricOpen(void) 
{
#if HAVE_ISOPEN
   ASSERT(!ofElectric.is_open());
#endif /* HAVE_ISOPEN */
   return Open(ELECTRIC);
}


int OutputHandler::AbstractOpen(void) 
{ 
#if HAVE_ISOPEN
   ASSERT(!ofAbstract.is_open());
#endif /* HAVE_ISOPEN */
   return Open(ABSTRACT);
}


int OutputHandler::InertiaOpen(void) 
{ 
#if HAVE_ISOPEN
   ASSERT(!ofInertia.is_open());
#endif /* HAVE_ISOPEN */
   return Open(INERTIA);
}


int OutputHandler::JointsOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofJoints.is_open());
#endif /* HAVE_ISOPEN */
   return Open(JOINTS);
}


int OutputHandler::ForcesOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofForces.is_open());
#endif /* HAVE_ISOPEN */
   return Open(FORCES);
}


int OutputHandler::BeamsOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofBeams.is_open());
#endif /* HAVE_ISOPEN */
   return Open(BEAMS);
}


int OutputHandler::RotorsOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofRotors.is_open());
#endif /* HAVE_ISOPEN */
   return Open(ROTORS);
}
#endif /* 0 */

int OutputHandler::RestartOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofRestart.is_open());
#endif /* HAVE_ISOPEN */
   return Open(RESTART);
}

#if 0
int OutputHandler::AerodynamicOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofAerodynamic.is_open());
#endif /* HAVE_ISOPEN */
   return Open(AERODYNAMIC);
}


int OutputHandler::HydraulicOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofHydraulic.is_open());
#endif /* HAVE_ISOPEN */
   return Open(HYDRAULIC);
}


int OutputHandler::PresNodesOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofPresNodes.is_open());
#endif /* HAVE_ISOPEN */
   return Open(PRESNODES);
}


int OutputHandler::LoadableOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofLoadable.is_open());
#endif /* HAVE_ISOPEN */
   return Open(LOADABLE);
}


int OutputHandler::GenelsOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofGenels.is_open());
#endif /* HAVE_ISOPEN */
   return Open(GENELS);
}
#endif /* 0 */

int OutputHandler::PartitionOpen(void)
{ 
#if HAVE_ISOPEN
  ASSERT(!ofPartition.is_open());
#endif /* HAVE_ISOPEN */
  return Open(PARTITION);
}

int OutputHandler::AdamsResOpen(void)
{ 
#if HAVE_ISOPEN
  ASSERT(!ofAdamsRes.is_open());
#endif /* HAVE_ISOPEN */
  return Open(ADAMSRES);
}

int OutputHandler::AdamsCmdOpen(void)
{
#if HAVE_ISOPEN
   ASSERT(!ofAdamsCmd.is_open());
#endif /* HAVE_ISOPEN */
   return Open(ADAMSCMD);
}

int OutputHandler::LogOpen(void) 
{
#if HAVE_ISOPEN
   ASSERT(!ofLog.is_open());
#endif /* HAVE_ISOPEN */
   return Open(LOG);
}


#if 0
int OutputHandler::AirPropsOpen(void)
{ 
#if HAVE_ISOPEN
   ASSERT(!ofAirProps.is_open());
#endif /* HAVE_ISOPEN */
   return Open(AIRPROPS);
}
#endif /* 0 */

/* Setta precisione e dimensioni campo */
const int iWidth = 7; /* Caratteri richiesti dalla notazione esponenziale */

void OutputHandler::SetWidth(int iNewWidth)
{
   ASSERT(iNewWidth > iWidth);
   if (iNewWidth > iWidth) {
      iCurrWidth = iNewWidth;
      iCurrPrecision = iCurrWidth-iWidth;
      for (int iCnt = 0; iCnt < LASTFILE; iCnt++) {
	 if(OutData[iCnt].fToUseDefaultPrecision && *OutData[iCnt].pof) {
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
	 if (OutData[iCnt].fToUseDefaultPrecision && *OutData[iCnt].pof) {
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
