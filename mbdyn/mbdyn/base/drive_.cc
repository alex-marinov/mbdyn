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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <mymath.h>
#include <dataman.h>
#include <drive_.h>
#include <dofdrive.h>
#include <filedrv.h>

#ifdef USE_MPI
#include <mpi++.h>
#endif /* USE_MPI */

/* StringDriveCaller - begin */

StringDriveCaller::StringDriveCaller(const DriveHandler* pDH,
				     const char* const sTmpStr)
: DriveCaller(pDH), sEvalStr(NULL), iEvalStrLen(0) 
{
   ASSERT(sTmpStr != NULL);
   iEvalStrLen = strlen(sTmpStr);
   SAFESTRDUP(sEvalStr, sTmpStr, DMmm);
}


StringDriveCaller::~StringDriveCaller(void)
{
   ASSERT(sEvalStr != NULL);
   if (sEvalStr != NULL) {
      SAFEDELETEARR(sEvalStr, DMmm);
   }
}


/* Copia */
DriveCaller* StringDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  StringDriveCaller, 
			  StringDriveCaller(pDrvHdl, sEvalStr), 
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& StringDriveCaller::Restart(ostream& out) const
{
   return out << "string, \"" << sEvalStr << "\"";
}

/* StringDriveCaller - end */


/* TimeDriveCaller - begin */

TimeDriveCaller::TimeDriveCaller(const DriveHandler* pDH)
: DriveCaller(pDH)
{
   NO_OP;
}


TimeDriveCaller::~TimeDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* TimeDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDrvHdl), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& TimeDriveCaller::Restart(ostream& out) const
{
   return out << "time";
}
 
/* TimeDriveCaller - end */


/* ConstDriveCaller - begin */

ConstDriveCaller::ConstDriveCaller(const DriveHandler* pDH, doublereal d)
: DriveCaller(pDH), dConst(d)
{
   NO_OP;
}


ConstDriveCaller::~ConstDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* ConstDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(pDrvHdl, dConst), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& ConstDriveCaller::Restart(ostream& out) const
{      
   return out << " const, " << dConst;
}
 
/* ConstDriveCaller - end */


/* LinearDriveCaller - begin */

LinearDriveCaller::LinearDriveCaller(const DriveHandler* pDH, 
				     doublereal d0, doublereal d1)
: DriveCaller(pDH), dC0(d0), dC1(d1)
{
   NO_OP;
}


LinearDriveCaller::~LinearDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* LinearDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, LinearDriveCaller, LinearDriveCaller(pDrvHdl, dC0, dC1), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& LinearDriveCaller::Restart(ostream& out) const
{
   return out << " linear, " << dC0 << ", " << dC1;
}
 
/* LinearDriveCaller - end */


/* ParabolicDriveCaller - begin */

ParabolicDriveCaller::ParabolicDriveCaller(const DriveHandler* pDH,
					   doublereal d0, 
					   doublereal d1,
					   doublereal d2)
: DriveCaller(pDH), dC0(d0), dC1(d1), dC2(d2)
{
   NO_OP;
}


ParabolicDriveCaller::~ParabolicDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* ParabolicDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, ParabolicDriveCaller, ParabolicDriveCaller(pDrvHdl, dC0, dC1, dC2), DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
ostream& ParabolicDriveCaller::Restart(ostream& out) const
{
   return out <<  " parabolic, " << dC0 << ", " << dC1 << ", " << dC2;
}
 
/* ParabolicDriveCaller - end */


/* CubicDriveCaller - begin */

CubicDriveCaller::CubicDriveCaller(const DriveHandler* pDH, 
				   doublereal d0, doublereal d1, 
				   doublereal d2, doublereal d3)
: DriveCaller(pDH), dC0(d0), dC1(d1), dC2(d2), dC3(d3) 
{
   NO_OP;
}


/* Copia */
DriveCaller* CubicDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, CubicDriveCaller, CubicDriveCaller(pDrvHdl, dC0, dC1, dC2, dC3), DMmm);
   
   return pDC;
}


CubicDriveCaller::~CubicDriveCaller(void)
{
   NO_OP;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& CubicDriveCaller::Restart(ostream& out) const
{
   return out << " cubic, " 
     << dC0 << ", " << dC1 << ", " << dC2 << ", " << dC3;
}

/* CubicDriveCaller - end */


/* StepDriveCaller - begin */

StepDriveCaller::StepDriveCaller(const DriveHandler* pDH, 
				 doublereal d1, doublereal d2, doublereal d3) 
: DriveCaller(pDH),
dStepTime(d1), dStepValue(d2), dInitialValue(d3)
{
   NO_OP;
}


StepDriveCaller::~StepDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* StepDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  StepDriveCaller, 
			  StepDriveCaller(pDrvHdl, 
					  dStepTime, 
					  dStepValue, 
					  dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& StepDriveCaller::Restart(ostream& out) const
{
   return out
     << " step, " 
     << dStepTime << ", " 
     << dStepValue << ", " 
     << dInitialValue;
}

/* StepDriveCaller - end */


/* DoubleStepDriveCaller - begin */

DoubleStepDriveCaller::DoubleStepDriveCaller(const DriveHandler* pDH, 
					     doublereal d1, doublereal d2,
					     doublereal d3, doublereal d4)
: DriveCaller(pDH), dStepTime(d1), dStepValue(d2),
dEndStepTime(d3), dInitialValue(d4) {
   NO_OP;
}


DoubleStepDriveCaller::~DoubleStepDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* DoubleStepDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  DoubleStepDriveCaller, 
			  DoubleStepDriveCaller(pDrvHdl, 
						dStepTime, 
						dStepValue,
						dEndStepTime,
						dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& DoubleStepDriveCaller::Restart(ostream& out) const
{
   return out 
     << " double step, "
     << dStepTime << ", " 
     << dStepValue << ", "
     << dEndStepTime << ", "
     << dInitialValue;
}

/* DoubleStepDriveCaller - end */


/* RampDriveCaller - begin */

RampDriveCaller::RampDriveCaller(const DriveHandler* pDH, 
				 doublereal d1, doublereal d2, 
				 doublereal d3, doublereal d4)
: DriveCaller(pDH), dSlope(d1), dStartTime(d2), dEndTime(d3), dInitialValue(d4)
{
   NO_OP;
}


RampDriveCaller::~RampDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* RampDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  RampDriveCaller, 
			  RampDriveCaller(pDrvHdl, 
					  dSlope,
					  dStartTime,
					  dEndTime,
					  dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& RampDriveCaller::Restart(ostream& out) const
{
   return out
     << "ramp, "
     << dSlope << ", " 
     << dStartTime << ", "
     << dEndTime << ", "
     << dInitialValue;
}

/* RampDriveCaller - end */


/* DoubleRampDriveCaller - begin */

DoubleRampDriveCaller::DoubleRampDriveCaller(const DriveHandler* pDH,
					     doublereal d1, doublereal d2,
					     doublereal d3, 
					     doublereal d4, doublereal d5, 
					     doublereal d6,
					     doublereal d7)
: DriveCaller(pDH),
dAscendingSlope(d1), dAscendingStartTime(d2), 
dAscendingEndTime(d3),
dDescendingSlope(d4), dDescendingStartTime(d5), 
dDescendingEndTime(d6), 
dInitialValue(d7)
{
   NO_OP;
}


DoubleRampDriveCaller::~DoubleRampDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* DoubleRampDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  DoubleRampDriveCaller, 
			  DoubleRampDriveCaller(pDrvHdl, 
						dAscendingSlope,
						dAscendingStartTime,
						dAscendingEndTime,
						dDescendingSlope,
						dDescendingStartTime,
						dDescendingEndTime, 
						dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& DoubleRampDriveCaller::Restart(ostream& out) const
{
   return out
     << " double ramp, "
     << dAscendingSlope << ", " 
     << dAscendingStartTime << ", "
     << dAscendingEndTime << ", "
     << dDescendingSlope << ", "
     << dDescendingStartTime << ", "
     << dDescendingEndTime << ", "
     << dInitialValue;
}

/* DoubleRampDriveCaller - end */


/* SineDriveCaller - begin */

SineDriveCaller::SineDriveCaller(const DriveHandler* pDH,
				 doublereal d1, doublereal d2, doublereal d3, 
				 integer iNumCyc, doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), dOmega(d2), dAmplitude(d3),
iNumCycles(iNumCyc), dInitialValue(d4), fNeverEnd(0)
{
   /* Onde di seno che partono da zero ed arrivano a0 */
   if (iNumCycles > 0) {
      dEndTime = dStartTime
	+2.*M_PI/dOmega*(doublereal(iNumCycles)-.5);
      dFinalValue = dInitialValue;		
      /* Onde di seno che continuano all'infinito */
   } else if(iNumCycles == 0) {	     
      dEndTime = 0.;
      fNeverEnd = flag(1);	
      
      /* Onde di seno che partono da 0 ed arrivano ad 1 
       * con tangente orizzontale */
   } else {
      dEndTime = dStartTime
	+2.*M_PI/dOmega*(doublereal(-iNumCycles)-.75);
      dFinalValue = dInitialValue+dAmplitude;
   }
}


SineDriveCaller::~SineDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* SineDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  SineDriveCaller, 
			  SineDriveCaller(pDrvHdl, 
					  dStartTime,
					  dOmega,
					  dAmplitude,
					  iNumCycles,
					  dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& SineDriveCaller::Restart(ostream& out) const
{
   return out 
     << " sine, "
     << dStartTime << ", "
     << dOmega << ", "
     << dAmplitude << ", "
     << iNumCycles << ", "
     << dInitialValue;
}

/* SineDriveCaller - end */


/* CosineDriveCaller - begin */

CosineDriveCaller::CosineDriveCaller(const DriveHandler* pDH,
				     doublereal d1, doublereal d2, 
				     doublereal d3, 
				     integer iNumCyc, doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), dOmega(d2), dAmplitude(d3),
iNumCycles(iNumCyc), dInitialValue(d4), fNeverEnd(0)
{
   /* Onde di coseno che partono da zero ed arrivano a 0 */
   if (iNumCycles > 0) {
      dEndTime = dStartTime
	+2.*M_PI/dOmega*doublereal(iNumCycles);
      dFinalValue = dInitialValue;			
      /* Onde di coseno che continuano all'infinito */
   } else if(iNumCycles == 0) {	     
      dEndTime = 0.;
      fNeverEnd = flag(1);		
      /* Onde di coseno che partono da 0 ed arrivano ad 1 
       * con tangente orizzontale */
   } else {
      dEndTime = dStartTime
	+2.*M_PI/dOmega*(doublereal(-iNumCycles)-.5);
      dFinalValue = dInitialValue+2.*dAmplitude;
   }
}


CosineDriveCaller::~CosineDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* CosineDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  CosineDriveCaller, 
			  CosineDriveCaller(pDrvHdl,
					    dStartTime,
					    dOmega,
					    dAmplitude,
					    iNumCycles,
					    dInitialValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& CosineDriveCaller::Restart(ostream& out) const
{
   return out
     << " cosine, "
     << dStartTime << ", "
     << dOmega << ", "
     << dAmplitude << ", "
     << iNumCycles << ", "
     << dInitialValue;
}

/* CosineDriveCaller - end */


/* FreqSweepDriveCaller - begin */

FreqSweepDriveCaller::FreqSweepDriveCaller(const DriveHandler* pDH,
					   doublereal d1, 
					   const DriveCaller* pOmega, 
					   const DriveCaller* pAmplitude,
					   doublereal d2,
					   doublereal d3,
					   doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), pOmega(pOmega), pAmplitude(pAmplitude),
dInitialValue(d2), dEndTime(d3), dFinalValue(d4), 
fNeverEnd(0)
{
   ASSERT(pOmega != NULL);
   ASSERT(pAmplitude != NULL);
   
   if (dEndTime <= dStartTime) {
      fNeverEnd = flag(1);
   }
}


FreqSweepDriveCaller::~FreqSweepDriveCaller(void)
{
   NO_OP;
}


/* Copia */
DriveCaller* FreqSweepDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  FreqSweepDriveCaller, 
			  FreqSweepDriveCaller(pDrvHdl, 
					       dStartTime, 
					       pOmega->pCopy(),
					       pAmplitude->pCopy(),
					       dInitialValue,
					       dEndTime, 
					       dFinalValue),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& FreqSweepDriveCaller::Restart(ostream& out) const
{
   return out 
     << " frequency sweep, "
     << dStartTime << ", ",
     pOmega->Restart(out) << ", ",
     pAmplitude->Restart(out) << ", "     
     << dInitialValue << ", "
     << dEndTime << ", "
     << dFinalValue;
}

/* FreqSweepDriveCaller - end */


/* ExpDriveCaller - begin */

ExpDriveCaller::ExpDriveCaller(const DriveHandler* pDH, 
			       doublereal dA, doublereal dT, 
			       doublereal dS, doublereal dI)
: DriveCaller(pDH), 
dAmplitude(dA), dTimeConst(dT), dStartTime(dS), dInitialValue(dI)
{
   NO_OP;
}


ExpDriveCaller::~ExpDriveCaller(void)
{
   NO_OP; 
}
 

/* Copia */
DriveCaller* ExpDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  ExpDriveCaller, 
			  ExpDriveCaller(pDrvHdl, 
					 dAmplitude,
					 dTimeConst,
					 dStartTime,
					 dInitialValue), 
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& ExpDriveCaller::Restart(ostream& out) const
{
   return out
     << " exponential, " 
     << dAmplitude << ", "
     << dTimeConst << ", "
     << dStartTime << ", "
     << dInitialValue;
}

/* ExpDriveCaller - end */


/* RandDriveCaller - begin */

RandDriveCaller::RandDriveCaller(const DriveHandler* pDH, 
				 doublereal dA, doublereal dR,
				 doublereal dS, doublereal dE, integer iS)
: DriveCaller(pDH), dAmplitude(dA), dRefVal(dR),
dStartTime(dS), dEndTime(dE), iSteps(iS)
{
   iBase = rand();
   iRandDriveNumber = ((DriveHandler*)pDrvHdl)->iRandInit(iSteps);
}


RandDriveCaller::~RandDriveCaller(void)
{
   NO_OP; 
}


/* Copia */
DriveCaller* RandDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  RandDriveCaller, 
			  RandDriveCaller(pDrvHdl,
					  dAmplitude,
					  dRefVal,
					  dStartTime, 
					  dEndTime,
					  iSteps),
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& RandDriveCaller::Restart(ostream& out) const
{
   out
     << " random, " 
     << dAmplitude << ", "
     << dRefVal << ", "
     << dStartTime << ", "
     << dEndTime;
   if(iSteps > 1) {
      out << ", steps, " << iSteps;
   }
   return out;
}

/* RandDriveCaller - end */


/* PiecewiseLinearDriveCaller - begin */

PiecewiseLinearDriveCaller::PiecewiseLinearDriveCaller(const DriveHandler* pDH,
		unsigned int i,
		doublereal *p)
: DriveCaller(pDH), iNumPoints(i), pPoints(p), pVals(p+i)
{
	ASSERT(i >= 2);
	ASSERT(p != NULL);
}

PiecewiseLinearDriveCaller::~PiecewiseLinearDriveCaller(void)
{
	if (pPoints != NULL) {
		SAFEDELETEARR(pPoints, DMmm);
	}
}

/* Copia */
DriveCaller* PiecewiseLinearDriveCaller::pCopy(void) const
{
	doublereal *p = NULL;
	SAFENEWARR(p, doublereal, 2*iNumPoints, DMmm);
	for (unsigned int i = 0; i < 2*iNumPoints; i++) {
		p[i] = pPoints[i];
	}
	
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC,
			PiecewiseLinearDriveCaller,
			PiecewiseLinearDriveCaller(pDrvHdl, iNumPoints, p),
			DMmm);

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
ostream& PiecewiseLinearDriveCaller::Restart(ostream& out) const
{
	out << "piecewise linear, " << iNumPoints;

	for (unsigned int i = 0; i < iNumPoints; i++) {
		out << ", " << pPoints[i] << ", " << pVals[i];
	}

	return out;
}

/* PiecewiseLinearDriveCaller - end */


/* DriveArrayCaller - begin */

DriveArrayCaller::DriveArrayCaller(const DriveHandler* pDH, 
				   unsigned short int i, 
				   const DriveCaller** ppDC)
: DriveCaller(pDH), iNumDrivers(i), ppDriveCallers(ppDC)
{
#ifdef DEBUG
   ASSERT(iNumDrivers > 0);
   ASSERT(ppDriveCallers != NULL);
   for (int i = 0; i < iNumDrivers; i++) {
      ASSERT(ppDC[i] != NULL);
   }
#endif /* DEBUG */
}


DriveArrayCaller::~DriveArrayCaller(void)
{
   ASSERT(ppDriveCallers != NULL);
   ASSERT(iNumDrivers > 0);
   for (int i = 0; i < iNumDrivers; i++) {	
      SAFEDELETE(ppDriveCallers[i], DMmm);
   }
   SAFEDELETEARR(ppDriveCallers, DMmm);
}


/* Copia */
DriveCaller* DriveArrayCaller::pCopy(void) const
{
   ASSERT(ppDriveCallers != NULL);
   ASSERT(iNumDrivers > 0);
   DriveCaller** ppDC = NULL;
   SAFENEWARR(ppDC, DriveCaller*, iNumDrivers, DMmm);
   for (int i = 0; i < iNumDrivers; i++) {
      ppDC[i] = ppDriveCallers[i]->pCopy();
   }
   
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  DriveArrayCaller, 
			  DriveArrayCaller(pDrvHdl, 
					   iNumDrivers,
					   (const DriveCaller**)ppDC), 
			  DMmm);
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
ostream& DriveArrayCaller::Restart(ostream& out) const
{
   out << " array, " << iNumDrivers;
   for (int i = 0; i < iNumDrivers; i++) {
      ASSERT(ppDriveCallers[i] != NULL);
      out << ", ", ppDriveCallers[i]->Restart(out);
   }
   return out;
}

/* DriveArrayCaller - end */


/* Legge i dati dei drivers built-in */

DriveCaller* ReadDriveData(const DataManager* pDM,
			   MBDynParser& HP, 
			   const DriveHandler* pDrvHdl)
{
   DEBUGCOUTFNAME("ReadDriveData()");
   
   const char* sKeyWords[] = {
      "time",
	"null",
	"one",
	"const",
	"linear",
	"parabolic",
	"cubic",
	"function",
	"step",
	"double" "step",
	"ramp",
	"double" "ramp",
	"sine",
	"cosine",
	"frequency" "sweep",
	"exponential",
	"random",
	"piecewise" "linear",
	"file", 
	"string",
	"dof",
	"array"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,

	TIME = 0,
	NULLDRIVE,
	ONEDRIVE,
	CONST,
	LINEAR,
	PARABOLIC,
	CUBIC,
	FUNCTION,
	STEP,
	DOUBLESTEP,
	RAMP,
	DOUBLERAMP,
	SINE,
	COSINE,
	FREQUENCYSWEEP,
	EXPONENTIAL,
	RANDOM,
	PIECEWISELINEAR,
	FILEDRIVE,      
	STRING,
	DOF,
	ARRAY,

	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);   
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord;
   if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
      CurrKeyWord = CONST;
   } else if (CurrKeyWord == FUNCTION) {
      CurrKeyWord = KeyWords(HP.GetWord());
   }   
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      cout << "drive type: " << sKeyWords[CurrKeyWord] << endl;
   }   
#endif   

   DriveCaller* pDC = NULL;
   
   switch (CurrKeyWord) {

      /* time */
    case TIME: {
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      TimeDriveCaller,
			      TimeDriveCaller(pDrvHdl),
			      DMmm);
     break;  
    }
      
      /* driver nullo */
    case NULLDRIVE: {
       /* allocazione e creazione */
        SAFENEWWITHCONSTRUCTOR(pDC,
			      NullDriveCaller,
			      NullDriveCaller(pDrvHdl),
			      DMmm);
      
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* driver unitario*/
    case ONEDRIVE: {
       /* allocazione e creazione */
        SAFENEWWITHCONSTRUCTOR(pDC,
			      OneDriveCaller,
			      OneDriveCaller(pDrvHdl),
			      DMmm);
      
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* driver costante */
    case CONST: {
       /* lettura dei dati specifici */
       doublereal dConst = HP.GetReal();
       DEBUGCOUT("Const value: " << dConst << endl);

       /* allocazione e creazione */
       if (dConst == 0.) {
          SAFENEWWITHCONSTRUCTOR(pDC,
			         NullDriveCaller,
			         NullDriveCaller(pDrvHdl),
			         DMmm);
       } else if (dConst == 1.) {
          SAFENEWWITHCONSTRUCTOR(pDC,
			         OneDriveCaller,
			         OneDriveCaller(pDrvHdl),
			         DMmm);
      
       } else {
          SAFENEWWITHCONSTRUCTOR(pDC,
	   		         ConstDriveCaller,
			         ConstDriveCaller(pDrvHdl, dConst), 
			         DMmm);
       }
       
       break;
    }
      
      /* retta */
    case LINEAR: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      LinearDriveCaller,
			      LinearDriveCaller(pDrvHdl, dC0, dC1), 
			      DMmm);
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* parabola */
    case PARABOLIC: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << endl);
       
       doublereal dC2 = HP.GetReal();
       DEBUGCOUT("Parabolic slope: " << dC2 << endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      ParabolicDriveCaller,
			      ParabolicDriveCaller(pDrvHdl, dC0, dC1, dC2), 
			      DMmm);
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* cubica */
    case CUBIC: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << endl);
       
       doublereal dC2 = HP.GetReal();
       DEBUGCOUT("Parabolic slope: " << dC2 << endl);
       
       doublereal dC3 = HP.GetReal();
       DEBUGCOUT("Cubic slope: " << dC3 << endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      CubicDriveCaller,
			      CubicDriveCaller(pDrvHdl, dC0, dC1, dC2, dC3), 
			      DMmm);
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* Scalino semplice */
    case STEP: {	
       doublereal dStepTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dStepTime << endl);
       
       doublereal dStepValue = HP.GetReal(1.);
       DEBUGCOUT("Step Value: " << dStepValue << endl);
       
       doublereal dInitialValue = HP.GetReal();
       DEBUGCOUT("InitialValue: " << dInitialValue << endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      StepDriveCaller,
			      StepDriveCaller(pDrvHdl, dStepTime,
					      dStepValue, dInitialValue), 
			      DMmm);
       break;
    }
      
      /* Scalino doppio */
    case DOUBLESTEP: {	
       doublereal dStepTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dStepTime << endl);
       
       doublereal dEndStepTime = HP.GetReal();
       DEBUGCOUT("Final time: " << dEndStepTime << endl);
       
       if(dEndStepTime <= dStepTime) {	      
	  cerr << "Warning at line " 
	    << HP.GetLineData() 
	    << ": final time " << dEndStepTime
	    << " is less than or equal to initial time " << dStepTime 
	    << " in double step func drive" << endl;
       }
       
       doublereal dStepValue = HP.GetReal(1.);
       DEBUGCOUT("Step Value: " << dStepValue << endl);
       
       doublereal dInitialValue = HP.GetReal();     
       DEBUGCOUT("InitialValue: " << dInitialValue << endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DoubleStepDriveCaller,
			      DoubleStepDriveCaller(pDrvHdl, 
						    dStepTime, 
						    dStepValue,
						    dEndStepTime,
						    dInitialValue), 
			      DMmm);
       break;
    }
      
      /* Rampa saturata */
    case RAMP: {	
       doublereal dSlope = HP.GetReal(1.);
       DEBUGCOUT("Slope Value: " << dSlope << endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << endl);
       
       doublereal dFinalTime = HP.GetReal();
       DEBUGCOUT("Final time: " << dFinalTime << endl);
       
       if(dFinalTime <= dInitialTime) {	      
	  cerr << "Warning at line "
	    << HP.GetLineData() 
	    << ": final time " << dFinalTime
	    << " is less than or equal to initial time " << dInitialTime 
	    << " in ramp func drive" << endl;
       }	   
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC, 
			      RampDriveCaller,
			      RampDriveCaller(pDrvHdl, 
					      dSlope, dInitialTime,
					      dFinalTime, dInitialValue), 
			      DMmm);
       break;
    }
      
      /* Rampa doppia */
    case DOUBLERAMP: {	
       doublereal dAscendingSlope = HP.GetReal(1.);
       DEBUGCOUT("Ascending Slope Value: " << dAscendingSlope << endl);
       
       doublereal dAscendingInitialTime = HP.GetReal();
       DEBUGCOUT("Ascending Initial time: " << dAscendingInitialTime << endl);
       
       doublereal dAscendingFinalTime = HP.GetReal();
       DEBUGCOUT("Ascending Final time: " << dAscendingFinalTime << endl);
       
       if (dAscendingFinalTime <= dAscendingInitialTime) {	      
	  cerr << "Warning at line " 
	    << HP.GetLineData() << ": ascending final time " 
	    << dAscendingFinalTime
	    << " is less than or equal to ascending initial time " 
	    << dAscendingInitialTime 
	    << " in double ramp func drive" << endl;
	  THROW(ErrGeneric());
       }
       
       doublereal dDescendingSlope = HP.GetReal(-1.);
       DEBUGCOUT("Descending Slope Value: " << dDescendingSlope << endl);
       
       doublereal dDescendingInitialTime = HP.GetReal();
       DEBUGCOUT("Descending Initial time: " << dDescendingInitialTime << endl);
       
       if (dDescendingInitialTime < dAscendingFinalTime) {
	  cerr << "Warning at line " 
	    << HP.GetLineData() << ": descending initial time " 
	    << dDescendingInitialTime
	    << " is less than ascending final time " 
	    << dAscendingFinalTime 
	    << " in double ramp func drive" << endl;
	  THROW(ErrGeneric());
       }	   
       
       doublereal dDescendingFinalTime = HP.GetReal();
       DEBUGCOUT("Descending Final time: " << dDescendingFinalTime << endl);
       
       if (dDescendingFinalTime <= dDescendingInitialTime) {	      
	  cerr << "Warning at line " 
	    << HP.GetLineData() << ": descending final time " 
	    << dDescendingFinalTime
	    << " is less than descending initial time " 
	    << dDescendingInitialTime 
	    << " in double ramp func drive" << endl;
	  THROW(ErrGeneric());
       }	   
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DoubleRampDriveCaller,
			      DoubleRampDriveCaller(pDrvHdl,
						    dAscendingSlope,
						    dAscendingInitialTime,
						    dAscendingFinalTime,
						    dDescendingSlope,
						    dDescendingInitialTime,
						    dDescendingFinalTime,
						    dInitialValue), 
			      DMmm);
       break;
    }
      
      /* Seno e coseno (limitati, illimitati, saturati) */
    case SINE:
    case COSINE: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << endl);
       
       doublereal dOmega = HP.GetReal(1.);
       DEBUGCOUT("Omega: " << dOmega << endl);
       
       doublereal dAmplitude = HP.GetReal();
       DEBUGCOUT("Amplitude: " << dAmplitude << endl);
       
       integer iNumCycles = HP.GetInt();
       DEBUGCOUT("Number of cycles: " << iNumCycles << endl);
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << endl);
       
       if (CurrKeyWord == SINE) {			    
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 SineDriveCaller,
				 SineDriveCaller(pDrvHdl,
						 dInitialTime,
						 dOmega,
						 dAmplitude,
						 iNumCycles,
						 dInitialValue), 
				 DMmm);
       } else if(CurrKeyWord == COSINE) {			    
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 CosineDriveCaller,
				 CosineDriveCaller(pDrvHdl,
						   dInitialTime, 
						   dOmega,
						   dAmplitude,
						   iNumCycles,
						   dInitialValue), 
				 DMmm);
       }
       
       break;
    }
      
    case FREQUENCYSWEEP: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << endl);
       
       DriveCaller* pOmega = ReadDriveData(pDM, HP, pDrvHdl);
       // DEBUGCOUT("Omega: " << dOmega << endl);
       
       DriveCaller* pAmplitude = ReadDriveData(pDM, HP, pDrvHdl);
       // DEBUGCOUT("Amplitude: " << dAmplitude << endl);
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("Initial value: " << dInitialValue << endl);
       
       doublereal dFinalTime = HP.GetReal();
       DEBUGCOUT("Final time: " << dFinalTime << endl);
       
       doublereal dFinalValue = HP.GetReal();	   
       DEBUGCOUT("Final value: " << dFinalValue << endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      FreqSweepDriveCaller,
			      FreqSweepDriveCaller(pDrvHdl,
						   dInitialTime,
						   pOmega,
						   pAmplitude,
						   dInitialValue,
						   dFinalTime,
						   dFinalValue), 
			      DMmm);       
       
       break;
    }
      
      /* esponenziale */
    case EXPONENTIAL: {	
       doublereal dAmplitude = HP.GetReal(1.);
       DEBUGCOUT("Amplitude value: " << dAmplitude << endl);
       
       doublereal dTimeConst = HP.GetReal();
       DEBUGCOUT("Time constant value: " << dTimeConst << endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << endl);
       
       doublereal dInitialValue = HP.GetReal();
       DEBUGCOUT("Initial value: " << dInitialValue << endl);
       
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      ExpDriveCaller,
			      ExpDriveCaller(pDrvHdl, dAmplitude, dTimeConst,
					     dInitialTime, dInitialValue), 
			      DMmm);
       break;
    }
      
      /* Numero casuale */
    case RANDOM: {	
       doublereal dAmplitude = HP.GetReal(1.);
       DEBUGCOUT("Amplitude value: " << dAmplitude << endl);
       
       doublereal dRefVal = HP.GetReal();
       DEBUGCOUT("Mean value: " << dRefVal << endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << endl);
       
       doublereal dFinalTime = HP.GetReal();
       DEBUGCOUT("Final time: " << dFinalTime << endl);
       
       /* Type of random number (additional data) */
       integer iSteps = 1;
       while (1) {
	  if (HP.IsKeyWord("steps")) {
	     iSteps = HP.GetInt();
	     if (iSteps <= 0) {		    
		cerr << "Warning: Steps number " << iSteps 
		  << " is illegal; resorting to default value" << endl;
		iSteps = 1;
	     }
	     DEBUGCOUT("Force changes every " << iSteps 
		       << " steps" << endl);
	  } else if (HP.IsKeyWord("seed")) {
	     if (HP.IsKeyWord("time")) {
		DEBUGCOUT("(Re)Seeding random numbers with current time ...");
		srand(time(NULL));
	     } else {
		DEBUGCOUT("(Re)Seeding random numbers with given value ...");
		srand(HP.GetInt());
	     }
	  } else {
	     break;
	  }	      
       }
       
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      RandDriveCaller,
			      RandDriveCaller(pDrvHdl, 
					      dAmplitude, 
					      dRefVal,
					      dInitialTime, 
					      dFinalTime,
					      iSteps),
			      DMmm);
       break;
    }

      /* Lineare a pezzi */
    case PIECEWISELINEAR: {
       unsigned int n = HP.GetInt();
       DEBUGCOUT("number of points: " << n << endl);

       if (n < 2) {
	  cerr << "Need at least two points for piecewise linear drive at line "
	    << HP.GetLineData() << endl;
	  THROW(DataManager::ErrGeneric());
       }

       doublereal *p = NULL;
       SAFENEWARR(p, doublereal, 2*n, DMmm);
       p[0] = HP.GetReal();
       p[n] = HP.GetReal();
       for (unsigned int i = 1; i < n; i++) {
	  p[i] = HP.GetReal();
	  if (p[i] <= p[i-1]) {
	     cerr << "point " << p[i]
	       << " is smaller than or equal to preceding point " << p[i-1]
	       << " at line " << HP.GetLineData();
	     THROW(DataManager::ErrGeneric());
	  }
	  p[i+n] = HP.GetReal();
       }

       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
		       PiecewiseLinearDriveCaller,
		       PiecewiseLinearDriveCaller(pDrvHdl, n, p),
		       DMmm);
       
       break;
    }
      
      /* driver stringa da valutare */
    case STRING: {	     
       /* lettura dei dati specifici */
       const char* sTmp = HP.GetStringWithDelims();
       DEBUGCOUT("String to evaluate: \"" << sTmp << '\"' << endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      StringDriveCaller,
			      StringDriveCaller(pDrvHdl, sTmp), 
			      DMmm);

       break;
    }
      
      /* driver legato ad un grado di liberta' nodale */
    case DOF: {
       if (pDM == NULL) {
	  cerr << "sorry, since the driver is not owned by a DataManager" << endl
	    << "no DOF dependent drivers are allowed;" << endl
	    << "aborting ...";	  
	  THROW(DataManager::ErrGeneric());
       }

       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       
#ifdef USE_MPI
       if (MPI::COMM_WORLD.Get_size() > 1) {
          cerr << "warning: add explicit connection entry for the "
            << psNodeNames[SD.pNode->GetNodeType()] 
	    << "(" << SD.pNode->GetLabel() << ") dof drive"
            " at line " << HP.GetLineData() << endl;
       }
#endif /* USE_MPI */
       
       /* Chiamata ricorsiva a leggere il drive supplementare */
       DriveCaller* pTmp = ReadDriveData(pDM, HP, pDrvHdl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DofDriveCaller,
			      DofDriveCaller(pDrvHdl, pTmp, SD),
			      DMmm);
       
       HP.PutKeyTable(K);
       break;
    }

      
      /* driver legato ad un grado di liberta' nodale */
    case ARRAY: {
       unsigned short int iNumDr = (unsigned short int)HP.GetInt();
       if (iNumDr == 0) {	      
	  cerr << "Sorry, at least one driver is required" << endl;
	  THROW(ErrGeneric());
       } else if (iNumDr == 1) {
	  /* creazione di un driver normale mediante chiamata ricorsiva */
	  pDC = ReadDriveData(pDM, HP, pDrvHdl);
	  HP.PutKeyTable(K);
       } else {
	  DriveCaller** ppDC = NULL;
	  SAFENEWARR(ppDC, DriveCaller*, iNumDr, DMmm);
	  for (int i = 0; i < iNumDr; i++) {
	     ppDC[i] = ReadDriveData(pDM, HP, pDrvHdl);
	     HP.PutKeyTable(K);
	  }
	  
	  /* allocazione e creazione array */
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 DriveArrayCaller,
				 DriveArrayCaller(pDrvHdl, iNumDr, 
						  (const DriveCaller**)ppDC),
				 DMmm);
       }
       break;
    }
      
      
      /* drive file */
    case FILEDRIVE: {	     
       /* lettura dei dati specifici */
       unsigned int uL = HP.GetInt();
       FileDrive* pDrv = (FileDrive*)pDM->pFindDrive(Drive::FILEDRIVE, uL);
       if (pDrv == NULL) {
	  cerr << "line " << HP.GetLineData() 
	    << ": can't find FileDrive(" << uL << ")" << endl;
	  THROW(ErrGeneric());
       }
              
       integer id = 1;
       if (HP.fIsArg()) {
	  id = HP.GetInt(id);
       }
       
       /* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
			       FileDriveCaller,
			       FileDriveCaller(pDM->pGetDrvHdl(), pDrv, id),
			       DMmm);     
       
       break;
    }
      
    default: {
       cerr << "unknown drive type at line " << HP.GetLineData() << endl;       
       THROW(DataManager::ErrGeneric());
    }	
   }
   
   ASSERT(pDC != NULL);
   return pDC;   
} /* ReadDriveData */
