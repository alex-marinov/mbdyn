/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "ac/math.h"
#include "ac/float.h"

#ifdef USE_MPI
#include "mbcomm.h"
#endif /* USE_MPI */

#include "dataman.h"
#include "drive_.h"
#include "dofdrive.h"
#include "privdrive.h"
#include "filedrv.h"
#include "ddrive.h"

/* StringDriveCaller - begin */

StringDriveCaller::StringDriveCaller(const DriveHandler* pDH,
				     const char* const sTmpStr)
: DriveCaller(pDH), sEvalStr(NULL), iEvalStrLen(0) 
{
   ASSERT(sTmpStr != NULL);
   iEvalStrLen = strlen(sTmpStr);
   SAFESTRDUP(sEvalStr, sTmpStr);
}


StringDriveCaller::~StringDriveCaller(void)
{
   ASSERT(sEvalStr != NULL);
   if (sEvalStr != NULL) {
      SAFEDELETEARR(sEvalStr);
   }
}


/* Copia */
DriveCaller* StringDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  StringDriveCaller, 
			  StringDriveCaller(pDrvHdl, sEvalStr));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& StringDriveCaller::Restart(std::ostream& out) const
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
   SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDrvHdl));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& TimeDriveCaller::Restart(std::ostream& out) const
{
   return out << "time";
}
 
/* TimeDriveCaller - end */


/* ConstDriveCaller - begin */

ConstDriveCaller::ConstDriveCaller(doublereal d)
: DriveCaller(0), dConst(d)
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
   SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(dConst));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& ConstDriveCaller::Restart(std::ostream& out) const
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
   SAFENEWWITHCONSTRUCTOR(pDC, LinearDriveCaller, LinearDriveCaller(pDrvHdl, dC0, dC1));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& LinearDriveCaller::Restart(std::ostream& out) const
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
   SAFENEWWITHCONSTRUCTOR(pDC, ParabolicDriveCaller, ParabolicDriveCaller(pDrvHdl, dC0, dC1, dC2));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& ParabolicDriveCaller::Restart(std::ostream& out) const
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
   SAFENEWWITHCONSTRUCTOR(pDC, CubicDriveCaller, CubicDriveCaller(pDrvHdl, dC0, dC1, dC2, dC3));
   
   return pDC;
}


CubicDriveCaller::~CubicDriveCaller(void)
{
   NO_OP;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& CubicDriveCaller::Restart(std::ostream& out) const
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
					  dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& StepDriveCaller::Restart(std::ostream& out) const
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
						dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& DoubleStepDriveCaller::Restart(std::ostream& out) const
{
   return out 
     << " double step, "
     << dStepTime << ", " 
     << dEndStepTime << ", "
     << dStepValue << ", "
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
					  dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& RampDriveCaller::Restart(std::ostream& out) const
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
						dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& DoubleRampDriveCaller::Restart(std::ostream& out) const
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
					  dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& SineDriveCaller::Restart(std::ostream& out) const
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


/* FourierSeriesDriveCaller - begin */

FourierSeriesDriveCaller::FourierSeriesDriveCaller(const DriveHandler* pDH,
	doublereal dStartTime,
	doublereal dOmega,
	std::vector<doublereal>& a, 
	integer iNumCyc,
	doublereal dInitialValue)
: DriveCaller(pDH),
dStartTime(dStartTime),
dOmega(dOmega),
iNumCycles(iNumCyc),
dInitialValue(dInitialValue),
bNeverEnd(false)
{
	ASSERT(iNumCycles >= 0);

	if (iNumCycles > 0) {
		dEndTime = dStartTime
			+ 2.*M_PI/dOmega*doublereal(iNumCycles);
	
	/* Onde di coseno che continuano all'infinito */
	} else if (iNumCycles == 0) {	     
		dEndTime = 0.;
		bNeverEnd = true;
	}

	amplitudes.resize(a.size());
	for (unsigned i = 0; i < a.size(); i++) {
		amplitudes[i] = a[i];
	}
}


FourierSeriesDriveCaller::~FourierSeriesDriveCaller(void)
{
	NO_OP;
}


/* Copia */
DriveCaller*
FourierSeriesDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC, 
		FourierSeriesDriveCaller, 
		FourierSeriesDriveCaller(pDrvHdl,
			dStartTime,
			dOmega,
			amplitudes,
			iNumCycles,
			dInitialValue));
 
	return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream&
FourierSeriesDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " fourier series, "
		<< dStartTime << ", "
		<< dOmega << ", "
		<< amplitudes.size()/2 << ", "
		<< 2.*amplitudes[0] << ", ";
	for (unsigned i = 1; i < amplitudes.size(); i++) {
		out << amplitudes[i] << ", ";
	}
	out
		<< iNumCycles << ", "
		<< dInitialValue;
}

/* FourierSeriesDriveCaller - end */


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
					    dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& CosineDriveCaller::Restart(std::ostream& out) const
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
dStartTime(d1), pOmega((DriveCaller*)pOmega),
pAmplitude((DriveCaller*)pAmplitude),
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
					       dFinalValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& FreqSweepDriveCaller::Restart(std::ostream& out) const
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
					 dInitialValue));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& ExpDriveCaller::Restart(std::ostream& out) const
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
					  iSteps));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& RandDriveCaller::Restart(std::ostream& out) const
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


/* MeterDriveCaller - begin */

MeterDriveCaller::MeterDriveCaller(const DriveHandler* pDH, 
				 doublereal dS, doublereal dE, integer iS)
: DriveCaller(pDH),
dStartTime(dS),
dEndTime(dE),
iSteps(iS)
{
   iMeterDriveNumber = ((DriveHandler*)pDrvHdl)->iMeterInit(iSteps);
}


MeterDriveCaller::~MeterDriveCaller(void)
{
   NO_OP; 
}


/* Copia */
DriveCaller* MeterDriveCaller::pCopy(void) const
{
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  MeterDriveCaller, 
			  MeterDriveCaller(pDrvHdl,
					  dStartTime, 
					  dEndTime,
					  iSteps));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& MeterDriveCaller::Restart(std::ostream& out) const
{
   out
     << " meter, " 
     << dStartTime << ", "
     << dEndTime;
   if(iSteps > 1) {
      out << ", steps, " << iSteps;
   }
   return out;
}

/* MeterDriveCaller - end */


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
		SAFEDELETEARR(pPoints);
	}
}

/* Copia */
DriveCaller* PiecewiseLinearDriveCaller::pCopy(void) const
{
	doublereal *p = NULL;
	SAFENEWARR(p, doublereal, 2*iNumPoints);
	for (unsigned int i = 0; i < 2*iNumPoints; i++) {
		p[i] = pPoints[i];
	}
	
	DriveCaller* pDC = NULL;
	SAFENEWWITHCONSTRUCTOR(pDC,
			PiecewiseLinearDriveCaller,
			PiecewiseLinearDriveCaller(pDrvHdl, iNumPoints, p));

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream& PiecewiseLinearDriveCaller::Restart(std::ostream& out) const
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
      SAFEDELETE(ppDriveCallers[i]);
   }
   SAFEDELETEARR(ppDriveCallers);
}


/* Copia */
DriveCaller* DriveArrayCaller::pCopy(void) const
{
   ASSERT(ppDriveCallers != NULL);
   ASSERT(iNumDrivers > 0);
   DriveCaller** ppDC = NULL;
   SAFENEWARR(ppDC, DriveCaller*, iNumDrivers);
   for (int i = 0; i < iNumDrivers; i++) {
      ppDC[i] = ppDriveCallers[i]->pCopy();
   }
   
   DriveCaller* pDC = NULL;
   SAFENEWWITHCONSTRUCTOR(pDC, 
			  DriveArrayCaller, 
			  DriveArrayCaller(pDrvHdl, 
					   iNumDrivers,
					   (const DriveCaller**)ppDC));
   
   return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */   
std::ostream& DriveArrayCaller::Restart(std::ostream& out) const
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

DriveCaller *
ReadDriveData(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
   DEBUGCOUTFNAME("ReadDriveData()");
   
   const DriveHandler* pDrvHdl = 0;
   if (pDM != 0) {
      pDrvHdl = pDM->pGetDrvHdl();
   }

   const char* sKeyWords[] = {
      "time",
	"null",
	"one",	/* deprecated */
	"unit",
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
	"fourier" "series",
	"frequency" "sweep",
	"exponential",
	"random",
	"meter",
	"piecewise" "linear",
	"file", 
	"string",
	"dof",
	"element",
	"drive",
	"array",
	NULL
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,

	TIME = 0,
	NULLDRIVE,
	ONEDRIVE,
	UNITDRIVE,
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
	FOURIERSERIES,
	FREQUENCYSWEEP,
	EXPONENTIAL,
	RANDOM,
	METER,
	PIECEWISELINEAR,
	FILEDRIVE,      
	STRING,
	DOF,
	ELEMENT,
	DRIVE,
	ARRAY,

	LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);
   
   /* lettura del tipo di drive */   
   KeyWords CurrKeyWord;
   if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
      CurrKeyWord = CONST;
   } else if (CurrKeyWord == FUNCTION) {
      CurrKeyWord = KeyWords(HP.GetWord());
   }   
   
#ifdef DEBUG   
   if (CurrKeyWord >= 0) {      
      std::cout<< "drive type: " << sKeyWords[CurrKeyWord] << std::endl;
   }   
#endif   

   if (pDrvHdl == 0 && !bDeferred) { 
      switch (CurrKeyWord) {
       case NULLDRIVE:
       case ONEDRIVE:
       case UNITDRIVE:
       case CONST:
       case DRIVE:
       case ARRAY:

       /* leave the checks for later... */
       case UNKNOWN:
	 /* drives that don't need the data manager... */
	 break;

       default:
	 silent_cerr(sKeyWords[CurrKeyWord] << " drive needs data manager "
			 "at line " << HP.GetLineData() << std::endl);
	 throw DataManager::ErrNeedDataManager();
      }
   }
   
   DriveCaller* pDC = NULL;
   
   switch (CurrKeyWord) {

      /* time */
    case TIME:
     /* allocazione e creazione */
     SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDrvHdl));
     break;  
      
      /* driver nullo */
    case NULLDRIVE:
       /* allocazione e creazione */
        SAFENEW(pDC, NullDriveCaller);
      
       /* scrittura dei dati specifici */	     
       
       break;
      
      /* driver unitario*/
    case ONEDRIVE:
    case UNITDRIVE:
       /* allocazione e creazione */
        SAFENEW(pDC, OneDriveCaller);
      
       /* scrittura dei dati specifici */	     
       
       break;
 
      /* driver costante */
    case CONST: {
       /* lettura dei dati specifici */
       doublereal dConst = HP.GetReal();
       DEBUGCOUT("Const value: " << dConst << std::endl);

       /* allocazione e creazione */
       if (dConst == 0.) {
          SAFENEW(pDC, NullDriveCaller);
       } else if (dConst == 1.) {
          SAFENEW(pDC, OneDriveCaller);
      
       } else {
	  SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller,
  			  ConstDriveCaller(dConst));
       }
       
       break;
    }
      
      /* retta */
    case LINEAR: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << std::endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      LinearDriveCaller,
			      LinearDriveCaller(pDrvHdl, dC0, dC1));
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* parabola */
    case PARABOLIC: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << std::endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << std::endl);
       
       doublereal dC2 = HP.GetReal();
       DEBUGCOUT("Parabolic slope: " << dC2 << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      ParabolicDriveCaller,
			      ParabolicDriveCaller(pDrvHdl, dC0, dC1, dC2));
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* cubica */
    case CUBIC: {
       /* lettura dei dati specifici */
       doublereal dC0 = HP.GetReal();
       DEBUGCOUT("Offset: " << dC0 << std::endl);
       
       doublereal dC1 = HP.GetReal();
       DEBUGCOUT("Slope: " << dC1 << std::endl);
       
       doublereal dC2 = HP.GetReal();
       DEBUGCOUT("Parabolic slope: " << dC2 << std::endl);
       
       doublereal dC3 = HP.GetReal();
       DEBUGCOUT("Cubic slope: " << dC3 << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      CubicDriveCaller,
			      CubicDriveCaller(pDrvHdl, dC0, dC1, dC2, dC3));
       
       /* scrittura dei dati specifici */	     
       
       break;
    }
      
      /* Scalino semplice */
    case STEP: {	
       doublereal dStepTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dStepTime << std::endl);
       
       doublereal dStepValue = HP.GetReal(1.);
       DEBUGCOUT("Step Value: " << dStepValue << std::endl);
       
       doublereal dInitialValue = HP.GetReal();
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      StepDriveCaller,
			      StepDriveCaller(pDrvHdl, dStepTime,
					      dStepValue, dInitialValue));
       break;
    }
      
      /* Scalino doppio */
    case DOUBLESTEP: {	
       doublereal dStepTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dStepTime << std::endl);
       
       doublereal dEndStepTime = HP.GetReal();
       DEBUGCOUT("Final time: " << dEndStepTime << std::endl);
       
       if(dEndStepTime <= dStepTime) {	      
	  silent_cerr("Warning at line " 
	    << HP.GetLineData() 
	    << ": final time " << dEndStepTime
	    << " is less than or equal to initial time " << dStepTime 
	    << " in double step func drive" << std::endl);
       }
       
       doublereal dStepValue = HP.GetReal(1.);
       DEBUGCOUT("Step Value: " << dStepValue << std::endl);
       
       doublereal dInitialValue = HP.GetReal();     
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DoubleStepDriveCaller,
			      DoubleStepDriveCaller(pDrvHdl, 
						    dStepTime, 
						    dStepValue,
						    dEndStepTime,
						    dInitialValue));
       break;
    }
      
      /* Rampa saturata */
    case RAMP: {	
       doublereal dSlope = HP.GetReal(1.);
       DEBUGCOUT("Slope Value: " << dSlope << std::endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

       doublereal dFinalTime = DBL_MAX;
       if (!HP.IsKeyWord("forever")) {
	  dFinalTime = HP.GetReal();
       }
       DEBUGCOUT("Final time: " << dFinalTime << std::endl);
       
       if(dFinalTime <= dInitialTime) {	      
	  silent_cerr("Warning at line "
	    << HP.GetLineData() 
	    << ": final time " << dFinalTime
	    << " is less than or equal to initial time " << dInitialTime 
	    << " in ramp func drive" << std::endl);
       }	   
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC, 
			      RampDriveCaller,
			      RampDriveCaller(pDrvHdl, 
					      dSlope, dInitialTime,
					      dFinalTime, dInitialValue));
       break;
    }
      
      /* Rampa doppia */
    case DOUBLERAMP: {	
       doublereal dAscendingSlope = HP.GetReal(1.);
       DEBUGCOUT("Ascending Slope Value: " << dAscendingSlope << std::endl);
       
       doublereal dAscendingInitialTime = HP.GetReal();
       DEBUGCOUT("Ascending Initial time: " << dAscendingInitialTime << std::endl);
       
       doublereal dAscendingFinalTime = HP.GetReal();
       DEBUGCOUT("Ascending Final time: " << dAscendingFinalTime << std::endl);
       
       if (dAscendingFinalTime <= dAscendingInitialTime) {	      
	  silent_cerr("Warning at line " 
	    << HP.GetLineData() << ": ascending final time " 
	    << dAscendingFinalTime
	    << " is less than or equal to ascending initial time " 
	    << dAscendingInitialTime 
	    << " in double ramp func drive" << std::endl);
	  throw ErrGeneric();
       }
       
       doublereal dDescendingSlope = HP.GetReal(-1.);
       DEBUGCOUT("Descending Slope Value: " << dDescendingSlope << std::endl);
       
       doublereal dDescendingInitialTime = DBL_MAX;
       if (!HP.IsKeyWord("forever")) {
          dDescendingInitialTime = HP.GetReal();
       }
       DEBUGCOUT("Descending Initial time: " << dDescendingInitialTime << std::endl);
       
       if (dDescendingInitialTime < dAscendingFinalTime) {
	  silent_cerr("Warning at line " 
	    << HP.GetLineData() << ": descending initial time " 
	    << dDescendingInitialTime
	    << " is less than ascending final time " 
	    << dAscendingFinalTime 
	    << " in double ramp func drive" << std::endl);
	  throw ErrGeneric();
       }	   
       
       doublereal dDescendingFinalTime = HP.GetReal();
       DEBUGCOUT("Descending Final time: " << dDescendingFinalTime << std::endl);
       
       if (dDescendingFinalTime <= dDescendingInitialTime) {	      
	  silent_cerr("Warning at line " 
	    << HP.GetLineData() << ": descending final time " 
	    << dDescendingFinalTime
	    << " is less than descending initial time " 
	    << dDescendingInitialTime 
	    << " in double ramp func drive" << std::endl);
	  throw ErrGeneric();
       }	   
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DoubleRampDriveCaller,
			      DoubleRampDriveCaller(pDrvHdl,
						    dAscendingSlope,
						    dAscendingInitialTime,
						    dAscendingFinalTime,
						    dDescendingSlope,
						    dDescendingInitialTime,
						    dDescendingFinalTime,
						    dInitialValue));
       break;
    }
      
      /* Seno e coseno (limitati, illimitati, saturati) */
    case SINE:
    case COSINE: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       doublereal dOmega = HP.GetReal(1.);
       DEBUGCOUT("Omega: " << dOmega << std::endl);
       
       doublereal dAmplitude = HP.GetReal();
       DEBUGCOUT("Amplitude: " << dAmplitude << std::endl);
       
       integer iNumCycles;
       if (HP.IsKeyWord("forever")) {
	  iNumCycles = 0;
       } else if (HP.IsKeyWord("one")) {
	  iNumCycles = 1;
       } else if (HP.IsKeyWord("half")) {
	  iNumCycles = -1;
       } else {
          iNumCycles = HP.GetInt();
       }
       DEBUGCOUT("Number of cycles: " << iNumCycles << std::endl);
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       if (CurrKeyWord == SINE) {			    
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 SineDriveCaller,
				 SineDriveCaller(pDrvHdl,
						 dInitialTime,
						 dOmega,
						 dAmplitude,
						 iNumCycles,
						 dInitialValue));
       } else if(CurrKeyWord == COSINE) {			    
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 CosineDriveCaller,
				 CosineDriveCaller(pDrvHdl,
						   dInitialTime, 
						   dOmega,
						   dAmplitude,
						   iNumCycles,
						   dInitialValue));
       }
       
       break;
    }
      
      /* Seno e coseno (limitati, illimitati, saturati) */
    case FOURIERSERIES: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       doublereal dOmega = HP.GetReal(1.);
       DEBUGCOUT("Omega: " << dOmega << std::endl);

       int n = HP.GetInt();
       if (n <= 0) {
	       silent_cerr("FourierSeriesDriveCaller: invalid order " << n
		       << " at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }

       std::vector<doublereal> a(1 + 2*n); 
       for (unsigned i = 0; i < 1 + 2*unsigned(n); i++) {
	       a[i] = HP.GetReal();
       }
       /* don't remember why, but the series starts with a_0/2 */
       a[0] /= 2.;
       
       integer iNumCycles;
       if (HP.IsKeyWord("forever")) {
	  iNumCycles = 0;
       } else if (HP.IsKeyWord("one")) {
	  iNumCycles = 1;
       } else {
          iNumCycles = HP.GetInt();
       }
       DEBUGCOUT("Number of cycles: " << iNumCycles << std::endl);
       if (iNumCycles < 0) {
	       silent_cerr("FourierSeriesDriveCaller: invalid number of cycles "
		       << iNumCycles << " at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
		FourierSeriesDriveCaller,
		FourierSeriesDriveCaller(pDrvHdl,
			dInitialTime,
			dOmega,
			a,
			iNumCycles,
			dInitialValue));
       
       break;
    }
      
    case FREQUENCYSWEEP: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       DriveCaller* pOmega = 0;
       pOmega = HP.GetDriveCaller();
       // DEBUGCOUT("Omega: " << dOmega << std::endl);
       
       DriveCaller* pAmplitude = 0;
       pAmplitude = HP.GetDriveCaller();
       // DEBUGCOUT("Amplitude: " << dAmplitude << std::endl);
       
       doublereal dInitialValue = HP.GetReal();	   
       DEBUGCOUT("Initial value: " << dInitialValue << std::endl);
       
       doublereal dFinalTime = DBL_MAX;
       if (!HP.IsKeyWord("forever")) {
          dFinalTime = HP.GetReal();
       }
       DEBUGCOUT("Final time: " << dFinalTime << std::endl);
       
       doublereal dFinalValue = HP.GetReal();	   
       DEBUGCOUT("Final value: " << dFinalValue << std::endl);
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      FreqSweepDriveCaller,
			      FreqSweepDriveCaller(pDrvHdl,
						   dInitialTime,
						   pOmega,
						   pAmplitude,
						   dInitialValue,
						   dFinalTime,
						   dFinalValue));
       
       break;
    }
      
      /* esponenziale */
    case EXPONENTIAL: {	
       doublereal dAmplitude = HP.GetReal(1.);
       DEBUGCOUT("Amplitude value: " << dAmplitude << std::endl);
       
       doublereal dTimeConst = HP.GetReal();
       DEBUGCOUT("Time constant value: " << dTimeConst << std::endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       doublereal dInitialValue = HP.GetReal();
       DEBUGCOUT("Initial value: " << dInitialValue << std::endl);
       
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      ExpDriveCaller,
			      ExpDriveCaller(pDrvHdl, dAmplitude, dTimeConst,
					     dInitialTime, dInitialValue));
       break;
    }
      
      /* Numero casuale */
    case RANDOM: {	
       doublereal dAmplitude = HP.GetReal(1.);
       DEBUGCOUT("Amplitude value: " << dAmplitude << std::endl);
       
       doublereal dRefVal = HP.GetReal();
       DEBUGCOUT("Mean value: " << dRefVal << std::endl);
       
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       doublereal dFinalTime = DBL_MAX;
       if (!HP.IsKeyWord("forever")) {
          dFinalTime = HP.GetReal();
       }
       DEBUGCOUT("Final time: " << dFinalTime << std::endl);
       
       /* Type of random number (additional data) */
       integer iSteps = 1;
       while (true) {
	  if (HP.IsKeyWord("steps")) {
	     iSteps = HP.GetInt();
	     if (iSteps <= 0) {		    
		silent_cerr("Warning: Steps number " << iSteps 
		  << " is illegal; resorting to default value" << std::endl);
		iSteps = 1;
	     }
	     DEBUGCOUT("Force changes every " << iSteps 
		       << " steps" << std::endl);
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
					      iSteps));
       break;
    }

      /* spike every N steps */
    case METER: {	
       doublereal dInitialTime = HP.GetReal();
       DEBUGCOUT("Initial time: " << dInitialTime << std::endl);
       
       doublereal dFinalTime = DBL_MAX;
       if (!HP.IsKeyWord("forever")) {
          dFinalTime = HP.GetReal();
       }
       DEBUGCOUT("Final time: " << dFinalTime << std::endl);
       
       /* Type of random number (additional data) */
       integer iSteps = 1;
       while (true) {
	  if (HP.IsKeyWord("steps")) {
	     iSteps = HP.GetInt();
	     if (iSteps <= 0) {		    
		silent_cerr("Warning: Steps number " << iSteps 
		  << " is illegal; resorting to default value" << std::endl);
		iSteps = 1;
	     }
	     DEBUGCOUT("Force changes every " << iSteps 
		       << " steps" << std::endl);
	  }else {
	     break;
	  }	      
       }
       
       
       SAFENEWWITHCONSTRUCTOR(pDC,
			      MeterDriveCaller,
			      MeterDriveCaller(pDrvHdl, 
					      dInitialTime, 
					      dFinalTime,
					      iSteps));
       break;
    }

      /* Lineare a pezzi */
    case PIECEWISELINEAR: {
       unsigned int n = HP.GetInt();
       DEBUGCOUT("number of points: " << n << std::endl);

       if (n < 2) {
	  silent_cerr("Need at least two points for piecewise linear drive at line "
	    << HP.GetLineData() << std::endl);
	  throw DataManager::ErrGeneric();
       }

       doublereal *p = NULL;
       SAFENEWARR(p, doublereal, 2*n);
       p[0] = HP.GetReal();
       p[n] = HP.GetReal();
       for (unsigned int i = 1; i < n; i++) {
	  p[i] = HP.GetReal();
	  if (p[i] <= p[i-1]) {
	     silent_cerr("point " << p[i]
	       << " is smaller than or equal to preceding point " << p[i-1]
	       << " at line " << HP.GetLineData() << std::endl);
	     throw DataManager::ErrGeneric();
	  }
	  p[i+n] = HP.GetReal();
       }

       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
		       PiecewiseLinearDriveCaller,
		       PiecewiseLinearDriveCaller(pDrvHdl, n, p));
       
       break;
    }
      
      /* driver stringa da valutare */
    case STRING: {	     
       /* lettura dei dati specifici */
       const char* sTmp = HP.GetStringWithDelims();
       DEBUGCOUT("String to evaluate: \"" << sTmp << '\"' << std::endl);
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      StringDriveCaller,
			      StringDriveCaller(pDrvHdl, sTmp));

       break;
    }
      
      /* driver legato ad un grado di liberta' nodale */
    case DOF: {
       if (pDM == NULL) {
	  silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
	    << "no DOF dependent drivers are allowed;" << std::endl
	    << "aborting ..." << std::endl);	  
	  throw DataManager::ErrGeneric();
       }

       ScalarDof SD = ReadScalarDof(pDM, HP, 1);
       
#ifdef USE_MPI
       if (MPI::Is_initialized() && MBDynComm.Get_size() > 1) {
          silent_cerr("warning: add explicit connection entry for "
            << psNodeNames[SD.pNode->GetNodeType()] 
	    << "(" << SD.pNode->GetLabel() << ") dof drive"
            " at line " << HP.GetLineData() << std::endl);
       }
#endif /* USE_MPI */
       
       /* Chiamata ricorsiva a leggere il drive supplementare */
       DriveCaller* pTmp = 0;
       pTmp = HP.GetDriveCaller();
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      DofDriveCaller,
			      DofDriveCaller(pDrvHdl, pTmp, SD));
       
       break;
    }

      /* driver legato ai dati privati di un elemento */
    case ELEMENT: {
       if (pDM == NULL) {
	  silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
	    << "no element dependent drivers are allowed;" << std::endl
	    << "aborting ..." << std::endl);
	  throw DataManager::ErrGeneric();
       }

       unsigned uLabel = HP.GetInt();
       KeyTable Kel(HP, psReadElemsElems);
       int k = HP.IsKeyWord();
       if (k == -1) {
	       const char *s = HP.GetString();
	       silent_cerr("unknown element type \"" << s
			       << "\" at line " << HP.GetLineData()
			       << std::endl);
	       throw ErrGeneric();
       }
       Elem *pElem = (Elem*)pDM->pFindElem(Elem::Type(k), uLabel);
       if (pElem ==  NULL) {
	       silent_cerr("unable to find " << psElemNames[Elem::Type(k)]
			       << "(" << uLabel << ") at line "
			       << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }
       unsigned int iMaxIndex = pElem->iGetNumPrivData();
       unsigned int iIndex = 0;
       const char *sIndexName = NULL;
       if (HP.IsKeyWord("string")) {
	       const char *s = HP.GetStringWithDelims();
	       iIndex = pElem->iGetPrivDataIdx(s);
	       SAFESTRDUP(sIndexName, s);

       } else if (HP.IsKeyWord("index")) {
	       iIndex = HP.GetInt();

       } else if (iMaxIndex == 1) {
	       iIndex = 1;

       } else {
	       silent_cerr("need a private data index for " 
			       << psElemNames[pElem->GetElemType()] 
			       << "(" << pElem->GetLabel() << ") "
			       "at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }

       if (iIndex < 1 || iIndex > iMaxIndex) {
	       silent_cerr("illegal index " << iIndex << " for "
			       << psElemNames[pElem->GetElemType()] 
			       << "(" << pElem->GetLabel() << ") "
			       "at line " << HP.GetLineData() << std::endl);
	       throw ErrGeneric();
       }

#ifdef USE_MPI
       /* FIXME: todo ... */
       if (MPI::Is_initialized() && MBDynComm.Get_size() > 1) {
          silent_cerr("warning: add explicit connection entry for "
            << psElemNames[pElem->GetElemType()] 
	    << "(" << pElem->GetLabel() << ") element drive"
            " at line " << HP.GetLineData() << std::endl);
       }
#endif /* USE_MPI */
       
       /* Chiamata ricorsiva a leggere il drive supplementare */
       DriveCaller* pTmp = 0;
       pTmp = HP.GetDriveCaller();
       
       /* allocazione e creazione */
       SAFENEWWITHCONSTRUCTOR(pDC,
			      PrivDriveCaller,
			      PrivDriveCaller(pDrvHdl, pTmp,
				      pElem, iIndex, sIndexName));
       if (sIndexName) {
	       SAFEDELETEARR(sIndexName);
       }
       break;
    }

    case DRIVE: {
       DriveCaller *pD1 = ReadDriveData(pDM, HP, bDeferred);
       DriveCaller *pD2 = ReadDriveData(pDM, HP, bDeferred);

       SAFENEWWITHCONSTRUCTOR(pDC,
		       DriveDriveCaller,
		       DriveDriveCaller(pDrvHdl, pD1, pD2));
       break;
    }

      
      /* driver legato ad un grado di liberta' nodale */
    case ARRAY: {
       unsigned short int iNumDr = (unsigned short int)HP.GetInt();
       if (iNumDr == 0) {	      
	  silent_cerr("Sorry, at least one driver is required" << std::endl);
	  throw ErrGeneric();
       } else if (iNumDr == 1) {
	  /* creazione di un driver normale mediante chiamata ricorsiva */
	  pDC = HP.GetDriveCaller();
	  // HP.PutKeyTable(K);
       } else {
	  DriveCaller** ppDC = NULL;
	  SAFENEWARR(ppDC, DriveCaller*, iNumDr);
	  for (int i = 0; i < iNumDr; i++) {
	     ppDC[i] = HP.GetDriveCaller();
	  }
	  // HP.PutKeyTable(K);
	  
	  /* allocazione e creazione array */
	  SAFENEWWITHCONSTRUCTOR(pDC,
				 DriveArrayCaller,
				 DriveArrayCaller(pDrvHdl, iNumDr, 
						  (const DriveCaller**)ppDC));
       }
       break;
    }
      
      
      /* drive file */
    case FILEDRIVE: {	     
       /* lettura dei dati specifici */
       unsigned int uL = HP.GetInt();
       FileDrive* pDrv = (FileDrive*)pDM->pFindDrive(Drive::FILEDRIVE, uL);
       if (pDrv == NULL) {
	  silent_cerr("line " << HP.GetLineData() 
	    << ": can't find FileDrive(" << uL << ")" << std::endl);
	  throw ErrGeneric();
       }
              
       integer id = 1;
       if (HP.IsArg()) {
	  id = HP.GetInt(id);
       }

       doublereal da = 1.;
       if (HP.IsKeyWord("amplitude")) {
	       da = HP.GetReal();
       }
       
       /* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
			       FileDriveCaller,
			       FileDriveCaller(pDM->pGetDrvHdl(), pDrv, id, da));
       
       break;
    }
      
    default: {
       silent_cerr("unknown drive type at line " << HP.GetLineData() << std::endl);
       throw DataManager::ErrGeneric();
    }	
   }
   
   ASSERT(pDC != NULL);
   return pDC;   
} /* ReadDriveData */
