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

#ifndef DRIVE__H
#define DRIVE__H

/* include generali */
#include <ac/sstream>

/* include per il debug */
#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

/* include del programma */
#include <mathp.h>
#include <output.h>
#include <withlab.h>

#include <drive.h>

/* StringDriveCaller - begin */

class StringDriveCaller : public DriveCaller {
 private:
   char* sEvalStr;
   int iEvalStrLen;

 public:
   StringDriveCaller(const DriveHandler* pDH, const char* const sTmpStr);
   ~StringDriveCaller(void);
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   inline doublereal dGet(void) const;
};


inline doublereal StringDriveCaller::dGet(const doublereal& dVar) const
{
   ((DriveHandler*)DriveCaller::pDrvHdl)->SetVar(dVar);
   
#ifdef DEBUG
   
   /* La variabile temporanea serve solo per il debug; 
    * il codice in effetti avrebbe dovuto essere:
    * 
    *     return DriveCaller::pDrvHdl->dGet(InStr);
    * 
    * Tuttavia si ritiene che la modifica non sia troppo onerosa */
   
   do {
#if defined(HAVE_SSTREAM)
      std::istringstream in(sEvalStr);
#else /* HAVE_STRSTREAM_H */
      istrstream in(sEvalStr);
#endif /* HAVE_STRSTREAM_H */
      InputStream In(in);
      silent_cout("StringDriveCaller::dGet(): "
	<< DriveCaller::pDrvHdl->dGet(In) << std::endl);
   } while (0);
#endif /* DEBUG */

#if defined(HAVE_SSTREAM)
   std::istringstream in(sEvalStr);
#else /* HAVE_STRSTREAM_H */ 
   istrstream in(sEvalStr);
#endif /* HAVE_STRSTREAM_H */
   InputStream In(in);

   return  DriveCaller::pDrvHdl->dGet(In);
}

inline doublereal StringDriveCaller::dGet(void) const
{
#if defined(HAVE_SSTREAM)
   std::istringstream in(sEvalStr);
#else /* HAVE_STRSTREAM_H */
   istrstream in(sEvalStr);
#endif /* HAVE_STRSTREAM_H */
   InputStream In(in);

   return DriveCaller::pDrvHdl->dGet(In);
}

/* StringDriveCaller - end */


/* TimeDriveCaller - begin */

class TimeDriveCaller : public DriveCaller {
 public:
   TimeDriveCaller(const DriveHandler* pDH);
   virtual ~TimeDriveCaller(void);
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   inline doublereal dGet(void) const;
};


inline doublereal TimeDriveCaller::dGet(const doublereal& dVar) const
{
   return dVar;
}


inline doublereal TimeDriveCaller::dGet(void) const
{
   return pDrvHdl->dGetTime();
}

/* TimeDriveCaller - end */


/* LinearDriveCaller - begin */

class LinearDriveCaller : public DriveCaller {
 private:
   doublereal dC0;
   doublereal dC1;
 public:
   LinearDriveCaller(const DriveHandler* pDH, doublereal d0, doublereal d1);
   virtual ~LinearDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal LinearDriveCaller::dGet(const doublereal& dVar) const
{
   return dC0+dC1*dVar;
}

/* LinearDriveCaller - end */


/* ParabolicDriveCaller - begin */

class ParabolicDriveCaller : public DriveCaller {
 private:
   doublereal dC0;
   doublereal dC1;
   doublereal dC2;
   
 public:
   ParabolicDriveCaller(const DriveHandler* pDH,
			doublereal d0, doublereal d1, doublereal d2);
   virtual ~ParabolicDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal ParabolicDriveCaller::dGet(const doublereal& dVar) const
{
   return dC0+dVar*(dC1+dC2*dVar);
}

/* ParabolicDriveCaller - end */


/* CubicDriveCaller - begin */

class CubicDriveCaller : public DriveCaller {
 private:
   doublereal dC0;
   doublereal dC1;
   doublereal dC2;
   doublereal dC3;
   
 public:
   CubicDriveCaller(const DriveHandler* pDH,
		    doublereal d0, doublereal d1,
		    doublereal d2, doublereal d3);
   virtual ~CubicDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal CubicDriveCaller::dGet(const doublereal& dVar) const
{
   return dC0+dVar*(dC1+dVar*(dC2+dC3*dVar));
}

/* CubicDriveCaller - end */


/* StepDriveCaller - begin */

class StepDriveCaller : public DriveCaller {
 private:
   doublereal dStepTime;
   doublereal dStepValue;
   doublereal dInitialValue;
   
 public:
   StepDriveCaller(const DriveHandler* pDH, 
		   doublereal d1, doublereal d2, doublereal d3);
   ~StepDriveCaller(void);

   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal StepDriveCaller::dGet(const doublereal& dVar) const
{
   if(dVar >= dStepTime) {
      return dStepValue;
   }
   /* else if dVar < dStepTime */
   return dInitialValue;
}

/* StepDriveCaller - end */


/* DoubleStepDriveCaller - begin */

class DoubleStepDriveCaller : public DriveCaller {
 private:
   doublereal dStepTime;
   doublereal dStepValue;
   doublereal dEndStepTime;
   doublereal dInitialValue;
   
 public:
   DoubleStepDriveCaller(const DriveHandler* pDH, 
			 doublereal d1, doublereal d2,
			 doublereal d3, doublereal d4);
   ~DoubleStepDriveCaller(void);

   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal DoubleStepDriveCaller::dGet(const doublereal& dVar) const
{
   if(dVar >= dStepTime && dVar < dEndStepTime) {
      return dStepValue;
   }
   /* else if dVar < dStepTime || dVar >= dEndStepTime */
   return dInitialValue;
}

/* DoubleStepDriveCaller - end */


/* RampDriveCaller - begin */

class RampDriveCaller : public DriveCaller {
 private:
   doublereal dSlope;
   doublereal dStartTime;
   doublereal dEndTime;
   doublereal dInitialValue;
   
 public:
   RampDriveCaller(const DriveHandler* pDH, 
		   doublereal d1, doublereal d2, doublereal d3, doublereal d4);
   ~RampDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal RampDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar < dStartTime) {
      return dInitialValue;
   }
   if (dVar > dEndTime) {
      return dInitialValue+dSlope*(dEndTime-dStartTime);
   } /* else if dVar >= dStartTime && dVar < dEndTime */
   return dInitialValue+dSlope*(dVar-dStartTime);
}

/* RampDriveCaller - end */


/* DoubleRampDriveCaller - begin */

class DoubleRampDriveCaller : public DriveCaller {
 private:
   doublereal dAscendingSlope;
   doublereal dAscendingStartTime;
   doublereal dAscendingEndTime;
   doublereal dDescendingSlope;
   doublereal dDescendingStartTime;
   doublereal dDescendingEndTime;
   doublereal dInitialValue;
   
 public:
   DoubleRampDriveCaller(const DriveHandler* pDH,
			 doublereal d1, doublereal d2, doublereal d3, 
			 doublereal d4, doublereal d5, doublereal d6,
			 doublereal d7);
   ~DoubleRampDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;

   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal DoubleRampDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar < dAscendingStartTime) {
      return dInitialValue;
   }
   if (dVar >= dAscendingStartTime && dVar < dAscendingEndTime) {
      return dInitialValue+dAscendingSlope*(dVar-dAscendingStartTime);
   }
   if (dVar >= dAscendingEndTime && dVar < dDescendingStartTime) {
      return dInitialValue
	      +dAscendingSlope*(dAscendingEndTime-dAscendingStartTime);
   }
   if (dVar >= dDescendingStartTime && dVar < dDescendingEndTime) {
      return dInitialValue
	      +dAscendingSlope*(dAscendingEndTime-dAscendingStartTime)
	      +dDescendingSlope*(dVar-dDescendingStartTime);
   } /* else if dVar >= dDescendingEndTime */
   return dInitialValue
	   +dAscendingSlope*(dAscendingEndTime-dAscendingStartTime)
	   +dDescendingSlope*(dDescendingEndTime-dDescendingStartTime);
}

/* DoubleRampDriveCaller - end */


/* SineDriveCaller - begin */

class SineDriveCaller : public DriveCaller {
 private:
   doublereal dStartTime;
   doublereal dOmega;
   doublereal dAmplitude;
   integer iNumCycles;
   doublereal dInitialValue;   
   doublereal dEndTime;
   doublereal dFinalValue;
   flag fNeverEnd;
   
 public:
   SineDriveCaller(const DriveHandler* pDH,
		   doublereal d1, doublereal d2, doublereal d3, 
		   integer iNumCyc, doublereal d4);
   ~SineDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal SineDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar <= dStartTime) {
      return dInitialValue;
   }
   if (fNeverEnd || dVar < dEndTime) {
      return dInitialValue+dAmplitude*sin(dOmega*(dVar-dStartTime));
   }
   /* else if dVar > dEndTime */
   return dFinalValue;
}

/* SineDriveCaller - end */


/* CosineDriveCaller - begin */

class CosineDriveCaller : public DriveCaller {
 private:
   doublereal dStartTime;
   doublereal dOmega;
   doublereal dAmplitude;
   integer iNumCycles;
   doublereal dInitialValue;
   doublereal dEndTime;
   doublereal dFinalValue;
   flag fNeverEnd;
   
 public:
   CosineDriveCaller(const DriveHandler* pDH,
		     doublereal d1, doublereal d2, doublereal d3, 
		     integer iNumCyc, doublereal d4);
   ~CosineDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;

   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal CosineDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar < dStartTime) {
      return dInitialValue;
   }
   if (fNeverEnd || dVar < dEndTime) {
      return dInitialValue+dAmplitude*(1.-cos(dOmega*(dVar-dStartTime)));
   } /* else if dTime > dEndTime */
   return dFinalValue;
}

/* CosineDriveCaller - end */


/* FreqSweepDriveCaller - begin */

class FreqSweepDriveCaller : public DriveCaller {
 private:
   doublereal dStartTime;       
   DriveCaller *pOmega;
   DriveCaller *pAmplitude;   
   doublereal dInitialValue;   
   doublereal dEndTime;
   doublereal dFinalValue;
   flag fNeverEnd;
   
 public:
   FreqSweepDriveCaller(const DriveHandler* pDH,
			doublereal d1, 
			const DriveCaller* pOmega, 
			const DriveCaller* pAmplitude,
			doublereal d2,
			doublereal d3,
			doublereal d4);
   ~FreqSweepDriveCaller(void);
  
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal FreqSweepDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar <= dStartTime) {
      return dInitialValue;
   }
   if (fNeverEnd || dVar < dEndTime) {
      return dInitialValue
	      +pAmplitude->dGet(dVar)*sin(pOmega->dGet(dVar)*(dVar-dStartTime));
   }
   /* else if dVar > dEndTime */
   return dFinalValue;
}

/* FreqSweepDriveCaller - end */


/* ExpDriveCaller - begin */

class ExpDriveCaller : public DriveCaller {
 private:
   doublereal dAmplitude;
   doublereal dTimeConst;
   doublereal dStartTime;
   doublereal dInitialValue;
   
 public:
   ExpDriveCaller(const DriveHandler* pDH, 
		  doublereal dA, doublereal dT, 
		  doublereal dS, doublereal dI);
   virtual ~ExpDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal ExpDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar <= dStartTime) {
      return dInitialValue;
   }
   return dInitialValue
	   +dAmplitude*(1.-exp((dStartTime-dVar)/dTimeConst));
}

/* ExpDriveCaller - end */


/* RandDriveCaller - begin */

class RandDriveCaller : public DriveCaller {
 private:
   long int iBase;
   doublereal dAmplitude;
   doublereal dRefVal;
   doublereal dStartTime;
   doublereal dEndTime;
   integer iSteps;
   integer iRandDriveNumber;
   
 public:
   RandDriveCaller(const DriveHandler* pDH, 
		   doublereal dA, doublereal dR,
		   doublereal dS, doublereal dE, integer iS);
   virtual ~RandDriveCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;
   
   inline doublereal dGet(const doublereal& dVar) const;
   /* inline doublereal dGet(void) const; */
};


inline doublereal RandDriveCaller::dGet(const doublereal& dVar) const
{
   if (dVar < dStartTime || dVar > dEndTime) {
      return dRefVal;
   }
   
   doublereal dRand = doublereal( (((unsigned long)pDrvHdl->iGetRand(iRandDriveNumber))
				   +iBase )%RAND_MAX);
   return dRefVal+dAmplitude*(2.*dRand/doublereal(RAND_MAX)-1.);
}

/* RandDriveCaller - end */


/* PiecewiseLinearDriveCaller - begin */

class PiecewiseLinearDriveCaller : public DriveCaller {
private:
	unsigned int iNumPoints;
	doublereal *pPoints;
	doublereal *pVals;

public:
	PiecewiseLinearDriveCaller(const DriveHandler* pDH,
			unsigned int i, doublereal *p);
	virtual ~PiecewiseLinearDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;
	
	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	inline doublereal dGet(const doublereal& dVar) const;
};

inline doublereal
PiecewiseLinearDriveCaller::dGet(const doublereal& dVar) const
{
	if (dVar <= pPoints[0]) {
		return pVals[0];
	}

	if (dVar >= pPoints[iNumPoints-1]) {
		return pVals[iNumPoints-1];
	}

	for (unsigned int i = 1; i < iNumPoints; i++) {
		if (dVar == pPoints[i]) {
			return pVals[i];
		}
		if (dVar < pPoints[i]) {
			doublereal dx = pPoints[i]-pPoints[i-1];
			return ((dVar-pPoints[i-1])*pVals[i]
				+(pPoints[i]-dVar)*pVals[i-1])/dx;
		}
	}

	throw ErrGeneric();
}

/* PiecewiseLinearDriveCaller - end */

/* DriveArrayCaller - begin */

class DriveArrayCaller : public DriveCaller {
 private:
   unsigned short int iNumDrivers;
   const DriveCaller** ppDriveCallers;
   
 public:
   DriveArrayCaller(const DriveHandler* pDH, 
		    unsigned short int i, 
		    const DriveCaller** ppDC);
   virtual ~DriveArrayCaller(void);
   
   /* Copia */
   virtual DriveCaller* pCopy(void) const;
   
   /* Scrive il contributo del DriveCaller al file di restart */   
   virtual std::ostream& Restart(std::ostream& out) const;

   inline doublereal dGet(const doublereal& dVar) const;
   inline doublereal dGet(void) const;
};


inline doublereal DriveArrayCaller::dGet(const doublereal& dVar) const
{
   doublereal d = 0.;
   for (int i = 0; i < iNumDrivers; i++) {
      ASSERT(ppDriveCallers[i] != NULL);
      d += ppDriveCallers[i]->dGet(dVar);
   }
   return d;
}


inline doublereal DriveArrayCaller::dGet(void) const
{
   doublereal d = 0.;
   for (int i = 0; i < iNumDrivers; i++) {
      ASSERT(ppDriveCallers[i] != NULL);
      d += ppDriveCallers[i]->dGet();
   }
   return d;
}

/* DriveArrayCaller - end */

#endif /* DRIVE__H */

