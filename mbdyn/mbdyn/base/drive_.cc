/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>
#include <cfloat>

#ifdef USE_MPI
#include "mbcomm.h"
#endif /* USE_MPI */

#include "dataman.h"
#include "drive_.h"
#include "dofdrive.h"
#include "privdrive.h"
#include "ddrive.h"
#include "shdrive.h"

#include "filedrv.h"
#include "fixedstep.h"
#include "varstep.h"
#include "sockdrv.h"
#include "streamdrive.h"
#include "socketstreamdrive.h"
#include "bufferstreamdrive.h"

#ifdef USE_GINAC
#include "ginacdrive.h"
#endif // USE_GINAC

#ifdef STATIC_MODULES
#ifdef USE_OCTAVE
#include "module-octave/module-octave.h"
#endif // USE_OCTAVE
#include "module-nodedistdrive/module-nodedistdrive.h"
#include "module-minmaxdrive/module-minmaxdrive.h"
#endif // STATIC_MODULES

/* StringDriveCaller - begin */

#ifdef USE_EE
StringDriveCaller::StringDriveCaller(const DriveHandler* pDH,
	const std::string& sTmpStr, const ExpressionElement *expr)
: DriveCaller(pDH), sEvalStr(sTmpStr), m_expr(new SharedExpr(expr))
{
	ASSERT(!sEvalStr.empty());
	ASSERT(m_expr->Get() != 0);
}


StringDriveCaller::StringDriveCaller(const DriveHandler* pDH,
	const std::string& sTmpStr, SharedExprPtr_t expr)
: DriveCaller(pDH), sEvalStr(sTmpStr), m_expr(expr)
{
	ASSERT(!sEvalStr.empty());
	ASSERT(m_expr->Get() != 0);
}


StringDriveCaller::~StringDriveCaller(void)
{
	ASSERT(m_expr->Get() != 0);
}


/* Copia */
DriveCaller *
StringDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
		StringDriveCaller,
		StringDriveCaller(pDrvHdl, sEvalStr, m_expr->pCopy()));

	return pDC;
}

#else // ! USE_EE

StringDriveCaller::StringDriveCaller(const DriveHandler* pDH,
	const std::string& sTmpStr)
: DriveCaller(pDH), sEvalStr(sTmpStr)
{
	ASSERT(!sEvalStr.empty());
}


StringDriveCaller::~StringDriveCaller(void)
{
	NO_OP;
}

DriveCaller * 
StringDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		StringDriveCaller,
		StringDriveCaller(pDrvHdl, sEvalStr));

	return pDC;
}

#endif // ! USE_EE


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
StringDriveCaller::Restart(std::ostream& out) const
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
DriveCaller *
TimeDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDrvHdl));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
TimeDriveCaller::Restart(std::ostream& out) const
{
	return out << "time";
}

/* TimeDriveCaller - end */


/* TimeStepDriveCaller - begin */

TimeStepDriveCaller::TimeStepDriveCaller(const DriveHandler* pDH)
: DriveCaller(pDH)
{
	NO_OP;
}

TimeStepDriveCaller::~TimeStepDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
TimeStepDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, TimeStepDriveCaller, TimeStepDriveCaller(pDrvHdl));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
TimeStepDriveCaller::Restart(std::ostream& out) const
{
	return out << "timestep";
}

/* TimeStepDriveCaller - end */


/* MultDriveCaller - begin */

MultDriveCaller::MultDriveCaller(const DriveHandler *pDH,
	const DriveCaller* pDC1, const DriveCaller *pDC2)
: DriveCaller(pDH),
DO1(pDC1), DO2(pDC2)
{
	NO_OP;
}

MultDriveCaller::~MultDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
MultDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, MultDriveCaller,
		MultDriveCaller(pDrvHdl,
			DO1.pGetDriveCaller()->pCopy(), DO2.pGetDriveCaller()->pCopy()));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
MultDriveCaller::Restart(std::ostream& out) const
{
	return out << "mult, ",
		DO1.pGetDriveCaller()->Restart(out) << ", ";
		DO2.pGetDriveCaller()->Restart(out);
}

/* MultDriveCaller - end */


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
DriveCaller*
ConstDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(dConst));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
ConstDriveCaller::Restart(std::ostream& out) const
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
DriveCaller *
LinearDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, LinearDriveCaller, LinearDriveCaller(pDrvHdl, dC0, dC1));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
LinearDriveCaller::Restart(std::ostream& out) const
{
	return out << " linear, " << dC0 << ", " << dC1;
}

/* LinearDriveCaller - end */


/* ParabolicDriveCaller - begin */

ParabolicDriveCaller::ParabolicDriveCaller(const DriveHandler* pDH,
	doublereal d0, doublereal d1, doublereal d2)
: DriveCaller(pDH), dC0(d0), dC1(d1), dC2(d2)
{
	NO_OP;
}

ParabolicDriveCaller::~ParabolicDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
ParabolicDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, ParabolicDriveCaller, ParabolicDriveCaller(pDrvHdl, dC0, dC1, dC2));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
ParabolicDriveCaller::Restart(std::ostream& out) const
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

CubicDriveCaller::~CubicDriveCaller(void)
{
   NO_OP;
}

/* Copia */
DriveCaller *
CubicDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, CubicDriveCaller, CubicDriveCaller(pDrvHdl, dC0, dC1, dC2, dC3));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
CubicDriveCaller::Restart(std::ostream& out) const
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
DriveCaller *
StepDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		StepDriveCaller,
		StepDriveCaller(pDrvHdl, dStepTime, dStepValue, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
StepDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " step, " << dStepTime
		<< ", " << dStepValue
		<< ", " << dInitialValue;
}

/* StepDriveCaller - end */


/* DoubleStepDriveCaller - begin */

DoubleStepDriveCaller::DoubleStepDriveCaller(const DriveHandler* pDH,
	doublereal d1, doublereal d2, doublereal d3, doublereal d4)
: DriveCaller(pDH), dStepTime(d1), dStepValue(d2),
dEndStepTime(d3), dInitialValue(d4)
{
	NO_OP;
}

DoubleStepDriveCaller::~DoubleStepDriveCaller(void)
{
	NO_OP;
}


/* Copia */
DriveCaller *
DoubleStepDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		DoubleStepDriveCaller,
		DoubleStepDriveCaller(pDrvHdl, dStepTime, dStepValue, dEndStepTime, dInitialValue));
	return pDC;
}


/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
DoubleStepDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " double step, " << dStepTime
		<< ", " << dEndStepTime
		<< ", " << dStepValue
		<< ", " << dInitialValue;
}

/* DoubleStepDriveCaller - end */


/* RampDriveCaller - begin */

RampDriveCaller::RampDriveCaller(const DriveHandler* pDH,
	doublereal d1, doublereal d2, doublereal d3, doublereal d4)
: DriveCaller(pDH), dSlope(d1), dStartTime(d2), dEndTime(d3), dInitialValue(d4)
{
	NO_OP;
}

RampDriveCaller::~RampDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
RampDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		RampDriveCaller,
		RampDriveCaller(pDrvHdl, dSlope, dStartTime, dEndTime, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
RampDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< "ramp, " << dSlope
		<< ", " << dStartTime
		<< ", " << dEndTime
		<< ", " << dInitialValue;
}

/* RampDriveCaller - end */


/* DoubleRampDriveCaller - begin */

DoubleRampDriveCaller::DoubleRampDriveCaller(const DriveHandler* pDH,
	doublereal d1, doublereal d2, doublereal d3,
	doublereal d4, doublereal d5, doublereal d6,
	doublereal d7)
: DriveCaller(pDH),
dAscendingSlope(d1), dAscendingStartTime(d2), dAscendingEndTime(d3),
dDescendingSlope(d4), dDescendingStartTime(d5), dDescendingEndTime(d6),
dInitialValue(d7)
{
	NO_OP;
}

DoubleRampDriveCaller::~DoubleRampDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
DoubleRampDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		DoubleRampDriveCaller,
		DoubleRampDriveCaller(pDrvHdl, dAscendingSlope, dAscendingStartTime, dAscendingEndTime,
			dDescendingSlope, dDescendingStartTime, dDescendingEndTime, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
DoubleRampDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " double ramp, " << dAscendingSlope
		<< ", " << dAscendingStartTime
		<< ", " << dAscendingEndTime
		<< ", " << dDescendingSlope
		<< ", " << dDescendingStartTime
		<< ", " << dDescendingEndTime
		<< ", " << dInitialValue;
}

/* DoubleRampDriveCaller - end */


/* SineDriveCaller - begin */

SineDriveCaller::SineDriveCaller(const DriveHandler* pDH,
	doublereal d1, doublereal d2, doublereal d3,
	integer iNumCyc, doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), dOmega(d2), dAmplitude(d3),
iNumCycles(iNumCyc), dInitialValue(d4), bNeverEnd(false)
{
	/* Onde di seno che partono da zero ed arrivano a0 */
	if (iNumCycles > 0) {
		dEndTime = dStartTime + 2.*M_PI/dOmega*(doublereal(iNumCycles) - .5);
		dFinalValue = dInitialValue;

	/* Onde di seno che continuano all'infinito */
	} else if (iNumCycles == 0) {
		dEndTime = 0.;
		bNeverEnd = true;

	/* Onde di seno che partono da 0 ed arrivano ad 1
	 * con tangente orizzontale */
	} else {
		dEndTime = dStartTime + 2.*M_PI/dOmega*(doublereal(-iNumCycles) - .75);
		dFinalValue = dInitialValue + dAmplitude;
	}
}

SineDriveCaller::~SineDriveCaller(void)
{
   NO_OP;
}

/* Copia */
DriveCaller *
SineDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		SineDriveCaller,
		SineDriveCaller(pDrvHdl, dStartTime, dOmega, dAmplitude, iNumCycles, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
SineDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " sine, " << dStartTime
		<< ", " << dOmega
		<< ", " << dAmplitude
		<< ", " << iNumCycles
		<< ", " << dInitialValue;
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
DriveCaller *
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
	out
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

	return out;
}

/* FourierSeriesDriveCaller - end */


/* CosineDriveCaller - begin */

CosineDriveCaller::CosineDriveCaller(const DriveHandler* pDH,
	doublereal d1, doublereal d2,
	doublereal d3,
	integer iNumCyc, doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), dOmega(d2), dAmplitude(d3),
iNumCycles(iNumCyc), dInitialValue(d4), bNeverEnd(false)
{
	/* Onde di coseno che partono da zero ed arrivano a 0 */
	if (iNumCycles > 0) {
		dEndTime = dStartTime + 2.*M_PI/dOmega*doublereal(iNumCycles);
		dFinalValue = dInitialValue;

	/* Onde di coseno che continuano all'infinito */
	} else if (iNumCycles == 0) {
		dEndTime = 0.;
		bNeverEnd = true;

	/* Onde di coseno che partono da 0 ed arrivano ad 1
	 * con tangente orizzontale */
	} else {
		dEndTime = dStartTime + 2.*M_PI/dOmega*(doublereal(-iNumCycles) - .5);
		dFinalValue = dInitialValue + 2.*dAmplitude;
	}
}

CosineDriveCaller::~CosineDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller*
CosineDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		CosineDriveCaller,
		CosineDriveCaller(pDrvHdl, dStartTime, dOmega, dAmplitude, iNumCycles, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
CosineDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " cosine, " << dStartTime
		<< ", " << dOmega
		<< ", " << dAmplitude
		<< ", " << iNumCycles
		<< ", " << dInitialValue;
}

/* CosineDriveCaller - end */


/* TanhDriveCaller - begin */

TanhDriveCaller::TanhDriveCaller(const DriveHandler* pDH,
	doublereal ds, doublereal da,
	doublereal db, doublereal di)
: DriveCaller(pDH),
dStart(ds), dA(da), dB(db), dInitialValue(di)
{
	NO_OP;
}

TanhDriveCaller::~TanhDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller*
TanhDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		TanhDriveCaller,
		TanhDriveCaller(pDrvHdl, dStart, dA, dB, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
TanhDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " tanh, " << dStart
		<< ", " << dA
		<< ", " << dB
		<< ", " << dInitialValue;
}

/* TanhDriveCaller - end */


/* FreqSweepDriveCaller - begin */

FreqSweepDriveCaller::FreqSweepDriveCaller(const DriveHandler* pDH,
	doublereal d1,
	const DriveCaller* pOmega,
	const DriveCaller* pAmplitude,
	doublereal d2,
	doublereal d3,
	doublereal d4)
: DriveCaller(pDH),
dStartTime(d1), pOmega(pOmega),
pAmplitude(pAmplitude),
dInitialValue(d2), dEndTime(d3), dFinalValue(d4),
bNeverEnd(false)
{
	ASSERT(pOmega != 0);
	ASSERT(pAmplitude != 0);

	if (dEndTime <= dStartTime) {
		bNeverEnd = true;
	}
}

FreqSweepDriveCaller::~FreqSweepDriveCaller(void)
{
	if (pOmega) {
		SAFEDELETE(pOmega);
	}

	if (pAmplitude) {
		SAFEDELETE(pAmplitude);
	}
}

/* Copia */
DriveCaller *
FreqSweepDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		FreqSweepDriveCaller,
		FreqSweepDriveCaller(pDrvHdl, dStartTime, pOmega->pCopy(),
			pAmplitude->pCopy(), dInitialValue, dEndTime, dFinalValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
FreqSweepDriveCaller::Restart(std::ostream& out) const
{
	out
		<< " frequency sweep, " << dStartTime
		<< ", ", pOmega->Restart(out)
		<< ", ", pAmplitude->Restart(out)
		<< ", " << dInitialValue
		<< ", " << dEndTime
		<< ", " << dFinalValue;
	return out;
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
DriveCaller*
ExpDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		ExpDriveCaller,
		ExpDriveCaller(pDrvHdl, dAmplitude, dTimeConst, dStartTime, dInitialValue));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
ExpDriveCaller::Restart(std::ostream& out) const
{
	return out
		<< " exponential, " << dAmplitude
		<< ", " << dTimeConst
		<< ", " << dStartTime
		<< ", " << dInitialValue;
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
	iRandDriveNumber = pDrvHdl->iRandInit(iSteps);
}

RandDriveCaller::~RandDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
RandDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		RandDriveCaller,
		RandDriveCaller(pDrvHdl, dAmplitude, dRefVal, dStartTime, dEndTime, iSteps));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
RandDriveCaller::Restart(std::ostream& out) const
{
	out
		<< " random, " << dAmplitude
		<< ", " << dRefVal
		<< ", " << dStartTime
		<< ", " << dEndTime;
	if (iSteps > 1) {
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
	iMeterDriveNumber = pDrvHdl->iMeterInit(iSteps);
}

MeterDriveCaller::~MeterDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
MeterDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		MeterDriveCaller,
		MeterDriveCaller(pDrvHdl, dStartTime, dEndTime, iSteps));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
MeterDriveCaller::Restart(std::ostream& out) const
{
	out
		<< " meter, " << dStartTime
		<< ", " << dEndTime;
	if (iSteps > 1) {
		out << ", steps, " << iSteps;
	}
	return out;
}

/* MeterDriveCaller - end */


/* ClosestNextDriveCaller - begin */

ClosestNextDriveCaller::ClosestNextDriveCaller(const DriveHandler* pDH,
	doublereal dS, doublereal dE, const DriveCaller *pIncrement)
: DriveCaller(pDH),
dStartTime(dS), dEndTime(dE),
pIncrement(pIncrement)
{
	iDriveNumber = pDrvHdl->iClosestNextInit(pIncrement, dStartTime);
}

ClosestNextDriveCaller::~ClosestNextDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
ClosestNextDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		ClosestNextDriveCaller,
		ClosestNextDriveCaller(pDrvHdl, dStartTime, dEndTime,
		pIncrement->pCopy()));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
ClosestNextDriveCaller::Restart(std::ostream& out) const
{
	out
		<< " closest next, " << dStartTime
		<< ", " << dEndTime
		<< ", ", pIncrement->Restart(out);

	return out;
}

/* ClosestNextDriveCaller - end */


/* DirectDriveCaller - begin */

DirectDriveCaller::DirectDriveCaller(const DriveHandler* pDH)
: DriveCaller(pDH)
{
	NO_OP;
}

DirectDriveCaller::~DirectDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
DirectDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		DirectDriveCaller,
		DirectDriveCaller(pDrvHdl));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
DirectDriveCaller::Restart(std::ostream& out) const
{
	return out << " direct";
}

/* DirectDriveCaller - end */


/* PiecewiseLinearDriveCaller - begin */

PiecewiseLinearDriveCaller::PiecewiseLinearDriveCaller(const DriveHandler* pDH,
	unsigned int i, doublereal *p)
: DriveCaller(pDH), iNumPoints(i), pPoints(p), pVals(p + i)
{
	ASSERT(i >= 2);
	ASSERT(p != 0);
}

PiecewiseLinearDriveCaller::~PiecewiseLinearDriveCaller(void)
{
	if (pPoints != 0) {
		SAFEDELETEARR(pPoints);
	}
}

/* Copia */
DriveCaller *
PiecewiseLinearDriveCaller::pCopy(void) const
{
	doublereal *p = 0;
	SAFENEWARR(p, doublereal, 2*iNumPoints);
	for (unsigned int i = 0; i < 2*iNumPoints; i++) {
		p[i] = pPoints[i];
	}

	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
			PiecewiseLinearDriveCaller,
			PiecewiseLinearDriveCaller(pDrvHdl, iNumPoints, p));

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
PiecewiseLinearDriveCaller::Restart(std::ostream& out) const
{
	out << "piecewise linear, " << iNumPoints;

	for (unsigned int i = 0; i < iNumPoints; i++) {
		out << ", " << pPoints[i] << ", " << pVals[i];
	}

	return out;
}

/* PiecewiseLinearDriveCaller - end */


/* DriveArrayCaller - begin */

DriveArrayCaller::DriveArrayCaller(const DriveHandler* pDH, dcv_t& DC)
: DriveCaller(pDH), m_dc(DC)
{
#ifdef DEBUG
	ASSERT(!m_dc.empty());
	for (dcv_t::const_iterator i = m_dc.begin(); i != m_dc.end(); ++i) {
		ASSERT((*i) != 0);
	}
#endif /* DEBUG */
}

DriveArrayCaller::~DriveArrayCaller(void)
{
	ASSERT(!m_dc.empty());

	for (dcv_t::iterator i = m_dc.begin(); i != m_dc.end(); ++i) {
		SAFEDELETE(*i);
	}
}

/* Copia */
DriveCaller *
DriveArrayCaller::pCopy(void) const
{
	ASSERT(!m_dc.empty());

	dcv_t DC(m_dc.size());

	for (unsigned i = 0; i < m_dc.size(); i++) {
		DC[i] = m_dc[i]->pCopy();
	}

	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		DriveArrayCaller,
		DriveArrayCaller(pDrvHdl, DC));

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
DriveArrayCaller::Restart(std::ostream& out) const
{
	out << " array, " << m_dc.size();

	for (dcv_t::const_iterator i = m_dc.begin(); i != m_dc.end(); ++i) {
		ASSERT((*i) != 0);

		out << ", ", (*i)->Restart(out);
	}

	return out;
}

/* DriveArrayCaller - end */

/* PeriodicDriveCaller - begin */

PeriodicDriveCaller::PeriodicDriveCaller(const DriveHandler *pDH,
	const DriveCaller* pDC, doublereal dT0, doublereal dPeriod)
: DriveCaller(pDH),
DO(pDC), dT0(dT0), dPeriod(dPeriod)
{
	NO_OP;
}

PeriodicDriveCaller::~PeriodicDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller *
PeriodicDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, PeriodicDriveCaller,
		PeriodicDriveCaller(pDrvHdl,
			DO.pGetDriveCaller()->pCopy(), dT0, dPeriod));
	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
PeriodicDriveCaller::Restart(std::ostream& out) const
{
	return out << "periodic, " << dT0 << ", " << dPeriod << ", ",
		DO.pGetDriveCaller()->Restart(out);
}

/* PeriodicDriveCaller - end */

/* PostponedDriveCaller - begin */

void
PostponedDriveCaller::Check(void) const
{
	if (!DO.pGetDriveCaller()) {
		const DriveCaller *pDC = HP.GetDrive(uDriveLabel);
		if (pDC == 0) {
			silent_cerr("PostponedDriveCaller: unable to resolve postponed drive caller \"" << uDriveLabel << "\"" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		DO.Set(pDC->pCopy());
		const_cast<PostponedDriveCaller *>(this)->SetDrvHdl(DO.pGetDriveCaller()->pGetDrvHdl());
	}
}

PostponedDriveCaller::PostponedDriveCaller(MBDynParser& HP, unsigned uLabel)
: DriveCaller(0),
HP(HP), uDriveLabel(uLabel), DO(0)
{
	NO_OP;
}

PostponedDriveCaller::~PostponedDriveCaller(void)
{
	NO_OP;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
PostponedDriveCaller::Restart(std::ostream& out) const
{
	Check();

	out << ", deferred, " << uDriveLabel << " /* actual drive: ",
	    DO.pGetDriveCaller()->Restart(out) << " */";

	return out;
}

/* Copia */
DriveCaller *
PostponedDriveCaller::pCopy(void) const
{
	DriveCaller *pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC, PostponedDriveCaller,
		PostponedDriveCaller(HP, uDriveLabel));
	return pDC;
}

/* PostponedDriveCaller - end */


/* bag that contains functions to parse drives */

typedef std::map<std::string, DriveRead *, ltstrcase> DriveFuncMapType;
static DriveFuncMapType DriveFuncMap;

struct DriveWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return DriveFuncMap.find(s) != DriveFuncMap.end();
	};
};

static DriveWordSetType DriveWordSet;

bool
SetDriveData(const std::string& name, DriveRead *rf)
{
	pedantic_cout("registering drive \"" << name << "\"" << std::endl);
	return DriveFuncMap.insert(DriveFuncMapType::value_type(name, rf)).second;
}

/* Reads drives */

Drive *
ReadDriveData(unsigned uLabel, const DataManager* pDM, MBDynParser& HP)
{
	DEBUGCOUTFNAME("ReadDriveData()");

	const char *s = HP.IsWord(DriveWordSet);
	if (s == 0) {
		silent_cerr("ReadDriveData(" << uLabel << "): "
			"unknown drive type "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveFuncMapType::iterator func = DriveFuncMap.find(std::string(s));
	if (func == DriveFuncMap.end()) {
		silent_cerr("unknown drive type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(uLabel, pDM, HP);
}

DriveRead::~DriveRead(void)
{
	NO_OP;
}

/* bag that contains functions to parse drive callers */

typedef std::map<std::string, DriveCallerRead *, ltstrcase> DriveCallerFuncMapType;
static DriveCallerFuncMapType DriveCallerFuncMap;

struct DriveCallerWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return DriveCallerFuncMap.find(s) != DriveCallerFuncMap.end();
	};
};

static DriveCallerWordSetType DriveCallerWordSet;

bool
SetDriveCallerData(const char *name, DriveCallerRead *rf)
{
	pedantic_cout("registering drive caller \"" << name << "\"" << std::endl);
	return DriveCallerFuncMap.insert(DriveCallerFuncMapType::value_type(name, rf)).second;
}

/* Reads drive callers */

DriveCaller *
ReadDriveCallerData(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	DEBUGCOUTFNAME("ReadDriveCallerData()");

	const char *s = HP.IsWord(DriveCallerWordSet);
	if (s == 0) {
		pedantic_cerr("ReadDriveCallerData: warning, assuming \"const\" drive caller at line " << HP.GetLineData() << std::endl);
		s = "const";
	}

	DriveCallerFuncMapType::iterator func = DriveCallerFuncMap.find(std::string(s));
	if (func == DriveCallerFuncMap.end()) {
		silent_cerr("unknown drive caller type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP, bDeferred);
}

void
DriveCallerRead::NeedDM(const DataManager* pDM, MBDynParser& HP, bool bDeferred,
	const char *const name)
{
	if (pDM == 0 && !bDeferred) {
		silent_cerr("\"" << name << "\" drive caller needs data manager "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrNeedDataManager(MBDYN_EXCEPT_ARGS);
	}
}

DriveCallerRead::~DriveCallerRead(void)
{
	NO_OP;
}

void DriveCallerRead::ReadOutput(DriveCaller* pDC, const DataManager* pDM, MBDynParser& HP)
{
	flag fOutput = pDM->bOutputDriveCallers()
					? DriveCaller::OUTPUT_VALUE
					: 0;

	if (HP.IsKeyWord("output")) {
		bool bGotKeyWord = false;

		if (HP.IsKeyWord("value")) {
			bGotKeyWord = true;

			if (HP.GetYesNoOrBool()) {
				fOutput |= DriveCaller::OUTPUT_VALUE;
			} else {
				fOutput &= ~DriveCaller::OUTPUT_VALUE;
       			}
		}

		if (HP.IsKeyWord("derivative")) {
			bGotKeyWord = true;
            
			if (HP.GetYesNoOrBool()) {
				fOutput |= DriveCaller::OUTPUT_DERIVATIVE;
			} else {
				fOutput &= ~DriveCaller::OUTPUT_DERIVATIVE;
			}
		}

		if (!bGotKeyWord) {
			if (HP.GetYesNoOrBool()) {
				fOutput |= DriveCaller::OUTPUT_VALUE;

				if (pDC->bIsDifferentiable()) {
					fOutput |= DriveCaller::OUTPUT_DERIVATIVE;
				}
			} else {
				fOutput = flag(0);
			}
		}
	}

	if (fOutput & DriveCaller::OUTPUT_DERIVATIVE) {
		if (!pDC->bIsDifferentiable()) {
			silent_cerr("warning invalid output flag: drive caller(" << pDC->GetLabel()
					<< ") is not differentiable at line "
					<< HP.GetLineData() << std::endl);

			fOutput &= ~DriveCaller::OUTPUT_DERIVATIVE;
		}
	}

	if (fOutput) {
		fOutput |= DriveCaller::OUTPUT;
	}
    
	pDC->SetOutputFlag(fOutput);

	flag fTrace = 0;

	if (HP.IsKeyWord("output" "trace")) {
		bool bGotKeyWord = false;

		if (HP.IsKeyWord("value")) {
			bGotKeyWord = true;

			if (HP.GetYesNoOrBool()) {
				fTrace |= DriveCaller::TRACE_VALUE;
       			} else {
				fTrace &= ~DriveCaller::TRACE_VALUE;
			}
		}

		if (HP.IsKeyWord("derivative")) {
			bGotKeyWord = true;

			if (HP.GetYesNoOrBool()) {
				fTrace |= DriveCaller::TRACE_DERIVATIVE;
			} else {
				fTrace &= ~DriveCaller::TRACE_DERIVATIVE;
			}
		}

		if (!bGotKeyWord) {
			if (HP.GetYesNoOrBool()) {
				fTrace |= DriveCaller::TRACE_VALUE;

				if (pDC->bIsDifferentiable()) {
					fTrace |= DriveCaller::TRACE_DERIVATIVE;
				}
			} else {
				fTrace = flag(0);
			}
		}
	}

	if (fTrace & DriveCaller::TRACE_DERIVATIVE) {
		if (!pDC->bIsDifferentiable()) {
			silent_cerr("warning invalid trace flag: drive caller (" << pDC->GetLabel()
					<< ") is not differentiable at line "
					<< HP.GetLineData() << std::endl);

			fTrace &= ~DriveCaller::TRACE_DERIVATIVE;
		}
	}

	if (fTrace) {
		fTrace |= DriveCaller::TRACE;
	}
    
	pDC->SetTraceFlag(fTrace);
}

struct TimeDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
TimeDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "time");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC, TimeDriveCaller, TimeDriveCaller(pDrvHdl));

	return pDC;
}

struct TimeStepDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
TimeStepDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "timestep");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC, TimeStepDriveCaller, TimeStepDriveCaller(pDrvHdl));

	return pDC;
}

struct MultDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
MultDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "mult");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller* pDC1 = HP.GetDriveCaller();
	DriveCaller* pDC2 = HP.GetDriveCaller();

	DriveCaller *pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		MultDriveCaller,
		MultDriveCaller(pDrvHdl, pDC1, pDC2));

	return pDC;
}

struct NullDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* /* pDM */ , MBDynParser& /* HP */ , bool /* bDeferred */ );
};

DriveCaller *
NullDCR::Read(const DataManager* /* pDM */ , MBDynParser& /* HP */ , bool /* bDeferred */ )
{
	DriveCaller *pDC = 0;

        SAFENEW(pDC, NullDriveCaller);

	return pDC;
}

struct OneDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* /* pDM */ , MBDynParser& /* HP */ , bool /* bDeferred */ );
};

DriveCaller *
OneDCR::Read(const DataManager* /* pDM */ , MBDynParser& /* HP */ , bool /* bDeferred */ )
{
	DriveCaller *pDC = 0;

        SAFENEW(pDC, OneDriveCaller);

	return pDC;
}

struct ConstDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool /* bDeferred */ );
};

DriveCaller *
ConstDCR::Read(const DataManager* pDM, MBDynParser& HP, bool /* bDeferred */ )
{
	DriveCaller *pDC = 0;

	/* lettura dei dati specifici */
	doublereal dConst = HP.GetReal();
	DEBUGCOUT("Const value: " << dConst << std::endl);

	/* allocazione e creazione */
	if (dConst == 0.) {
		SAFENEW(pDC, NullDriveCaller);

	} else if (dConst == 1.) {
		SAFENEW(pDC, OneDriveCaller);

	} else {
		SAFENEWWITHCONSTRUCTOR(pDC, ConstDriveCaller, ConstDriveCaller(dConst));
	}

	return pDC;
}

struct LinearDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
LinearDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "linear");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* lettura dei dati specifici */
	doublereal dC0 = HP.GetReal();
	DEBUGCOUT("Offset: " << dC0 << std::endl);

	doublereal dC1 = HP.GetReal();
	DEBUGCOUT("Slope: " << dC1 << std::endl);

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		LinearDriveCaller,
		LinearDriveCaller(pDrvHdl, dC0, dC1));

	return pDC;
}

struct ParabolicDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ParabolicDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "parabolic");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

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

	return pDC;
}

struct CubicDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
CubicDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "cubic");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

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

	return pDC;
}

/*
 * this allows each drive to be preceded by the keyword "function"
 */
struct FunctionDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
FunctionDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	return HP.GetDriveCaller(bDeferred);
}

struct StepDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
StepDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "step");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	doublereal dStepTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dStepTime << std::endl);

	doublereal dStepValue = HP.GetReal(1.);
	DEBUGCOUT("Step Value: " << dStepValue << std::endl);

	doublereal dInitialValue = HP.GetReal();
	DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);

	SAFENEWWITHCONSTRUCTOR(pDC,
		StepDriveCaller,
		StepDriveCaller(pDrvHdl, dStepTime, dStepValue, dInitialValue));

	return pDC;
}

struct DoubleStepDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
DoubleStepDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "double step");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	doublereal dStepTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dStepTime << std::endl);

	doublereal dEndStepTime = HP.GetReal();
	DEBUGCOUT("Final time: " << dEndStepTime << std::endl);

	if (dEndStepTime <= dStepTime) {
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

	return pDC;
}

struct RampDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
RampDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "ramp");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Rampa saturata */
	doublereal dSlope = HP.GetReal(1.);
	DEBUGCOUT("Slope Value: " << dSlope << std::endl);

	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	doublereal dFinalTime = std::numeric_limits<double>::max();
	if (!HP.IsKeyWord("forever")) {
		dFinalTime = HP.GetReal();
	}
	DEBUGCOUT("Final time: " << dFinalTime << std::endl);

	if (dFinalTime <= dInitialTime) {
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
		RampDriveCaller(pDrvHdl, dSlope, dInitialTime,
			dFinalTime, dInitialValue));

	return pDC;
}

struct DoubleRampDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
DoubleRampDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "double ramp");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Rampa doppia */
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
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dDescendingSlope = HP.GetReal(-1.);
	DEBUGCOUT("Descending Slope Value: " << dDescendingSlope << std::endl);

	doublereal dDescendingInitialTime = std::numeric_limits<doublereal>::max();
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
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dDescendingFinalTime;
	if (HP.IsKeyWord("forever")) {
		dDescendingFinalTime = std::numeric_limits<doublereal>::max();

	} else {
		dDescendingFinalTime = HP.GetReal();
	}

	DEBUGCOUT("Descending Final time: " << dDescendingFinalTime << std::endl);

	if (dDescendingFinalTime <= dDescendingInitialTime) {
		silent_cerr("Warning at line "
			<< HP.GetLineData() << ": descending final time "
			<< dDescendingFinalTime
			<< " is less than descending initial time "
			<< dDescendingInitialTime
			<< " in double ramp func drive" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dInitialValue = HP.GetReal();
	DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);

	SAFENEWWITHCONSTRUCTOR(pDC,
		DoubleRampDriveCaller,
		DoubleRampDriveCaller(pDrvHdl,
			dAscendingSlope, dAscendingInitialTime, dAscendingFinalTime,
			dDescendingSlope, dDescendingInitialTime, dDescendingFinalTime,
			dInitialValue));

	return pDC;
}

struct SineCosineDCR {
protected:
	virtual ~SineCosineDCR(void) { NO_OP; };
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred, bool bSine);
};

DriveCaller *
SineCosineDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred, bool bSine)
{
	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Seno e coseno (limitati, illimitati, saturati) */
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

	if (bSine) {
		SAFENEWWITHCONSTRUCTOR(pDC,
			SineDriveCaller,
			SineDriveCaller(pDrvHdl, dInitialTime, dOmega, dAmplitude, iNumCycles, dInitialValue));

	} else {
		SAFENEWWITHCONSTRUCTOR(pDC,
			CosineDriveCaller,
			CosineDriveCaller(pDrvHdl, dInitialTime, dOmega, dAmplitude, iNumCycles, dInitialValue));
	}

	return pDC;
}

struct SineDCR : public DriveCallerRead, protected SineCosineDCR {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
SineDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "sine");

	return SineCosineDCR::Read(pDM, HP, bDeferred, true);
}

struct CosineDCR : public DriveCallerRead, protected SineCosineDCR {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
CosineDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "cosine");

	return SineCosineDCR::Read(pDM, HP, bDeferred, false);
}

struct TanhDCR : public DriveCallerRead {
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
TanhDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "tanh");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Seno e coseno (limitati, illimitati, saturati) */
	doublereal dStart = HP.GetReal();
	DEBUGCOUT("Initial time: " << dStart << std::endl);

	doublereal dA = HP.GetReal(1.);
	DEBUGCOUT("Amplitude: " << dA << std::endl);

	doublereal dB = HP.GetReal();
	DEBUGCOUT("Slope: " << dB << std::endl);

	doublereal dInitialValue = HP.GetReal();
	DEBUGCOUT("InitialValue: " << dInitialValue << std::endl);

	SAFENEWWITHCONSTRUCTOR(pDC,
		TanhDriveCaller,
		TanhDriveCaller(pDrvHdl, dStart, dA, dB, dInitialValue));

	return pDC;
}

struct FourierSeriesDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
FourierSeriesDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "Fourier series");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	doublereal dOmega = HP.GetReal(1.);
	DEBUGCOUT("Omega: " << dOmega << std::endl);

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("FourierSeriesDriveCaller: invalid order " << n
			<< " at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	return pDC;
}

struct FrequencySweepDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
FrequencySweepDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "frequency sweep");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	DriveCaller* pOmega = 0;
	pOmega = HP.GetDriveCaller();

	DriveCaller* pAmplitude = 0;
	pAmplitude = HP.GetDriveCaller();

	doublereal dInitialValue = HP.GetReal();
	DEBUGCOUT("Initial value: " << dInitialValue << std::endl);

	doublereal dFinalTime = std::numeric_limits<double>::max();
	if (!HP.IsKeyWord("forever")) {
		dFinalTime = HP.GetReal();
	}
	DEBUGCOUT("Final time: " << dFinalTime << std::endl);

	doublereal dFinalValue = HP.GetReal();
	DEBUGCOUT("Final value: " << dFinalValue << std::endl);

	SAFENEWWITHCONSTRUCTOR(pDC,
		FreqSweepDriveCaller,
		FreqSweepDriveCaller(pDrvHdl, dInitialTime, pOmega, pAmplitude,
			dInitialValue, dFinalTime, dFinalValue));

	return pDC;
}

struct ExponentialDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ExponentialDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "exponential");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

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

	return pDC;
}

struct RandomDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
RandomDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "random");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Numero casuale */
	doublereal dAmplitude = HP.GetReal(1.);
	DEBUGCOUT("Amplitude value: " << dAmplitude << std::endl);

	doublereal dRefVal = HP.GetReal();
	DEBUGCOUT("Mean value: " << dRefVal << std::endl);

	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	doublereal dFinalTime = std::numeric_limits<double>::max();
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
				srand(time(0));
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
		RandDriveCaller(pDrvHdl, dAmplitude, dRefVal, dInitialTime, dFinalTime, iSteps));

	return pDC;
}

struct MeterDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
MeterDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "meter");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* spike every N steps */
	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	doublereal dFinalTime = std::numeric_limits<double>::max();
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
		} else {
			break;
		}
	}

	SAFENEWWITHCONSTRUCTOR(pDC,
		MeterDriveCaller,
		MeterDriveCaller(pDrvHdl, dInitialTime, dFinalTime, iSteps));

	return pDC;
}

struct ClosestNextDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ClosestNextDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "closest next");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	doublereal dInitialTime = HP.GetReal();
	DEBUGCOUT("Initial time: " << dInitialTime << std::endl);

	doublereal dFinalTime = std::numeric_limits<double>::max();
	if (!HP.IsKeyWord("forever")) {
		dFinalTime = HP.GetReal();
	}
	DEBUGCOUT("Final time: " << dFinalTime << std::endl);

	const DriveCaller *pIncrement = HP.GetDriveCaller(bDeferred);

	SAFENEWWITHCONSTRUCTOR(pDC,
		ClosestNextDriveCaller,
		ClosestNextDriveCaller(pDrvHdl, dInitialTime, dFinalTime,
			pIncrement));

	return pDC;
}

struct DirectDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
DirectDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "direct");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
		DirectDriveCaller,
		DirectDriveCaller(pDrvHdl));

	return pDC;
}

struct PiecewiseLinearDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
PiecewiseLinearDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "piecewise linear");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* Lineare a pezzi */
	unsigned int n = HP.GetInt();
	DEBUGCOUT("number of points: " << n << std::endl);

	if (n < 2) {
		silent_cerr("Need at least two points for piecewise linear drive at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal *p = 0;
	SAFENEWARR(p, doublereal, 2*n);
	p[0] = HP.GetReal();
	p[n] = HP.GetReal();
	for (unsigned int i = 1; i < n; i++) {
		p[i] = HP.GetReal();
		if (p[i] <= p[i-1]) {
			silent_cerr("point " << p[i]
				<< " is smaller than or equal to preceding point " << p[i-1]
				<< " at line " << HP.GetLineData() << std::endl);
			SAFEDELETE(p);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		p[i+n] = HP.GetReal();
	}

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		PiecewiseLinearDriveCaller,
		PiecewiseLinearDriveCaller(pDrvHdl, n, p));

	return pDC;
}

struct StringDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
StringDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "string");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* driver stringa da valutare */
	/* lettura dei dati specifici */
	std::string s(HP.GetStringWithDelims());

#ifdef USE_EE
	std::istringstream in(s);
	InputStream In(in);

	ExpressionElement *expr = HP.GetMathParser().GetExpr(In);

	if (pedantic_out) {
		std::string out = EEStrOut(expr);
		pedantic_cout("StringDriveCaller: \"" << s << "\" => \"" << out << "\"" << std::endl);
	}

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		StringDriveCaller,
		StringDriveCaller(pDrvHdl, s, expr));

#else // ! USE_EE

// #define TRIM_ALL_SPACES
#define TRIM_ALL_SPACES_BUT_ONE

#if defined(TRIM_ALL_SPACES)
	for (std::string::iterator i = s.begin(); i != s.end();) {
		if (isspace(*i)) {
			i = s.erase(i);

		} else {
			++i;
		}
	}
#elif defined(TRIM_ALL_SPACES_BUT_ONE)
	bool bString(false);
	for (std::string::iterator i = s.begin(); i != s.end();) {
		if (isspace(*i)) {
			if (!bString) {
				bString = true;
				++i;

			} else {
				i = s.erase(i);
			}

		} else {
			if (bString) {
				bString = false;
			}
			++i;
		}
	}
#endif // TRIM_ALL_SPACES_BUT_ONE

	DEBUGCOUT("String to evaluate: \"" << s << '\"' << std::endl);

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		StringDriveCaller,
		StringDriveCaller(pDrvHdl, s));

#endif // ! USE_EE

	return pDC;
}

struct DofDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
DofDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "dof");

	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* driver legato ad un grado di liberta' nodale */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no DOF dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ScalarDof SD = ReadScalarDof(pDM, HP, true, true);

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

	return pDC;
}

/*
 * shared by "node" and "element" private data drives
 */
struct SimulationEntityDCR : public DriveCallerRead {
protected:
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP,
		const SimulationEntity *pSE, const std::string& msg);
};

DriveCaller *
SimulationEntityDCR::Read(const DataManager* pDM,
	MBDynParser& HP, const SimulationEntity *pSE, const std::string& msg)
{
	const DriveHandler* pDrvHdl = pDM->pGetDrvHdl();
	DriveCaller *pDC = 0;

	unsigned int iIndex = 0;
	std::string sIndexName;
	if (HP.IsKeyWord("string")) {
		sIndexName = HP.GetStringWithDelims();
		if (sIndexName.empty()) {
			silent_cerr("empty string for " << msg
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		iIndex = pSE->iGetPrivDataIdx(sIndexName.c_str());
		if (iIndex == 0) {
			silent_cerr("illegal string \"" << sIndexName << "\""
				" for " << msg
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (HP.IsKeyWord("index")) {
		silent_cout("\"index\" deprecated for " << msg
			<< "; use \"string\" instead at line " << HP.GetLineData() << std::endl);

		iIndex = HP.GetInt();

	} else if (pSE->iGetNumPrivData() == 1) {
		iIndex = 1;

	} else {
		silent_cerr("need a private data index for " << msg
			<< " at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (iIndex < 1 || iIndex > pSE->iGetNumPrivData()) {
		silent_cerr("illegal index " << iIndex << " for " << msg
			<< " at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef USE_MPI
	/* FIXME: todo ... */
	if (MPI::Is_initialized() && MBDynComm.Get_size() > 1) {
		silent_cerr("warning: add explicit connection entry for " << msg
			<< " at line " << HP.GetLineData() << std::endl);
	}
#endif /* USE_MPI */

	/* Chiamata ricorsiva a leggere il drive supplementare */
	DriveCaller* pTmp = 0;
	pTmp = HP.GetDriveCaller();

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		PrivDriveCaller,
		PrivDriveCaller(pDrvHdl, pTmp, pSE, iIndex, sIndexName.c_str()));

	return pDC;
}

struct ElementDCR : public SimulationEntityDCR {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ElementDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "element");

	/* driver legato ai dati privati di un elemento */
	if (pDM == 0) {
		silent_cerr("since the driver is not owned by a DataManager" << std::endl
			<< "no element dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	unsigned uLabel = HP.GetInt();
	KeyTable Kel(HP, psReadElemsElems);
	int k = HP.IsKeyWord();
	if (k == -1) {
		const char *s = HP.GetString();
		silent_cerr("unknown element type \"" << s
			<< "\" at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Type(Label) */
	std::ostringstream os;
	os << psElemNames[Elem::Type(k)] << "(" << uLabel << ")";

	const Elem *pElem = pDM->pFindElem(Elem::Type(k), uLabel);
	if (pElem == 0) {
		silent_cerr("unable to find " << os.str() << " at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return SimulationEntityDCR::Read(pDM, HP, pElem, os.str());
}

struct NodeDCR : public SimulationEntityDCR {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
NodeDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "node");

	/* driver legato ai dati privati di un nodo */
	if (pDM == 0) {
		silent_cerr("since the driver is not owned by a DataManager" << std::endl
			<< "no node dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	unsigned uLabel = HP.GetInt();
	KeyTable Kel(HP, psReadNodesNodes);
	int k = HP.IsKeyWord();
	if (k == -1) {
		const char *s = HP.GetString();
		silent_cerr("unknown node type \"" << s
			<< "\" at line " << HP.GetLineData()
			<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Type(Label) */
	std::ostringstream os;
	os << psNodeNames[Elem::Type(k)] << "(" << uLabel << ")";

	const Node *pNode = pDM->pFindNode(Node::Type(k), uLabel);
	if (pNode == 0) {
		silent_cerr("unable to find " << os.str() << " at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return SimulationEntityDCR::Read(pDM, HP, pNode, os.str());
}

struct DriveDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
DriveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	DriveCaller *pD1 = HP.GetDriveCaller(bDeferred);
	DriveCaller *pD2 = HP.GetDriveCaller(bDeferred);

	SAFENEWWITHCONSTRUCTOR(pDC,
		DriveDriveCaller,
		DriveDriveCaller(pDrvHdl, pD1, pD2));

	return pDC;
}

struct SHDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
SHDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	DriveCaller *pFunc = HP.GetDriveCaller(bDeferred);
	DriveCaller *pTrigger = HP.GetDriveCaller(bDeferred);

	doublereal dVal0 = 0.;
	if (HP.IsKeyWord("initial" "value")) {
		dVal0 = HP.GetReal();
	}

	SAFENEWWITHCONSTRUCTOR(pDC,
		SHDriveCaller,
		SHDriveCaller(pDrvHdl, pFunc, pTrigger, dVal0));

	return pDC;
}

struct ArrayDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
ArrayDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	DriveCaller *pDC = 0;

	/* driver legato ad un grado di liberta' nodale */
	unsigned short int iNumDr = (unsigned short int)HP.GetInt();
	if (iNumDr == 0) {
		silent_cerr("Sorry, at least one driver is required" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);

	/* creazione di un driver normale mediante chiamata ricorsiva */
	} else if (iNumDr == 1) {
		pDC = HP.GetDriveCaller();

	} else {
		DriveArrayCaller::dcv_t DC(iNumDr);
		for (int i = 0; i < iNumDr; i++) {
			DC[i] = HP.GetDriveCaller();
		}

		/* allocazione e creazione array */
		SAFENEWWITHCONSTRUCTOR(pDC,
			DriveArrayCaller,
			DriveArrayCaller(pDrvHdl, DC));
	}

	return pDC;
}

struct FileDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
FileDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "file");

	/* driver legato ai driver */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no driver dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const DriveHandler* pDrvHdl = pDM->pGetDrvHdl();
	DriveCaller *pDC = 0;

	/* drive file */
	/* lettura dei dati specifici */
	unsigned int uL = HP.GetInt();
	FileDrive* pDrv = dynamic_cast<FileDrive*>(pDM->pFindDrive(Drive::FILEDRIVE, uL));
	if (pDrv == 0) {
		silent_cerr("line " << HP.GetLineData()
			<< ": can't find FileDrive(" << uL << ")" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer id = 1;
	if (HP.IsArg()) {
		id = HP.GetInt(id);
		if (id < 1 || id > pDrv->iGetNumDrives()) {
			silent_cerr("line " << HP.GetLineData()
				<< ": invalid column number " << id
				<< " (must be between 1 and " 
				<< pDrv->iGetNumDrives() << ")" << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	doublereal da = 1.;
	if (HP.IsKeyWord("amplitude")) {
		da = HP.GetReal();
	}

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		FileDriveCaller,
		FileDriveCaller(pDrvHdl, pDrv, id, da));

	return pDC;
}

struct PeriodicDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
PeriodicDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	NeedDM(pDM, HP, bDeferred, "periodic");

	/* driver legato ai driver */
	if (pDM == 0) {
		silent_cerr("sorry, since the driver is not owned by a DataManager" << std::endl
			<< "no driver dependent drivers are allowed;" << std::endl
			<< "aborting..." << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const DriveHandler* pDrvHdl = pDM->pGetDrvHdl();
	DriveCaller *pDC = 0;

	doublereal dT0 = HP.GetReal();
	doublereal dPeriod = HP.GetReal();
	if (dPeriod <= 0.) {
		silent_cerr("PeriodicDriveCaller: invalid negative or null period at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const DriveCaller *pDC1 = HP.GetDriveCaller();

	/* allocazione e creazione */
	SAFENEWWITHCONSTRUCTOR(pDC,
		PeriodicDriveCaller,
		PeriodicDriveCaller(pDrvHdl, pDC1, dT0, dPeriod));

	return pDC;
}


struct PostponedDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

DriveCaller *
PostponedDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	integer label = HP.GetInt();
	if (label < 0) {
		silent_cerr("PostponedDriveCaller: invalid negative label \"" << label << "\" at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* allocazione e creazione */
	DriveCaller *pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		PostponedDriveCaller,
		PostponedDriveCaller(HP, unsigned(label)));

	return pDC;
}



static unsigned d_done;

void
InitDriveData(void)
{
	if (::d_done++ > 0) {
		return;
	}

	SetDriveData("fixed" "step", new FixedStepDR);
	SetDriveData("variable" "step", new VariableStepDR);
#ifdef USE_SOCKET
	SetDriveData("socket", new SocketDR);
	SetDriveData("socket" "stream", new StreamDR("socket stream"));
	SetDriveData("rtai" "input", new StreamDR("RTAI input"));
	SetDriveData("stream", new StreamDR);
#endif // USE_SOCKET
	SetDriveData("buffer" "stream", new BufferStreamDR);

	/* NOTE: add here initialization of new built-in drives;
	 * alternative ways to register new custom drives are:
	 * - call SetDriveData() from anywhere in the code
	 * - write a module that calls SetDriveData() from inside a function
	 *   called module_init(), and run-time load it using "module load"
	 *   in the input file.
	 */
}

void
DestroyDriveData(void)
{
	if (::d_done == 0) {
		silent_cerr("DestroyDriveData() called once too many" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--::d_done > 0) {
		return;
	}

	/* free stuff */
	for (DriveFuncMapType::iterator i = DriveFuncMap.begin(); i != DriveFuncMap.end(); ++i) {
		delete i->second;
	}
	DriveFuncMap.clear();
}


static unsigned dc_done;

void
InitDriveCallerData(void)
{
	if (::dc_done++ > 0) {
		return;
	}

	SetDriveCallerData("array", new ArrayDCR);
	SetDriveCallerData("closest" "next", new ClosestNextDCR);
	SetDriveCallerData("const", new ConstDCR);
	SetDriveCallerData("cosine", new CosineDCR);
	SetDriveCallerData("cubic", new CubicDCR);
	SetDriveCallerData("direct", new DirectDCR);
	SetDriveCallerData("dof", new DofDCR);
	SetDriveCallerData("double" "ramp", new DoubleRampDCR);
	SetDriveCallerData("double" "step", new DoubleStepDCR);
	SetDriveCallerData("drive", new DriveDCR);
	SetDriveCallerData("element", new ElementDCR);
	SetDriveCallerData("exponential", new ExponentialDCR);
	SetDriveCallerData("file", new FileDCR);
	SetDriveCallerData("fourier" "series", new FourierSeriesDCR);
	SetDriveCallerData("frequency" "sweep", new FrequencySweepDCR);
	SetDriveCallerData("function", new FunctionDCR);
#ifdef USE_GINAC
	SetDriveCallerData("ginac", new GiNaCDCR);
#endif // USE_GINAC
	SetDriveCallerData("linear", new LinearDCR);
	SetDriveCallerData("meter", new MeterDCR);
	SetDriveCallerData("mult", new MultDCR);
	SetDriveCallerData("node", new NodeDCR);
	SetDriveCallerData("null", new NullDCR);
	SetDriveCallerData("one", new OneDCR);	/* deprecated */
	SetDriveCallerData("parabolic", new ParabolicDCR);
	SetDriveCallerData("periodic", new PeriodicDCR);
	SetDriveCallerData("piecewise" "linear", new PiecewiseLinearDCR);
	SetDriveCallerData("postponed", new PostponedDCR);
	SetDriveCallerData("ramp", new RampDCR);
	SetDriveCallerData("random", new RandomDCR);
	SetDriveCallerData("sample" "and" "hold", new SHDCR);
	SetDriveCallerData("sine", new SineDCR);
	SetDriveCallerData("step", new StepDCR);
	SetDriveCallerData("string", new StringDCR);
	SetDriveCallerData("tanh", new TanhDCR);
	SetDriveCallerData("time", new TimeDCR);
	SetDriveCallerData("timestep", new TimeStepDCR);
	SetDriveCallerData("unit", new OneDCR);

	/* NOTE: add here initialization of new built-in drive callers;
	 * alternative ways to register new custom drive callers are:
	 * - call SetDriveCallerData() from anywhere in the code
	 * - write a module that calls SetDriveCallerData() from inside a function
	 *   called module_init(), and run-time load it using "module load"
	 *   in the input file.
	 */
#ifdef STATIC_MODULES
#ifdef USE_OCTAVE
	mbdyn_octave_set();
#endif // USE_OCTAVE
	nodedistdrive_set();
	minmaxdrive_set();
#endif // STATIC_MODULES
}

void
DestroyDriveCallerData(void)
{
	if (::dc_done == 0) {
		silent_cerr("DestroyDriveCallerData() called once too many" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--::dc_done > 0) {
		return;
	}

	/* free stuff */
	for (DriveCallerFuncMapType::iterator i = DriveCallerFuncMap.begin(); i != DriveCallerFuncMap.end(); ++i) {
		delete i->second;
	}
	DriveCallerFuncMap.clear();
}
