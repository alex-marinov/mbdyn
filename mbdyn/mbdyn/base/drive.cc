/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "drive.h"

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
: Parser(mp),
pTime(0),
pTimeStep(0),
pStep(0),
pVar(0),
pXCurr(0),
pXPrimeCurr(0),
iCurrStep(0),
Meter(0),
Rand(0),
ClosestNext(0),
SH(0)
{
#ifdef USE_MULTITHREAD
	if (pthread_mutex_init(&parser_mutex, NULL)) {
		silent_cerr("DriveHandler::DriveHandler(): mutex init failed"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif /* USE_MULTITHREAD */

	NamedValue *v;

	/* Inserisce la variabile Time nella tabella dei simboli; sara'
	 * mantenuta aggiornata dal DriveHandler */
	v = Parser.GetSymbolTable().Get("Time");
	if (v == 0) {
		pTime = Parser.GetSymbolTable().Put("Time", Real(0));
		if (pTime == 0) {
			silent_cerr("DriveHandler::DriveHandler(): "
				"error while inserting symbol 'Time'"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (!v->IsVar()) {
			silent_cerr("Symbol 'Time' must be a variable"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pTime = dynamic_cast<Var *>(v);
	}

	/* Inserisce la variabile TimeStep nella tabella dei simboli; sara'
	 * mantenuta aggiornata dal DriveHandler */
	v = Parser.GetSymbolTable().Get("TimeStep");
	if (v == 0) {
		pTimeStep = Parser.GetSymbolTable().Put("TimeStep", Real(-1.));
		if (pTimeStep == 0) {
			silent_cerr("DriveHandler::DriveHandler(): "
				"error while inserting symbol 'TimeStep'"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (!v->IsVar()) {
			silent_cerr("Symbol 'TimeStep' must be a variable"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pTimeStep = dynamic_cast<Var *>(v);
	}

	/* Inserisce la variabile Step nella tabella dei simboli; sara'
	 * mantenuta aggiornata dal DriveHandler */
	v = Parser.GetSymbolTable().Get("Step");
	if (v == 0) {
		pStep = Parser.GetSymbolTable().Put("Step", Int(-1));
		if (pStep == 0) {
			silent_cerr("DriveHandler::DriveHandler(): "
				"error while inserting symbol 'Step'"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (!v->IsVar()) {
			silent_cerr("Symbol 'Step' must be a variable"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pStep = dynamic_cast<Var *>(v);
	}

	/* Inserisce la variabile Var nella tabella dei simboli; sara'
	 * mantenuta aggiornata dai DriveCaller attraverso il DriveHandler */
	v = Parser.GetSymbolTable().Get("Var");
	if (v == 0) {
		pVar = Parser.GetSymbolTable().Put("Var", Real(0));
		if (pVar == 0) {
			silent_cerr("DriveHandler::DriveHandler(): "
				"error while insterting symbol 'Var'"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		if (!v->IsVar()) {
			silent_cerr("Symbol 'Var' must be a variable"
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pVar = dynamic_cast<Var *>(v);
	}

	/* Calcola il seed di riferimento per i numeri random */
	srand(time(NULL));
}

DriveHandler::~DriveHandler(void)
{
#ifdef USE_MULTITHREAD
	pthread_mutex_destroy(&parser_mutex);
#endif /* USE_MULTITHREAD */

	for (std::vector<MyMeter *>::iterator i = Meter.begin();
		i != Meter.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	for (std::vector<MyRand *>::iterator i = Rand.begin();
		i != Rand.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	for (std::vector<MyClosestNext*>::iterator i = ClosestNext.begin();
		i != ClosestNext.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	for (std::vector<MySH *>::iterator i = SH.begin();
		i != SH.end(); ++i)
	{
		SAFEDELETE(*i);
	}
}

void
DriveHandler::SetTime(const doublereal& dt, const doublereal& dts,
	const integer& s)
{
	/* Setta la variabile Time nella tabella dei simboli */
	ASSERT(pTime != 0);
	pTime->SetVal(dt);

	/* Setta la variabile TimeStep nella tabella dei simboli */
	if (dts >= 0.) {
		ASSERT(pTimeStep != 0);
		pTimeStep->SetVal(dts);
	}

	/* Setta la variabile Step nella tabella dei simboli */
	if (s >= 0) {
		ASSERT(pStep != 0);
		pStep->SetVal(s);

		/* in case of new step */
		if (s != iCurrStep) {
			ASSERT(iCurrStep + 1 == s);
			iCurrStep = s;

			for (std::vector<MyMeter *>::iterator i = Meter.begin();
				i != Meter.end(); ++i)
			{
				(*i)->Set();
			}

			for (std::vector<MyRand *>::iterator i = Rand.begin();
				i != Rand.end(); ++i)
			{
				(*i)->Set();
			}

			for (std::vector<MyClosestNext*>::iterator i = ClosestNext.begin();
				i != ClosestNext.end(); ++i)
			{
				(*i)->Set();
			}

			for (std::vector<MySH *>::iterator i = SH.begin();
				i != SH.end(); ++i)
			{
				(*i)->Set();
			}
		}
	}
}


void
DriveHandler::LinkToSolution(const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	pXCurr = const_cast<VectorHandler *>(&XCurr);
	pXPrimeCurr = const_cast<VectorHandler *>(&XPrimeCurr);
}

integer
DriveHandler::iRandInit(integer iSteps)
{
	MyRand* pmr = 0;
	integer iNumber = Rand.size();
	SAFENEWWITHCONSTRUCTOR(pmr,
		MyRand,
		MyRand((unsigned int)iNumber, iSteps, rand()));
	Rand.push_back(pmr);

	return iNumber;
}

integer
DriveHandler::iMeterInit(integer iSteps)
{
	MyMeter* pmm = 0;
	integer iNumber = Meter.size();
	SAFENEWWITHCONSTRUCTOR(pmm,
		MyMeter,
		MyMeter((unsigned int)iNumber, iSteps));
	Meter.push_back(pmm);

	return iNumber;
}

integer
DriveHandler::iClosestNextInit(const DriveCaller *pIncrement,
	doublereal dStartTime)
{
	MyClosestNext* pmc = 0;
	integer iNumber = ClosestNext.size();
	SAFENEWWITHCONSTRUCTOR(pmc,
		MyClosestNext,
		MyClosestNext((unsigned int)iNumber, this,
			pIncrement, dStartTime));
	ClosestNext.push_back(pmc);

	return iNumber;
}


integer
DriveHandler::iSHInit(const DriveCaller *pFunc, const DriveCaller *pTrigger,
	const doublereal dVal0)
{
	MySH* pms = 0;
	integer iNumber = SH.size();
	SAFENEWWITHCONSTRUCTOR(pms,
		MySH,
		MySH((unsigned int)iNumber, pFunc, pTrigger, dVal0));
	SH.push_back(pms);

	return iNumber;
}

const DriveCaller *
DriveHandler::pGetSHFunc(integer iNumber) const
{
	return SH[iNumber]->pGetFunc();
}

const DriveCaller *
DriveHandler::pGetSHTrigger(integer iNumber) const
{
	return SH[iNumber]->pGetTrigger();
}

const doublereal
DriveHandler::dGetSHVal0(integer iNumber) const
{
	return SH[iNumber]->dGetVal0();
}

void
DriveHandler::PutSymbolTable(Table& T)
{
	Parser.PutSymbolTable(T);
}


void
DriveHandler::SetVar(const doublereal& dVar)
{
	ASSERT(pVar != 0);
	pVar->SetVal(dVar);
}


doublereal
DriveHandler::dGet(InputStream& InStr) const
{
	doublereal d;

#ifdef USE_MULTITHREAD
	// FIXME: risk recursive lock
	pthread_mutex_lock(&parser_mutex);
#endif /* USE_MULTITHREAD */

	try {
		d = Parser.GetLastStmt(InStr);

	} catch (MBDynErrBase e) {
		silent_cerr("StringDrive: " << e.what() << std::endl);
		throw e;

#if 0
	} catch (ErrGeneric e) {
		silent_cerr("StringDrive: " << e.what() << std::endl);
		throw e;
#endif

	} catch (...) {
		silent_cerr("StringDrive error" << std::endl);
		throw;
	}

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

DriveHandler::MyClosestNext::MyClosestNext(unsigned int uLabel,
	const DriveHandler *pDH,
	const DriveCaller *pIncrement,
	doublereal dStartTime)
: WithLabel(uLabel), pDH(pDH), Increment(pIncrement), bMustSetNext(false), dNext(dStartTime)
{
	NO_OP;
}

DriveHandler::MyClosestNext::~MyClosestNext(void)
{
	NO_OP;
}

DriveHandler::MySH::MySH(unsigned int uLabel,
	const DriveCaller *pFunc,
	const DriveCaller *pTrigger,
	const doublereal dVal0)
: WithLabel(uLabel), dVal0(dVal0), dVal(dVal0), Func(pFunc), Trigger(pTrigger)
{
	NO_OP;
}

DriveHandler::MySH::~MySH(void)
{
	NO_OP;
}

const DriveCaller *
DriveHandler::MySH::pGetFunc(void) const
{
	return Func.pGetDriveCaller();
}

const DriveCaller *
DriveHandler::MySH::pGetTrigger(void) const
{
	return Trigger.pGetDriveCaller();
}

const doublereal
DriveHandler::MySH::dGetVal0(void) const
{
	return dVal0;
}

/* DriveHandler - end */

/* DriveCaller - begin */

DriveCaller::DriveCaller(const DriveHandler* pDH)
: pDrvHdl((DriveHandler *)pDH)
{
	NO_OP;
}

DriveCaller::~DriveCaller(void)
{
	NO_OP;
}

void
DriveCaller::SetDrvHdl(const DriveHandler* pDH)
{
	pDrvHdl = const_cast<DriveHandler *>(pDH);
}

const DriveHandler *
DriveCaller::pGetDrvHdl(void) const
{
	return pDrvHdl;
}

doublereal
DriveCaller::dGetP(const doublereal& dVar) const
{
	/* shouldn't get called if not differentiable,
	 * or should be overridden if differentiable */
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
DriveCaller *
NullDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEW(pDC, NullDriveCaller);

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
NullDriveCaller::Restart(std::ostream& out) const
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
DriveCaller *
OneDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;
	SAFENEW(pDC, OneDriveCaller);

	return pDC;
}

/* Scrive il contributo del DriveCaller al file di restart */
std::ostream&
OneDriveCaller::Restart(std::ostream& out) const
{
	return out << "one";
}

/* OneDriveCaller - end */


/* DriveOwner - begin */

DriveOwner::DriveOwner(const DriveCaller* pDC)
: pDriveCaller(const_cast<DriveCaller*>(pDC))
{
	NO_OP;
}

DriveOwner::DriveOwner(const DriveOwner& drive)
: pDriveCaller(drive.pDriveCaller ? drive.pDriveCaller->pCopy() : 0)
{
	NO_OP;
}

DriveOwner::~DriveOwner(void)
{
	if (pDriveCaller != 0) {
		SAFEDELETE(pDriveCaller);
	}
}

void
DriveOwner::Set(const DriveCaller* pDC)
{
	ASSERT(pDC != 0);
	if (pDriveCaller != 0) {
		DEBUGCOUT("warning: the original pointer to a drive caller is not null!" << std::endl);
		SAFEDELETE(pDriveCaller);
	}
	pDriveCaller = const_cast<DriveCaller*>(pDC);
}

DriveCaller *
DriveOwner::pGetDriveCaller(void) const
{
	return pDriveCaller;
}

doublereal
DriveOwner::dGet(const doublereal& dVal) const
{
	return pDriveCaller->dGet(dVal);
}


doublereal
DriveOwner::dGet(void) const
{
	return pDriveCaller->dGet();
}

bool
DriveOwner::bIsDifferentiable(void) const
{
	return pDriveCaller->bIsDifferentiable();
}

doublereal
DriveOwner::dGetP(const doublereal& dVal) const
{
	return pDriveCaller->dGetP(dVal);
}


doublereal
DriveOwner::dGetP(void) const
{
	return pDriveCaller->dGetP();
}

/* DriveOwner - end */

