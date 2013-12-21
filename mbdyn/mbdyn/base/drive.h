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

#ifndef DRIVE_H
#define DRIVE_H


/* include generali */
#include <time.h>
#include "ac/f2c.h"
#include "ac/pthread.h"

/* include per il debug */
#include "myassert.h"
#include "mynewmem.h"

/* include del programma */
#include "mathp.h"
#include "output.h"
#include "solman.h"
#include "withlab.h"

extern const char* psDriveNames[];
extern const char* psReadControlDrivers[];

/* Drive - begin */

/* Classe dei drivers, ovvero oggetti che sono in grado di restituire
 * il valore di certi parametri utilizzati da elementi in funzione
 * di altri parametri indipendenti, principalmente del tempo.
 * Per ora la principale classe derivata da questi e' costituita
 * da funzioni ottenute per interpolazione di dati forniti dall'utente.
 * Esistono altri due tipi di drivers che, per varie ragioni, non
 * sono stati derivati da questa classe, ma incorporati
 * direttamente nel DriveHandler (vedi sotto). Questi due tipi restituiscono
 * il valore mediante valutazione di funzioni prestabilite in dipendenza da
 * parametri forniti dal valutatore, dato dalla classe DriveCaller.
 * Nel primo caso la funzione e' "built-in" nel DriveHandler, nel secondo
 * si ottiene mediante valutazione nel MathParser di una stringa fornita
 * dall'utente. La tabella dei simboli del MathParser contiene di default
 * la variabile Time, che viene mantenuta aggiornata e puo' essere usata
 * per valutare il valore delle espressioni fornite dall'utente.
 *
 * Quindi:
 * - ogni elemento che usa drivers contiene gli opportuni DriveCaller.
 * - il DriveCaller contiene i dati necessari per chiamare il DriveHandler.
 * - il DriveHandler restituisce il valore richiesto:
 *   - valutando le funzioni "built-in"
 *   - interrogando i singoli Drive.
 * */

class DriveHandler;
class DriveCaller;

class Drive : public WithLabel {
public:
	// Tipi di drive
	enum Type {
		UNKNOWN = -1,

		FILEDRIVE = 0,

		LASTDRIVETYPE
	};

public:
	enum Bailout {
		BO_NONE		= 0x0,
		BO_UPPER	= 0x1,
		BO_LOWER	= 0x2,
		BO_ANY		= (BO_UPPER | BO_LOWER)
	};

protected:
	const DriveHandler* pDrvHdl;
	static doublereal dReturnValue;

public:
	Drive(unsigned int uL, const DriveHandler* pDH);
	virtual ~Drive(void);

	/* Tipo del drive (usato solo per debug ecc.) */
	virtual Drive::Type GetDriveType(void) const = 0;

	virtual std::ostream& Restart(std::ostream& out) const = 0;

	virtual void ServePending(const doublereal& t) = 0;
};

/* Drive - end */

/* DriveOwner - begin */

/* Possessore di DriveCaller, ne garantisce la corretta distruzione */

class DriveOwner {
protected:
	DriveCaller* pDriveCaller;

public:
	DriveOwner(const DriveCaller* pDC = 0);
	DriveOwner(const DriveOwner& drive);
	virtual ~DriveOwner(void);

	void Set(const DriveCaller* pDC);
	DriveCaller* pGetDriveCaller(void) const;

	doublereal dGet(const doublereal& dVar) const;
	doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	bool bIsDifferentiable(void) const;
	doublereal dGetP(const doublereal& dVar) const;
	doublereal dGetP(void) const;
};

/* DriveOwner - end */




/* DriveHandler - begin */

/* Classe che gestisce la valutazione del valore dei Drive. L'oggetto e'
 * posseduto dal DataManager, che lo costruisce e lo mette a disposizione
 * dei DriveCaller. Contiene i Drive definiti dall'utente piu' i due
 * "built-in" che valutano funzioni prestabilite o fornite dall'utente
 * sotto forma di stringa */

class DriveHandler {
	friend class DataManager;
	friend class RandDriveCaller;
	friend class MeterDriveCaller;
	friend class ClosestNextDriveCaller;
	friend class SHDriveCaller;

private:
#ifdef USE_MULTITHREAD
	mutable pthread_mutex_t parser_mutex;
#endif /* USE_MULTITHREAD */
	MathParser& Parser;

	/* variabili predefinite: tempo e variabile generica */
	Var* pTime;
	Var* pTimeStep;
	Var* pStep;
	Var* pVar;

	static doublereal dDriveHandlerReturnValue; /* Usato per ritornare un reference */

	mutable VectorHandler* pXCurr;
	mutable VectorHandler* pXPrimeCurr;

	integer iCurrStep;

	/* For meters */
	class MyMeter : public WithLabel {
	protected:
		integer iSteps;
		mutable integer iCurr;

	public:
		MyMeter(unsigned int uLabel, integer iS = 1);
		virtual ~MyMeter(void);

		inline integer iGetSteps(void) const;
		inline bool bGetMeter(void) const;
		virtual inline void Set(void);
	};

	/* For random drivers */
	class MyRand : public MyMeter {
	protected:
		integer iRand;

	public:
		MyRand(unsigned int uLabel, integer iS = 1, integer iR = 0);
		virtual ~MyRand(void);

		inline integer iGetSteps(void) const;
		inline integer iGetRand(void) const;
		virtual inline void Set(void);
	};

	/* For closest next drivers */
	class MyClosestNext : public WithLabel {
	protected:
		const DriveHandler *pDH;
		const DriveOwner Increment;
		bool bMustSetNext;
		doublereal dNext;

	public:
		MyClosestNext(unsigned int uLabel, const DriveHandler *pDH,
			const DriveCaller *pIncrementDC, doublereal dStartTime);
		virtual ~MyClosestNext(void);

		inline bool bGetClosestNext(void) const;
		virtual inline void Set(void);
	};

	/* For sample'n'hold */
	class MySH : public WithLabel {
	protected:
		const doublereal dVal0;
		mutable doublereal dVal;
		const DriveOwner Func;
		const DriveOwner Trigger;

	public:
		MySH(unsigned int uLabel,
			const DriveCaller *pFunc,
			const DriveCaller *pTrigger,
			const doublereal dVal0);
		virtual ~MySH(void);

		inline doublereal dGetSH(void) const;
		virtual inline void Set(void);
		const DriveCaller *pGetFunc(void) const;
		const DriveCaller *pGetTrigger(void) const;
		const doublereal dGetVal0(void) const;
	};

	std::vector<MyMeter *> Meter;
	std::vector<MyRand *> Rand;
	std::vector<MyClosestNext *> ClosestNext;
	std::vector<MySH *> SH;

protected:
	void SetTime(const doublereal& dt, const doublereal& dts,
		const integer& s);
	void LinkToSolution(const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	integer iRandInit(integer iSteps);
	integer iMeterInit(integer iSteps);
	integer iClosestNextInit(const DriveCaller *pIncrement, doublereal dStartTime);
	integer iSHInit(const DriveCaller *pFunc, const DriveCaller *pTrigger,
		const doublereal dVal0);

public:
	DriveHandler(MathParser &mp);
	~DriveHandler(void);

	void PutSymbolTable(Table& T);
	void SetVar(const doublereal& dVar);

	doublereal dGet(InputStream& InStr) const;

	inline doublereal dGetTime(void) const;
	inline doublereal dGetTimeStep(void) const;
	inline integer iGetStep(void) const;
	inline long int iGetRand(integer iNumber) const;
	inline bool bGetMeter(integer iNumber) const;
	inline bool bGetClosestNext(integer iNumber) const;
	inline doublereal dGetSH(integer iNumber) const;

	const DriveCaller *pGetSHFunc(integer iNumber) const;
	const DriveCaller *pGetSHTrigger(integer iNumber) const;
	const doublereal dGetSHVal0(integer iNumber) const;
};


inline integer
DriveHandler::MyRand::iGetRand(void) const
{
	return iRand;
}


inline void
DriveHandler::MyRand::Set(void)
{
	MyMeter::Set();
	if (iCurr == 0) {
		iRand = rand();
	}
}


inline integer
DriveHandler::MyMeter::iGetSteps(void) const
{
	return iSteps;
}


inline bool
DriveHandler::MyMeter::bGetMeter(void) const
{
	return iCurr == 0;
}


inline void
DriveHandler::MyMeter::Set(void)
{
	if (++iCurr == iSteps) {
		iCurr = 0;
	}
}


inline bool
DriveHandler::MyClosestNext::bGetClosestNext(void) const
{
	return pDH->dGetTime() >= dNext;
}


inline void
DriveHandler::MyClosestNext::Set(void)
{
	if (bMustSetNext) {
		do {
			dNext += Increment.dGet();
		} while (bGetClosestNext());

		bMustSetNext = false;

	} else if (bGetClosestNext()) {
		bMustSetNext = true;
	}
}


inline doublereal
DriveHandler::MySH::dGetSH(void) const
{
	return dVal;
}


inline void
DriveHandler::MySH::Set(void)
{
	if (Trigger.dGet()) {
		dVal = Func.dGet();
	}
}


inline doublereal
DriveHandler::dGetTime(void) const
{
	ASSERT(pTime != 0);
	return pTime->GetVal().GetReal();
}

inline doublereal
DriveHandler::dGetTimeStep(void) const
{
	ASSERT(pTimeStep != 0);
	return pTimeStep->GetVal().GetReal();
}

inline integer
DriveHandler::iGetStep(void) const
{
	ASSERT(pStep != 0);
	return pStep->GetVal().GetInt();
}


inline long int
DriveHandler::iGetRand(integer iNumber) const
{
	return Rand[iNumber]->iGetRand();
}


inline bool
DriveHandler::bGetMeter(integer iNumber) const
{
	return Meter[iNumber]->bGetMeter();
}

inline bool
DriveHandler::bGetClosestNext(integer iNumber) const
{
	return ClosestNext[iNumber]->bGetClosestNext();
}

inline doublereal
DriveHandler::dGetSH(integer iNumber) const
{
	return SH[iNumber]->dGetSH();
}

/* DriveHandler - end */



/* DriveCaller - begin */

/* Classe che chiama il drive handler per ottenere il valore del driver ad
 * un dato istante di tempo. Gli oggetti sono posseduti dai singoli elementi.
 * Ogni oggetto contiene i parametri che gli occorrono per la chiamata. */

class DriveCaller : public WithLabel {
protected:
	mutable DriveHandler* pDrvHdl;

public:
	DriveCaller(const DriveHandler* pDH);
	virtual ~DriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const = 0;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;

	/* Restituisce il valore del driver */
	virtual doublereal dGet(const doublereal& dVar) const = 0;
	virtual inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;

	/* allows to set the drive handler */
	virtual void SetDrvHdl(const DriveHandler* pDH);
	virtual const DriveHandler *pGetDrvHdl(void) const;
};

inline doublereal
DriveCaller::dGet(void) const
{
	return dGet(pDrvHdl->dGetTime());
}

inline bool
DriveCaller::bIsDifferentiable(void) const
{
	return false;
}

inline doublereal
DriveCaller::dGetP(void) const
{
	return dGetP(pDrvHdl->dGetTime());
}

/* DriveCaller - end */


/* NullDriveCaller - begin */

class NullDriveCaller : public DriveCaller {
public:
	NullDriveCaller(void);
	virtual ~NullDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Restituisce il valore del driver */
	virtual inline doublereal dGet(const doublereal& dVar) const;
	virtual inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
};

inline doublereal
NullDriveCaller::dGet(const doublereal& /* dVar */ ) const
{
	return 0.;
}

inline doublereal
NullDriveCaller::dGet(void) const
{
	return 0.;
}

inline bool
NullDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

inline doublereal
NullDriveCaller::dGetP(const doublereal& /* dVar */ ) const
{
	return 0.;
}

inline doublereal
NullDriveCaller::dGetP(void) const
{
	return 0.;
}

/* NullDriveCaller - end */


/* OneDriveCaller - begin */

class OneDriveCaller : public DriveCaller {
public:
	OneDriveCaller(void);
	virtual ~OneDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Restituisce il valore del driver */
	virtual inline doublereal dGet(const doublereal& dVar) const;
	virtual inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
};

inline doublereal
OneDriveCaller::dGet(const doublereal& /* dVar */ ) const
{
	return 1.;
}

inline doublereal
OneDriveCaller::dGet(void) const
{
	return 1.;
}

inline bool
OneDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

inline doublereal
OneDriveCaller::dGetP(const doublereal& /* dVar */ ) const
{
	return 0.;
}

inline doublereal
OneDriveCaller::dGetP(void) const
{
	return 0.;
}

/* OneDriveCaller - end */


/* ConstDriveCaller - begin */

class ConstDriveCaller : public DriveCaller {
private:
	doublereal dConst;

public:
	ConstDriveCaller(doublereal d);
	virtual ~ConstDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	inline doublereal dGet(const doublereal& /* dVar */ ) const;
	inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
};

inline doublereal
ConstDriveCaller::dGet(const doublereal& /* dVar */ ) const
{
	return dConst;
}

inline doublereal
ConstDriveCaller::dGet(void) const
{
	return dConst;
}

inline bool
ConstDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

inline doublereal
ConstDriveCaller::dGetP(const doublereal& /* dVar */ ) const
{
	return 0.;
}

inline doublereal
ConstDriveCaller::dGetP(void) const
{
	return 0.;
}

/* ConstDriveCaller - end */


class DataManager;
class MBDynParser;

/* prototype of the functional object: reads a drive caller */
struct DriveCallerRead {
protected:
	/* Helper */
	void
	NeedDM(const DataManager* pDM, MBDynParser& HP, bool bDeferred,
		const char *const name);

public:
	virtual ~DriveCallerRead(void);
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) = 0;
};

/* drive caller registration function: call to register one */
extern bool
SetDriveData(const char *name, DriveCallerRead* rf);

/* function that reads a drive caller */
extern DriveCaller*
ReadDriveData(const DataManager* pDM, MBDynParser& HP, bool bDeferred);

/* create/destroy */
extern void InitDriveData(void);
extern void DestroyDriveData(void);

#endif /* DRIVE_H */

