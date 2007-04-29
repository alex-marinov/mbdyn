/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifndef TPLDRIVE__H
#define TPLDRIVE__H

#include <drive.h>
#include <tpldrive.h>
#include <mbpar.h>

/* ZeroTplDriveCaller - begin */

template <class T>
class ZeroTplDriveCaller : public TplDriveCaller<T> {
public:
	ZeroTplDriveCaller(void) {
		NO_OP;
	};

	~ZeroTplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		typedef ZeroTplDriveCaller<T> dc;
		TplDriveCaller<T>* pDC = 0;

		SAFENEW(pDC, dc);

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "zero";
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		return out;
	};

	inline T Get(void) const {
		return T(0.);
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		return true;
	};

	inline T GetP(void) const {
		return T(0.);
	};

	inline int getNDrives(void) const {
		return 0;
	};
};

/* ZeroTplDriveCaller - end */


/* SingleTplDriveCaller - begin */

template <class T>
class SingleTplDriveCaller : public TplDriveCaller<T>, public DriveOwner {
protected:
	T t;

public:
	SingleTplDriveCaller(const DriveCaller* pDC, const T& x)
	: DriveOwner(pDC), t((T&)x) {
		NO_OP;
	};

	~SingleTplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		typedef SingleTplDriveCaller<T> dc;
		TplDriveCaller<T>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy(), t));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "single, ",
		Write(out, t, ", ") << ", ";
		return pGetDriveCaller()->Restart(out);
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		Write(out, t, ", ") << ", ";
		return pGetDriveCaller()->Restart(out);
	};

	inline T Get(void) const {
		return t*dGet();
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		return DriveOwner::bIsDifferentiable();
	};

	inline T GetP(void) const {
		return t*dGetP();
	};

	inline int getNDrives(void) const {
		return 1;
	};
};

/* Nota: in caso scalare, viene semplificata la classe in modo da
 *       usare solo il drive senza pesatura che viene assunta unitaria */

template<>
class SingleTplDriveCaller<doublereal>
: public TplDriveCaller<doublereal>, public DriveOwner {
public:
	SingleTplDriveCaller(const DriveCaller* pDC, const doublereal& = 0.)
	: DriveOwner(pDC) {
		NO_OP;
	};

	~SingleTplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<doublereal>* pCopy(void) const {
		TplDriveCaller<doublereal>* pDC = 0;

		typedef SingleTplDriveCaller<doublereal> dc;
		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy()));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "single, ";
		return pGetDriveCaller()->Restart(out);
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		return pGetDriveCaller()->Restart(out);
	};

	inline doublereal Get(void) const {
		return dGet();
	};

	inline bool bIsDifferentiable(void) const {
		return DriveOwner::bIsDifferentiable();
	};

	inline doublereal GetP(void) const {
		return dGetP();
	};

	inline int getNDrives(void) const {
		return 1;
	};
};

/* SingleTplDriveCaller - end */


/* ArrayTplDriveCaller - begin */

template <class T>
struct DrivesArray {
	DriveCaller* pDriveCaller;
	T t;
};


template <class T>
class ArrayTplDriveCaller : public TplDriveCaller<T> {
protected:
	DrivesArray<T>* pDrivesArray;
	unsigned short int iNumDrives;

public:
	ArrayTplDriveCaller(unsigned short int i, DrivesArray<T>* pDA)
	: pDrivesArray(pDA), iNumDrives(i) {
		ASSERT(i > 0);
		ASSERT(pDA != 0);
	};

	~ArrayTplDriveCaller(void) {
		for (int i = 0; i < iNumDrives; i++) {
			SAFEDELETE(pDrivesArray[i].pDriveCaller);
		}
		SAFEDELETEARR(pDrivesArray);
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		typedef DrivesArray<T> da;
		da* pDA = 0;

		SAFENEWARR(pDA, da, iNumDrives);

		for (int i = 0; i < iNumDrives; i++) {
			pDA[i].pDriveCaller = pDrivesArray[i].pDriveCaller->pCopy();
			pDA[i].t = pDrivesArray[i].t;
		}

		typedef ArrayTplDriveCaller<T> dc;
		TplDriveCaller<T>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(iNumDrives, pDA));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "array, " << iNumDrives;
		for (int i = 0; i < iNumDrives; i++) {
			out << ", ",
				Write(out, pDrivesArray[i].t, ", ") << ", ",
				pDrivesArray[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (int i = 0; i < iNumDrives; i++) {
			out << ", ",
				Write(out, pDrivesArray[i].t, ", ") << ", ",
				pDrivesArray[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	inline T Get(void) const {
		T v = 0.;
		for (int i = 0; i < iNumDrives; i++) {
			v += (pDrivesArray[i].t)*(pDrivesArray[i].pDriveCaller->dGet());
		}
		return v;
	};

	inline bool bIsDifferentiable(void) const {
		for (int i = 0; i < iNumDrives; i++) {
			if (!pDrivesArray[i].pDriveCaller->bIsDifferentiable()) {
				return false;
			}
		}
		return true;
	};

	inline T GetP(void) const {
		T v = 0.;
		for (int i = 0; i < iNumDrives; i++) {
			v += (pDrivesArray[i].t)*(pDrivesArray[i].pDriveCaller->dGetP());
		}
		return v;
	};

	inline int getNDrives(void) const {
		return iNumDrives;
	};
};

template<>
class ArrayTplDriveCaller<doublereal> : public TplDriveCaller<doublereal> {
protected:
	DrivesArray<doublereal>* pDrivesArray;
	unsigned short int iNumDrives;

public:
	ArrayTplDriveCaller(unsigned short int i,
		DrivesArray<doublereal>* pDA)
	: pDrivesArray(pDA), iNumDrives(i) {
		ASSERT(i > 0);
		ASSERT(pDA != 0);
	};

	virtual ~ArrayTplDriveCaller(void) {
		for (int i = 0; i < iNumDrives; i++) {
			SAFEDELETE(pDrivesArray[i].pDriveCaller);
		}
		SAFEDELETEARR(pDrivesArray);
	};

	/* copia */
	virtual TplDriveCaller<doublereal>* pCopy(void) const {
		typedef DrivesArray<doublereal> da;
		da* pDA = 0;

		SAFENEWARR(pDA, da, iNumDrives);

		for (int i = 0; i < iNumDrives; i++) {
			pDA[i].pDriveCaller = pDrivesArray[i].pDriveCaller->pCopy();
			pDA[i].t = pDrivesArray[i].t;
		}

		typedef ArrayTplDriveCaller<doublereal> dc;
		TplDriveCaller<doublereal>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(iNumDrives, pDA));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "array, " << iNumDrives;
		for (int i = 0; i < iNumDrives; i++) {
			out << ", ", pDrivesArray[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (int i = 0; i < iNumDrives; i++) {
			out << ", ", pDrivesArray[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	inline doublereal Get(void) const {
		doublereal v = 0.;
		for (int i = 0; i < iNumDrives; i++) {
			v += pDrivesArray[i].pDriveCaller->dGet();
		}
		return v;
	};

	inline int getNDrives(void) const {
		return iNumDrives;
	};
};

/* Nota: di questa classe non viene scritta esplicitamente la versione
 *       per reali in quanto il coefficiente moltiplicativo
 *       puo' essere usato per un'ulteriore pesatura del drive */

/* ArrayTplDriveCaller - end */


extern doublereal GetT(MBDynParser& HP, const doublereal& t);

template <class T> T GetT(MBDynParser& HP, const T& t)
{
	return HP.Get(t);
}


/* enum delle parole chiave */
template <class T>
class ReadTplDriveKeyWords {
public:
	enum KeyWords {
		UNKNOWNTPLDRIVE = -1,
		NULLTPLDRIVE = 0,
			ZERO,		/* deprecated */
		SINGLE,
		ARRAY,
		LASTKEYWORD
	};
};


template <class T>
TplDriveCaller<T>* ReadTplDrive(const DataManager* pDM,
	MBDynParser& HP,
	const T& t)
{
	DEBUGCOUT("Entering ReadTplDrive" << std::endl);

	const char* sKeyWords[] = {
		"null",
			"zero",		/* deprecated */
		"single",
		"array",
		0
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	TplDriveCaller<T>* pTplDC = 0;

	/* Nota: siccome capita che non ci sia la deformazione imposta,
	 * faccio chiedere se e' desiderata */
	if (!HP.IsArg()) {
		DriveCaller* pDC = 0;
		SAFENEW(pDC, NullDriveCaller);

		T t(0.);
		SAFENEWWITHCONSTRUCTOR(pTplDC,
			SingleTplDriveCaller<T>,
			SingleTplDriveCaller<T>(pDC, t));

	} else { /* c'e' la deformazione */
		int CurrKeyWord = HP.IsKeyWord();

		/* valore di default (pericoloso!) */
		if (CurrKeyWord == ReadTplDriveKeyWords<T>::UNKNOWNTPLDRIVE) {
			CurrKeyWord = ReadTplDriveKeyWords<T>::SINGLE;
		}

restart:
		switch (CurrKeyWord) {
		case ReadTplDriveKeyWords<T>::ZERO:
			silent_cerr("\"zero\" is deprecated "
				"for null template drive "
				"at line " << HP.GetLineData() << ";"
				" use \"null\" instead." << std::endl);
			/* fallthru */
		case ReadTplDriveKeyWords<T>::NULLTPLDRIVE:
			SAFENEW(pTplDC, ZeroTplDriveCaller<T>);
			break;

		case ReadTplDriveKeyWords<T>::SINGLE: {
			T t(0.);

			t = GetT(HP, t);

			DriveCaller* pDC = HP.GetDriveCaller();

			SAFENEWWITHCONSTRUCTOR(pTplDC,
				SingleTplDriveCaller<T>,
				SingleTplDriveCaller<T>(pDC, t));
			} break;

		case ReadTplDriveKeyWords<T>::ARRAY: {
			unsigned short int iNumDr = HP.GetInt();
			if (iNumDr == 0) {
				silent_cerr("At least one drive is required "
					"in array template drive" << std::endl);
				throw ErrGeneric();

			} else if (iNumDr == 1) {
				CurrKeyWord = ReadTplDriveKeyWords<T>::SINGLE;
				goto restart;
			} /* else */


			DrivesArray<T>* pDA = 0;
			SAFENEWARR(pDA, DrivesArray<T>, iNumDr);

			for (unsigned short int i = 0; i < iNumDr; i++) {
				T t(0.);
				pDA[i].t = GetT(HP, t);
				pDA[i].pDriveCaller = HP.GetDriveCaller();
			}

			SAFENEWWITHCONSTRUCTOR(pTplDC,
				ArrayTplDriveCaller<T>,
				ArrayTplDriveCaller<T>(iNumDr, pDA));

			} break;

		default:
			silent_cerr("Unknown template drive type" << std::endl);
			throw ErrGeneric();
		}
	}

	if (pTplDC == 0) {
		silent_cerr("Error in allocation of template drive" << std::endl);
		throw ErrMemory();
	}

	return pTplDC;
}

#endif /* TPLDRIVE__H */

