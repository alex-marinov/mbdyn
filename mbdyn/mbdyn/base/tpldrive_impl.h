/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef TPLDRIVE_IMPL_H
#define TPLDRIVE_IMPL_H

#include <typeinfo>
#include "tpldrive.h"
#include "mbpar.h"

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

	inline T Get(const doublereal& dVar) const {
		return mb_zero<T>();
	};

	inline T Get(void) const {
		return mb_zero<T>();
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		return true;
	};

	inline T GetP(void) const {
		return mb_zero<T>();
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
	: DriveOwner(pDC), t(const_cast<T&>(x)) {
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

	inline T Get(const doublereal& dVar) const {
		return t*dGet(dVar);
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

	inline doublereal Get(const doublereal& dVar) const {
		return dGet(dVar);
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

template <class T>
TplDriveCaller<T> *
DC2TDC(DriveCaller *pDC, T& t)
{
	typedef SingleTplDriveCaller<T> STDC;

	TplDriveCaller<T> *p = 0;
	SAFENEWWITHCONSTRUCTOR(p, STDC, STDC(pDC, t));
	return p;
}


/* SingleTplDriveCaller - end */

#endif // TPLDRIVE_IMPL_H

