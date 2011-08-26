/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#ifndef TPLDRIVE_H
#define TPLDRIVE_H

#include "drive.h"
#include "matvec3.h"
#include "matvec6.h"

// forward declaration, since "reffrm.h" needs "tpldrive.h"
class ReferenceFrame;

/* TplDriveCaller - begin */

template <class T>
class TplDriveCaller {
public:
	virtual ~TplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const = 0;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;
	virtual std::ostream& Restart_int(std::ostream& out) const = 0;

	/* Restituisce il valore del driver */
	virtual T Get(const doublereal& dVar) const = 0;
	virtual T Get(void) const {
		return Get(0.);
	};

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const {
		return false;
	};
	virtual T GetP(void) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

	virtual int getNDrives(void) const = 0;
};

/* TplDriveCaller - end */


/* TplDriveOwner - begin */

template <class T>
class TplDriveOwner {
protected:
	mutable TplDriveCaller<T>* pTplDriveCaller;

public:
	TplDriveOwner(const TplDriveCaller<T>* pDC = 0)
	: pTplDriveCaller((TplDriveCaller<T>*)pDC) {
		NO_OP;
	};

	virtual ~TplDriveOwner(void) {
		if (pTplDriveCaller != 0) {
			SAFEDELETE(pTplDriveCaller);
		}
	};

	void Set(const TplDriveCaller<T>* pDC) {
		ASSERT(pDC != 0);
		if (pTplDriveCaller != 0) {
			SAFEDELETE(pTplDriveCaller);
		}
		pTplDriveCaller = const_cast<TplDriveCaller<T>*>(pDC);
	};

	TplDriveCaller<T>* pGetDriveCaller(void) const {
		return pTplDriveCaller;
	};

	T Get(const doublereal& dVar) const {
		return pTplDriveCaller->Get(dVar);
	};

	T Get(void) const {
		return pTplDriveCaller->Get();
	};

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const {
		return pTplDriveCaller->bIsDifferentiable();
	};
	virtual T GetP(void) const {
		return pTplDriveCaller->GetP();
	};
};

/* TplDriveOwner - end */

/* functions that read a template drive caller */
extern TplDriveCaller<doublereal> *
ReadDC1D(const DataManager* pDM, MBDynParser& HP);
extern TplDriveCaller<Vec3> *
ReadDC3D(const DataManager* pDM, MBDynParser& HP);
extern TplDriveCaller<Vec6> *
ReadDC6D(const DataManager* pDM, MBDynParser& HP);

extern TplDriveCaller<Vec3> *
ReadDCVecRel(const DataManager* pDM, MBDynParser& HP, const ReferenceFrame& rf);
extern TplDriveCaller<Vec3> *
ReadDCVecAbs(const DataManager* pDM, MBDynParser& HP, const ReferenceFrame& rf);

/* prototype of the template functional object: reads a template drive caller */
template <class T>
struct TplDriveCallerRead {
	virtual ~TplDriveCallerRead<T>( void ) { NO_OP; };
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) = 0;
};

/* template drive caller registration functions: call to register one */
extern bool
SetDC1D(const char *name, TplDriveCallerRead<doublereal> *rf);
extern bool
SetDC3D(const char *name, TplDriveCallerRead<Vec3> *rf);
extern bool
SetDC6D(const char *name, TplDriveCallerRead<Vec6> *rf);

/* create/destroy */
extern void InitTplDC(void);
extern void DestroyTplDC(void);

#endif /* TPLDRIVE_H */

