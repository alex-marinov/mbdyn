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

#ifndef TPLDRIVE_IMPL_H
#define TPLDRIVE_IMPL_H

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

template <class T>
TplDriveCaller<T> *
ReadTplDC(const DataManager* pDM, MBDynParser& HP, const T& t)
{
	// FIXME: ugly!
	if (typeid(T) == typeid(doublereal)) {
		return static_cast<TplDriveCaller<T> *>((void *)ReadDC1D(pDM, HP));

	} else if (typeid(T) == typeid(Vec3)) {
		return static_cast<TplDriveCaller<T> *>((void *)ReadDC3D(pDM, HP));

	} else if (typeid(T) == typeid(Vec6)) {
		return static_cast<TplDriveCaller<T> *>((void *)ReadDC6D(pDM, HP));
	}

	silent_cerr("unknown type in ReadTplDC "
		"at line " << HP.GetLineData() << std::endl);
	throw ErrGeneric();
}

#endif // TPLDRIVE_IMPL_H

