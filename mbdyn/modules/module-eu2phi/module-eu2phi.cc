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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "tpldrive.h"
#include "Rot.hh"

class Eu2PhiWrap : public TplDriveCaller<Vec3> {
private:
	TplDriveCaller<Vec3> *m_pDC;
	OrientationDescription m_od;

public:
	Eu2PhiWrap(TplDriveCaller<Vec3>* pDC, OrientationDescription od) : m_pDC(pDC), m_od(od) {
		NO_OP;
	};

	virtual ~Eu2PhiWrap(void) {
		delete m_pDC;
	};
 
	/* copia */
	virtual TplDriveCaller<Vec3>* pCopy(void) const {
		return new Eu2PhiWrap(m_pDC->pCopy(), m_od);
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		return out;
	};

	/* Restituisce il valore del driver */
	virtual Vec3 Get(const doublereal& dVar) const {
		Vec3 E(m_pDC->Get(dVar));

		Mat3x3 R;
		switch (m_od) {
		case EULER_123:
			R = EulerAngles123_2MatR(E);
			break;

		case EULER_313:
			R = EulerAngles313_2MatR(E);
			break;

		case EULER_321:
			R = EulerAngles321_2MatR(E);
			break;

		default:
			ASSERT(0);
		}

		return RotManip::VecRot(R);
	};

	virtual Vec3 Get(void) const {
		return Get(0.);
	};

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const {
		return false;
	};
	virtual Vec3 GetP(void) const {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

	virtual int getNDrives(void) const {
		return 1;
	};
};

/* prototype of the functional object: reads a drive caller */
struct Eu2PhiDCR : public TplDriveCallerRead<Vec3> {
	virtual TplDriveCaller<Vec3> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		if (HP.IsKeyWord("help")) {
			silent_cout(
"Eu2PhiWrap: converts a TplDriveCaller<Vec3> containing three Euler angles\n"
"into the corresponding Euler vector\n"
"\n"
"Syntax:\n"
"    eu2phi ,\n"
"        [ help , ]\n"
"        [ format , { euler123 | euler313 | euler321 } , ]\n"
"        (TplDriveCaller<Vec3>) <drive>\n");
		}

		OrientationDescription od(EULER_123);

		if (HP.IsKeyWord("format")) {
			od = ReadOrientationDescription(HP);
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
				break;

			default:
				silent_cerr("Eu2PhiWrap: unhandled format at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		return new Eu2PhiWrap(HP.GetTplDriveCaller<Vec3>(), od);
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	TplDriveCallerRead<Vec3> *rf = new Eu2PhiDCR;

	if (!SetDC3D("eu2phi", rf)) {
		delete rf;

		silent_cerr("Eu2PhiDCR: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

