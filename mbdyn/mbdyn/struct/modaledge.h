/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2007-2015
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

/* Forza */

#ifndef MODALEDGE_H
#define MODALEDGE_H

#include "modalext.h"
#include "aerodyn.h"

/* ExtForceEDGE - begin */

class ExtForceEDGE : public ExtModalForceBase {
protected:
	const AirProperties *pAP;

public:
	ExtForceEDGE(DataManager *pDM);
	virtual ~ExtForceEDGE(void);

	bool
	Prepare(ExtFileHandlerBase *pEFH, unsigned uLabel, bool bRigid, unsigned uModes);
};

/* ExtForceEDGE - end */

/* ExtRigidForceEDGE - begin */

class ExtRigidForceEDGE : public ExtForceEDGE {
public:
	ExtRigidForceEDGE(DataManager *pDM);

	unsigned
	Recv(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& a);

	void
	Send(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP);
};

/* ExtRigidForceEDGE - end */

/* ExtModalForceEDGE - begin */

class ExtModalForceEDGE : public ExtForceEDGE {
public:
	ExtModalForceEDGE(DataManager *pDM);

	unsigned
	Recv(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned& uLabel,
		Vec3& f, Vec3& m, std::vector<doublereal>& a);

	void
	Send(ExtFileHandlerBase *pEFH, unsigned uFlags, unsigned uLabel,
		const Vec3& x, const Mat3x3& R, const Vec3& v, const Vec3& w,
		const std::vector<doublereal>& q,
		const std::vector<doublereal>& qP);
};

/* ExtModalForceEDGE - end */

#endif // MODALEDGE_H

