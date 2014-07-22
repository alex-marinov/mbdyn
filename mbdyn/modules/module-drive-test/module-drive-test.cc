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
#include "drive.h"

class DriveTestCaller : public DriveCaller {
private:
	DriveCaller *m_pDC;
	std::string m_fname;
	doublereal m_dInitialTime;
	doublereal m_dFinalTime;
	doublereal m_dTimeStep;

public:
	DriveTestCaller(DriveCaller *pDC,
		const std::string& fname,
		doublereal dInitialTime,
		doublereal dFinalTime,
		doublereal dTimeStep);
	virtual ~DriveTestCaller(void);
 
	/* Copia */
	virtual DriveCaller* pCopy(void) const;
 
	/* Scrive il contributo del DriveCaller al file di restart */   
	virtual std::ostream& Restart(std::ostream& out) const;
 
	inline doublereal dGet(const doublereal& dVar) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
};

DriveTestCaller::DriveTestCaller(DriveCaller *pDC,
		const std::string& fname,
	doublereal dInitialTime,
	doublereal dFinalTime,
	doublereal dTimeStep)
: DriveCaller(0),
m_pDC(pDC),
m_fname(fname),
m_dInitialTime(dInitialTime),
m_dFinalTime(dFinalTime),
m_dTimeStep(dTimeStep)
{
	NO_OP;
}

DriveTestCaller::~DriveTestCaller(void)
{
	NO_OP;
}

DriveCaller *
DriveTestCaller::pCopy(void) const
{
	return new DriveTestCaller(m_pDC->pCopy(), m_fname,
		m_dInitialTime, m_dFinalTime, m_dTimeStep);
}

std::ostream&
DriveTestCaller::Restart(std::ostream& out) const
{
	return out << "drive test, ", m_pDC->Restart(out)
		<< ", \"" << m_fname << "\""
		<< ", " << m_dInitialTime
		<< ", " << m_dFinalTime
		<< ", " << m_dTimeStep;
}

inline doublereal 
DriveTestCaller::dGet(const doublereal& dVar) const 
{
	std::ofstream out(m_fname.c_str());

	for (doublereal d = m_dInitialTime; d <= m_dFinalTime; d += m_dTimeStep) {
		out << d << " " << m_pDC->dGet(d) << std::endl;
	}

	return m_pDC->dGet(dVar); 
}

inline bool
DriveTestCaller::bIsDifferentiable(void) const
{
	return m_pDC->bIsDifferentiable();
}

inline doublereal 
DriveTestCaller::dGetP(const doublereal& dVar) const
{
	ASSERT(bIsDifferentiable());

	std::ofstream out(m_fname.c_str());

	for (doublereal d = m_dInitialTime; d <= m_dFinalTime; d += m_dTimeStep) {
		out << d << " " << m_pDC->dGetP(d) << std::endl;
	}

	return m_pDC->dGetP(dVar);;
}

/* prototype of the functional object: reads a drive caller */
struct DTDCR : public DriveCallerRead {
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) {
		DriveCaller *pDC = HP.GetDriveCaller();
		std::string fname(HP.GetFileName());
		doublereal dInitialTime(HP.GetReal());
		doublereal dFinalTime(HP.GetReal());
		doublereal dTimeStep(HP.GetReal());
		return new DriveTestCaller(pDC, fname, dInitialTime, dFinalTime, dTimeStep);
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	DriveCallerRead	*rf = new DTDCR;

	if (!SetDriveData("drive" "test", rf)) {
		delete rf;

		silent_cerr("DriveTest: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

