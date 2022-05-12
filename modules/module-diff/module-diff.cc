/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include <iostream>
#include <string>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

class ModuleDiffDrive 
: virtual public Elem, public UserDefinedElem, public DriveOwner {
public:
	enum Method 
	{
		X2MV,
		X2PVD3,
		ALPHAXBETAV,
		IMPLICIT_EULER
	};
private:
	// add private data
	integer iNumDrives;
	std::vector<DriveCaller*> m_pDCs;
	const DataManager* m_pDM;
	Method m_diffMethod;
	doublereal m_dAlpha;
	doublereal m_dBeta;
	doublereal m_dTimePrev;
	std::vector<doublereal> m_dX;
	std::vector<doublereal> m_dXP;
	std::vector<doublereal> m_dXPrev;
	std::vector<doublereal> m_dXPPrev;
	std::vector<doublereal> m_dXPInitial;
public:
	ModuleDiffDrive(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleDiffDrive(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);
	virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
	unsigned int iGetNumPrivData(void) const;
	unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

ModuleDiffDrive::ModuleDiffDrive(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
DriveOwner(0),
m_pDM(pDM),
iNumDrives(0),
m_dAlpha(0),
m_dBeta(0),
m_dTimePrev(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	diff drive						\n"
"Author: 	Andrea Zanoni <andrea.zanoni@polimi.it>			\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
"									\n"
"USAGE:									\n"
"diff drive,\n"
"	(integer)<num_drives>,\n"
"	(DriveCaller)<drive_1>,\n"
"	[(DriveCaller)<drive_2>,]\n"
"	[...]\n"
"	[(DriveCaller)<drive_n>,]\n"
"	method, {implicit euler | x2mv| x2pvd3 | alphaxbetav, alpha, beta,}\n"
"	[, initial value, (real)<initial_value_1>]\n"
"	[, initial value, (real)<initial_value_1>]\n"
"	[, ...]\n"
"	[, initial value, (real)<initial_value_n>]\n"
"									\n"
" Private data:								\n"
"	d<N>	echo of input drive N					\n"
"	dp<N>	computed derivative of drive N				\n"
"									\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	iNumDrives = HP.GetInt();
	if (iNumDrives <= 0) {
		silent_cerr("ModuleDiffDrive: invalid number of drives provided" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	m_pDCs.resize(iNumDrives);
	m_dX.resize(iNumDrives, 0);
	m_dXP.resize(iNumDrives, 0);
	m_dXPrev.resize(iNumDrives, 0);
	m_dXPPrev.resize(iNumDrives, 0);
	for (std::vector<DriveCaller*>::iterator it = m_pDCs.begin(); it != m_pDCs.end(); ++it) {
		*it = HP.GetDriveCaller();
	}

	if (HP.IsKeyWord("method")) {
		if (HP.IsKeyWord("x2mv")) {
			m_diffMethod = ModuleDiffDrive::X2MV;
		} else if (HP.IsKeyWord("x2pvd3")) {
			m_diffMethod = ModuleDiffDrive::X2PVD3;
		} else if (HP.IsKeyWord("alphaxbetav")) {
			m_diffMethod = ModuleDiffDrive::ALPHAXBETAV;
			m_dAlpha = HP.GetReal();
			m_dBeta = HP.GetReal();
			if (m_dAlpha < 0) {
				silent_cerr("ModuleDiffDrive: invalid value for Alpha." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			if (m_dBeta < 0) {
				silent_cerr("ModuleDiffDrive: invalid value for Beta." << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			
		} else if (HP.IsKeyWord("implicit" "euler")) {
			m_diffMethod = ModuleDiffDrive::IMPLICIT_EULER;
		} else {
			silent_cerr("ModuleDiffDrive: unrecognised method type." << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	} else {
		silent_cerr("ModuleDiffDrive: method specification expected." << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_dXPInitial.resize(iNumDrives);
	integer jj = 0;
	while (jj < iNumDrives) {
		if (HP.IsKeyWord("initial" "value")) {
			m_dXPInitial[jj] = HP.GetReal();
			jj = jj + 1;
		} else {
			break;
		}
	}

	// initialize previous time with the time of element instantiation
	m_dTimePrev = pDM->dGetTime();
}

void
ModuleDiffDrive::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleDiffDrive::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleDiffDrive::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// NO_OP: nothing to be done
	WorkVec.ResizeReset(0);

	return WorkVec;
}

unsigned int
ModuleDiffDrive::iGetNumPrivData(void) const
{
	return 2*iNumDrives;
}

unsigned int
ModuleDiffDrive::iGetPrivDataIdx(const char *s) const
{
	std::string ss(s);
	unsigned i = ss.find("d");
	if (i != std::string::npos) {
		if ((ss.at(i + 1)) == 'p') {
			std::string snum = ss.substr(i + 2);
			return (std::stoi(snum) + iNumDrives);
		} else {
			std::string snum = ss.substr(i + 1);
			return std::stoi(snum);
		}
	}
	silent_cerr("ModuleDiffDrive: " << s << "is an invalid name for a private datum." << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
ModuleDiffDrive::dGetPrivData(unsigned int i) const
{
	if (i > 0 && i <= iNumDrives) {
		return m_dX[i - 1]; 
	} else if (i > iNumDrives && i <= 2*iNumDrives) {
		return m_dXP[i - iNumDrives - 1];
	}
	silent_cerr("ModuleDiffDrive: invalid private data index." << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
ModuleDiffDrive::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ModuleDiffDrive::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ModuleDiffDrive::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleDiffDrive::Restart(std::ostream& out) const
{
	return out << "# ModuleDiffDrive: not implemented" << std::endl;
}

unsigned int
ModuleDiffDrive::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleDiffDrive::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleDiffDrive::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// NO_OP: nothing to be done
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleDiffDrive::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// NO_OP: nothing to be done
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}


ModuleDiffDrive::~ModuleDiffDrive(void)
{
	// destroy private data
	for (std::vector<DriveCaller*>::iterator it = m_pDCs.begin(); it != m_pDCs.end(); ++it) {
		delete *it;
	}
}

void
ModuleDiffDrive::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		if (OH.UseText(OutputHandler::LOADABLE)) {
			std::ostream& out = OH.Loadable();
			out << std::setw(8) << GetLabel();	// 1 label
			for (integer ii = 0; ii < iNumDrives; ii++)
			{
				out << " " << m_dX[ii];		// even columns: drives echos
				out << " " << m_dXP[ii];	// odd columns: drive diffs
			}
			out << std::endl;
		}
		if (OH.UseBinary(OutputHandler::LOADABLE)) {
		// TODO
		}
	}
}

void
ModuleDiffDrive::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	// get elapsed time since last call: should be the timestep
	// FIXME: can we integrate backward in time? CHECK!
	doublereal h = m_pDM->dGetTime() - m_dTimePrev;
	DEBUGLCOUT(MYDEBUG_SOL, "ModuleDiffDrive: h = " << m_pDM->dGetTime() << " - " << m_dTimePrev << " = " << h << std::endl);
	m_dTimePrev = m_pDM->dGetTime();
	ASSERT(h > 0);
	if (h <= 0) {
		silent_cerr("ModuleDiffDrive: null or negative computed timestep." << std::endl);
		return;
	}

	// compute derivative according to selected method
	switch (m_diffMethod) {
		case X2MV:
			for (integer ii = 0; ii < iNumDrives; ii++) {
				m_dX[ii] = m_pDCs[ii]->dGet();
				m_dXP[ii] = (m_dX[ii] - m_dXPrev[ii])*(2./h) - m_dXPPrev[ii];
				m_dXPrev[ii] = m_dX[ii];
				m_dXPPrev[ii] = m_dXP[ii];
			}
			break;
		case X2PVD3:
			for (integer ii = 0; ii < iNumDrives; ii++) {
				m_dX[ii] = m_pDCs[ii]->dGet();
				m_dXP[ii] = (m_dX[ii] - m_dXPrev[ii])*(2./3./h) - m_dXPPrev[ii]/3.;
				m_dXPrev[ii] = m_dX[ii];
				m_dXPPrev[ii] = m_dXP[ii];
			}
			break;
		case ALPHAXBETAV:
			for (integer ii = 0; ii < iNumDrives; ii++) {
				m_dX[ii] = m_pDCs[ii]->dGet();
				m_dXP[ii] = (m_dX[ii] - m_dXPrev[ii])*(m_dAlpha/h) - m_dXPPrev[ii]*m_dBeta;
				m_dXPrev[ii] = m_dX[ii];
				m_dXPPrev[ii] = m_dXP[ii];
			}
			break;
		case IMPLICIT_EULER:
			for (integer ii = 0; ii < iNumDrives; ii++) {
				m_dX[ii] = m_pDCs[ii]->dGet();
				m_dXP[ii] = (m_dX[ii] - m_dXPrev[ii])/h;
				m_dXPrev[ii] = m_dX[ii];
				m_dXPPrev[ii] = m_dXP[ii];
			}
			break;
	}
}

bool
diffdrive_set(void)
{
	UserDefinedElemRead *rf = new UDERead<ModuleDiffDrive>;

	if (!SetUDE("diffdrive", rf)) {
		delete rf;

		silent_cerr("ModuleDiffDrive: diffdrive_set() failed" << std::endl);

		return false;
	}

	return true;
}

#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ModuleDiffDrive>;

	if (!diffdrive_set()) {
		delete rf;

		silent_cerr("ModuleDiffDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}
#endif // ! STATIC_MODULES
