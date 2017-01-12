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
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

class ModuleMDS
: virtual public Elem, public UserDefinedElem {
private:
	// add private data
	doublereal dM, dD, dK;
	doublereal dX, dV, dVP;
	DriveOwner F;

public:
	ModuleMDS(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleMDS(void);

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
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;
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

ModuleMDS::ModuleMDS(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	msd							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"	All rights reserved						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	dM = HP.GetReal();
	dD = HP.GetReal();
	dK = HP.GetReal();

	F.Set(HP.GetDriveCaller());

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleMDS::~ModuleMDS(void)
{
	// destroy private data
	NO_OP;
}

void
ModuleMDS::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel()	// 1:	label
			<< " " << dX			// 
			<< " " << dV			// 
			<< " " << dVP			// 
			<< std::endl;
	}
}

void
ModuleMDS::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 2;
}

VariableSubMatrixHandler& 
ModuleMDS::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	WM.Resize(2, 2);

	integer iFirstIndex = iGetFirstIndex();

	WM.PutRowIndex(1, iFirstIndex + 1);
	WM.PutColIndex(1, iFirstIndex + 1);
	WM.PutRowIndex(2, iFirstIndex + 2);
	WM.PutColIndex(2, iFirstIndex + 2);

	/*
	 *	        [ -m   0 ]
	 *	r_/y' = [        ]
	 *	        [  0  -1 ]
	 *
	 *	       [ -d  -k ]
	 *	r_/y = [        ]
	 *	       [  1   0 ]
	 *
	 *	                            [ (m + dCoef*d) dCoef*k ]
	 *	J = -(r_/y' + dCoef*r_/y) = [                       ]
	 *	                            [    -dCoef        1.   ]
	 */

#if 1
	// old style
	WM.PutCoef(1, 1, dM + dCoef*dD);
	WM.PutCoef(1, 2, dCoef*dK);
	WM.PutCoef(2, 1, -dCoef);
	WM.PutCoef(2, 2, 1.);
#endif

#if 0
	// cleaner (less efficient?)
	WM(1, 1) = dM + dCoef*dD;
	WM(1, 2) = dCoef*dK;
	WM(2, 1) = -dCoef;
	WM(2, 2) = 1.;
#endif

	return WorkMat;
}

SubVectorHandler& 
ModuleMDS::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkVec.ResizeReset(2);

	integer iFirstIndex = iGetFirstIndex();

	WorkVec.PutRowIndex(1, iFirstIndex + 1);
	WorkVec.PutRowIndex(2, iFirstIndex + 2);

	dVP = XPrimeCurr(iFirstIndex + 1);
	doublereal dXP = XPrimeCurr(iFirstIndex + 2);
	dV = XCurr(iFirstIndex + 1);
	dX = XCurr(iFirstIndex + 2);

	/*
	 *	    { f(t) - m*v' - d*v - d*k }
	 *	r = {                         }
	 *	    {         v - x'          }
	 */

#if 0
	WorkVec.PutCoef(1, F.dGet() - dM*dVP - dD*dV - dK*dX);
	WorkVec.PutCoef(2, dV - dXP);
#endif

	WorkVec(1) = F.dGet() - dM*dVP - dD*dV - dK*dX;
	WorkVec(2) = dV - dXP;

	return WorkVec;
}

unsigned int
ModuleMDS::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleMDS::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ModuleMDS::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ModuleMDS::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

unsigned int
ModuleMDS::iGetNumDof(void) const
{
	return 2;
}

DofOrder::Order
ModuleMDS::GetDofType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

DofOrder::Order
ModuleMDS::GetEqType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

std::ostream&
ModuleMDS::Restart(std::ostream& out) const
{
	return out << "# ModuleMDS: not implemented" << std::endl;
}

unsigned int
ModuleMDS::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleMDS::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleMDS::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleMDS::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf = new UDERead<ModuleMDS>;

	if (!SetUDE("mds", rf)) {
		delete rf;

		silent_cerr("module-mds: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

