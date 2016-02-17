/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
#include <cstring>
#include <cstdlib>
#include <vector>

#include "dataman.h"
#include "userelem.h"

class ModuleController
: virtual public Elem, public UserDefinedElem {
private:
	std::vector<DriveCaller *> measures;

	// outputs = A * measures
	FullMatrixHandler A;

public:
	ModuleController(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleController(void);

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
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
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

ModuleController::ModuleController(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	controller						\n"
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

	if (!HP.IsKeyWord("measures")) {
		silent_cerr("ModuleController(" << uLabel << "): "
			"\"measures\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer nm = HP.GetInt();
	if (nm <= 0) {
		silent_cerr("ModuleController(" << uLabel << "): "
			"invalid measures number at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	measures.resize(nm);
	for (std::vector<DriveCaller *>::iterator i = measures.begin();
		i != measures.end(); ++i)
	{
		*i = HP.GetDriveCaller(false);
	}

	if (!HP.IsKeyWord("outputs")) {
		silent_cerr("ModuleController(" << uLabel << "): "
			"\"outputs\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer no = HP.GetInt();
	if (no <= 0) {
		silent_cerr("ModuleController(" << uLabel << "): "
			"invalid outputs number at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	A.Resize(no, nm);

	for (integer ir = 1; ir <= no; ir++) {
		for (integer ic = 1; ic <= nm; ic++) {
			A(ir, ic) = HP.GetReal();
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleController::~ModuleController(void)
{
	// destroy private data
	for (std::vector<DriveCaller *>::iterator i = measures.begin();
		i != measures.end(); ++i)
	{
		delete *i;
	}
}

void
ModuleController::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		std::ostream& out = OH.Loadable();

		out << std::setw(8) << GetLabel();

		for (std::vector<DriveCaller *>::const_iterator i = measures.begin();
			i != measures.end(); ++i)
		{
			out << " " << (*i)->dGet();
		}

		for (integer i = 1; i <= A.iGetNumRows(); i++) {
			out << " " << dGetPrivData(i);
		}

		out << std::endl;
	}
}

void
ModuleController::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleController::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleController::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

unsigned int
ModuleController::iGetNumPrivData(void) const
{
	return A.iGetNumRows();
}

unsigned int
ModuleController::iGetPrivDataIdx(const char *s) const
{
	if (s[0] != 'x') {
		return 0;
	}
	s++;

	if (s[0] != '[') {
		return 0;
	}

	s++;
	if (!std::isdigit(s[0])) {
		return 0;
	}

	char *next;
	unsigned long idx = std::strtoul(s, &next, 10);
	if (next == s || next[0] != ']') {
		return 0;
	}
	s = next + 1;

	if (s[0] != '\0') {
		return 0;
	}

	return idx;
}

doublereal
ModuleController::dGetPrivData(unsigned int i) const
{
	doublereal d = 0.;
	for (integer ic = 1; ic <= A.iGetNumCols(); ic++) {
		d += A(i, ic)*measures[ic - 1]->dGet();
	}

	return d;
}

int
ModuleController::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
ModuleController::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
ModuleController::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleController::Restart(std::ostream& out) const
{
	return out << "# ModuleController: not implemented" << std::endl;
}

unsigned int
ModuleController::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleController::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleController::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleController::InitialAssRes(
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
	UserDefinedElemRead *rf = new UDERead<ModuleController>;

	if (!SetUDE("controller", rf)) {
		delete rf;

		silent_cerr("module-controller: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

