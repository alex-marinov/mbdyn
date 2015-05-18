/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

// ModuleBullet: begin

class ModuleBullet
: virtual public Elem, public UserDefinedElem {
private:
	// add private data

public:
	ModuleBullet(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleBullet(void);

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

ModuleBullet::ModuleBullet(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	bullet							\n"
"Author: 	Vivek Kumar <vivekkumar0893@gmail.com>			\n"
"		Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
"									\n"
"All rights reserved							\n"
"									\n"
"Syntax:								\n"
"	user defined : <label> , bullet , ... ;				\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleBullet::~ModuleBullet(void)
{
	NO_OP;
}

void
ModuleBullet::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		// std::ostream& out = OH.Loadable();

		// write output as appropriate
	}
}

void
ModuleBullet::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	// set *piNumRows to 6 * number_of_nodes
	// set *piNumCols to 1 for residual only
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleBullet::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleBullet::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);

	// compute contact forces and write them in contribution to residual

	return WorkVec;
}

unsigned int
ModuleBullet::iGetNumPrivData(void) const
{
	// return number of private data
	return 0;
}

unsigned int
ModuleBullet::iGetPrivDataIdx(const char *s) const
{
	// parse string and compute index of requested private data

	// shouldn't get here until private data are defined
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
ModuleBullet::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 1 && i <= iGetNumPrivData());

	// compute requested private data

	// shouldn't get here until private data are defined
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
ModuleBullet::iGetNumConnectedNodes(void) const
{
	// return number of connected nodes
	return 0;
}

void
ModuleBullet::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	// resize to number of connected nodes
	// fill array with connected nodes
	connectedNodes.resize(0);
	// connectedNodes[0] = m_pNode;
}

void
ModuleBullet::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleBullet::Restart(std::ostream& out) const
{
	return out << "# ModuleBullet: not implemented" << std::endl;
}

unsigned int
ModuleBullet::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleBullet::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleBullet::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleBullet::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

// ModuleBullet: end

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<ModuleBullet>;

	if (!SetUDE("bullet", rf1)) {
		delete rf1;

		silent_cerr("ModuleBullet: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

