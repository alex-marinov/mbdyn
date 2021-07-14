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
#include "indvel.h"

class ModuleIndVel
: virtual public Elem, public UserDefinedElem, public InducedVelocity {
private:
	// TODO: define per-point structure
	struct PointData {
		const Elem *pEl;
		unsigned counter;
		Vec3 F;
		Vec3 M;
		doublereal dW;
		Vec3 X;

	public:
		PointData(void) { NO_OP; };
		~PointData(void) { NO_OP; };
	};

	std::vector<PointData> m_data;
	typedef std::vector<PointData> PD;
	PD::iterator m_data_iter;

	// add private data
	int iFirstAssembly;

public:
	ModuleIndVel(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~ModuleIndVel(void);

	virtual Elem::Type GetElemType(void) const;

	// induced velocity specific calls
	virtual InducedVelocity::Type GetInducedVelocityType(void) const;
	virtual bool bSectionalForces(void) const;
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3&) const;
	virtual void AddSectionalForce(Elem::Type type,
		const Elem *pEl, unsigned uPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W);

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

ModuleIndVel::ModuleIndVel(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
InducedVelocity(uLabel, 0, 0, flag(0)),
iFirstAssembly(2)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	IndVel							\n"
"Author: 	Pierangelo Masarati <masarati@aero.polimi.it>		\n"
"Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements a dummy induced velocity model.	\n"
"									\n"
"	All rights reserved.						\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// do something useful
	pCraft = pDM->ReadNode<StructNode, Node::STRUCTURAL>(HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleIndVel::~ModuleIndVel(void)
{
	// destroy private data
	NO_OP;
}

Elem::Type
ModuleIndVel::GetElemType(void) const
{
	return Elem::LOADABLE;
}

InducedVelocity::Type
ModuleIndVel::GetInducedVelocityType(void) const
{
	return InducedVelocity::USER_DEFINED;
}

bool
ModuleIndVel::bSectionalForces(void) const
{
	return true;
}

Vec3
ModuleIndVel::GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return Zero3;
}

void
ModuleIndVel::AddSectionalForce(Elem::Type type,
	const Elem *pEl, unsigned uPnt,
	const Vec3& F, const Vec3& M, doublereal dW,
	const Vec3& X, const Mat3x3& R,
	const Vec3& V, const Vec3& W)
{
	std::cerr << "ModuleIndVel(" << GetLabel() << ")::AddSectionalForce: "
		<< psElemNames[type] << "(" << pEl->GetLabel() << "):" << uPnt << std::endl;

	if (iFirstAssembly == 1) {
		unsigned idx = m_data.size();
		m_data.resize(idx + 1);
		m_data[idx].pEl = pEl;
		if (idx == 0 || (idx > 0 && pEl != m_data[idx - 1].pEl)) {
			m_data[idx].counter = 0;
		} else {
			m_data[idx].counter = m_data[idx - 1].counter + 1;
		}

		m_data[idx].F = Zero3;
		m_data[idx].M = Zero3;
		m_data[idx].dW = 0.;
		m_data[idx].X = Zero3;

	} else {
		const Mat3x3& R(pCraft->GetRCurr());

		// resolve force, moment and point in craft's reference frame
		m_data_iter->F = R.MulTV(F);
		m_data_iter->M = R.MulTV(M);
		m_data_iter->dW = dW;
		m_data_iter->X = R.MulTV(X - pCraft->GetXCurr());

		++m_data_iter;
	}
}

void
ModuleIndVel::Output(OutputHandler& OH) const
{
	std::cerr << "ModuleIndVel(" << GetLabel() << ")::Output: size=" << m_data.size() << std::endl;

	// should do something useful
	if (bToBeOutput()) {
		std::ostream& out = OH.Loadable();

		for (PD::const_iterator i = m_data.begin(); i != m_data.end(); ++i) {
			out << GetLabel() << "#" << i->pEl->GetLabel() << "#" << i->counter
				<< " " << i->F
				<< " " << i->M
				<< " " << i->dW
				<< " " << i->X
				<< std::endl;
		}
	}
}

void
ModuleIndVel::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler& 
ModuleIndVel::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	// should do something useful
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleIndVel::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	std::cerr << "ModuleIndVel(" << GetLabel() << ")::AssRes: iFirstAssembly=" << iFirstAssembly << std::endl;

	if (iFirstAssembly) {
		iFirstAssembly--;
	}

	// re-initlize iterator to loop over point data
	// when AddSectionalForces is called
	m_data_iter = m_data.begin();

	// should do something useful
	WorkVec.ResizeReset(0);

	return WorkVec;
}

unsigned int
ModuleIndVel::iGetNumPrivData(void) const
{
	return 0;
}

int
ModuleIndVel::iGetNumConnectedNodes(void) const
{
	return 1;
}

void
ModuleIndVel::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

void
ModuleIndVel::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleIndVel::Restart(std::ostream& out) const
{
	return out << "# ModuleIndVel: not implemented" << std::endl;
}

unsigned int
ModuleIndVel::iGetInitialNumDof(void) const
{
	return 0;
}

void 
ModuleIndVel::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleIndVel::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
ModuleIndVel::InitialAssRes(
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
	UserDefinedElemRead *rf = new UDERead<ModuleIndVel>;

	if (!SetUDE("indvel", rf)) {
		delete rf;

		silent_cerr("module-indvel: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

