/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "dataman.h"
#include "elem.h"
#include "gravity.h"
#include "aerodyn.h"

/* Elem - begin */

Elem::Elem(unsigned int uL, flag fOut)
: WithLabel(uL), ToBeOutput(fOut),
m_uInverseDynamicsFlags(InverseDynamics::ERGONOMY|InverseDynamics::RIGHT_HAND_SIDE)
{
	NO_OP;
}

Elem::~Elem(void) 
{
	NO_OP;
}

/* assemblaggio matrici per autovalori */
void
Elem::AssMats(VariableSubMatrixHandler& /* WorkMatA */ ,
	VariableSubMatrixHandler& /* WorkMatB */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ ) 
{
	silent_cerr(psElemNames[GetElemType()] << "(" << GetLabel() << "): "
		"AssMats() not implemented yet" << std::endl);
}

bool
Elem::bInverseDynamics(void) const
{
	return false;
}

void
Elem::SetInverseDynamicsFlags(unsigned uIDF)
{
	m_uInverseDynamicsFlags = uIDF;
}

unsigned
Elem::GetInverseDynamicsFlags(void) const
{
	return m_uInverseDynamicsFlags;
}

bool
Elem::bIsErgonomy(void) const
{
	return (m_uInverseDynamicsFlags & InverseDynamics::ERGONOMY);
}

bool
Elem::bIsRightHandSide(void) const
{
	return (m_uInverseDynamicsFlags & InverseDynamics::RIGHT_HAND_SIDE);
}

/* inverse dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
Elem::AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) 
{ 
	silent_cerr(psElemNames[GetElemType()] << "(" << GetLabel() << "): "
		"Elem::AssJac() for inverse dynamics not implemented yet" << std::endl);
	return WorkMat;
};

/* inverse dynamics residual assembly */
SubVectorHandler&
Elem::AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder) 
{ 
	silent_cerr(psElemNames[GetElemType()] << "(" << GetLabel() << "): "
		"Elem::AssRes(" << invdyn2str(iOrder) << ") for inverse dynamics not implemented yet" << std::endl);
	return WorkVec;
};
	
unsigned int 
Elem::iGetNumDof(void) const
{
	return 0;
}

std::ostream&
Elem::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

void
Elem::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i <= 0);
	desc.resize(0);
}

std::ostream&
Elem::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out;
}

void
Elem::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i <= 0);
	desc.resize(0);
}

DofOrder::Order 
Elem::GetDofType(unsigned int) const
{
	silent_cerr(psElemNames[GetElemType()] << "(" << GetLabel() << "): "
		"GetDofType() is undefined because element "
		"has no degrees of freedom" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

Elem::Type
str2elemtype(const char *const s)
{
	for (int i = 0; i < Elem::LASTELEMTYPE; i++) {
		if (strcasecmp(s, psReadElemsElems[i]) == 0) {
			return Elem::Type(i);
		}
	}

	return Elem::UNKNOWN;
}

/* Elem - end */

/* database of registered element types */
typedef std::map<std::string, ElemRead *, ltstrcase> ElemFuncMapType;
static ElemFuncMapType ElemFuncMap;

/* element parsing checkers */
struct ElemWordSetType : public HighParser::WordSet {
	/*
	 * returns true if the string "s" is recognized as an element type
	 * hijacks the database of registered element types for consistency
	 */
	bool IsWord(const std::string& s) const {
		return ElemFuncMap.find(std::string(s)) != ElemFuncMap.end();
	};
};

static ElemWordSetType ElemWordSet;

/* element type registration functions: call to register one */
bool
SetElem(const char *name, ElemRead *rf)
{
	pedantic_cout("registering element \"" << name << "\""
		<< std::endl );
	return ElemFuncMap.insert(ElemFuncMapType::value_type(name, rf)).second;
}

/* ElemWithDofs - begin */

ElemWithDofs::ElemWithDofs(unsigned int uL,
	const DofOwner* pDO, flag fOut)
: Elem(uL, fOut), DofOwnerOwner(const_cast<DofOwner *>(pDO))
{
	NO_OP;
}

ElemWithDofs::~ElemWithDofs(void)
{
	NO_OP;
}

/* ElemWithDofs - end */


/* SubjectToInitialAssembly - begin */

SubjectToInitialAssembly::SubjectToInitialAssembly(void) 
{
	NO_OP;
}


SubjectToInitialAssembly::~SubjectToInitialAssembly(void)
{
	NO_OP;
}

/* SubjectToInitialAssembly - end */


/* InitialAssemblyElem - begin */

InitialAssemblyElem::InitialAssemblyElem(unsigned int uL, flag fOut)
: Elem(uL, fOut)
{ 
	NO_OP;
}

InitialAssemblyElem::~InitialAssemblyElem(void)
{
	NO_OP;
}

/* InitialAssemblyElem - end */
