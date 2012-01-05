/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Nested elements
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "nestedelem.h"

NestedElem::NestedElem(const Elem* pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
InitialAssemblyElem(pE->GetLabel(), pE->fToBeOutput()),
AerodynamicElem(pE->GetLabel(), dynamic_cast<const ElemWithDofs *>(pE) ? dynamic_cast<const ElemWithDofs *>(pE)->pGetDofOwner() : 0, pE->fToBeOutput()),
ElemGravityOwner(pE->GetLabel(), pE->fToBeOutput()),
pElem(const_cast<Elem *>(pE))
{
	ASSERT(pE != NULL);
}


NestedElem::~NestedElem(void)
{
	ASSERT(pElem != 0);

  	if (pElem != 0) {
		SAFEDELETE(pElem);
	}
}

Elem *
NestedElem::pGetElem(void) const
{
	NestedElem *pNE = dynamic_cast<NestedElem *>(pElem);
	if (pNE) {
		return pNE->pGetElem();
	}

	return pElem;
}

void
NestedElem::OutputPrepare(OutputHandler& OH)
{
	ASSERT(pElem != NULL);
	pElem->OutputPrepare(OH);
}

void
NestedElem::Output(OutputHandler& OH) const
{
	ASSERT(pElem != NULL);
	pElem->Output(OH);
}

void
NestedElem::SetOutputFlag(flag f)
{
	ASSERT(pElem != NULL);
	pElem->SetOutputFlag(f);
}

Elem::Type
NestedElem::GetElemType(void) const
{
	ASSERT(pElem != NULL);
	return pElem->GetElemType(); 
}

/* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
unsigned int
NestedElem::iGetNumDof(void) const
{ 
	ASSERT(pElem != NULL);
	return pElem->iGetNumDof();
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
NestedElem::GetDofType(unsigned int i) const
{
	ASSERT(pElem != NULL);
	return pElem->GetDofType(i);
}

DofOrder::Order
NestedElem::GetEqType(unsigned int i) const
{
	ASSERT(pElem != NULL);
	return pElem->GetEqType(i);
}

std::ostream&
NestedElem::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(pElem != NULL);
	return pElem->DescribeDof(out, prefix, bInitial);
}

void
NestedElem::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	ASSERT(pElem != NULL);
	return pElem->DescribeDof(desc, bInitial, i);
}

std::ostream&
NestedElem::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(pElem != NULL);
	return pElem->DescribeEq(out, prefix, bInitial);
}

void
NestedElem::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	ASSERT(pElem != NULL);
	return pElem->DescribeEq(desc, bInitial, i);
}


/* Dimensioni del workspace */
void
NestedElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	ASSERT(pElem != NULL);
	pElem->WorkSpaceDim(piNumRows, piNumCols);
}

/*
 * Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator
 */
void
NestedElem::BeforePredict(VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev) const
{
	ASSERT(pElem != NULL);
     	pElem->BeforePredict(X, XP, XPrev, XPPrev);
}

void
NestedElem::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	ASSERT(pElem != NULL);
	pElem->AfterPredict(X, XP);
}

void
NestedElem::SetValue(DataManager *pdm,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	ASSERT(pElem != NULL);
	pElem->SetValue(pdm, X, XP, ph);
}

/* Aggiorna dati in base alla soluzione */
void
NestedElem::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);
	pElem->Update(XCurr, XPrimeCurr);
}

  
void
NestedElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	ASSERT(pElem != NULL);
	pElem->AfterConvergence(X, XP);
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
NestedElem::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);
	return pElem->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
}
 
void
NestedElem::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);

	pElem->AssMats(WorkMatA, WorkMatB, XCurr, XPrimeCurr);
}
 
/* assemblaggio residuo */
SubVectorHandler&
NestedElem::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);
	return pElem->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
}

/*
 * Metodi per l'estrazione di dati "privati".
 * Si suppone che l'estrattore li sappia interpretare.
 * Come default non ci sono dati privati estraibili
 */
unsigned int
NestedElem::iGetNumPrivData(void) const
{
	return pElem->iGetNumPrivData();
}

/*
 * Maps a string (possibly with substrings) to a private data;
 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
 * in case of unrecognized data; error must be handled by caller
 */
unsigned int
NestedElem::iGetPrivDataIdx(const char *s) const
{
	return pElem->iGetPrivDataIdx(s);
}

/*
 * Returns the current value of a private data
 * with 0 < i <= iGetNumPrivData()
 */
doublereal
NestedElem::dGetPrivData(unsigned int i) const
{
	return pElem->dGetPrivData(i);
}

/* *******PER IL SOLUTORE PARALLELO******** *
 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs
 */
int
NestedElem::GetNumConnectedNodes(void) const
{
	ASSERT(pElem != NULL);
	return pElem->GetNumConnectedNodes();
}

void
NestedElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
	ASSERT(pElem != NULL);
	pElem->GetConnectedNodes(connectedNodes);
}

/* inverse dynamics capable element */
bool
NestedElem::bInverseDynamics(void) const
{
	ASSERT(pElem != NULL);
	return pElem->bInverseDynamics();
}

/* Inverse Dynamics: */
void
NestedElem::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	ASSERT(pElem != NULL);
	pElem->Update(XCurr, iOrder);
}
 
/* Inverse Dynamics: */
void
NestedElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP, const VectorHandler& XPP)
{
	ASSERT(pElem != NULL);
	pElem->AfterConvergence(X, XP, XPP);
}

/* inverse dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
NestedElem::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	ASSERT(pElem != NULL);
	return pElem->AssJac(WorkMat, XCurr);
}

/* inverse dynamics residual assembly */
SubVectorHandler&
NestedElem::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr,
	const VectorHandler& XPrimePrimeCurr,
	InverseDynamics::Order iOrder)
{
	ASSERT(pElem != NULL);
	return pElem->AssRes(WorkVec, XCurr, XPrimeCurr, XPrimePrimeCurr, iOrder);
}

/* InitialAssemblyElem */
unsigned int
NestedElem::iGetInitialNumDof(void) const
{
	ASSERT(pElem != NULL);
	InitialAssemblyElem *pIAE = dynamic_cast<InitialAssemblyElem*>(pElem);
	if (pIAE) {
		return pIAE->iGetInitialNumDof();
	}
	return 0;
}

void
NestedElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	ASSERT(pElem != NULL);
	InitialAssemblyElem *pIAE = dynamic_cast<InitialAssemblyElem*>(pElem);
	if (pIAE) {
		pIAE->InitialWorkSpaceDim(piNumRows, piNumCols);
	}
}

VariableSubMatrixHandler& 
NestedElem::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	ASSERT(pElem != NULL);
	InitialAssemblyElem *pIAE = dynamic_cast<InitialAssemblyElem*>(pElem);
	if (pIAE) {
		return pIAE->InitialAssJac(WorkMat, XCurr);
	}

	WorkMat.SetNullMatrix();
	return WorkMat;
}

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
NestedElem::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	ASSERT(pElem != NULL);
	InitialAssemblyElem *pIAE = dynamic_cast<InitialAssemblyElem*>(pElem);
	if (pIAE) {
		return pIAE->InitialAssRes(WorkVec, XCurr);
	}

	WorkVec.Resize(0);
	return WorkVec;
}

/* AerodynamicElem */
AerodynamicElem::Type
NestedElem::GetAerodynamicElemType(void) const
{
	ASSERT(pElem != NULL);
	AerodynamicElem *pAE = dynamic_cast<AerodynamicElem *>(pElem);
	if (pAE) {
		return pAE->GetAerodynamicElemType();
	}

	return AerodynamicElem::UNKNOWN;
}

bool
NestedElem::NeedsAirProperties(void) const
{
	ASSERT(pElem != NULL);
	AerodynamicElem *pAE = dynamic_cast<AerodynamicElem *>(pElem);
	if (pAE) {
		return pAE->NeedsAirProperties();
	}

	return false;
}

void
NestedElem::PutAirProperties(const AirProperties* pAP)
{
	ASSERT(pElem != NULL);
	ASSERT(NeedsAirProperties());
	AerodynamicElem *pAE = dynamic_cast<AerodynamicElem *>(pElem);
	if (pAE) {
		return pAE->PutAirProperties(pAP);
	}
}

const InducedVelocity *
NestedElem::pGetInducedVelocity(void) const
{
	ASSERT(pElem != NULL);
	AerodynamicElem *pAE = dynamic_cast<AerodynamicElem *>(pElem);
	if (pAE) {
		return pAE->pGetInducedVelocity();
	}

	return 0;
}

/* ElemGravityOwner */
Vec3
NestedElem::GetS_int(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetS_int();
	}

	return Zero3;
}

Mat3x3
NestedElem::GetJ_int(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetJ_int();
	}

	return Zero3x3;
}

Vec3
NestedElem::GetB_int(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetB_int();
	}

	return Zero3;
}

// NOTE: gravity owners must provide the momenta moment
// with respect to the origin of the global reference frame!
Vec3
NestedElem::GetG_int(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetG_int();
	}

	return Zero3;
}

doublereal
NestedElem::dGetM(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->dGetM();
	}

	return 0.;
}

Vec3
NestedElem::GetS(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetS();
	}

	return Zero3;
}

Mat3x3
NestedElem::GetJ(void) const
{
	ASSERT(pElem != NULL);
	ElemGravityOwner *pEGO = dynamic_cast<ElemGravityOwner *>(pElem);
	if (pEGO) {
		return pEGO->GetJ();
	}

	return Zero3x3;
}

/* ElemWithDofs */
const DofOwner*
NestedElem::pGetDofOwner(void) const
{
	ASSERT(pElem != NULL);
	ElemWithDofs *pEwD = dynamic_cast<ElemWithDofs *>(pElem);
	if (pEwD) {
		return pEwD->pGetDofOwner();
	}

	return 0;

}

integer
NestedElem::iGetFirstIndex(void) const
{
	ASSERT(pElem != NULL);
	ElemWithDofs *pEwD = dynamic_cast<ElemWithDofs *>(pElem);
	if (pEwD) {
		return pEwD->iGetFirstIndex();
	}

	return -1;
}

void
NestedElem::SetInitialValue(VectorHandler& X)
{
	ASSERT(pElem != NULL);
	ElemWithDofs *pEwD = dynamic_cast<ElemWithDofs *>(pElem);
	if (pEwD) {
		pEwD->SetInitialValue(X);
	}
}

