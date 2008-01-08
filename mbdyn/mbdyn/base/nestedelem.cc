/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "nestedelem.h"

NestedElem::NestedElem(DataManager *pdm, const Elem* pE)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
pDM(pdm),
pElem(const_cast<Elem *>(pE))
{
	ASSERT(pDC != NULL);
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
NestedElem::SetInitialValue(VectorHandler& X) const
{
	ASSERT(pElem != NULL);
	ElemWithDofs*	pEwD = dynamic_cast<ElemWithDofs *>(pElem);
	if (pEwD) {
		pEwD->SetInitialValue(X);
	}
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

/* Inverse Dynamics: */
void
NestedElem::Update(const VectorHandler& XCurr, int iOrder)
{
	ASSERT(pElem != NULL);
	pElem->Update(XCurr, iOrder);
}
 
   
void
NestedElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	ASSERT(pElem != NULL);
	pElem->AfterConvergence(X, XP);
}

/* Inverse Dynamics: */
void
NestedElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP, const VectorHandler& XPP)
{
	ASSERT(pElem != NULL);
	pElem->AfterConvergence(X, XP, XPP);
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
NestedElem::GetConnectedNodes(std::vector<const Node *>& connectedNodes)
{
	ASSERT(pElem != NULL);
	pElem->GetConnectedNodes(connectedNodes);
}

