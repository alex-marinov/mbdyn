/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* Driven elements:
 * elements that are used depending on the (boolean) value
 * of a driver. Example: a driven joint is assembled only 
 * if the driver is true, otherwise there is no joint and
 * the reaction unknowns are set to zero
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "driven.h"

DrivenElem::DrivenElem(const DriveCaller* pDC, const Elem* pE)
: Elem(pE->GetLabel(), Elem::DRIVEN, pE->fToBeOutput()),
DriveOwner(pDC),
pElem((Elem*)pE)
{
	ASSERT(pDC != NULL);
	ASSERT(pE != NULL);
}


DrivenElem::~DrivenElem(void)
{
	ASSERT(pElem != NULL);
  	if (pElem != NULL) {
		SAFEDELETE(pElem);
	}
}

void
DrivenElem::Output(OutputHandler& OH) const
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
		pElem->Output(OH);
	}	 
}

void
DrivenElem::SetOutputFlag(flag f)
{
	ASSERT(pElem != NULL);
	pElem->SetOutputFlag(f);
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
DrivenElem::Restart(std::ostream& out) const
{
	ASSERT(pElem != NULL);
	out << "driven: " << GetLabel() << ", ",
		pGetDriveCaller()->Restart(out) << ", ",
		pElem->Restart(out);
	return out;
}

/* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
unsigned int
DrivenElem::iGetNumDof(void) const
{ 
	ASSERT(pElem != NULL);
	return pElem->iGetNumDof();
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
DrivenElem::GetDofType(unsigned int i) const
{
	ASSERT(pElem != NULL);
	return pElem->GetDofType(i);
}

/* Dimensioni del workspace */
void
DrivenElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	ASSERT(pElem != NULL);
	pElem->WorkSpaceDim(piNumRows, piNumCols);
}

/*
 * Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator
 */
void
DrivenElem::BeforePredict(VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev) const
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
     		pElem->BeforePredict(X, XP, XPrev, XPPrev);
	}
}

void
DrivenElem::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
		pElem->AfterPredict(X, XP);
	}
}

/* Aggiorna dati in base alla soluzione */
void
DrivenElem::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
		pElem->Update(XCurr, XPrimeCurr);
	}
}
   
void
DrivenElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
		pElem->AfterConvergence(X, XP);
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
DrivenElem::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);

	if (dGet() != 0.) {
		return pElem->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
	}

	unsigned int iNumDofs = pElem->iGetNumDof();
	if (iNumDofs == 0) {
		WorkMat.SetNullMatrix();

	} else {
		SparseSubMatrixHandler& WM = WorkMat.SetSparse();
		WM.ResizeReset(iNumDofs, 0);

		integer iFirstIndex = pElem->pGetElemWithDofs()->iGetFirstIndex();
  		
  		for (unsigned int iCnt = 1; iCnt <= iNumDofs; iCnt++) {
			WM.PutItem(iCnt, iFirstIndex+iCnt,
	   				iFirstIndex+iCnt, 1.);
 		}
	}
	return WorkMat;
}
 
void
DrivenElem::AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);

	if (dGet() != 0.) {
		pElem->AssMats(WorkMatA, WorkMatB, XCurr, XPrimeCurr);
		return;
	}

	WorkMatA.SetNullMatrix();

	unsigned int iNumDofs = pElem->iGetNumDof();
	if (iNumDofs == 0) {
		WorkMatB.SetNullMatrix();
	} else {
		SparseSubMatrixHandler& WM = WorkMatB.SetSparse();
		WM.ResizeReset(iNumDofs, 0);

		integer iFirstIndex = pElem->pGetElemWithDofs()->iGetFirstIndex();
  		
  		for (unsigned int iCnt = 1; iCnt <= iNumDofs; iCnt++) {
			WM.PutItem(iCnt, iFirstIndex+iCnt,
	   				iFirstIndex+iCnt, 1.);
 		}
	}
};
 
/* assemblaggio residuo */
SubVectorHandler&
DrivenElem::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != NULL);
	if (dGet() != 0.) {
		return pElem->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
	}

	unsigned int iNumDofs = pElem->iGetNumDof();
	if (iNumDofs == 0) {
		WorkVec.Resize(0);
	} else {
		WorkVec.ResizeReset(iNumDofs);

		integer iFirstIndex = pElem->pGetElemWithDofs()->iGetFirstIndex();

		for (unsigned int iCnt = 1; iCnt <= iNumDofs; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex+iCnt);
			WorkVec.PutCoef(iCnt, -XCurr.dGetCoef(iFirstIndex+iCnt));
		}
	}
	return WorkVec;
}

/* *******PER IL SOLUTORE PARALLELO******** *
 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs
 */
void
DrivenElem::GetConnectedNodes(int& NumNodes,
		Node::Type* NdTyps,
		unsigned int* NdLabels)
{
	ASSERT(pElem != NULL);
	return pElem->GetConnectedNodes(NumNodes, NdTyps, NdLabels);
}

