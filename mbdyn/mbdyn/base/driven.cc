/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "driven.h"
#include "joint.h"

DrivenElem::DrivenElem(DataManager *pdm,
		const DriveCaller* pDC, const Elem* pE,
		SimulationEntity::Hints *ph)
: Elem(pE->GetLabel(), pE->fToBeOutput()),
NestedElem(pE),
DriveOwner(pDC),
pDM(pdm),
pHints(ph),
bActive(false)
{
	ASSERT(pDC != 0);

	bActive = (pDC->dGet() != 0.);
}


DrivenElem::~DrivenElem(void)
{
	ASSERT(pElem != 0);

	if (pHints != 0) {
		for (unsigned i = 0; i < pHints->size(); i++) {
			delete (*pHints)[i];
		}
		delete pHints;
		pHints = 0;
	}
}

bool
DrivenElem::bIsActive(void) const
{
	return (dGet() != 0.);
}

void
DrivenElem::Output(OutputHandler& OH) const
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		pElem->Output(OH);
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
DrivenElem::Restart(std::ostream& out) const
{
	ASSERT(pElem != 0);

	out << "driven: " << GetLabel() << ", ",
		pGetDriveCaller()->Restart(out) << ", ";

	if (pHints) {
		for (unsigned i = 0; i < pHints->size(); i++) {
			out << "hint, \"";

			// TODO

			out << "\", ";
		}
	}

	pElem->Restart(out);

	return out;
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
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
     		pElem->BeforePredict(X, XP, XPrev, XPPrev);
	}
}

void
DrivenElem::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		if (!bActive) {
			bActive = true;
			pElem->SetValue(pDM, X, XP, pHints);
		}
		pElem->AfterPredict(X, XP);

	} else {
		if (bActive) {
			bActive = false;
		}
	}
}

void
DrivenElem::SetValue(DataManager *pdm,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	ASSERT(pElem != 0);

	if (dGet() != 0.) {
		pElem->SetValue(pdm, X, XP, ph);
	}
}

#if 0
void
DrivenElem::SetInitialValue(VectorHandler& X)
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		ElemWithDofs*	pEwD = dynamic_cast<ElemWithDofs *>(pElem);
		if (pEwD) {
			pEwD->SetInitialValue(X);
		}
	}
}
#endif

/* Aggiorna dati in base alla soluzione */
void
DrivenElem::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		pElem->Update(XCurr, XPrimeCurr);
	}
}

/* Inverse Dynamics: */
void
DrivenElem::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		pElem->Update(XCurr, iOrder);
	}
}

/* inverse dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
DrivenElem::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	ASSERT(pElem != 0);

	// must be a joint
	ASSERT(dynamic_cast<const Joint *>(pElem) != 0);
	// must either be torque, prescribed motion or ergonomy
	ASSERT(dynamic_cast<const Joint *>(pElem)->bIsTorque()
		|| dynamic_cast<const Joint *>(pElem)->bIsPrescribedMotion()
		|| dynamic_cast<const Joint *>(pElem)->bIsErgonomy());
	// dGet() must not be zero, otherwise AssJac() would not be called
	ASSERT(dGet() != 0.);

	return pElem->AssJac(WorkMat, XCurr);
}

/* inverse dynamics residual assembly */
SubVectorHandler&
DrivenElem::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr,
	const VectorHandler& XPrimePrimeCurr,
	InverseDynamics::Order iOrder)
{
	ASSERT(pElem != 0);

	if (dGet() != 0.) {
		return pElem->AssRes(WorkVec, XCurr, XPrimeCurr, XPrimePrimeCurr, iOrder);
	}

	WorkVec.Resize(0);

	return WorkVec;
}

/* Inverse Dynamics: */
void
DrivenElem::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP, const VectorHandler& XPP)
{
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		pElem->AfterConvergence(X, XP, XPP);
	}
}

void
DrivenElem::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	ASSERT(pElem != 0);
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
	ASSERT(pElem != 0);

	if (dGet() != 0.) {
		return pElem->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr);
	}

	unsigned int iNumDofs = pElem->iGetNumDof();
	if (iNumDofs == 0) {
		WorkMat.SetNullMatrix();

	} else {
		SparseSubMatrixHandler& WM = WorkMat.SetSparse();
		WM.ResizeReset(iNumDofs, 0);

		/* NOTE: must not fail, since iNumDofs != 0 */
		integer iFirstIndex = dynamic_cast<ElemWithDofs *>(pElem)->iGetFirstIndex();

  		for (unsigned int iCnt = 1; iCnt <= iNumDofs; iCnt++) {
  			doublereal J;

  			switch (pElem->GetDofType(iCnt - 1)) {
  			case DofOrder::ALGEBRAIC:
  				J = 1.;
  				break;
  			case DofOrder::DIFFERENTIAL:
  				J = dCoef;
  				break;
  			default:
  				ASSERT(0);
  				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
  			}

			WM.PutItem(iCnt, iFirstIndex+iCnt,
	   				iFirstIndex+iCnt, J);
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
	ASSERT(pElem != 0);

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

		/* NOTE: must not fail, since iNumDofs != 0 */
		integer iFirstIndex = dynamic_cast<ElemWithDofs *>(pElem)->iGetFirstIndex();

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
	ASSERT(pElem != 0);
	if (dGet() != 0.) {
		return pElem->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr);
	}

	unsigned int iNumDofs = pElem->iGetNumDof();
	if (iNumDofs == 0) {
		WorkVec.Resize(0);
	} else {
		WorkVec.ResizeReset(iNumDofs);

		/* NOTE: must not fail, since iNumDofs != 0 */
		integer iFirstIndex = dynamic_cast<ElemWithDofs *>(pElem)->iGetFirstIndex();

		for (unsigned int iCnt = 1; iCnt <= iNumDofs; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iFirstIndex+iCnt);
			WorkVec.PutCoef(iCnt, -XCurr(iFirstIndex+iCnt));
		}
	}

	return WorkVec;
}

/*
 * Returns the current value of a private data
 * with 0 < i <= iGetNumPrivData()
 */
doublereal
DrivenElem::dGetPrivData(unsigned int i) const
{
	if (dGet() != 0.) {
		return pElem->dGetPrivData(i);
	}

	/* safe default */
	return 0.;
}

/* InitialAssemblyElem */
unsigned int
DrivenElem::iGetInitialNumDof(void) const
{
	if (dGet() != 0.) {
		return NestedElem::iGetInitialNumDof();
	}

	return 0;
}

void
DrivenElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	if (dGet() != 0.) {
		NestedElem::InitialWorkSpaceDim(piNumRows, piNumCols);
	}
}

VariableSubMatrixHandler&
DrivenElem::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	if (dGet() != 0.) {
		return NestedElem::InitialAssJac(WorkMat, XCurr);
	}

	WorkMat.SetNullMatrix();
	return WorkMat;
}

SubVectorHandler&
DrivenElem::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	if (dGet() != 0.) {
		return NestedElem::InitialAssRes(WorkVec, XCurr);
	}

	WorkVec.Resize(0);
	return WorkVec;
}

/* ElemGravityOwner */
Vec3
DrivenElem::GetS_int(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetS_int();
	}

	return Zero3;
}

Mat3x3
DrivenElem::GetJ_int(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetJ_int();
	}

	return Zero3x3;
}

Vec3
DrivenElem::GetB_int(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetB_int();
	}

	return Zero3;
}

// NOTE: gravity owners must provide the momenta moment
// with respect to the origin of the global reference frame!
Vec3
DrivenElem::GetG_int(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetG_int();
	}

	return Zero3;
}

doublereal
DrivenElem::dGetM(void) const
{
	if (dGet() != 0.) {
		return NestedElem::dGetM();
	}

	return 0.;
}

Vec3
DrivenElem::GetS(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetS();
	}

	return Zero3;
}

Mat3x3
DrivenElem::GetJ(void) const
{
	if (dGet() != 0.) {
		return NestedElem::GetJ();
	}

	return Zero3x3;
}

/* ElemDofOwner */
void
DrivenElem::SetInitialValue(VectorHandler& X)
{
	if (dGet() != 0.) {
		return NestedElem::SetInitialValue(X);
	}
}

