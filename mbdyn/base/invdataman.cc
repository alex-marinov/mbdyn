/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/* Inverse Dynamics DataManager */

/*
 * Copyright 2008 Alessandro Fumagalli <alessandro.fumagalli@polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */


#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <strings.h>
#include <time.h>

/* for inertia element */
#include <set>

#include "dataman.h"
#include "friction.h"

#include "solver.h"
#include "invdyn.h"
#include "invsolver.h"
#include "constltp.h"
#include "dataman_.h"
/* add-ons for math parser */
#include "dofpgin.h"
#include "privpgin.h"
#include "dummypgin.h"
#include "modelns.h"

/* temporary? */
#include "beam.h"
#include "inertia.h"

/* To allow direct loading of modules */
#include "modules.h"

/* To handle  of Elem2Param */
#include "j2p.h"

/* deformable joint */
#include "beam2.h"
#include "body.h"
#include "joint.h"
#include "jointreg.h"
#include "Rot.hh"


void 
DataManager::LinkToSolution(VectorHandler& XCurr,
	VectorHandler& XPrimeCurr,
	VectorHandler& XPrimePrimeCurr,
	VectorHandler& LambdaCurr)
{
	pXCurr = &XCurr;
	pXPrimeCurr = &XPrimeCurr;
	pXPrimePrimeCurr = &XPrimePrimeCurr;
	pLambdaCurr = &LambdaCurr;

	DrvHdl.LinkToSolution(XCurr, XPrimeCurr);
}
void
DataManager::AssConstrJac(MatrixHandler& JacHdl)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	ASSERT(pWorkMat != NULL);
	ASSERT(Elems.begin() != Elems.end());

	AssConstrJac(JacHdl, ElemIter, *pWorkMat);
}

void
DataManager::AssConstrJac(MatrixHandler& JacHdl,
		VecIter<Elem *> &Iter,
		VariableSubMatrixHandler& WorkMat)
{
	DEBUGCOUT("Entering DataManager::AssJac()" << std::endl);

	JacHdl.Reset();

	InverseSolver *pIDS = dynamic_cast<InverseSolver *>(pSolver);

	switch (pIDS->GetProblemType()) {
	case InverseDynamics::FULLY_ACTUATED_COLLOCATED:
		for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
			j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
		{
			Joint *pJ = Cast<Joint>(j->second);
			if (pJ->bIsPrescribedMotion()) {
				ASSERT(pJ->bIsTorque());
				JacHdl += j->second->AssJac(WorkMat, *pXCurr);
			}
		}
		break;

	case InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED:
	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
		silent_cerr("DataManager::AssConstrJac(" << pIDS->GetProblemType() << ") not implemented yet" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED: {
		// NOTE: in this case, the name of the function is misleading,
		// since it assembles the entire problem's Jacobian matrix
		// and not only the Jacobian matrix of the constraints
		InverseDynamicsStepSolver *pIDSS = dynamic_cast<InverseDynamicsStepSolver *>(pSolver->pGetStepIntegrator());
		if (pIDSS->GetOrder() == InverseDynamics::INVERSE_DYNAMICS) {
			for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
				j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Joint *pJ = Cast<Joint>(j->second, true);
				if (pJ == 0) {
					bActive = false;
					pJ = Cast<Joint>(j->second, false);
				}

				if (pJ->bIsTorque() && bActive) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
					WorkMat.AddToT(JacHdl);

				} else if (pJ->bIsPrescribedMotion() || (pJ->bIsTorque() && !bActive)) {
					integer iNumDof = pJ->iGetNumDof();
					integer iFirstIndex = pJ->iGetFirstIndex();
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) = 1.;
					}
				}
			}

			// FIXME: regularization could be needed also in other phases
			// may need to be related to the state of the regularized joint
			for (ElemContainerType::iterator j = ElemData[Elem::JOINT_REGULARIZATION].ElemContainer.begin();
				j != ElemData[Elem::JOINT_REGULARIZATION].ElemContainer.end(); ++j)
			{
				JointRegularization *pJ = Cast<JointRegularization>(j->second, true);
				if (pJ) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
				}
			}

#if 0
			// this _should_ be harmless...
			for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
				n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
			{
				ASSERT(n->second->iGetNumDof() == 6);
				integer iFirstIndex = n->second->iGetFirstIndex();

				for (integer iCnt = 1; iCnt <= 6; iCnt++) {
					JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) += 1.;
				}
			}
#endif

		} else {
			doublereal dw1, dw2;
			pIDS->GetWeight(pIDSS->GetOrder(), dw1, dw2);
			if (dw1 > 0.) {
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					ASSERT(n->second->iGetNumDof() == 6);
					integer iFirstIndex = n->second->iGetFirstIndex();

					for (integer iCnt = 1; iCnt <= 6; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) += dw1;
					}
				}
			}

			if (dw2 > 0.) {
// DO NOT ENABLE, BROKEN
//#define TORQUE_OPTIMIZATION
#ifdef TORQUE_OPTIMIZATION
				if (pIDSS->GetOrder() == InverseDynamics::POSITION) {
					doublereal h = DrvHdl.dGetTimeStep();
					dw2 /= h*h;
				}
#endif // TORQUE_OPTIMIZATION

				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					for (int iRow = 1; iRow <= 3; iRow++) {
						JacHdl(iFirstIndex + iRow, iFirstIndex + iRow) += dm;

						for (int iCol = 1; iCol <= 3; iCol++) {
							JacHdl(iFirstIndex + 3 + iRow, iFirstIndex + 3 + iCol) += J(iRow, iCol);
						}
					}

					JacHdl(iFirstIndex + 3 + 1, iFirstIndex + 2) = -S(3);	// 1, 2
					JacHdl(iFirstIndex + 3 + 1, iFirstIndex + 3) = S(2);	// 1, 3
					JacHdl(iFirstIndex + 3 + 2, iFirstIndex + 3) = -S(1);	// 2, 3
					JacHdl(iFirstIndex + 3 + 2, iFirstIndex + 1) = S(3);	// 2, 1
					JacHdl(iFirstIndex + 3 + 3, iFirstIndex + 1) = -S(2);	// 3, 1
					JacHdl(iFirstIndex + 3 + 3, iFirstIndex + 2) = S(1);	// 3, 2

					JacHdl(iFirstIndex + 2, iFirstIndex + 3 + 1) = -S(3);	// 2, 1
					JacHdl(iFirstIndex + 3, iFirstIndex + 3 + 1) = S(2);	// 3, 1
					JacHdl(iFirstIndex + 3, iFirstIndex + 3 + 2) = -S(1);	// 3, 2
					JacHdl(iFirstIndex + 1, iFirstIndex + 3 + 2) = S(3);	// 1, 2
					JacHdl(iFirstIndex + 1, iFirstIndex + 3 + 3) = -S(2);	// 1, 3
					JacHdl(iFirstIndex + 2, iFirstIndex + 3 + 3) = S(1);	// 2, 3
				}
			}

			for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
				j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Joint *pJ = Cast<Joint>(j->second, true);
				if (pJ == 0) {
					bActive = false;
					pJ = Cast<Joint>(j->second, false);
				}

				if (pJ->bIsPrescribedMotion() && bActive) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
					WorkMat.AddToT(JacHdl);

				} else if (pJ->bIsErgonomy() && bActive && pIDSS->GetOrder() == InverseDynamics::POSITION) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);

				} else if (pJ->bIsTorque() || (pJ->bIsPrescribedMotion() && !bActive)) {
					integer iNumDof = pJ->iGetNumDof();
					integer iFirstIndex = pJ->iGetFirstIndex();
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) = 1.;
					}
				}
			}

			for (ElemContainerType::iterator j = ElemData[Elem::BEAM].ElemContainer.begin();
				j != ElemData[Elem::BEAM].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Beam2 *pB = Cast<Beam2>(j->second, true);
				if (pB == 0) {
					bActive = false;
					pB = Cast<Beam2>(j->second, false);
				}

				if (pB && pB->bIsErgonomy() && bActive && pIDSS->GetOrder() == InverseDynamics::POSITION) {
					JacHdl += pB->AssJac(WorkMat, *pXCurr);
				}
			}
		}
		} break;

	case InverseDynamics::UNDERDETERMINED_OVERACTUATED:
	case InverseDynamics::FULLY_DETERMINED_OVERACTUATED:{
		InverseDynamicsStepSolver *pIDSS = dynamic_cast<InverseDynamicsStepSolver *>(pSolver->pGetStepIntegrator());
		if (pIDSS->GetOrder() == InverseDynamics::INVERSE_DYNAMICS) {
			for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
				j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Joint *pJ = Cast<Joint>(j->second, true);
				if (pJ == 0) {
					bActive = false;
					pJ = Cast<Joint>(j->second, false);
				}

				if (pJ->bIsTorque() && bActive) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
					WorkMat.AddToT(JacHdl);

				} else if (pJ->bIsPrescribedMotion() || (pJ->bIsTorque() && !bActive)) {
					integer iNumDof = pJ->iGetNumDof();
					integer iFirstIndex = pJ->iGetFirstIndex();
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) = 1.;
					}
				}
			}

			// FIXME: regularization could be needed also in other phases
			// may need to be related to the state of the regularized joint
			for (ElemContainerType::iterator j = ElemData[Elem::JOINT_REGULARIZATION].ElemContainer.begin();
				j != ElemData[Elem::JOINT_REGULARIZATION].ElemContainer.end(); ++j)
			{
				JointRegularization *pJ = Cast<JointRegularization>(j->second, true);
				if (pJ) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
				}
			}

#if 0
			// this _should_ be harmless...
			for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
				n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
			{
				ASSERT(n->second->iGetNumDof() == 6);
				integer iFirstIndex = n->second->iGetFirstIndex();

				for (integer iCnt = 1; iCnt <= 6; iCnt++) {
					JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) += 1.;
				}
			}
#endif

		} else {
			doublereal dw1, dw2;
			pIDS->GetWeight(pIDSS->GetOrder(), dw1, dw2);
			if (dw1 > 0.) {
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					ASSERT(n->second->iGetNumDof() == 6);
					integer iFirstIndex = n->second->iGetFirstIndex();

					for (integer iCnt = 1; iCnt <= 6; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) += dw1;
					}
				}
			}

			if (dw2 > 0.) {
// DO NOT ENABLE, BROKEN
//#define TORQUE_OPTIMIZATION
#ifdef TORQUE_OPTIMIZATION
				if (pIDSS->GetOrder() == InverseDynamics::POSITION) {
					doublereal h = DrvHdl.dGetTimeStep();
					dw2 /= h*h;
				}
#endif // TORQUE_OPTIMIZATION

				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					for (int iRow = 1; iRow <= 3; iRow++) {
						JacHdl(iFirstIndex + iRow, iFirstIndex + iRow) += dm;

						for (int iCol = 1; iCol <= 3; iCol++) {
							JacHdl(iFirstIndex + 3 + iRow, iFirstIndex + 3 + iCol) += J(iRow, iCol);
						}
					}

					JacHdl(iFirstIndex + 3 + 1, iFirstIndex + 2) = -S(3);	// 1, 2
					JacHdl(iFirstIndex + 3 + 1, iFirstIndex + 3) = S(2);	// 1, 3
					JacHdl(iFirstIndex + 3 + 2, iFirstIndex + 3) = -S(1);	// 2, 3
					JacHdl(iFirstIndex + 3 + 2, iFirstIndex + 1) = S(3);	// 2, 1
					JacHdl(iFirstIndex + 3 + 3, iFirstIndex + 1) = -S(2);	// 3, 1
					JacHdl(iFirstIndex + 3 + 3, iFirstIndex + 2) = S(1);	// 3, 2

					JacHdl(iFirstIndex + 2, iFirstIndex + 3 + 1) = -S(3);	// 2, 1
					JacHdl(iFirstIndex + 3, iFirstIndex + 3 + 1) = S(2);	// 3, 1
					JacHdl(iFirstIndex + 3, iFirstIndex + 3 + 2) = -S(1);	// 3, 2
					JacHdl(iFirstIndex + 1, iFirstIndex + 3 + 2) = S(3);	// 1, 2
					JacHdl(iFirstIndex + 1, iFirstIndex + 3 + 3) = -S(2);	// 1, 3
					JacHdl(iFirstIndex + 2, iFirstIndex + 3 + 3) = S(1);	// 2, 3
				}
			}

			for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
				j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Joint *pJ = Cast<Joint>(j->second, true);
				if (pJ == 0) {
					bActive = false;
					pJ = Cast<Joint>(j->second, false);
				}

				if (pJ->bIsPrescribedMotion() && bActive) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);
					WorkMat.AddToT(JacHdl);

				} else if (pJ->bIsErgonomy() && bActive && pIDSS->GetOrder() == InverseDynamics::POSITION) {
					JacHdl += pJ->AssJac(WorkMat, *pXCurr);

				} else if (pJ->bIsTorque() || (pJ->bIsPrescribedMotion() && !bActive)) {
					integer iNumDof = pJ->iGetNumDof();
					integer iFirstIndex = pJ->iGetFirstIndex();
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						JacHdl(iFirstIndex + iCnt, iFirstIndex + iCnt) = 1.;
					}
				}
			}

			for (ElemContainerType::iterator j = ElemData[Elem::BEAM].ElemContainer.begin();
				j != ElemData[Elem::BEAM].ElemContainer.end(); ++j)
			{
				bool bActive(true);
				Beam2 *pB = Cast<Beam2>(j->second, true);
				if (pB == 0) {
					bActive = false;
					pB = Cast<Beam2>(j->second, false);
				}

				if (pB && pB->bIsErgonomy() && bActive && pIDSS->GetOrder() == InverseDynamics::POSITION) {
					JacHdl += pB->AssJac(WorkMat, *pXCurr);
				}
			}
		}
	}


	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* Constraint residual assembly: */
void
DataManager::AssConstrRes(VectorHandler& ResHdl, InverseDynamics::Order iOrder)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	AssConstrRes(ResHdl, ElemIter, *pWorkVec, iOrder);
}

void
DataManager::AssConstrRes(VectorHandler& ResHdl,
	VecIter<Elem *> &Iter,
	SubVectorHandler& WorkVec,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	// TODO

	bool ChangedEqStructure(false);
	InverseSolver *pIDS = dynamic_cast<InverseSolver *>(pSolver);

	switch (pIDS->GetProblemType()) {
	case InverseDynamics::FULLY_ACTUATED_COLLOCATED:
		for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
			j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
		{
			Joint *pJ = Cast<Joint>(j->second);
			if (pJ->bIsPrescribedMotion()) {
				ASSERT(pJ->bIsTorque());
				try {
					ResHdl += pJ->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, iOrder);
				}
				catch (Elem::ChangedEquationStructure& e) {
					ResHdl += WorkVec;
					ChangedEqStructure = true;
				}
			}
		}
		break;

	case InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED:
		//break;

	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
		silent_cerr("DataManager::AssConstrRes(" << pIDS->GetProblemType() << ") not implemented yet" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED: {
		// NOTE: in this case, the name of the function is misleading,
		// since it assembles the entire problem's residual
		// and not only the residual of the constraints
		doublereal dw1, dw2;
		pIDS->GetWeight(iOrder, dw1, dw2);
		switch (iOrder) {
		case InverseDynamics::POSITION:
			if (dw1 > 0.) {
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DX(pNode->GetXPrev() - pNode->GetXCurr());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*DX(iCnt);
					}

					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DTheta(RotManip::VecRot(pNode->GetRPrev().MulMT(pNode->GetRCurr())));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*DTheta(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
#ifdef TORQUE_OPTIMIZATION
				doublereal h = DrvHdl.dGetTimeStep();
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsRightHandSide()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DXP((pNode->GetXCurr() - pNode->GetXPrev())/h);
					Vec3 DXPP((DXP - pNode->GetVPrev())/h);
					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DW(RotManip::VecRot(pNode->GetRCurr().MulMT(pNode->GetRPrev()))/h);
					Vec3 DWP((DW - pNode->GetWPrev())/h);

					Vec3 XRes(DXPP*dm + DWP.Cross(S) + DW.Cross(DW.Cross(S)));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) -= XRes(iCnt);
					}

					Vec3 ThetaRes(S.Cross(DXPP) + J*DWP + DW.Cross(J*DW));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) -= ThetaRes(iCnt);
					}
				}
#else // !TORQUE_OPTIMIZATION
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DX(pNode->GetXPrev() - pNode->GetXCurr());
					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DTheta(RotManip::VecRot(pNode->GetRPrev().MulMT(pNode->GetRCurr())));

					Vec3 XRes(DX*dm - S.Cross(DTheta));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += XRes(iCnt);
					}

					Vec3 ThetaRes(S.Cross(DX) + J*DTheta);
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += ThetaRes(iCnt);
					}
				}
#endif // !TORQUE_OPTIMIZATION
			}
			break;

		case InverseDynamics::VELOCITY:
			if (dw1 > 0.) {
				doublereal h = DrvHdl.dGetTimeStep();

				// xp_k = xp_km1/3 + 2/3*(x_k - x_km1)/h
				// xp_k = 2*(x_k - x_km1)/h - xp_km1
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					const Vec3& VPrev(pNode->GetVPrev());
						(void) VPrev; // silence set but not used warning
					const Vec3& XPrev(pNode->GetXPrev());
					const Vec3& XCurr(pNode->GetXCurr());

// #define USE_2XmV 1	// 2nd order, a-stable, oscillations
// #define USE_2XpVd3 1	// 1st order, a/l-stable, more accurate
#define USE_X 1		// 1st order, l-stable (implicit Euler), less accurate (alpha == 1.)
// #define USE_alphaX_betaV	(1.0)	// alpha = 1. + |rho|; beta = (1 - alpha) for 1st order accuracy

#if USE_2XmV
					Vec3 VRef((XCurr - XPrev)*(2./h) - VPrev);
#elif USE_2XpVd3
					Vec3 VRef((XCurr - XPrev)*(2./3./h) + VPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 VRef((XCurr - XPrev)*(USE_alphaX_betaV/h) + VPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 VRef((XCurr - XPrev)/h);
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*VRef(iCnt);
					}

					const Vec3& WPrev(pNode->GetWPrev());
						(void) WPrev; // silence set but not used warning
					const Mat3x3& RPrev(pNode->GetRPrev());
					const Mat3x3& RCurr(pNode->GetRCurr());

#if USE_2XmV
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(2./h) - WPrev);
#elif USE_2XpVd3
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(2./3./h) + WPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(USE_alphaX_betaV/h) + WPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))/h);
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*WRef(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 BPrev(pNode->GetVPrev()*dm - S.Cross(pNode->GetWPrev()));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += BPrev(iCnt);
					}

					Vec3 GPrev(S.Cross(pNode->GetVPrev()) + J*pNode->GetWPrev());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += GPrev(iCnt);
					}
				}
			}
			break;

		case InverseDynamics::ACCELERATION:
			if (dw1 > 0.) {
				doublereal h = DrvHdl.dGetTimeStep();

				// xpp_k = xpp_km1/3 + 2/3*(xp_k - xp_km1)/h
				// xpp_k = 2*(xp_k - xp_km1)/h - xpp_km1
				// xpp_k = xpp_km1 + 6*(xp_k + xp_km1)/h - 12*(x_k - x_km1)/h^2
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					const Vec3& XPPPrev(pNode->GetXPPPrev());
						(void) XPPPrev; // silence set but not used warning
					const Vec3& VPrev(pNode->GetVPrev());
					const Vec3& VCurr(pNode->GetVCurr());

#if USE_2XmV
					Vec3 XPPRef((VCurr - VPrev)*(2./h) - XPPPrev);
#elif USE_2XpVd3
					Vec3 XPPRef((VCurr - VPrev)*(2./3./h) + XPPPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 XPPRef((VCurr - VPrev)*(USE_alphaX_betaV/h) + XPPPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 XPPRef((VCurr - VPrev)/h);
#endif
#if 0
					const Vec3& XPrev(pNode->GetXPrev());
					const Vec3& XCurr(pNode->GetXCurr());

					Vec3 XPPRef(XPPPrev + (VCurr + VPrev)*(6./h) - (XCurr - XPrev)*(12./h/h));
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*XPPRef(iCnt);
					}

					const Vec3& WPPrev(pNode->GetWPPrev());
						(void) WPPrev; // silence set but not used warning
					const Vec3& WPrev(pNode->GetWPrev());
					const Vec3& WCurr(pNode->GetWCurr());

#if USE_2XmV
					Vec3 WPRef((WCurr - WPrev)*(2./h) - WPPrev);
#elif USE_2XpVd3
					Vec3 WPRef((WCurr - WPrev)*(2./3./h) + WPPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 WPRef((WCurr - WPrev)*(USE_alphaX_betaV/h) + WPPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 WPRef((WCurr - WPrev)/h);
#endif
#if 0
					const Mat3x3& RPrev(pNode->GetRPrev());
					const Mat3x3& RCurr(pNode->GetRCurr());

					Vec3 WPRef(WPPrev + (WCurr + WPrev)*(6./h) - RotManip::VecRot(RCurr.MulMT(RPrev))*(12./h/h));
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*WPRef(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 BPPrev(pNode->GetXPPPrev()*dm - S.Cross(pNode->GetWPPrev()));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += BPPrev(iCnt);
					}

					Vec3 GPPrev(S.Cross(pNode->GetXPPPrev()) + J*pNode->GetWPPrev());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += GPPrev(iCnt);
					}
				}
			}
			break;

		default:
			ASSERT(0);
		}

		for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
			j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
		{
			bool bActive(true);
			Joint *pJ = Cast<Joint>(j->second, true);
			if (pJ == 0) {
				bActive = false;
				pJ = Cast<Joint>(j->second, false);
			}

			if ((pJ->bIsPrescribedMotion()
				|| (iOrder == InverseDynamics::POSITION && pJ->bIsErgonomy())) && bActive)
			{
				try {
					ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, iOrder);
				}
				catch (Elem::ChangedEquationStructure& e) {
					ResHdl += WorkVec;
					ChangedEqStructure = true;
				}

			} else if (pJ->bIsTorque() || (pJ->bIsPrescribedMotion() && !bActive)) {
				integer iNumDof = pJ->iGetNumDof();
				integer iFirstIndex = pJ->iGetFirstIndex();
				if (iOrder == InverseDynamics::POSITION) {
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						ResHdl(iFirstIndex + iCnt) = -(*pXCurr)(iFirstIndex + iCnt);
					}

				} else {
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						ResHdl(iFirstIndex + iCnt) = 0.;
					}
				}
			}
		}

		for (ElemContainerType::iterator j = ElemData[Elem::BEAM].ElemContainer.begin();
			j != ElemData[Elem::BEAM].ElemContainer.end(); ++j)
		{
			bool bActive(true);
			Beam2 *pB = Cast<Beam2>(j->second, true);
			if (pB == 0) {
				bActive = false;
				pB = Cast<Beam2>(j->second, false);
			}

			if (pB && (iOrder == InverseDynamics::POSITION && pB->bIsErgonomy()) && bActive)
			{
				try {
					ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, iOrder);
				}
				catch (Elem::ChangedEquationStructure& e) {
					ResHdl += WorkVec;
					ChangedEqStructure = true;
				}
			}
		}

		} break;
	case InverseDynamics::UNDERDETERMINED_OVERACTUATED:
	case InverseDynamics::FULLY_DETERMINED_OVERACTUATED:{
		doublereal dw1, dw2;
		pIDS->GetWeight(iOrder, dw1, dw2);
		switch (iOrder) {
		case InverseDynamics::POSITION:
			if (dw1 > 0.) {
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DX(pNode->GetXPrev() - pNode->GetXCurr());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*DX(iCnt);
					}

					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DTheta(RotManip::VecRot(pNode->GetRPrev().MulMT(pNode->GetRCurr())));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*DTheta(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
#ifdef TORQUE_OPTIMIZATION
				doublereal h = DrvHdl.dGetTimeStep();
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsRightHandSide()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DXP((pNode->GetXCurr() - pNode->GetXPrev())/h);
					Vec3 DXPP((DXP - pNode->GetVPrev())/h);
					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DW(RotManip::VecRot(pNode->GetRCurr().MulMT(pNode->GetRPrev()))/h);
					Vec3 DWP((DW - pNode->GetWPrev())/h);

					Vec3 XRes(DXPP*dm + DWP.Cross(S) + DW.Cross(DW.Cross(S)));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) -= XRes(iCnt);
					}

					Vec3 ThetaRes(S.Cross(DXPP) + J*DWP + DW.Cross(J*DW));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) -= ThetaRes(iCnt);
					}
				}
#else // !TORQUE_OPTIMIZATION
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 DX(pNode->GetXPrev() - pNode->GetXCurr());
					// VecRot(Rp*Rc^T) = -VecRot(Rc*Rp^T)
					Vec3 DTheta(RotManip::VecRot(pNode->GetRPrev().MulMT(pNode->GetRCurr())));

					Vec3 XRes(DX*dm - S.Cross(DTheta));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += XRes(iCnt);
					}

					Vec3 ThetaRes(S.Cross(DX) + J*DTheta);
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += ThetaRes(iCnt);
					}
				}
#endif // !TORQUE_OPTIMIZATION
			}
			break;

		case InverseDynamics::VELOCITY:
			if (dw1 > 0.) {
				doublereal h = DrvHdl.dGetTimeStep();

				// xp_k = xp_km1/3 + 2/3*(x_k - x_km1)/h
				// xp_k = 2*(x_k - x_km1)/h - xp_km1
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					const Vec3& VPrev(pNode->GetVPrev());
						(void) VPrev; // silence set but not used warning
					const Vec3& XPrev(pNode->GetXPrev());
					const Vec3& XCurr(pNode->GetXCurr());

// #define USE_2XmV 1	// 2nd order, a-stable, oscillations
// #define USE_2XpVd3 1	// 1st order, a/l-stable, more accurate
#define USE_X 1		// 1st order, l-stable (implicit Euler), less accurate (alpha == 1.)
// #define USE_alphaX_betaV	(1.0)	// alpha = 1. + |rho|; beta = (1 - alpha) for 1st order accuracy

#if USE_2XmV
					Vec3 VRef((XCurr - XPrev)*(2./h) - VPrev);
#elif USE_2XpVd3
					Vec3 VRef((XCurr - XPrev)*(2./3./h) + VPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 VRef((XCurr - XPrev)*(USE_alphaX_betaV/h) + VPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 VRef((XCurr - XPrev)/h);
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*VRef(iCnt);
					}

					const Vec3& WPrev(pNode->GetWPrev());
						(void) WPrev; // silence set but not used warning
					const Mat3x3& RPrev(pNode->GetRPrev());
					const Mat3x3& RCurr(pNode->GetRCurr());

#if USE_2XmV
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(2./h) - WPrev);
#elif USE_2XpVd3
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(2./3./h) + WPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))*(USE_alphaX_betaV/h) + WPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 WRef(RotManip::VecRot(RCurr.MulMT(RPrev))/h);
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*WRef(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 BPrev(pNode->GetVPrev()*dm - S.Cross(pNode->GetWPrev()));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += BPrev(iCnt);
					}

					Vec3 GPrev(S.Cross(pNode->GetVPrev()) + J*pNode->GetWPrev());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += GPrev(iCnt);
					}
				}
			}
			break;

		case InverseDynamics::ACCELERATION:
			if (dw1 > 0.) {
				doublereal h = DrvHdl.dGetTimeStep();

				// xpp_k = xpp_km1/3 + 2/3*(xp_k - xp_km1)/h
				// xpp_k = 2*(xp_k - xp_km1)/h - xpp_km1
				// xpp_k = xpp_km1 + 6*(xp_k + xp_km1)/h - 12*(x_k - x_km1)/h^2
				for (NodeContainerType::const_iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
					n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
				{
					const StructNode *pNode = dynamic_cast<const StructNode *>(n->second);

					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					const Vec3& XPPPrev(pNode->GetXPPPrev());
						(void) XPPPrev; // silence set but not used warning
					const Vec3& VPrev(pNode->GetVPrev());
					const Vec3& VCurr(pNode->GetVCurr());

#if USE_2XmV
					Vec3 XPPRef((VCurr - VPrev)*(2./h) - XPPPrev);
#elif USE_2XpVd3
					Vec3 XPPRef((VCurr - VPrev)*(2./3./h) + XPPPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 XPPRef((VCurr - VPrev)*(USE_alphaX_betaV/h) + XPPPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 XPPRef((VCurr - VPrev)/h);
#endif
#if 0
					const Vec3& XPrev(pNode->GetXPrev());
					const Vec3& XCurr(pNode->GetXCurr());

					Vec3 XPPRef(XPPPrev + (VCurr + VPrev)*(6./h) - (XCurr - XPrev)*(12./h/h));
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += dw1*XPPRef(iCnt);
					}

					const Vec3& WPPrev(pNode->GetWPPrev());
						(void) WPPrev; // silence set but not used warning
					const Vec3& WPrev(pNode->GetWPrev());
					const Vec3& WCurr(pNode->GetWCurr());

#if USE_2XmV
					Vec3 WPRef((WCurr - WPrev)*(2./h) - WPPrev);
#elif USE_2XpVd3
					Vec3 WPRef((WCurr - WPrev)*(2./3./h) + WPPrev/3.);
#elif defined(USE_alphaX_betaV)
					Vec3 WPRef((WCurr - WPrev)*(USE_alphaX_betaV/h) + WPPrev*(1 - USE_alphaX_betaV));
#elif USE_X
					Vec3 WPRef((WCurr - WPrev)/h);
#endif
#if 0
					const Mat3x3& RPrev(pNode->GetRPrev());
					const Mat3x3& RCurr(pNode->GetRCurr());

					Vec3 WPRef(WPPrev + (WCurr + WPrev)*(6./h) - RotManip::VecRot(RCurr.MulMT(RPrev))*(12./h/h));
#endif
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += dw1*WPRef(iCnt);
					}
				}
			}

			if (dw2 > 0.) {
				for (ElemContainerType::const_iterator b = ElemData[Elem::BODY].ElemContainer.begin();
					b != ElemData[Elem::BODY].ElemContainer.end(); ++b)
				{
					if (!b->second->bIsErgonomy()) {
						continue;
					}

					const Body *pBody(Cast<Body>(b->second));
					doublereal dm(pBody->dGetM()*dw2);
					Vec3 S(pBody->GetS()*dw2);
					Mat3x3 J = (pBody->GetJ()*dw2);
					const StructNode *pNode = pBody->pGetNode();
					ASSERT(pNode->iGetNumDof() == 6);
					integer iFirstIndex = pNode->iGetFirstIndex();

					Vec3 BPPrev(pNode->GetXPPPrev()*dm - S.Cross(pNode->GetWPPrev()));
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + iCnt) += BPPrev(iCnt);
					}

					Vec3 GPPrev(S.Cross(pNode->GetXPPPrev()) + J*pNode->GetWPPrev());
					for (integer iCnt = 1; iCnt <= 3; iCnt++) {
						ResHdl(iFirstIndex + 3 + iCnt) += GPPrev(iCnt);
					}
				}
			}
			break;

		default:
			ASSERT(0);
		}

		for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
			j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
		{
			bool bActive(true);
			Joint *pJ = Cast<Joint>(j->second, true);
			if (pJ == 0) {
				bActive = false;
				pJ = Cast<Joint>(j->second, false);
			}

			if ((pJ->bIsPrescribedMotion()
				|| (iOrder == InverseDynamics::POSITION && pJ->bIsErgonomy())) && bActive)
			{
				try {
					ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, iOrder);
				}
				catch (Elem::ChangedEquationStructure& e) {
					ResHdl += WorkVec;
					ChangedEqStructure = true;
				}

			} else if (pJ->bIsTorque() || (pJ->bIsPrescribedMotion() && !bActive)) {
				integer iNumDof = pJ->iGetNumDof();
				integer iFirstIndex = pJ->iGetFirstIndex();
				if (iOrder == InverseDynamics::POSITION) {
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						ResHdl(iFirstIndex + iCnt) = -(*pXCurr)(iFirstIndex + iCnt);
					}

				} else {
					for (int iCnt = 1; iCnt <= iNumDof; iCnt++) {
						ResHdl(iFirstIndex + iCnt) = 0.;
					}
				}
			}
		}

		for (ElemContainerType::iterator j = ElemData[Elem::BEAM].ElemContainer.begin();
			j != ElemData[Elem::BEAM].ElemContainer.end(); ++j)
		{
			bool bActive(true);
			Beam2 *pB = Cast<Beam2>(j->second, true);
			if (pB == 0) {
				bActive = false;
				pB = Cast<Beam2>(j->second, false);
			}

			if (pB && (iOrder == InverseDynamics::POSITION && pB->bIsErgonomy()) && bActive)
			{
				try {
					ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
						*pXPrimeCurr, *pXPrimePrimeCurr, iOrder);
				}
				catch (Elem::ChangedEquationStructure& e) {
					ResHdl += WorkVec;
					ChangedEqStructure = true;
				}
			}
		}
	}

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (ChangedEqStructure) {
		throw ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

/* Equilibrium residual assembly, no constraints */
void
DataManager::AssRes(VectorHandler& ResHdl)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);
	AssRes(ResHdl, ElemIter, *pWorkVec);

	for (ElemContainerType::iterator j = ElemData[Elem::INERTIA].ElemContainer.begin();
		j != ElemData[Elem::INERTIA].ElemContainer.end(); ++j)
	{
		bool bActive(true);
		Inertia *pI = Cast<Inertia>(j->second, true);
		if (pI == 0) {
			bActive = false;
			pI = Cast<Inertia>(j->second, false);
		}
			if (pI && bActive)
		{
			ResHdl += j->second->AssRes(*pWorkVec, *pXCurr,
				*pXPrimeCurr, *pXPrimePrimeCurr, 
				InverseDynamics::INVERSE_DYNAMICS);
		}
	}
}

void
DataManager::AssRes(VectorHandler& ResHdl,
	VecIter<Elem *> &Iter,
	SubVectorHandler& WorkVec)
{
	DEBUGCOUT("Entering AssRes()" << std::endl);

	// TODO
	// FIXME: could be as the rest?

	const Elem::Type ElemType[] = {
		Elem::BODY,
		Elem::BEAM,
		Elem::FORCE,
		
		Elem::LASTELEMTYPE
	};

	bool ChangedEqStructure(false);
	
	for (int et = 0; ElemType[et] != Elem::LASTELEMTYPE; et++) {
		for (ElemContainerType::iterator j = ElemData[ElemType[et]].ElemContainer.begin();
			j != ElemData[ElemType[et]].ElemContainer.end(); ++j)
		{
			if (!j->second->bIsRightHandSide()) {
				continue;
			}

			try {
				ResHdl += j->second->AssRes(WorkVec, *pXCurr, 
					*pXPrimeCurr, *pXPrimePrimeCurr, 
					InverseDynamics::INVERSE_DYNAMICS);
			}
			catch (Elem::ChangedEquationStructure& e) {
				ResHdl += WorkVec;
				ChangedEqStructure = true;
			}
		}
	}

	for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
		j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
	{
		bool bActive(true);
		Joint *pJ = Cast<Joint>(j->second, true);
		if (pJ == 0) {
			bActive = false;
			pJ = Cast<Joint>(j->second, false);
		}

		if (bActive && pJ->bIsRightHandSide()) {
			try {
				ResHdl += pJ->AssRes(WorkVec, *pXCurr, 
					*pXPrimeCurr, *pXPrimePrimeCurr, 
					InverseDynamics::INVERSE_DYNAMICS);
			}
			catch (Elem::ChangedEquationStructure& e) {
				ResHdl += WorkVec;
				ChangedEqStructure = true;
			}
		}
	}

	if (ChangedEqStructure) {
		throw ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

void
DataManager::Update(InverseDynamics::Order iOrder) const
{
	const VectorHandler* pV = 0;

	switch (iOrder) {
	case InverseDynamics::INVERSE_DYNAMICS:
		// Update constraints reactions (for output only...)
		for (ElemContainerType::const_iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
			j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
		{
			j->second->Update(*pLambdaCurr, iOrder);
		}
		return;

	// Nodes:
	case InverseDynamics::POSITION:
		// Update nodes positions
		pV = pXCurr;	
		break;

	case InverseDynamics::VELOCITY:
		// Update nodes velocities
		pV = pXPrimeCurr;
		break;

	case InverseDynamics::ACCELERATION:
		// Update nodes accelerations
		pV = pXPrimePrimeCurr;
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(pV != 0);
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->Update(*pV, iOrder);
	}
}

void
DataManager::IDAfterConvergence(void) const
{
	DEBUGCOUTFNAME("DataManager::IDAfterConvergence");

	// Nodes:
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		(*i)->AfterConvergence(*pXCurr, *pXPrimeCurr, *pXPrimePrimeCurr);
	}

	Elem* pEl = NULL;
	if (ElemIter.bGetFirst(pEl)) {
		do {
			pEl->AfterConvergence(*pXCurr, *pXPrimeCurr, *pXPrimePrimeCurr);
		} while (ElemIter.bGetNext(pEl));
	}
}

void
DataManager::IDDofOwnerSet(void)
{
	DEBUGCOUTFNAME("DataManager::IDDofOwnerSet");

	/* Setta i DofOwner dei nodi */
	for (NodeVecType::const_iterator i = Nodes.begin(); i != Nodes.end(); ++i) {
		DofOwner* pDO = const_cast<DofOwner *>((*i)->pGetDofOwner());
		pDO->iNumDofs = (*i)->iGetNumDof();
		iIDNodeTotNumDofs += pDO->iNumDofs;
	}

	for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
                j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
        {
		ElemWithDofs* pEWD = Cast<ElemWithDofs>(j->second);
		iIDJointTotNumDofs += pEWD->iGetNumDof();
        }
	
	/* Setta i DofOwner degli elementi (chi li possiede) */
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		DofOwner::Type DT = ElemData[iCnt].DofOwnerType;
		if (DT != DofOwner::UNKNOWN) {
			DEBUGLCOUT(MYDEBUG_INIT, "Elem type " << iCnt
					<< " (" << psElemNames[iCnt] << ")"
					<< std::endl);

			for (ElemContainerType::const_iterator p = ElemData[iCnt].ElemContainer.begin();
				p != ElemData[iCnt].ElemContainer.end(); ++p)
			{
				ElemWithDofs* pEWD = Cast<ElemWithDofs>(p->second);

				DEBUGLCOUT(MYDEBUG_INIT, "    " << psElemNames[pEWD->GetElemType()]
						<< "(" << pEWD->GetLabel() << ")" << std::endl);

				DofOwner* pDO = const_cast<DofOwner *>(pEWD->pGetDofOwner());
				pDO->iNumDofs = pEWD->iGetNumDof();
				DEBUGLCOUT(MYDEBUG_INIT, "    num dofs: " << pDO->iNumDofs << std::endl);
			}
		}
	}
}

void
DataManager::IDDofInit(void)
{  
	if (iTotDofOwners == 0) {	
		silent_cerr("DataManager::IDDofInit: no dof owners are defined" << std::endl);

		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
      
	/* Di ogni DofOwner setta il primo indice
	 * e calcola il numero totale di Dof:
	 * poichè per il problema inverso non si 
	 * possono aggiungere incognite diverse
	 * da posizione (velocità e accelerazione)
	 * dei nodi e reazioni vincolari, viene
	 * controllato che gli elementi che aggiungono 
	 * dof siano solo nodi e vincoli.
	 * 
	 * Per problemi mal posti (DoF nodi != Dof vincoli): 
	 * iTotDofs = (DoF nodi) + (Dof vincoli), altrimenti
	 * 
	 * Per problemi ben posti:
	 * iTotDofs = DoF nodi (= DoF vincoli)
	 */

	integer iRealTotDofs = 0;
	
	/* Mette gli indici ai DofOwner dei nodi strutturali: */
	/* contatore dei Dof dei nodi */
	integer iNodeIndex = 0;

	NodeContainerType::const_iterator i = NodeData[Node::STRUCTURAL].NodeContainer.begin();
	for (int iDOm1 = 0; iDOm1 < DofData[DofOwner::STRUCTURALNODE].iNum;
		++iDOm1, ++i)
	{
		DofOwner *pDO = const_cast<DofOwner *>(i->second->pGetDofOwner());
		unsigned iNumDofs = pDO->iNumDofs = i->second->iGetNumDof();
		if (iNumDofs > 0) {
			pDO->iFirstIndex = iNodeIndex;
			iNodeIndex += iNumDofs;
			iRealTotDofs += iNumDofs;

		} else {
			// note: this could only be possible for dummy nodes?
			pDO->iFirstIndex = -1;
			DEBUGCERR("warning, Structural(" << i->second->GetLabel() << ") "
				"(DofOwner #" << (iDOm1 + 1) << ") has 0 dofs" << std::endl);
		}
	}

	ASSERT(iNodeIndex == iIDNodeTotNumDofs);

	/* Gli indici dei nodi sono ok */
	
	if (bOutputAccels) {
		for (NodeContainerType::iterator n = NodeData[Node::STRUCTURAL].NodeContainer.begin();
			n != NodeData[Node::STRUCTURAL].NodeContainer.end(); ++n)
		{
			dynamic_cast<StructNode *>(n->second)->OutputAccelerations(true);
		}
	}

	/* Se il problema è ben posto, gli indici delle equazioni di vincolo
	 * hanno numerazione indipendente dai nodi. Altrimenti la numerazione 
	 * è a partire dagli indici dei nodi (per fare spazio alla matrice 
	 * peso nello jacobiano) */

	/* contatore dei Dof dei joint */
	integer iJointIndex;
	integer iPrescribedMotionJointIndex;
	integer iTorqueJointIndex;
	switch (dynamic_cast<InverseSolver *>(pSolver)->GetProblemType()) {
	case InverseDynamics::FULLY_ACTUATED_COLLOCATED:
		if (iIDNodeTotNumDofs != iIDJointTotNumDofs) {
			silent_cerr("DataManager::IDDofInit: nodeDoFs=" << iIDNodeTotNumDofs
				<< " jointDoFs=" << iIDJointTotNumDofs << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	case InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED:
		iJointIndex = 0;
		iPrescribedMotionJointIndex = 0;
		iTorqueJointIndex = 0;
		break;

	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED:
		iJointIndex = iNodeIndex;
		iPrescribedMotionJointIndex = iNodeIndex;
		iTorqueJointIndex = iNodeIndex;
		break;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (ElemContainerType::iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
		j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
	{
		Joint *pJ = Cast<Joint>(j->second);
		ASSERT(pJ != 0);

		DofOwner *pDO = const_cast<DofOwner *>(pJ->pGetDofOwner());
		ASSERT(pDO != 0);

		unsigned iNumDofs = pDO->iNumDofs;
		if (iNumDofs > 0) {
			pDO->iFirstIndex = iJointIndex;
			iRealTotDofs += iNumDofs;
			iJointIndex += iNumDofs;
			if (pJ->bIsPrescribedMotion()) {
				iPrescribedMotionJointIndex += iNumDofs;
			}
			if (pJ->bIsTorque()) {
				iTorqueJointIndex += iNumDofs;
			}

		} else {
			pDO->iFirstIndex = -1;
			DEBUGCERR("warning, Joint(" << j->second->GetLabel() << ") "
				"has 0 dofs" << std::endl);
		}

		const char *sTorque = pJ->bIsTorque() ? "T" : "_";
		const char *sPrescribedMotion = pJ->bIsPrescribedMotion() ? "P" : "_";
		const char *sErgonomy = pJ->bIsErgonomy() ? "E" : "_";
		const char *sRightHandSide = pJ->bIsRightHandSide() ? "R" : "_";

		silent_cout("Joint(" << pJ->GetLabel() << "): " << sTorque << sPrescribedMotion << sErgonomy << sRightHandSide << std::endl);
	}

	for (ElemContainerType::iterator j = ElemData[Elem::BEAM].ElemContainer.begin();
		j != ElemData[Elem::BEAM].ElemContainer.end(); ++j)
	{
		Beam2 *pB = Cast<Beam2>(j->second);
		ASSERT(pB != 0);

		const char *sTorque = "_";
		const char *sPrescribedMotion = "_";
		const char *sErgonomy = pB->bIsErgonomy() ? "E" : "_";
		const char *sRightHandSide = pB->bIsRightHandSide() ? "R" : "_";

		silent_cout("Beam2(" << pB->GetLabel() << "): " << sTorque << sPrescribedMotion << sErgonomy << sRightHandSide << std::endl);
	}

	switch (dynamic_cast<InverseSolver *>(pSolver)->GetProblemType()) {
	case InverseDynamics::FULLY_ACTUATED_COLLOCATED:
		ASSERT(iIDNodeTotNumDofs == iIDJointTotNumDofs);
		// fall thru

	case InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED:
		iTotDofs = iNodeIndex;
		break;

	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED:
		iTotDofs = iTorqueJointIndex;
		break;

	default:
		break;
	}

	iTotDofs = iJointIndex;

	DEBUGLCOUT(MYDEBUG_INIT, "iTotDofs = " << iTotDofs << std::endl);

 
	/* Crea la struttura dinamica dei Dof */
	if (iTotDofs == 0) {
		silent_cerr("DataManager::IDDofInit: no dofs defined" << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Dofs.resize(iRealTotDofs); /* Inizializza l'iteratore sui Dof */

	/* Inizializza la struttura dinamica dei Dof */
	integer iIndex = DofOwners[0].iFirstIndex;
	for (integer idx = 0; idx < iRealTotDofs; idx++) {
		Dofs[idx].iIndex = iIndex++;
		Dofs[idx].Order = DofOrder::DIFFERENTIAL;
	}
}  

int
DataManager::iIDGetNodeTotNumDofs(void) const
{
	return iIDNodeTotNumDofs;
}

int
DataManager::iIDGetJointTotNumDofs(void) const
{
	return iIDJointTotNumDofs;
}

int
DataManager::iIDGetTotNumDofs(void) const
{
	switch (dynamic_cast<InverseSolver *>(pSolver)->GetProblemType()) {
	case InverseDynamics::FULLY_ACTUATED_COLLOCATED:
	case InverseDynamics::FULLY_ACTUATED_NON_COLLOCATED:
		ASSERT(iIDNodeTotNumDofs == iIDJointTotNumDofs);
		return iIDNodeTotNumDofs;

	case InverseDynamics::UNDERDETERMINED_UNDERACTUATED_COLLOCATED:
	case InverseDynamics::UNDERDETERMINED_FULLY_ACTUATED:
		return iIDNodeTotNumDofs + iIDJointTotNumDofs;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
DataManager::IDSetTest(NonlinearSolverTestRange *pResTest, NonlinearSolverTestRange *pSolTest, bool bFullResTest)
{
	integer iFirstIndex = std::numeric_limits<integer>::max();
	integer iLastIndex = -1;

	for (NodeVecType::const_iterator n = Nodes.begin(); n != Nodes.end(); ++n) {
		integer iIndex = (*n)->iGetFirstIndex();
		integer iNumDofs = (*n)->iGetNumDof();

		if (iFirstIndex > iIndex) {
			iFirstIndex = iIndex;
		}

		if (iLastIndex < iIndex + iNumDofs) {
			iLastIndex = iIndex + iNumDofs;
		}
	}

	silent_cout("DataManager::IDSetTest(): SolTest range=[" << iFirstIndex + 1 << ":" << iLastIndex << "]" << std::endl);

	pSolTest->SetRange(iFirstIndex + 1, iLastIndex);

	iFirstIndex = std::numeric_limits<integer>::max();
	iLastIndex = -1;

	for (ElemContainerType::const_iterator j = ElemData[Elem::JOINT].ElemContainer.begin();
                j != ElemData[Elem::JOINT].ElemContainer.end(); ++j)
        {
		const Joint *pJ = Cast<Joint>(j->second);
		if (!pJ->bIsTorque() && !pJ->bIsPrescribedMotion()) {
			continue;
		}

		const ElemWithDofs* pEWD = Cast<ElemWithDofs>(j->second);
		integer iNumDofs = pEWD->iGetNumDof();
		if (iNumDofs > 0) {
			integer iIndex = pEWD->iGetFirstIndex();

			if (iFirstIndex > iIndex) {
				iFirstIndex = iIndex;
			}

			if (iLastIndex < iIndex + iNumDofs) {
				iLastIndex = iIndex + iNumDofs;
			}
		}
	}

	if (bFullResTest) {
		// NOTE: assumes *ALL* residual elements need to be evaluated!  Use with care...
		iFirstIndex = 0;
	}

	silent_cout("DataManager::IDSetTest(): ResTest range=[" << iFirstIndex + 1 << ":" << iLastIndex << "]" << std::endl);

	pResTest->SetRange(iFirstIndex + 1, iLastIndex);
}
