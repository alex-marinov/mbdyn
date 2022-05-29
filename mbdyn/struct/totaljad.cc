/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "totaljad.h"

TotalJointAd::TotalJointAd(unsigned int uL,
                           const DofOwner *pDO,
                           bool bPos[3],
                           bool bVel[3],
                           TplDriveCaller<Vec3> *const pDCPos[3],
                           bool bRot[3],
                           bool bAgv[3],	/* Agv stands for AnGular Velocity */
                           TplDriveCaller<Vec3> *const pDCRot[3],
                           const StructNodeAd *pN1,
                           const Vec3& f1Tmp,
                           const Mat3x3& R1hTmp,
                           const Mat3x3& R1hrTmp,
                           const StructNodeAd *pN2,
                           const Vec3& f2Tmp,
                           const Mat3x3& R2hTmp,
                           const Mat3x3& R2hrTmp,
                           flag fOut)
:Elem(uL, fOut),
 TotalJoint(uL, pDO, bPos, bVel, pDCPos, bRot, bAgv, pDCRot, pN1, f1Tmp, R1hTmp, R1hrTmp, pN2, f2Tmp, R2hTmp, R2hrTmp, fOut),
 pNode1(pN1),
 pNode2(pN2)
{
     ASSERT(pN1 != nullptr);
     ASSERT(pN2 != nullptr);
}

TotalJointAd::~TotalJointAd()
{

};

VariableSubMatrixHandler&
TotalJointAd::AssJac(VariableSubMatrixHandler& WorkMat,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUT("Entering TotalJointAd::AssJac()" << std::endl);

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
TotalJointAd::AssJac(VectorHandler& JacY,
                     const VectorHandler& Y,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr,
                     VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
TotalJointAd::AssRes(SubVectorHandler& WorkVec,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUT("Entering TotalJointAd::AssRes()" << std::endl);

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);
     return WorkVec;
}

template <typename T>
void
TotalJointAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                     doublereal dCoef,
                     const sp_grad::SpGradientVectorHandler<T>& XCurr,
                     const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                     sp_grad::SpFunctionCall func)
{
        using namespace sp_grad;

        if (iGetNumDof() == 0) {
             return;
        }

        const integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
        const integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
        const integer iFirstReactionIndex = iGetFirstIndex();

        SpColVectorA<T, 3> F, M;

        for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
                XCurr.dGetCoef(iFirstReactionIndex + 1 + iPosEqIndex[iCnt], F(iPosIncid[iCnt]), 1.);
        }

        for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
                XCurr.dGetCoef(iFirstReactionIndex + 1 + iRotEqIndex[iCnt], M(iRotIncid[iCnt]), 1.);
        }

        for (unsigned iCnt = 0; iCnt < nVelConstraints; iCnt++) {
                XCurr.dGetCoef(iFirstReactionIndex + 1 + iVelEqIndex[iCnt], F(iVelIncid[iCnt]), 1.);
        }

        for (unsigned iCnt = 0; iCnt < nAgvConstraints; iCnt++) {
                XCurr.dGetCoef(iFirstReactionIndex + 1 + iAgvEqIndex[iCnt], M(iAgvIncid[iCnt]), 1.);
        }

        SpColVectorA<T, 3> X1, X2, V1, V2, W1, W2;
        SpMatrixA<T, 3, 3> R1, R2;

        pNode1->GetXCurr(X1, dCoef, func);
        pNode1->GetRCurr(R1, dCoef, func);
        pNode2->GetXCurr(X2, dCoef, func);
        pNode2->GetRCurr(R2, dCoef, func);

        if (nVelConstraints) {
                pNode1->GetVCurr(V1, dCoef, func);
                pNode2->GetVCurr(V2, dCoef, func);
        }

        if (nVelConstraints || nAgvConstraints) {
                pNode1->GetWCurr(W1, dCoef, func);
                pNode2->GetWCurr(W2, dCoef, func);
        }

        SpColVector<T, 3> b2 = R2 * f2;
        SpColVector<T, 3> b1 = X2 + b2 - X1;

        SpMatrix<T, 3, 3> R1R1h = R1 * R1h;
        SpMatrix<T, 3, 3> R1r = R1 * R1hr;

        SpColVectorA<T, 3> XDelta;

        if (nPosConstraints) {
                XDelta = Transpose(R1R1h) * b1 - tilde_f1 - XDrv.Get();
        }

        SpColVectorA<T, 3> VDelta;

        if (nVelConstraints) {
                VDelta = Transpose(R1R1h) * (Cross(b1, W1) + V2 - Cross(b2, W2) - V1) - XDrv.Get();
        }

        SpColVectorA<T, 3> ThetaDelta;

        if (nRotConstraints) {
                SpMatrix<T, 3, 3> R2r = R2 * R2hr;

                SpColVector<doublereal, 3> ThetaDrvTmp = ThetaDrv.Get();
                if (nRotConstraints < 3) {
                        for (int i = nRotConstraints; i < 3; i++) {
                                // zero out unconstrained components of drive
                                ThetaDrvTmp(iRotIncid[i]) = 0.;
                        }
                        // add remnant to make ThetaDelta as close to zero
                        // as possible
                        ThetaDrvTmp += ThetaDeltaRemnant;
                }
                SpMatrix<doublereal, 3, 3> R0 = MatRotVec(ThetaDrvTmp);
                SpMatrix<T, 3, 3> RDelta = Transpose(R1r) * (R2r * Transpose(R0));
                ThetaDelta = VecRotMat(RDelta);
        }

        SpColVectorA<T, 3> WDelta;

        if (nAgvConstraints) {
                WDelta = Transpose(R1r) * (W2 - W1) - ThetaDrv.Get();
        }

        SpColVector<T, 3> FTmp = R1R1h * F;
        SpColVector<T, 3> MTmp = R1r * M;

        WorkVec.AddItem(iNode1FirstMomIndex + 1, FTmp);
        WorkVec.AddItem(iNode1FirstMomIndex + 4, SpColVector<T, 3>(MTmp + Cross(b1, FTmp)));

        WorkVec.AddItem(iNode2FirstMomIndex + 1, SpColVector<T, 3>(-FTmp));
        WorkVec.AddItem(iNode2FirstMomIndex + 4, SpColVector<T, 3>(-MTmp - Cross(b2, FTmp)));

        for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
                WorkVec.AddItem(iFirstReactionIndex + iPosEqIndex[iCnt] + 1, XDelta(iPosIncid[iCnt]) / -dCoef);
        }

        for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
                WorkVec.AddItem(iFirstReactionIndex + iRotEqIndex[iCnt] + 1, ThetaDelta(iRotIncid[iCnt]) / -dCoef);
        }

        for (unsigned iCnt = 0; iCnt < nVelConstraints; iCnt++) {
                WorkVec.AddItem(iFirstReactionIndex + iVelEqIndex[iCnt] + 1, -VDelta(iVelIncid[iCnt]));
        }

        for (unsigned iCnt = 0; iCnt < nAgvConstraints; iCnt++) {
                WorkVec.AddItem(iFirstReactionIndex + iAgvEqIndex[iCnt] + 1, -WDelta(iAgvIncid[iCnt]));
        }

        UpdateThetaDelta(ThetaDelta);
        UpdateF(F);
        UpdateM(M);
}

void
TotalJointAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12 + nConstraints ;
     *piNumCols = 0;
}

void
TotalJointAd::UpdateThetaDelta(const sp_grad::SpColVector<doublereal, 3>& ThetaDeltaCurr)
{
     for (sp_grad::index_type i = 1; i <= ThetaDeltaCurr.iGetNumRows(); ++i) {
          ThetaDelta(i) = ThetaDeltaCurr(i);
     }
}

void
TotalJointAd::UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalJointAd::UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}

void
TotalJointAd::UpdateF(const sp_grad::SpColVector<doublereal, 3>& FCurr)
{
        for (sp_grad::index_type i = 1; i <= FCurr.iGetNumRows(); ++i) {
                F(i) = FCurr(i);
        }
}

void
TotalJointAd::UpdateF(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalJointAd::UpdateF(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}

void
TotalJointAd::UpdateM(const sp_grad::SpColVector<doublereal, 3>& MCurr)
{
        for (sp_grad::index_type i = 1; i <= MCurr.iGetNumRows(); ++i) {
                M(i) = MCurr(i);
        }
}

void
TotalJointAd::UpdateM(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalJointAd::UpdateM(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}

TotalPinJointAd::TotalPinJointAd(unsigned int uL,
                                 const DofOwner *pDO,
                                 bool bPos[3],
                                 bool bVel[3],
                                 TplDriveCaller<Vec3> *const pDCPos[3],
                                 bool bRot[3],
                                 bool bAgv[3],
                                 TplDriveCaller<Vec3> *const pDCRot[3],
                                 const Vec3& XcTmp,
                                 const Mat3x3& RchTmp,
                                 const Mat3x3& RchrTmp,
                                 const StructNodeAd *pN,
                                 const Vec3& fnTmp,
                                 const Mat3x3& RnhTmp,
                                 const Mat3x3& RnhrTmp,
                                 flag fOut)
:Elem(uL, fOut),
 TotalPinJoint(uL, pDO, bPos, bVel, pDCPos, bRot, bAgv, pDCRot, XcTmp, RchrTmp, RchrTmp, pN, fnTmp, RnhTmp, RnhrTmp, fOut),
 pNode(pN)
{
     ASSERT(pN != nullptr);
}

TotalPinJointAd::~TotalPinJointAd(void)
{
}

VariableSubMatrixHandler&
TotalPinJointAd::AssJac(VariableSubMatrixHandler& WorkMat,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{

     DEBUGCOUT("Entering TotalPinJointAd::AssJac()" << std::endl);

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
TotalPinJointAd::AssJac(VectorHandler& JacY,
                        const VectorHandler& Y,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr,
                        VariableSubMatrixHandler& WorkMat)
{
     using namespace sp_grad;

     SpGradientAssVec<GpGradProd>::AssJac(this,
                                          JacY,
                                          Y,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
}

SubVectorHandler&
TotalPinJointAd::AssRes(SubVectorHandler& WorkVec,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
     DEBUGCOUT("Entering TotalPinJointAd::AssRes()" << std::endl);

     sp_grad::SpGradientAssVec<doublereal>::AssRes(this,
                                                   WorkVec,
                                                   dCoef,
                                                   XCurr,
                                                   XPrimeCurr,
                                                   sp_grad::SpFunctionCall::REGULAR_RES);
     return WorkVec;
}

template <typename T>
void
TotalPinJointAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const sp_grad::SpGradientVectorHandler<T>& XCurr,
                        const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                        sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     if (iGetNumDof() == 0) {
          return;
     }

     const integer iNode1FirstMomIndex = pNode->iGetFirstMomentumIndex();
     const integer iFirstReactionIndex = iGetFirstIndex();

     SpColVectorA<T, 3> F, M;

     for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
          XCurr.dGetCoef(iFirstReactionIndex + 1 + iPosEqIndex[iCnt], F(iPosIncid[iCnt]), 1.);
     }

     for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
          XCurr.dGetCoef(iFirstReactionIndex + 1 + iRotEqIndex[iCnt], M(iRotIncid[iCnt]), 1.);
     }

     for (unsigned iCnt = 0; iCnt < nVelConstraints; iCnt++) {
          XCurr.dGetCoef(iFirstReactionIndex + 1 + iVelEqIndex[iCnt], F(iVelIncid[iCnt]), 1.);
     }

     for (unsigned iCnt = 0; iCnt < nAgvConstraints; iCnt++) {
          XCurr.dGetCoef(iFirstReactionIndex + 1 + iAgvEqIndex[iCnt], M(iAgvIncid[iCnt]), 1.);
     }

     SpColVectorA<T, 3> X1, V1, W1;
     SpMatrixA<T, 3, 3> R1;

     pNode->GetXCurr(X1, dCoef, func);
     pNode->GetRCurr(R1, dCoef, func);

     if (nVelConstraints) {
          pNode->GetVCurr(V1, dCoef, func);
     }

     if (nVelConstraints || nAgvConstraints) {
          pNode->GetWCurr(W1, dCoef, func);
     }

     SpColVector<T, 3> fn = R1 * tilde_fn;
     SpColVectorA<T, 3> XDelta;

     if (nPosConstraints) {
          XDelta = RchT * (X1 + fn) - tilde_Xc - XDrv.Get();
     }

     SpColVectorA<T, 3> VDelta;

     if (nVelConstraints) {
          VDelta = RchT * (V1 - Cross(fn, W1)) - XDrv.Get();
     }

     SpColVectorA<T, 3> ThetaDelta;

     if (nRotConstraints) {
          SpMatrix<T, 3, 3> Rnhr = R1 * tilde_Rnhr;

          SpColVector<doublereal, 3> ThetaDrvTmp = ThetaDrv.Get();

          if (nRotConstraints < 3) {
               for (int i = nRotConstraints; i < 3; i++) {
                    // zero out unconstrained components of drive
                    ThetaDrvTmp(iRotIncid[i]) = 0.;
               }
               // add remnant to make ThetaDelta as close to zero
               // as possible
               ThetaDrvTmp += ThetaDeltaRemnant;
          }
          SpMatrix<doublereal, 3, 3> R0 = MatRotVec(ThetaDrvTmp);
          SpMatrix<T, 3, 3> RDelta = RchrT * (Rnhr * Transpose(R0));
          ThetaDelta = VecRotMat(RDelta);
     }

     SpColVectorA<T, 3> WDelta;

     if (nAgvConstraints) {
          WDelta = RchT * W1 - ThetaDrv.Get();
     }

     SpColVector<T, 3> FTmp = -Rch * F;
     SpColVector<T, 3> MTmp = -Rchr * M;

     WorkVec.AddItem(iNode1FirstMomIndex + 1, FTmp);
     WorkVec.AddItem(iNode1FirstMomIndex + 4, SpColVector<T, 3>(MTmp + Cross(fn, FTmp)));

     for (unsigned iCnt = 0; iCnt < nPosConstraints; iCnt++) {
          WorkVec.AddItem(iFirstReactionIndex + iPosEqIndex[iCnt] + 1, XDelta(iPosIncid[iCnt]) / -dCoef);
     }

     for (unsigned iCnt = 0; iCnt < nRotConstraints; iCnt++) {
          WorkVec.AddItem(iFirstReactionIndex + iRotEqIndex[iCnt] + 1, ThetaDelta(iRotIncid[iCnt]) / -dCoef);
     }

     for (unsigned iCnt = 0; iCnt < nVelConstraints; iCnt++) {
          WorkVec.AddItem(iFirstReactionIndex + iVelEqIndex[iCnt] + 1, -VDelta(iVelIncid[iCnt]));
     }

     for (unsigned iCnt = 0; iCnt < nAgvConstraints; iCnt++) {
          WorkVec.AddItem(iFirstReactionIndex + iAgvEqIndex[iCnt] + 1, -WDelta(iAgvIncid[iCnt]));
     }

     UpdateThetaDelta(ThetaDelta);
     UpdateF(F);
     UpdateM(M);
}

void
TotalPinJointAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6 + nConstraints ;
     *piNumCols = 0;
}

void
TotalPinJointAd::UpdateThetaDelta(const sp_grad::SpColVector<doublereal, 3>& ThetaDeltaCurr)
{
     for (sp_grad::index_type i = 1; i <= ThetaDeltaCurr.iGetNumRows(); ++i) {
          ThetaDelta(i) = ThetaDeltaCurr(i);
     }
}

void
TotalPinJointAd::UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalPinJointAd::UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}

void
TotalPinJointAd::UpdateF(const sp_grad::SpColVector<doublereal, 3>& FCurr)
{
        for (sp_grad::index_type i = 1; i <= FCurr.iGetNumRows(); ++i) {
                F(i) = FCurr(i);
        }
}

void
TotalPinJointAd::UpdateF(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalPinJointAd::UpdateF(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}

void
TotalPinJointAd::UpdateM(const sp_grad::SpColVector<doublereal, 3>& MCurr)
{
        for (sp_grad::index_type i = 1; i <= MCurr.iGetNumRows(); ++i) {
                M(i) = MCurr(i);
        }
}

void
TotalPinJointAd::UpdateM(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&)
{
}

void
TotalPinJointAd::UpdateM(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&)
{
}
