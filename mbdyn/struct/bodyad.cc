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
#include "bodyad.h"

BodyAd::BodyAd(unsigned int uL,
               const StructNodeAd *pNode,
               doublereal dMass,
               const Vec3& Xgc,
               const Mat3x3& J,
               flag fOut)
     :Body(uL, pNode, dMass, Xgc, J, fOut)
{
     ASSERT(pNode != NULL);
     ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
     ASSERT(dMass > 0.);
}

BodyAd::~BodyAd(void)
{
}

template <typename T>
void
BodyAd::AssVecRBK_int(const sp_grad::SpColVector<T, 3>& STmp,
                      const sp_grad::SpMatrix<T, 3, 3>& JTmp,
                      sp_grad::SpGradientAssVec<T>& WorkVec,
                      doublereal dCoef,
                      sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     SpColVector<T, 3> X(3, 1), V(3, 1), W(3, 0);

     pNode->GetXCurr(X, dCoef, func);
     pNode->GetVCurr(V, dCoef, func);
     pNode->GetWCurr(W, dCoef, func);

     SpColVector<T, 3> s0 = X * dMass + STmp;

     // force
     SpColVector<T, 3> F = pRBK->GetXPP() * -dMass
          - Cross(pRBK->GetWP(), s0)
          - Cross(pRBK->GetW(), Cross(pRBK->GetW(), s0));

     // moment
     SpColVector<T, 3> a = pRBK->GetXPP()
          + Cross(pRBK->GetWP(), X)
          + Cross(pRBK->GetW(), Cross(pRBK->GetW(), X))
          + Cross(pRBK->GetW(), V);

     const SpColVector<T, 3> JTmpWRBK = JTmp * pRBK->GetW();

     SpColVector<T, 3> M = -Cross(STmp, a)
          - Cross(pRBK->GetW(), JTmpWRBK)
          - JTmp * pRBK->GetWP()
          - Cross(W, JTmpWRBK)
          + JTmp * Cross(W, pRBK->GetW())
          - Cross(V, Cross(pRBK->GetW(), STmp));

     const integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

     WorkVec.AddItem(iFirstMomentumIndex + 1, F);
     WorkVec.AddItem(iFirstMomentumIndex + 4, M);
}

void BodyAd::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                           const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     for (integer i = 1; i <= 3; ++i) {
          STmp(i) = S(i);
     }

     for (integer j = 1; j <= 3; ++j) {
          for (integer i = 1; i <= 3; ++i) {
               JTmp(i, j) = J(i, j);
          }
     }
}

DynamicBodyAd::DynamicBodyAd(unsigned int uL,
                             const DynamicStructNodeAd* pN,
                             doublereal dMass,
                             const Vec3& Xgc,
                             const Mat3x3& J,
                             flag fOut)
     :Elem(uL, fOut),
      Body(uL, pN, dMass, Xgc, J, fOut),
      BodyAd(uL, pN, dMass, Xgc, J, fOut),
      DynamicBody(uL, pN, dMass, Xgc, J, fOut),
      pNode(pN)
{
}

DynamicBodyAd::~DynamicBodyAd()
{
}

void
DynamicBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
DynamicBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                      doublereal dCoef,
                      const VectorHandler& XCurr,
                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("DynamicBodyAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::AssJac(this,
                                          WorkMat.SetSparseGradient(),
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
DynamicBodyAd::AssJac(VectorHandler& JacY,
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
DynamicBodyAd::AssRes(SubVectorHandler& WorkVec,
                      doublereal dCoef,
                      const VectorHandler& XCurr ,
                      const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("DynamicBodyAd::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename T>
void
DynamicBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                      doublereal dCoef,
                      const sp_grad::SpGradientVectorHandler<T>& XCurr,
                      const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                      sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 GravityAcceleration;
     bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
                                        GravityAcceleration);

     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();

     SpColVectorA<T, 3> V;
     SpColVectorA<T, 3> W;
     SpMatrixA<T, 3, 3> R;

     pNode->GetVCurr(V, dCoef, func);
     pNode->GetWCurr(W, dCoef, func);
     pNode->GetRCurr(R, dCoef, func);

     SpColVector<T, 3> STmp = R * S0;
     SpMatrix<T, 3, 3> JTmp = R * J0 * Transpose(R);
     SpColVector<T, 3> f1 = V * -dMass - Cross(W, STmp);
     SpColVector<T, 3> f2 = -Cross(STmp, V) - JTmp * W;

     WorkVec.AddItem(iFirstPositionIndex + 1, f1);
     WorkVec.AddItem(iFirstPositionIndex + 4, f2);

     if (g) {
          SpColVector<T, 3> f3 = GravityAcceleration * dMass;
          SpColVector<T, 3> f4 = Cross(STmp, GravityAcceleration);

          WorkVec.AddItem(iFirstPositionIndex + 7, f3);
          WorkVec.AddItem(iFirstPositionIndex + 10, f4);
     }

     if (pRBK) {
          AssVecRBK_int(STmp, JTmp, WorkVec, dCoef, func);
     }

     UpdateInertia(STmp, JTmp);
}

void DynamicBodyAd::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                                  const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     BodyAd::UpdateInertia(S, J);

     pNode->AddInertia(dMass, STmp, JTmp);
}

StaticBodyAd::StaticBodyAd(unsigned int uL,
                           const StaticStructNodeAd* pN,
                           doublereal dMass,
                           const Vec3& Xgc,
                           const Mat3x3& J,
                           flag fOut)
     : Elem(uL, fOut),
       Body(uL, pN, dMass, Xgc, J, fOut),
       BodyAd(uL, pN, dMass, Xgc, J, fOut),
       StaticBody(uL, pN, dMass, Xgc, J, fOut),
       pNode(pN)
{
}


StaticBodyAd::~StaticBodyAd()
{
}

void
StaticBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
StaticBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("StaticBodyAd::AssJac");

     using namespace sp_grad;

     SpGradientAssVec<SpGradient>::AssJac(this,
                                          WorkMat.SetSparseGradient(),
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_JAC);
     return WorkMat;
}

void
StaticBodyAd::AssJac(VectorHandler& JacY,
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
StaticBodyAd::AssRes(SubVectorHandler& WorkVec,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("StaticBodyAd::AssRes");

     using namespace sp_grad;

     SpGradientAssVec<doublereal>::AssRes(this,
                                          WorkVec,
                                          dCoef,
                                          XCurr,
                                          XPrimeCurr,
                                          SpFunctionCall::REGULAR_RES);

     return WorkVec;
}

template <typename T>
void
StaticBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                     doublereal dCoef,
                     const sp_grad::SpGradientVectorHandler<T>& XCurr,
                     const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                     sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 Acceleration(Zero3);
     bool g = GravityOwner::bGetGravity(pNode->GetXCurr(), Acceleration);

     /* W is uninitialized because its use is conditioned by w */
     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     if (!g && !pRBK) {
          return;
     }

     SpMatrixA<T, 3, 3> R;

     pNode->GetRCurr(R, dCoef, func);

     SpColVector<T, 3> STmp = R * S0;
     SpMatrix<T, 3, 3> JTmp = R * J0 * Transpose(R);

     if (g) {
          integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
          SpColVector<T, 3> FTmp = Acceleration * dMass;
          SpColVector<T, 3> MTmp = Cross(STmp, Acceleration);

          WorkVec.AddItem(iFirstMomentumIndex + 1, FTmp);
          WorkVec.AddItem(iFirstMomentumIndex + 4, MTmp);
     }

     if (pRBK) {
          AssVecRBK_int(STmp, JTmp, WorkVec, dCoef, func);
     }

     UpdateInertia(STmp, JTmp);
}

ModalBodyAd::ModalBodyAd(unsigned int uL,
                         const ModalNodeAd* pNode,
                         doublereal dMass,
                         const Vec3& Xgc,
                         const Mat3x3& J,
                         flag fOut)
     : Elem(uL, fOut),
       Body(uL, pNode, dMass, Xgc, J, fOut),
       DynamicBody(uL, pNode, dMass, Xgc, J, fOut),
       BodyAd(uL, pNode, dMass, Xgc, J, fOut),
       ModalBody(uL, pNode, dMass, Xgc, J, fOut),
       pNode(pNode)
{
}

ModalBodyAd::~ModalBodyAd(void)
{
}

void
ModalBodyAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
ModalBodyAd::AssJac(VariableSubMatrixHandler& WorkMat,
                    doublereal dCoef,
                    const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ModalBodyAd::AssJac");

     sp_grad::SpGradientAssVec<sp_grad::SpGradient>::AssJac(this,
                                                            WorkMat.SetSparseGradient(),
                                                            dCoef,
                                                            XCurr,
                                                            XPrimeCurr,
                                                            sp_grad::SpFunctionCall::REGULAR_JAC);

     return WorkMat;
}

void
ModalBodyAd::AssJac(VectorHandler& JacY,
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
ModalBodyAd::AssRes(SubVectorHandler& WorkVec,
                    doublereal dCoef,
                    const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("ModalBodyAd::AssRes");

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
ModalBodyAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const sp_grad::SpGradientVectorHandler<T>& XCurr,
                    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                    sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     Vec3 GravityAcceleration;
     const bool g = GravityOwner::bGetGravity(pNode->GetXCurr(),
                                              GravityAcceleration);

     if (pNode->pGetRBK()) {
          throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
     }

     SpColVector<T, 3> W(3, 3);

     pNode->GetWCurr(W, dCoef, func);

     SpMatrix<T, 3, 3> R(3, 3, 3);

     pNode->GetRCurr(R, dCoef, func);

     SpColVector<T, 3> STmp = R * S0;
     SpMatrix<T, 3, 3> JTmp = R * J0 * Transpose(R);

     SpColVector<T, 3> XPP(3, 1), WP(3, 1);

     pNode->GetXPPCurr(XPP, dCoef, func);
     pNode->GetWPCurr(WP, dCoef, func);

     SpColVector<T, 3> F = XPP * -dMass - Cross(WP, STmp) - Cross(W, Cross(W, STmp));
     SpColVector<T, 3> M = -Cross(STmp, XPP) - Cross(W, JTmp * W) - JTmp * WP;

     if (g) {
          F += GravityAcceleration * dMass;
          M += Cross(STmp, GravityAcceleration);
     }

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();

     WorkVec.AddItem(iFirstPositionIndex + 7, F);
     WorkVec.AddItem(iFirstPositionIndex + 10, M);

     UpdateInertia(STmp, JTmp);
}

void ModalBodyAd::UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                                const sp_grad::SpMatrix<doublereal, 3, 3>& J) const
{
     BodyAd::UpdateInertia(S, J);

     pNode->AddInertia(dMass, STmp, JTmp);
}
