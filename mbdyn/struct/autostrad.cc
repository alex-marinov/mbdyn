/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "autostrad.h"

AutomaticStructDispElemAd::AutomaticStructDispElemAd(const DynamicStructDispNodeAd* pN)
     :Elem(pN->GetLabel(), pN->fToBeOutput()),
      AutomaticStructDispElem(pN),
      pNode(pN)
{
}

AutomaticStructDispElemAd::~AutomaticStructDispElemAd()
{
}

void AutomaticStructDispElemAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 6;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
AutomaticStructDispElemAd::AssJac(VariableSubMatrixHandler& WorkMat,
                                  doublereal dCoef,
                                  const VectorHandler& XCurr,
                                  const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticStructElemAd::AssJac");

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
AutomaticStructDispElemAd::AssJac(VectorHandler& JacY,
                                  const VectorHandler& Y,
                                  doublereal dCoef,
                                  const VectorHandler& XCurr,
                                  const VectorHandler& XPrimeCurr,
                                  VariableSubMatrixHandler& WorkMat)
{
     DEBUGCOUTFNAME("AutomaticStructElemAd::AssJac");

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
AutomaticStructDispElemAd::AssRes(SubVectorHandler& WorkVec,
                                  doublereal dCoef,
                                  const VectorHandler& XCurr,
                                  const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticStructElemAd::AssRes");

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
AutomaticStructDispElemAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                                  doublereal dCoef,
                                  const sp_grad::SpGradientVectorHandler<T>& XCurr,
                                  const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                                  sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
     const integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

     SpColVector<T, 3> B(3, 1), BP(3, 1), V(3, 1);

     pNode->GetVCurr(V, dCoef, func);

     XCurr.GetVec(iFirstMomentumIndex + 1, B, dCoef);
     XPrimeCurr.GetVec(iFirstMomentumIndex + 1, BP, 1.);

     WorkVec.AddItem(iFirstPositionIndex + 1, B);

     SpColVector<T, 3> F = -BP;

     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     if (pRBK) {
          const Vec3& W0 = pRBK->GetW();

          F -= 2. * Cross(W0, B);
     }

     WorkVec.AddItem(iFirstPositionIndex + 4, F);

     UpdateState(B, BP);
}

inline void
AutomaticStructDispElemAd::UpdateState(const sp_grad::SpColVector<doublereal, 3>& BTmp,
                                       const sp_grad::SpColVector<doublereal, 3>& BPTmp)
{
     using namespace sp_grad;

     for (integer i = 1; i <= 3; ++i) {
          B(i) = BTmp(i);
          BP(i) = BPTmp(i);
     }

     // reset instantaneous inertia properties
     m = 0.;
}

AutomaticStructElemAd::AutomaticStructElemAd(const DynamicStructNodeAd* pN)
     :Elem(pN->GetLabel(), pN->fToBeOutput()),
      AutomaticStructElem(pN),
      pNode(const_cast<DynamicStructNodeAd*>(pN))
{
     pNode->SetAutoStr(this);
}

AutomaticStructElemAd::~AutomaticStructElemAd()
{
}

void AutomaticStructElemAd::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
     *piNumRows = 12;
     *piNumCols = 0;
}

VariableSubMatrixHandler&
AutomaticStructElemAd::AssJac(VariableSubMatrixHandler& WorkMat,
                              doublereal dCoef,
                              const VectorHandler& XCurr,
                              const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticStructElemAd::AssJac");

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
AutomaticStructElemAd::AssJac(VectorHandler& JacY,
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
AutomaticStructElemAd::AssRes(SubVectorHandler& WorkVec,
                              doublereal dCoef,
                              const VectorHandler& XCurr,
                              const VectorHandler& XPrimeCurr)
{
     DEBUGCOUTFNAME("AutomaticStructElemAd::AssRes");

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
AutomaticStructElemAd::AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
                              doublereal dCoef,
                              const sp_grad::SpGradientVectorHandler<T>& XCurr,
                              const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
                              sp_grad::SpFunctionCall func)
{
     using namespace sp_grad;

     const integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
     const integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();

     SpColVector<T, 3> B(3, 1), G(3, 1), BP(3, 1), GP(3, 1), V(3, 1);

     pNode->GetVCurr(V, dCoef, func);

     XCurr.GetVec(iFirstMomentumIndex + 1, B, dCoef);
     XCurr.GetVec(iFirstMomentumIndex + 4, G, dCoef);
     XPrimeCurr.GetVec(iFirstMomentumIndex + 1, BP, 1.);
     XPrimeCurr.GetVec(iFirstMomentumIndex + 4, GP, 1.);

     /*
      * Momentum and momenta moment (about node):
      *
      * B = m V + W /\ S
      *
      * G = S /\ V + J W
      *
      * Bp = F
      *
      * Gp + V /\ B = M
      */
     WorkVec.AddItem(iFirstPositionIndex + 1, B);
     WorkVec.AddItem(iFirstPositionIndex + 4, G);

     SpColVector<T, 3> F = -BP;
     SpColVector<T, 3> M = -Cross(V, B) - GP;

     // relative frame dynamics contribution
     // (see tecman, "Dynamics in a Relative Reference Frame")
     const RigidBodyKinematics *pRBK = pNode->pGetRBK();

     if (pRBK) {
          const Vec3& W0 = pRBK->GetW();

          F -= 2. * Cross(W0, B);
          M -= Cross(W0, G);
     }

     WorkVec.AddItem(iFirstPositionIndex + 7, F);
     WorkVec.AddItem(iFirstPositionIndex + 10, M);

     UpdateState(B, BP, G, GP);
}

inline void
AutomaticStructElemAd::UpdateState(const sp_grad::SpColVector<doublereal, 3>& BTmp,
                                   const sp_grad::SpColVector<doublereal, 3>& BPTmp,
                                   const sp_grad::SpColVector<doublereal, 3>& GTmp,
                                   const sp_grad::SpColVector<doublereal, 3>& GPTmp)
{
     using namespace sp_grad;

     for (integer i = 1; i <= 3; ++i) {
          B(i) = BTmp(i);
          BP(i) = BPTmp(i);
          G(i) = GTmp(i);
          GP(i) = GPTmp(i);
     }

     // reset instantaneous inertia properties
     m = 0.;
     S = Zero3;
     J = Zero3x3;
}
