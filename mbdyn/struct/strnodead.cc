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
#include "strnodead.h"

StructDispNodeAd::StructDispNodeAd(unsigned int uL,
                                   const DofOwner* pDO,
                                   const Vec3& X0,
                                   const Vec3& V0,
                                   const StructNode *pRN,
                                   const RigidBodyKinematics *pRBK,
                                   doublereal dPosStiff,
                                   doublereal dVelStiff,
                                   OrientationDescription od,
                                   flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut),
 XY(::Zero3)
{
     DEBUGCERR("StructDispNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
}

StructDispNodeAd::~StructDispNodeAd()
{
}

void
StructDispNodeAd::UpdateJac(const VectorHandler& Y, doublereal dCoef)
{
     integer iFirstDofIndex = iGetFirstIndex();

     for (integer i = 1; i <= 3; ++i) {
          XY(i) = Y(iFirstDofIndex + i);
     }
}

DynamicStructDispNodeAd::DynamicStructDispNodeAd(unsigned int uL,
                                                 const DofOwner* pDO,
                                                 const Vec3& X0,
                                                 const Vec3& V0,
                                                 const StructNode *pRN,
                                                 const RigidBodyKinematics *pRBK,
                                                 doublereal dPosStiff,
                                                 doublereal dVelStiff,
                                                 OrientationDescription od,
                                                 flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut),
 DynamicStructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut),
 StructDispNodeAd(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut)
{
     DEBUGCERR("DynamicStructDispNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
}

DynamicStructDispNodeAd::~DynamicStructDispNodeAd()
{
}

inline integer
DynamicStructDispNodeAd::iGetInitialFirstIndexPrime() const
{
     // FIXME: Is it correct this way?
     return iGetFirstIndex() + 3;
}


StaticStructDispNodeAd::StaticStructDispNodeAd(unsigned int uL,
                                               const DofOwner* pDO,
                                               const Vec3& X0,
                                               const Vec3& V0,
                                               const StructNode *pRN,
                                               const RigidBodyKinematics *pRBK,
                                               doublereal dPosStiff,
                                               doublereal dVelStiff,
                                               OrientationDescription od,
                                               flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut),
 StaticStructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut),
 StructDispNodeAd(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, od, fOut)
{
     DEBUGCERR("StaticStructDispNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
}

StaticStructDispNodeAd::~StaticStructDispNodeAd()
{
}

inline integer
StaticStructDispNodeAd::iGetInitialFirstIndexPrime() const
{
     // FIXME: Is it correct this way?
     return iGetFirstIndex() + 3;
}

StructNodeAd::StructNodeAd(unsigned int uL,
                           const DofOwner* pDO,
                           const Vec3& X0,
                           const Mat3x3& R0,
                           const Vec3& V0,
                           const Vec3& W0,
                           const StructNode *pRN,
                           const RigidBodyKinematics *pRBK,
                           doublereal dPosStiff,
                           doublereal dVelStiff,
                           bool bOmRot,
                           OrientationDescription ood,
                           flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructDispNodeAd(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 eCurrFunc(sp_grad::UNKNOWN_FUNC),
 gY(::Zero3),
 bNeedRotation(false),
 bUpdateRotation(true),
 bUpdateRotationGradProd(true)
{
     DEBUGCERR("StructNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
     DEBUGCERR("R=" << GetRCurr() << "\n");
     DEBUGCERR("W=" << GetWCurr() << "\n");
}

StructNodeAd::~StructNodeAd()
{
}

inline void StructNodeAd::GetWCurrInitAss(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     W = WCurr;
}

inline void StructNodeAd::GetWCurrInitAss(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
          iFirstDofIndex = iGetFirstIndex() + 9;
          break;

     default:
          SP_GRAD_ASSERT(false);
     }

     W.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          W(i).Reset(WCurr(i), iFirstDofIndex + i, -1.);
     }
}

inline void StructNodeAd::GetWCurrInitAss(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     // FIXME: Should be eliminated
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

template <typename T>
void StructNodeAd::UpdateRotation(const Mat3x3& RRef, const Vec3& WRef, const sp_grad::SpColVector<T, 3>& g, const sp_grad::SpColVector<T, 3>& gP, sp_grad::SpMatrix<T, 3, 3>& R, sp_grad::SpColVector<T, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);

     using namespace sp_grad;

     const T d = 4. / (4. + Dot(g, g));

     SpMatrix<T, 3, 3> RDelta(3, 3, 8);

     const T tmp1 = -g(3) * g(3);
     const T tmp2 = -g(2) * g(2);
     const T tmp3 = -g(1) * g(1);
     const T tmp4 = g(1) * g(2) * 0.5;
     const T tmp5 = g(2) * g(3) * 0.5;
     const T tmp6 = g(1) * g(3) * 0.5;

     RDelta(1,1) = (tmp1 + tmp2) * d * 0.5 + 1;
     RDelta(1,2) = (tmp4 - g(3)) * d;
     RDelta(1,3) = (tmp6 + g(2)) * d;
     RDelta(2,1) = (g(3) + tmp4) * d;
     RDelta(2,2) = (tmp1 + tmp3) * d * 0.5 + 1.;
     RDelta(2,3) = (tmp5 - g(1)) * d;
     RDelta(3,1) = (tmp6 - g(2)) * d;
     RDelta(3,2) = (tmp5 + g(1)) * d;
     RDelta(3,3) = (tmp2 + tmp3) * d * 0.5 + 1.;

     R = EvalUnique(RDelta * RRef);

     switch (func) {
     case SpFunctionCall::INITIAL_ASS_JAC:
          GetWCurrInitAss(W, dCoef, func);
          break;

     case SpFunctionCall::INITIAL_DER_JAC:
     case SpFunctionCall::REGULAR_JAC:
     {
          SpMatrix<T, 3, 3> G(3, 3, 8);

          const T tmp7 = 0.5 * g(1) * d;
          const T tmp8 = 0.5 * g(2) * d;
          const T tmp9 = 0.5 * g(3) * d;

          G(1,1) = d;
          G(1,2) = -tmp9;
          G(1,3) = tmp8;
          G(2,1) = tmp9;
          G(2,2) = d;
          G(2,3) = -tmp7;
          G(3,1) = -tmp8;
          G(3,2) = tmp7;
          G(3,3) = d;

          W = EvalUnique(G * gP + RDelta * WRef); // Note that the first index of gP and g must be the same in order to work!
     }
     break;

     default:
          SP_GRAD_ASSERT(false);
     }
}

void StructNodeAd::UpdateRotation(doublereal dCoef) const
{
     SP_GRAD_ASSERT(bNeedRotation);

     if (bUpdateRotation) {
          sp_grad::SpColVectorA<sp_grad::SpGradient, 3, 1> gCurr_grad, gPCurr_grad;

          GetgCurr(gCurr_grad, dCoef, eCurrFunc);
          GetgPCurr(gPCurr_grad, dCoef, eCurrFunc);

          UpdateRotation(RRef, WRef, gCurr_grad, gPCurr_grad, RCurr_grad, WCurr_grad, dCoef, eCurrFunc);

#if defined(DEBUG)
          {
               using sp_grad::index_type;

               const double dTol = sqrt(std::numeric_limits<doublereal>::epsilon());

               Mat3x3 RDelta(CGR_Rot::MatR, gCurr);

               Mat3x3 RCurr_tmp = RDelta * RRef;

               bool bErr = false;

               for (index_type i = 1; i <= 3; ++i) {
                    for (index_type j = 1; j <= 3; ++j) {
                         if (std::abs(RCurr_grad(i, j).dGetValue() - RCurr_tmp(i, j)) > dTol) {
                              bErr = true;
                         }
                    }
               }

               for (index_type i = 1; i <= 3; ++i) {
                    for (index_type j = 1; j <= 3; ++j) {
                         if (std::abs(RCurr_grad(i, j).dGetValue() - RCurr(i, j)) > dTol) {
                              bErr = true;
                         }
                    }

                    if (std::abs(WCurr_grad(i).dGetValue() - WCurr(i)) > dTol) {
                         bErr = true;
                    }
               }

               if (bErr) {
                    std::cerr << "gCurr=" << gCurr << std::endl;
                    std::cerr << "RCurr=" << std::endl;
                    for (integer i = 1; i <= 3; ++i) {
                         for (integer j = 1; j <= 3; ++j) {
                              std::cerr << RCurr(i, j) << " ";
                         }
                         std::cerr << std::endl;
                    }

                    std::cerr << "RRef=" << std::endl;
                    for (integer i = 1; i <= 3; ++i) {
                         for (integer j = 1; j <= 3; ++j) {
                              std::cerr << RRef(i, j) << " ";
                         }
                         std::cerr << std::endl;
                    }

                    std::cerr << "RPrev=" << std::endl;
                    for (integer i = 1; i <= 3; ++i) {
                         for (integer j = 1; j <= 3; ++j) {
                              std::cerr << RPrev(i, j) << " ";
                         }
                         std::cerr << std::endl;
                    }

                    std::cerr << "RCurr_grad=" << std::endl;

                    for (integer i = 1; i <= 3; ++i) {
                         for (integer j = 1; j <= 3; ++j) {
                              std::cerr << RCurr_grad(i, j).dGetValue() << " ";
                         }
                         std::cerr << std::endl;
                    }

                    std::cerr << "WCurr=" << WCurr << std::endl;
                    std::cerr << "WCurr_grad=";
                    for (integer i = 1; i <= 3; ++i) {
                         std::cerr << WCurr_grad(i).dGetValue() << " ";
                    }
                    std::cerr << std::endl;
                    SP_GRAD_ASSERT(false);
               }
          }
#endif
          bUpdateRotation = false;
     }
}

void StructNodeAd::UpdateRotation(const VectorHandler& Y, doublereal dCoef) const
{
     SP_GRAD_ASSERT(bNeedRotation);

     if (bUpdateRotationGradProd) {
          sp_grad::SpColVectorA<sp_grad::GpGradProd, 3> gCurr_gradp, gPCurr_gradp;

          GetgCurr(gCurr_gradp, dCoef, eCurrFunc);
          GetgPCurr(gPCurr_gradp, dCoef, eCurrFunc);

          UpdateRotation(RRef, WRef, gCurr_gradp, gPCurr_gradp, RCurr_gradp, WCurr_gradp, dCoef, eCurrFunc);
          bUpdateRotationGradProd = false;
     }
}

void
StructNodeAd::Update(const VectorHandler& X, const VectorHandler& XP)
{
     InvalidateGradients();

     StructNode::Update(X, XP);
}

void
StructNodeAd::DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP)
{
     InvalidateGradients();

     StructNode::DerivativesUpdate(X, XP);
}

void
StructNodeAd::InitialUpdate(const VectorHandler& X)
{
     InvalidateGradients();

     StructNode::InitialUpdate(X);
}

void StructNodeAd::UpdateJac(doublereal dCoef)
{
     if (bNeedRotation) {
          UpdateRotation(dCoef);
     }
}

void StructNodeAd::UpdateJac(const VectorHandler& Y, doublereal dCoef)
{
     StructDispNodeAd::UpdateJac(Y, dCoef);

     UpdateJacRotation(Y, dCoef);
}

void StructNodeAd::UpdateJacRotation(const VectorHandler& Y, doublereal dCoef)
{
     bUpdateRotationGradProd = true;

     if (bNeedRotation) {
          integer iFirstIndex = iGetFirstIndex();

          for (integer i = 1; i <= 3; ++i) {
               gY(i) = Y(iFirstIndex + i + 3);
          }

          UpdateRotation(Y, dCoef);
     }
}

void
StructNodeAd::SetInitialValue(VectorHandler& X)
{
     StructNode::SetInitialValue(X);

     eCurrFunc = sp_grad::INITIAL_ASS_JAC;
}

void
StructNodeAd::SetValue(DataManager *pDM,
                       VectorHandler& X, VectorHandler& XP,
                       SimulationEntity::Hints *ph)
{
     InvalidateGradients();

     eCurrFunc = sp_grad::REGULAR_JAC;

     StructNode::SetValue(pDM, X, XP, ph);
}

void
StructNodeAd::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
     InvalidateGradients();

     StructNode::AfterPredict(X, XP);
}

void
StructNodeAd::AfterConvergence(const VectorHandler& X,
                               const VectorHandler& XP,
                               const VectorHandler& XPP)
{
     InvalidateGradients();
}

void StructNodeAd::InvalidateGradients() const
{
     bUpdateRotationGradProd = bUpdateRotation = true;
}

DynamicStructNodeAd::DynamicStructNodeAd(unsigned int uL,
                                         const DofOwner* pDO,
                                         const Vec3& X0,
                                         const Mat3x3& R0,
                                         const Vec3& V0,
                                         const Vec3& W0,
                                         const StructNode *pRN,
                                         const RigidBodyKinematics *pRBK,
                                         doublereal dPosStiff,
                                         doublereal dVelStiff,
                                         bool bOmRot,
                                         OrientationDescription ood,
                                         flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 DynamicStructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructDispNodeAd(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut), 
 DynamicStructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructNodeAd(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut)
{
     DEBUGCERR("DynamicStructNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
     DEBUGCERR("R=" << GetRCurr() << "\n");
     DEBUGCERR("W=" << GetWCurr() << "\n");
}

DynamicStructNodeAd::~DynamicStructNodeAd()
{
}

void DynamicStructNodeAd::Update(const VectorHandler& X, const VectorHandler& XP)
{
     InvalidateGradients();

     DynamicStructNode::Update(X, XP);
}

integer DynamicStructNodeAd::iGetInitialFirstIndexPrime() const
{
     return iGetFirstIndex() + 6;
}

StaticStructNodeAd::StaticStructNodeAd(unsigned int uL,
                                       const DofOwner* pDO,
                                       const Vec3& X0,
                                       const Mat3x3& R0,
                                       const Vec3& V0,
                                       const Vec3& W0,
                                       const StructNode *pRN,
                                       const RigidBodyKinematics *pRBK,
                                       doublereal dPosStiff,
                                       doublereal dVelStiff,
                                       bool bOmRot,
                                       OrientationDescription ood,
                                       flag fOut)
:StructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StaticStructDispNode(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut), 
 StructDispNodeAd(uL, pDO, X0, V0, pRN, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StaticStructNode(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructNodeAd(uL, pDO, X0, R0, V0, W0, pRN, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut)
{
     DEBUGCERR("StaticStructNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
     DEBUGCERR("R=" << GetRCurr() << "\n");
     DEBUGCERR("W=" << GetWCurr() << "\n");
}

StaticStructNodeAd::~StaticStructNodeAd()
{
}

integer StaticStructNodeAd::iGetInitialFirstIndexPrime() const
{
     return iGetFirstIndex() + 6;
}

void StaticStructNodeAd::Update(const VectorHandler& X, const VectorHandler& XP)
{
     InvalidateGradients();

     StaticStructNode::Update(X, XP);
}

ModalNodeAd::ModalNodeAd(unsigned int uL,
                         const DofOwner* pDO,
                         const Vec3& X0,
                         const Mat3x3& R0,
                         const Vec3& V0,
                         const Vec3& W0,
                         const RigidBodyKinematics *pRBK,
                         doublereal dPosStiff,
                         doublereal dVelStiff,
                         bool bOmRot,
                         OrientationDescription ood,
                         flag fOut)
:StructDispNode(uL, pDO, X0, V0, 0, pRBK, dPosStiff, dVelStiff, ood, fOut),
 DynamicStructDispNode(uL, pDO, X0, V0, 0, pRBK, dPosStiff, dVelStiff, ood, fOut),
 StructNode(uL, pDO, X0, R0, V0, W0, 0, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructDispNodeAd(uL, pDO, X0, V0, 0, pRBK, dPosStiff, dVelStiff, ood, fOut), 
 ModalNode(uL, pDO, X0, R0, V0, W0, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 StructNodeAd(uL, pDO, X0, R0, V0, W0, 0, pRBK, dPosStiff, dVelStiff, bOmRot, ood, fOut),
 XPPY(::Zero3), WPY(::Zero3)
{
     DEBUGCERR("StaticStructNodeAd(" << GetLabel() << "\n");
     DEBUGCERR("X=" << GetXCurr() << "\n");
     DEBUGCERR("V=" << GetVCurr() << "\n");
     DEBUGCERR("R=" << GetRCurr() << "\n");
     DEBUGCERR("W=" << GetWCurr() << "\n");
}

ModalNodeAd::~ModalNodeAd()
{
}

integer
ModalNodeAd::iGetInitialFirstIndexPrime() const
{
        // FIXME: Don't know how it should be implemented!
        throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

void
ModalNodeAd::Update(const VectorHandler& X, const VectorHandler& XP)
{
     InvalidateGradients();

     ModalNode::Update(X, XP);
}

void
ModalNodeAd::DerivativesUpdate(const VectorHandler& X,
                               const VectorHandler& XP)
{
     InvalidateGradients();

     ModalNode::DerivativesUpdate(X, XP);
}

void
ModalNodeAd::AfterConvergence(const VectorHandler& X,
                              const VectorHandler& XP)
{
     // Must override DynamicStructNode's function!
}

void ModalNodeAd::UpdateJac(const VectorHandler& Y, doublereal dCoef)
{
        StructNodeAd::UpdateJac(Y, dCoef);

        const integer iFirstIndex = iGetFirstIndex();

        for (sp_grad::index_type i = 1; i <= 3; ++i) {
                XPPY(i) = Y(iFirstIndex + i + 6);
                WPY(i) = Y(iFirstIndex + i + 9);
        }
}
