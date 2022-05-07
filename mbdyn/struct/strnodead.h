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

#ifndef STRNODEAD_H
#define STRNODEAD_H

#include "strnode.h"
#include "sp_gradient.h"
#include "sp_gradient_op.h"
#include "sp_matrix_base.h"

class StructDispNodeAd: virtual public StructDispNode {
protected:
     StructDispNodeAd(unsigned int uL,
                      const DofOwner* pDO,
                      const Vec3& X0,
                      const Vec3& V0,
                      const StructNode *pRN,
                      const RigidBodyKinematics *pRBK,
                      doublereal dPosStiff,
                      doublereal dVelStiff,
                      OrientationDescription od,
                      flag fOut);
     virtual ~StructDispNodeAd();

public:
     using StructDispNode::GetXCurr;
     using StructDispNode::GetVCurr;
     
     inline void
     GetXCurr(sp_grad::SpColVector<doublereal, 3>& X,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     inline void
     GetXCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& X,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     inline void
     GetVCurr(sp_grad::SpColVector<doublereal, 3>& V,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     inline void
     GetVCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& V,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     inline void
     GetXCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& X,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     inline void
     GetVCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& V,
              doublereal dCoef,
              sp_grad::SpFunctionCall func) const;

     virtual void
     UpdateJac(const VectorHandler& Y, doublereal dCoef) override;

protected:
     // Returns the dof index -1 of VCurr during initial assembly
     virtual integer iGetInitialFirstIndexPrime() const=0;

private:
     Vec3 XY;
};

class DynamicStructDispNodeAd: virtual public DynamicStructDispNode, public StructDispNodeAd {
public:
     DynamicStructDispNodeAd(unsigned int uL,
                             const DofOwner* pDO,
                             const Vec3& X0,
                             const Vec3& V0,
                             const StructNode *pRN,
                             const RigidBodyKinematics *pRBK,
                             doublereal dPosStiff,
                             doublereal dVelStiff,
                             OrientationDescription od,
                             flag fOut);

     virtual ~DynamicStructDispNodeAd();

protected:
     virtual inline integer iGetInitialFirstIndexPrime() const override;
};

class StaticStructDispNodeAd: virtual public StaticStructDispNode, public StructDispNodeAd {
public:
     StaticStructDispNodeAd(unsigned int uL,
                            const DofOwner* pDO,
                            const Vec3& X0,
                            const Vec3& V0,
                            const StructNode *pRN,
                            const RigidBodyKinematics *pRBK,
                            doublereal dPosStiff,
                            doublereal dVelStiff,
                            OrientationDescription od,
                            flag fOut);

     virtual ~StaticStructDispNodeAd();

protected:
     virtual inline integer iGetInitialFirstIndexPrime() const override;
};

class StructNodeAd: virtual public StructNode, virtual public StructDispNodeAd {
protected:
     StructNodeAd(unsigned int uL,
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
                  flag fOut);

     virtual ~StructNodeAd();

public:
     using StructNode::GetgCurr;
     using StructNode::GetgPCurr;
     using StructNode::GetRCurr;
     using StructNode::GetWCurr;
     
     inline void GetgCurr(sp_grad::SpColVector<doublereal, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetgCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetgPCurr(sp_grad::SpColVector<doublereal, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetgPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetRCurr(sp_grad::SpMatrix<doublereal, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetRCurr(sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetWCurr(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetWCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetgCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetgPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetRCurr(sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void GetWCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     virtual void Update(const VectorHandler& X, const VectorHandler& XP) override;
     virtual void InitialUpdate(const VectorHandler& X) override;
     virtual void DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP) override;
     virtual void UpdateJac(doublereal dCoef) override;
     virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;
     virtual void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP, SimulationEntity::Hints *ph) override;
     virtual void SetInitialValue(VectorHandler& X) override;
     virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) override;
     virtual void AfterConvergence(const VectorHandler& X, 
                                   const VectorHandler& XP, 
                                   const VectorHandler& XPP) override;

protected:
     inline void InvalidateGradients() const;
     inline void UpdateJacRotation(const VectorHandler& Y, doublereal dCoef);
     
private:
     void UpdateRotation(doublereal dCoef) const;
     void UpdateRotation(const VectorHandler& Y, doublereal dCoef) const;
     
     template <typename T>
     inline void UpdateRotation(const Mat3x3& RRef, const Vec3& WRef, const sp_grad::SpColVector<T, 3>& g, const sp_grad::SpColVector<T, 3>& gP, sp_grad::SpMatrix<T, 3, 3>& RCurr, sp_grad::SpColVector<T, 3>& WCurr, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetWCurrInitAss(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetWCurrInitAss(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetWCurrInitAss(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     sp_grad::SpFunctionCall eCurrFunc;
     mutable sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3, 3> RCurr_grad;
     mutable sp_grad::SpColVectorA<sp_grad::SpGradient, 3, 3> WCurr_grad;
     mutable sp_grad::SpMatrixA<sp_grad::GpGradProd, 3, 3> RCurr_gradp;
     mutable sp_grad::SpColVectorA<sp_grad::GpGradProd, 3> WCurr_gradp;
     mutable Vec3 gY;
     mutable bool bNeedRotation, bUpdateRotation, bUpdateRotationGradProd;
};

class DynamicStructNodeAd: public DynamicStructNode, public StructNodeAd
{
public:
     DynamicStructNodeAd(unsigned int uL,
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
                         flag fOut);

     virtual ~DynamicStructNodeAd();
     virtual void Update(const VectorHandler&, const VectorHandler&) override;
     
protected:
     virtual integer iGetInitialFirstIndexPrime() const override;
};

class StaticStructNodeAd: public StaticStructNode, public StructNodeAd
{
public:
     StaticStructNodeAd(unsigned int uL,
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
                        flag fOut);

     virtual ~StaticStructNodeAd();
     virtual void Update(const VectorHandler&, const VectorHandler&) override;
     
protected:
     virtual integer iGetInitialFirstIndexPrime() const override;
};

class ModalNodeAd: public ModalNode, public StructNodeAd {
public:
     ModalNodeAd(unsigned int uL,
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
                 flag fOut);

     virtual ~ModalNodeAd();

     virtual void Update(const VectorHandler& X, const VectorHandler& XP) override;

     virtual void DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP) override;
     
     virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

     using StructNode::GetXPPCurr;
     using StructNode::GetWPCurr;

     inline void
     GetXPPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;
                
     inline void
     GetXPPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void
     GetXPPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void
     GetWPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;
                
     inline void
     GetWPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     inline void
     GetWPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

     virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;

protected:
     virtual inline integer iGetInitialFirstIndexPrime() const override;
     
private:
     Vec3 XPPY, WPY;
};

inline void StructDispNodeAd::GetXCurr(sp_grad::SpColVector<doublereal, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     X = StructDispNode::GetXCurr();
}

inline void StructDispNodeAd::GetXCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
          SP_GRAD_ASSERT(dCoef == 1.);
     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
          iFirstDofIndex = StructDispNode::iGetFirstIndex();
          break;

     default:
          SP_GRAD_ASSERT(false);
     }

     X.ResizeReset(3, 1);

     const Vec3& XCurr = StructDispNode::GetXCurr();

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          X(i).Reset(XCurr(i), iFirstDofIndex + i, -dCoef);
     }
}

inline void StructDispNodeAd::GetXCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     X.ResizeReset(3, 1);

     const Vec3& XCurr = StructDispNode::GetXCurr();

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          X(i).Reset(XCurr(i), -dCoef * XY(i));
     }
}

inline void StructDispNodeAd::GetVCurr(sp_grad::SpColVector<doublereal, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     V = StructDispNode::GetVCurr();
}

inline void StructDispNodeAd::GetVCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
          SP_GRAD_ASSERT(dCoef == 1.);
          iFirstDofIndex = StructDispNode::iGetInitialFirstIndexPrime();
          break;

     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
          iFirstDofIndex = StructDispNode::iGetFirstIndex();
          break;

     default:
          SP_GRAD_ASSERT(false);
     }

     V.ResizeReset(3, 1);

     const Vec3& VCurr = StructDispNode::GetVCurr();

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          V(i).Reset(VCurr(i), iFirstDofIndex + i, -1.);
     }
}

inline void StructDispNodeAd::GetVCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     V.ResizeReset(3, 1);

     const Vec3& VCurr = StructDispNode::GetVCurr();

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          V(i).Reset(VCurr(i), -XY(i));
     }
}

inline void StructNodeAd::GetgCurr(sp_grad::SpColVector<doublereal, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     g = gCurr;
}
     
inline void StructNodeAd::GetgCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  SP_GRAD_ASSERT(dCoef == 1.);

     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex();
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }

     g.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  g(i).Reset(gCurr(i), iFirstDofIndex + i + 3, -dCoef);
     }
}

inline void StructNodeAd::GetgPCurr(sp_grad::SpColVector<doublereal, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     gP = gPCurr;
}
                
inline void StructNodeAd::GetgPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  gP = gPCurr;
	  return;
	  
     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex() + 3;
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }

     gP.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  gP(i).Reset(gPCurr(i), iFirstDofIndex + i, -1.);
     }
}


inline void StructNodeAd::GetgCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     g.ResizeReset(3, 1);
     
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          g(i).Reset(gCurr(i), -dCoef * gY(i));
     }
}

inline void StructNodeAd::GetgPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     gP.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          gP(i).Reset(gPCurr(i), -gY(i));
     }
}

inline void StructNodeAd::GetRCurr(sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotationGradProd);
     
     R = RCurr_gradp;
}

inline void StructNodeAd::GetWCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotationGradProd);
     
     W = WCurr_gradp;
}

inline void StructNodeAd::GetRCurr(sp_grad::SpMatrix<doublereal, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     bNeedRotation = true;
     R = RCurr;
}

inline void StructNodeAd::GetRCurr(sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotation);
     
     R = RCurr_grad;
}

inline void StructNodeAd::GetWCurr(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     bNeedRotation = true;
     W = WCurr;
}

inline void StructNodeAd::GetWCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotation);

     W = WCurr_grad;
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
     // FIXME: Not needed and should be eliminated
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

inline integer
ModalNodeAd::iGetInitialFirstIndexPrime() const
{
	// FIXME: Don't know how it should be implemented!
	throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}

inline void
ModalNodeAd::GetXPPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	XPP = XPPCurr;
}
                
inline void
ModalNodeAd::GetXPPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	sp_grad::index_type iFirstDofIndex = -1;

	switch (func) {
	case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
	case sp_grad::SpFunctionCall::REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		SP_GRAD_ASSERT(false);
	}

	XPP.ResizeReset(3, 1);
	
	for (sp_grad::index_type i = 1; i <= 3; ++i) {
	     XPP(i).Reset(XPPCurr(i), iFirstDofIndex + i + 6, -1.);
	}		
}

inline void
ModalNodeAd::GetXPPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          XPP(i).Reset(XPPCurr(i), -XPPY(i));
     }
}
                
inline void
ModalNodeAd::GetWPCurr(sp_grad::SpColVector<doublereal, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	WP = WPCurr;
}
                
inline void
ModalNodeAd::GetWPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	sp_grad::index_type iFirstDofIndex = -1;

	switch (func) {
	case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
	case sp_grad::SpFunctionCall::REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		SP_GRAD_ASSERT(false);
	}

	WP.ResizeReset(3, 1);
	
	for (sp_grad::index_type i = 1; i <= 3; ++i) {
	     WP(i).Reset(WPCurr(i), iFirstDofIndex + i + 9, -1.);
	}
}

inline void
ModalNodeAd::GetWPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          WP(i).Reset(WPCurr(i), -WPY(i));
     }
}

#endif
