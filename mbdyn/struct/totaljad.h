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

#ifndef TOTALJAD_H
#define TOTALJAD_H

#include "totalj.h"
#include <sp_gradient.h>
#include <sp_matrix_base.h>
#include <sp_matvecass.h>
#include "strnodead.h"

class TotalJointAd: public TotalJoint {
public:
     TotalJointAd(unsigned int uL,
                  const DofOwner *pDO,
                  bool bPos[3],
                  bool bVel[3],
                  TplDriveCaller<Vec3> *const pDCPos[3],
                  bool bRot[3],
                  bool bAgv[3],
                  TplDriveCaller<Vec3> *const pDCRot[3],
                  const StructNodeAd* pN1,
                  const Vec3& f1Tmp,
                  const Mat3x3& R1hTmp,
                  const Mat3x3& R1hrTmp,
                  const StructNodeAd* pN2,
                  const Vec3& f2Tmp,
                  const Mat3x3& R2hTmp,
                  const Mat3x3& R2hrTmp,
                  flag fOut);

     ~TotalJointAd();

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     template <typename T>
     void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            sp_grad::SpFunctionCall func);

     virtual void
     WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     
private:
     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<doublereal, 3>& ThetaDelta);

     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& ThetaDelta);

     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& ThetaDelta);

     inline void
     UpdateF(const sp_grad::SpColVector<doublereal, 3>& FCurr);

     inline void
     UpdateF(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&);

     inline void
     UpdateF(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&);

     inline void
     UpdateM(const sp_grad::SpColVector<doublereal, 3>& MCurr);

     inline void
     UpdateM(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&);

     inline void
     UpdateM(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&);

private:
     const StructNodeAd* pNode1;
     const StructNodeAd* pNode2;
};

class TotalPinJointAd: public TotalPinJoint {
public:
     TotalPinJointAd(unsigned int uL,
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
                     const StructNodeAd* pN,
                     const Vec3& fnTmp,
                     const Mat3x3& RnhTmp,
                     const Mat3x3& RnhrTmp,
                     flag fOut);

     ~TotalPinJointAd();

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     template <typename T>
     void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            sp_grad::SpFunctionCall func);

     virtual void
     WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
     
private:
     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<doublereal, 3>& ThetaDelta);

     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& ThetaDelta);

     inline void
     UpdateThetaDelta(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& ThetaDelta);

     inline void
     UpdateF(const sp_grad::SpColVector<doublereal, 3>& FCurr);

     inline void
     UpdateF(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&);

     inline void
     UpdateF(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&);

     inline void
     UpdateM(const sp_grad::SpColVector<doublereal, 3>& MCurr);

     inline void
     UpdateM(const sp_grad::SpColVector<sp_grad::SpGradient, 3>&);

     inline void
     UpdateM(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>&);

private:
     const StructNodeAd* pNode;
};

#endif // TOTALJAD_H
