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

#ifndef BODYAD_H
#define BODYAD_H

#include "body.h"
#include "strnodead.h"
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_matvecass.h"

class BodyAd: virtual public Body {
public:
     BodyAd(unsigned int uL, const StructNodeAd *pNode,
            doublereal dMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp,
            flag fOut);

     virtual ~BodyAd();

protected:
     template <typename T>
     void
     AssVecRBK_int(const sp_grad::SpColVector<T, 3>& STmp,
                   const sp_grad::SpMatrix<T, 3, 3>& JTmp,
                   sp_grad::SpGradientAssVec<T>& WorkVec,
                   doublereal dCoef,
                   sp_grad::SpFunctionCall func);
     
     void
     UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& STmp,
                   const sp_grad::SpMatrix<doublereal, 3, 3>& JTmp) const;

     void
     UpdateInertia(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& STmp,
                   const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& JTmp) const {}

     void
     UpdateInertia(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& STmp,
                   const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& JTmp) const {}     
};

class DynamicBodyAd: public DynamicBody, public BodyAd {
public:
     DynamicBodyAd(unsigned int uL, const DynamicStructNodeAd* pNodeTmp, 
                   doublereal dMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp, 
                   flag fOut);

     virtual ~DynamicBodyAd();
 
     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;
 
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
protected:
     using BodyAd::UpdateInertia;
     
     void
     UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& STmp,
                   const sp_grad::SpMatrix<doublereal, 3, 3>& JTmp) const;

private:
     const DynamicStructNodeAd* const pNode;
};

class StaticBodyAd: public StaticBody, public BodyAd {
public:
     StaticBodyAd(unsigned int uL, const StaticStructNodeAd* pNode,
                  doublereal dMass, const Vec3& Xgc, const Mat3x3& J,
                  flag fOut);

     virtual ~StaticBodyAd();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

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
private:
     const StaticStructNodeAd* const pNode;
};

class ModalBodyAd: public ModalBody, public BodyAd {
public:
     ModalBodyAd(unsigned int uL, const ModalNodeAd* pNodeTmp,
                 doublereal dMassTmp, const Vec3& XgcTmp, const Mat3x3& JTmp,
                 flag fOut);

     virtual ~ModalBodyAd();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

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
     
private:
     using BodyAd::UpdateInertia;
     
     void UpdateInertia(const sp_grad::SpColVector<doublereal, 3>& S,
                        const sp_grad::SpMatrix<doublereal, 3, 3>& J) const;
     
     const ModalNodeAd* const pNode;
};

#endif /* BODYAD_H */
