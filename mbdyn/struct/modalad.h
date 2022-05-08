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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef MODALAD_H
#define MODALAD_H

#include "modal.h"
#include "strnodead.h"
#include <sp_gradient.h>
#include <sp_matrix_base.h>
#include <sp_matvecass.h>

class ModalAd: public Modal {
public:
     ModalAd(unsigned int uL,
             const ModalNodeAd* pModalNodeTmp,
             const Vec3& x0,
             const Mat3x3& R0,
             const DofOwner* pDO,
             unsigned int N,
             unsigned int NS,
             unsigned int NFN,
             doublereal dMass,
             const Vec3& STmp,
             const Mat3x3& JTmp,
             const std::vector<unsigned int>& uModeNumber,
             MatNxN *pGenMass,
             MatNxN *pGenStiff,
             MatNxN *pGenDamp,
             const std::vector<std::string>& IdFEMNodes,
             Mat3xN *pN,
             const std::vector<Modal::StrNodeData>& snd,
             Mat3xN *pPHIt,
             Mat3xN *pPHIr,
             Mat3xN *pModeShapest,
             Mat3xN *pModeShapesr,
             Mat3xN *pInv3,
             Mat3xN *pInv4,
             Mat3xN *pInv5,
             Mat3xN *pInv8,
             Mat3xN *pInv9,
             Mat3xN *pInv10,
             Mat3xN *pInv11,
             VecN *a,
             VecN *aP,
             flag fOut);

     virtual ~ModalAd();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
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
     AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     inline void GetACurr(sp_grad::index_type iMode, sp_grad::SpGradient& ai, doublereal dCoef, sp_grad::SpFunctionCall func) const {
          using namespace sp_grad;

          SP_GRAD_ASSERT(iMode >= 1);
          SP_GRAD_ASSERT(iMode <= a.iGetNumRows());
          SP_GRAD_ASSERT(func != SpFunctionCall::INITIAL_ASS_JAC || dCoef == 1.);

          const index_type iFirstIndex = iGetFirstIndex();

          ai.Reset(a(iMode), iFirstIndex + iMode, -dCoef);
     }

     inline void GetACurr(sp_grad::index_type iMode, doublereal& ai, doublereal, sp_grad::SpFunctionCall) const {
          SP_GRAD_ASSERT(iMode >= 1);
          SP_GRAD_ASSERT(iMode <= a.iGetNumRows());
          ai = a(iMode);
     }

     template <typename T>
     void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            sp_grad::SpFunctionCall func);

private:
     inline void
     UpdateStrNodeData(StrNodeData& oNode,
                       const sp_grad::SpColVector<doublereal, 3>& d1tot,
                       const sp_grad::SpMatrix<doublereal, 3, 3>& R1tot,
                       const sp_grad::SpColVector<doublereal, 3>& F,
                       const sp_grad::SpColVector<doublereal, 3>& M,
                       const sp_grad::SpMatrix<doublereal, 3, 3>& R2);
     inline void
     UpdateStrNodeData(StrNodeData& oNode,
                       const sp_grad::SpColVector<sp_grad::SpGradient, 3>& d1tot,
                       const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R1tot,
                       const sp_grad::SpColVector<sp_grad::SpGradient, 3>& F,
                       const sp_grad::SpColVector<sp_grad::SpGradient, 3>& M,
                       const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R2) {}

     inline void
     UpdateStrNodeData(StrNodeData& oNode,
                       const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& d1tot,
                       const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R1tot,
                       const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& F,
                       const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& M,
                       const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R2) {}

     inline void
     UpdateModalNode(const sp_grad::SpColVector<doublereal, 3>& x,
                     const sp_grad::SpMatrix<doublereal, 3, 3>& R);

     inline void
     UpdateModalNode(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& x,
                     const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R) {}

     inline void
     UpdateModalNode(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& x,
                     const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R) {}

     inline void
     UpdateState(const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& a,
                 const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& aP,
                 const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& b,
                 const sp_grad::SpColVector<doublereal, sp_grad::SpMatrixSize::DYNAMIC>& bP);

     inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::SpGradient, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::SpGradient, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::SpGradient, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::SpGradient, sp_grad::SpMatrixSize::DYNAMIC>&) { }

     inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::GpGradProd, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, sp_grad::SpMatrixSize::DYNAMIC>&,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, sp_grad::SpMatrixSize::DYNAMIC>&) { }

     inline void
     UpdateInvariants(const sp_grad::SpColVector<doublereal, 3>& Inv3jaj,
                      const sp_grad::SpMatrix<doublereal, 3, 3>& Inv8jaj,
                      const sp_grad::SpMatrix<doublereal, 3, 3>& Inv9jkajak);

     inline void
     UpdateInvariants(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& Inv3jaj,
                      const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& Inv8jaj,
                      const sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& Inv9jkajak) {}

     inline void
     UpdateInvariants(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& Inv3jaj,
                      const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& Inv8jaj,
                      const sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& Inv9jkajak) {}

private:
     const ModalNodeAd* const pModalNode;
};

#endif /* MODALAD_H */
