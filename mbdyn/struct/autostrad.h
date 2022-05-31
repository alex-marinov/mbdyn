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

#ifndef AUTOSTRAD_H
#define AUTOSTRAD_H

#include "autostr.h"
#include "strnodead.h"
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_matvecass.h"

class AutomaticStructDispElemAd: public AutomaticStructDispElem {
public:
     explicit
     AutomaticStructDispElemAd(const DynamicStructDispNodeAd* pN);

     virtual ~AutomaticStructDispElemAd();

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
     inline void
     UpdateState(const sp_grad::SpColVector<doublereal, 3>& B,
                 const sp_grad::SpColVector<doublereal, 3>& BP);

     static inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& B,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& BP) {
     }

     static inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& B,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& BP) {
     }

     const DynamicStructDispNodeAd* const pNode;
};

class AutomaticStructElemAd: public AutomaticStructElem {
public:
     explicit
     AutomaticStructElemAd(const DynamicStructNodeAd* pN);

     virtual ~AutomaticStructElemAd();

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

     inline void
     UpdateState(const sp_grad::SpColVector<doublereal, 3>& B,
                 const sp_grad::SpColVector<doublereal, 3>& BP,
                 const sp_grad::SpColVector<doublereal, 3>& G,
                 const sp_grad::SpColVector<doublereal, 3>& GP);

     static inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& B,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& BP,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& G,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& GP) {
     }

     static inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& B,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& BP,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& G,
                 const sp_grad::SpColVector<sp_grad::GpGradProd, 3>& GP) {
     }

private:
     DynamicStructNodeAd* const pNode;
};

#endif /* AUTOSTRAD_H */
