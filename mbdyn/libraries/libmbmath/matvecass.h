/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2013(-2013) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef ___MATVECASS_H__INCLUDED___
#define ___MATVECASS_H__INCLUDED___

#include "gradient.h"
#include "matvec.h"
#include "submat.h"

namespace grad {

template <typename T>
class GradientVectorHandler;

template <>
class GradientVectorHandler<doublereal> {
public:
	GradientVectorHandler(const VectorHandler& vh)
		:vh(vh) {

	}

	void dGetCoef(integer iRow, doublereal& dVal, doublereal dCoef, LocalDofMap* pDofMap) const {
		dVal = vh.dGetCoef(iRow);
	}

	template <index_type N_rows>
	void GetVec(integer iRow, Vector<double, N_rows>& v, doublereal dCoef, LocalDofMap* pDofMap) const {
		for (integer i = 0; i < v.iGetNumRows(); ++i) {
			v(i + 1) = vh.dGetCoef(iRow + i);
		}
	}

private:
	const VectorHandler& vh;
};

template <index_type N_SIZE>
class GradientVectorHandler<Gradient<N_SIZE> > {
public:
	GradientVectorHandler(const VectorHandler& vh)
		:vh(vh) {

	}

	void dGetCoef(integer iRow, Gradient<N_SIZE>& gVal, doublereal dCoef, LocalDofMap* pDofMap) const {
		gVal.SetValuePreserve(vh.dGetCoef(iRow));
		gVal.DerivativeResizeReset(pDofMap, iRow, MapVectorBase::GLOBAL, -dCoef);
	}

	template <index_type N_rows>
	void GetVec(integer iRow, Vector<Gradient<N_SIZE>, N_rows>& v, doublereal dCoef, LocalDofMap* pDofMap) const {
		for (integer i = 0; i < v.iGetNumRows(); ++i) {
			Gradient<N_SIZE>& v_i = v(i + 1);
			v_i.SetValuePreserve(vh.dGetCoef(iRow + i));
			v_i.DerivativeResizeReset(pDofMap, iRow, iRow + v.iGetNumRows(), MapVectorBase::GLOBAL, 0.);
			v_i.SetDerivativeGlobal(iRow + i, -dCoef);
		}
	}
private:
	const VectorHandler& vh;
};

class GradientAssVecBase {
public:
	enum mode_t { RESET, APPEND };
};

template <typename T>
class GradientAssVec;

template <>
class GradientAssVec<doublereal>: GradientAssVecBase {
public:
	GradientAssVec(SubVectorHandler& vh, enum mode_t mode = RESET)
		:WorkVec(vh) {

		switch (mode) {
		case APPEND:
			iSubRow = WorkVec.iGetSize();
			break;

		case RESET:
			iSubRow = 0;
			WorkVec.Resize(iSubRow);
			break;

		default:
			GRADIENT_ASSERT(0);
		}
	}

	template <typename T>
	static void
	AssRes(T* pElem,
		   SubVectorHandler& WorkVec,
		   doublereal dCoef,
		   const VectorHandler& XCurr,
		   const VectorHandler& XPrimeCurr,
		   enum FunctionCall func,
		   enum mode_t mode = RESET) {
		const GradientVectorHandler<doublereal> XCurr_grad(XCurr);
		const GradientVectorHandler<doublereal> XPrimeCurr_grad(XPrimeCurr);
		GradientAssVec WorkVec_grad(WorkVec, mode);

		pElem->AssRes(WorkVec_grad, dCoef, XCurr_grad, XPrimeCurr_grad, func);
	}

	template <typename T>
	static void
	InitialAssRes(T* pElem,
		   SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr,
		   enum FunctionCall func,
		   enum mode_t mode = RESET) {
		const GradientVectorHandler<doublereal> XCurr_grad(XCurr);
		GradientAssVec WorkVec_grad(WorkVec, mode);

		pElem->InitialAssRes(WorkVec_grad, XCurr_grad, func);
	}

	void AddItem(const integer iRow, const double dCoef) {
		WorkVec.Resize(++iSubRow);

		GRADIENT_ASSERT(iSubRow <= WorkVec.iGetSize());

		WorkVec.PutItem(iSubRow, iRow, dCoef);
	}

	template <index_type N_rows>
	void AddItem(const integer iFirstRow, const Vector<doublereal, N_rows>& v) {
		// zero based index according to VectorHandler::Put(integer iRow, const Vec3& v)
		WorkVec.Resize(iSubRow + v.iGetNumRows());

		for (integer i = 0; i < v.iGetNumRows(); ++i) {
			GRADIENT_ASSERT(iSubRow + 1 <= WorkVec.iGetSize());
			WorkVec.PutItem(++iSubRow, iFirstRow + i, v(i + 1));
		}
	}

private:
	SubVectorHandler& WorkVec;
	integer iSubRow;
};

template <index_type N_SIZE>
class GradientAssVec<Gradient<N_SIZE> >: public GradientAssVecBase {
public:
	GradientAssVec(SparseSubMatrixHandler& mh, enum mode_t mode = RESET)
		:WorkMat(mh) {

		switch (mode) {
		case RESET:
			iNextItem = 1;
			WorkMat.ResizeReset(iNextItem, 0);
			WorkMat.PutItem(1, 1, 1, 0.); //FIXME: avoid SIGSEGV if the matrix is empty
			break;

		case APPEND:
			iNextItem = WorkMat.iGetNumRows() + 1;
			break;

		default:
			GRADIENT_ASSERT(0);
		}

	}

	template <typename T>
	static void AssJac(T* pElem,
					   SparseSubMatrixHandler& WorkMat,
					   doublereal dCoef,
					   const VectorHandler& XCurr,
					   const VectorHandler& XPrimeCurr,
					   enum FunctionCall func,
					   LocalDofMap* pDofMap,
					   enum mode_t mode = RESET) {

		const GradientVectorHandler<Gradient<N_SIZE> > XCurr_grad(XCurr);
		const GradientVectorHandler<Gradient<N_SIZE> > XPrimeCurr_grad(XPrimeCurr);

		GradientAssVec WorkMat_grad(WorkMat, mode);

		if (pDofMap) {
			pDofMap->Reset(func);
		}

		pElem->AssRes(WorkMat_grad, dCoef, XCurr_grad, XPrimeCurr_grad, func);
	}

	template <typename T>
	static void InitialAssJac(T* pElem,
					   	   	  SparseSubMatrixHandler& WorkMat,
					   	   	  const VectorHandler& XCurr,
					   	   	  enum FunctionCall func,
					   	   	  LocalDofMap* pDofMap,
					   	   	  enum mode_t mode = RESET) {

		const GradientVectorHandler<Gradient<N_SIZE> > XCurr_grad(XCurr);

		GradientAssVec WorkMat_grad(WorkMat, mode);

		if (pDofMap) {
			pDofMap->Reset(func);
		}

		pElem->InitialAssRes(WorkMat_grad, XCurr_grad, func);
	}

	void AddItem(integer iRow, const Gradient<N_SIZE>& g) {
		const index_type iNextItemNew = iNextItem + g.iGetLocalSize();

		if (iNextItemNew > WorkMat.iGetNumRows()) {
			WorkMat.Resize(iNextItemNew - 1, 0);
		}

		for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
			const double dCoef = g.dGetDerivativeLocal(i);
			const index_type iDofIndex = g.iGetGlobalDof(i);

			GRADIENT_ASSERT(iNextItem < iNextItemNew);
			GRADIENT_ASSERT(iNextItem <= WorkMat.iGetNumRows());

			WorkMat.PutItem(iNextItem++, iRow, iDofIndex, dCoef);
		}
	}

	template <index_type N_rows>
	void AddItem(integer iFirstRow, const Vector<Gradient<N_SIZE>, N_rows>& v) {
		// zero based index according to VectorHandler::Put(integer iRow, const Vec3& v)

		for (integer i = 0; i < N_rows; ++i) {
			AddItem(iFirstRow + i, v(i + 1));
		}
	}
private:
	SparseSubMatrixHandler& WorkMat;
	integer iNextItem;
};

} // namespace

#endif /* ___MATVECASS_H__INCLUDED___ */
