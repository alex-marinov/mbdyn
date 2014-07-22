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

#ifndef DGEEQU_H
#define DGEEQU_H

#include <algorithm>
#include <cassert>
#include <limits>
#include <cmath>
#include <vector>
#include <limits>
#include "mh.h"
#include "fullmh.h"
#include "solman.h"

#ifdef USE_LAPACK
#include "ac/lapack.h"
#endif // USE_LAPACK

class MatrixScaleBase
{
public:
	explicit MatrixScaleBase(const SolutionManager::ScaleOpt& scale);
	virtual ~MatrixScaleBase();
	inline VectorHandler& ScaleRightHandSide(VectorHandler& bVH) const;
	inline VectorHandler& ScaleSolution(VectorHandler& xVH) const;
	std::ostream& Report(std::ostream& os) const;
	const std::vector<doublereal>& GetRowScale()const{ return rowScale; }
	const std::vector<doublereal>& GetColScale()const{ return colScale; }
	bool bGetInitialized()const{ return !(rowScale.empty() && colScale.empty()); } // Allow row only or column only scaling

protected:
	inline MatrixHandler::Norm_t GetCondNumNorm()const;
	virtual std::ostream& vReport(std::ostream& os) const=0;
	inline void Prepare(const MatrixHandler& mh, integer& nrows, integer& ncols);
	void PrepareRows(const MatrixHandler& mh, integer& nrows);
	void PrepareCols(const MatrixHandler& mh, integer& ncols);
	inline bool bReport() const;
	std::vector<doublereal> rowScale, colScale;
	mutable doublereal dCondBefore, dCondAfter;
	const unsigned uFlags;
	bool bOK;

private:
	static VectorHandler&
	ScaleVector(VectorHandler& v, const std::vector<doublereal>& s);
};

template <typename T>
class MatrixScale: public MatrixScaleBase
{
public:
	inline explicit MatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~MatrixScale();
	inline T& ScaleMatrix(T& mh) const;
	inline bool ComputeScaleFactors(const T& mh);
	static MatrixScale<T>* Allocate(const SolutionManager::ScaleOpt& scale);

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)=0;
};

template <typename T>
class RowSumMatrixScale: public MatrixScale<T>
{
public:
	inline RowSumMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~RowSumMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	std::vector<doublereal> normRow;
};

template <typename T>
class RowMaxMatrixScale: public MatrixScale<T>
{
public:
	inline RowMaxMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~RowMaxMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	std::vector<doublereal> normRow;
};

template <typename T>
class ColSumMatrixScale: public MatrixScale<T>
{
public:
	inline ColSumMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~ColSumMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	std::vector<doublereal> normCol;
};

template <typename T>
class ColMaxMatrixScale: public MatrixScale<T>
{
public:
	inline ColMaxMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~ColMaxMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	std::vector<doublereal> normCol;
};

// computes scaling factors for a matrix handler that has an iterator
// based on lapack's dgeequ

template <typename T>
class LapackMatrixScale: public MatrixScale<T>
{
public:
	inline LapackMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~LapackMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	doublereal SMLNUM, BIGNUM;
	doublereal rowcnd, colcnd, amax;
};

// computes scaling factors for a matrix handler that has an iterator
// based on `A parallel Matrix Scaling Algorithm'
// from Patrick R. Amestoy, Iain S. Duff, Daniel Ruiz and Bora Ucar

template <typename T>
class IterativeMatrixScale: public MatrixScale<T>
{
public:
	inline IterativeMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~IterativeMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	std::vector<doublereal> DR, DC, normR, normC;
	const integer iMaxIter;
	const doublereal dTol;
	integer iIterTaken;
	doublereal maxNormR, maxNormC;
};

// computes scaling factors for a matrix handler that has an iterator
// based on
// R. Sinkhorn and P. Knopp. Concerning nonnegative matrices and doubly stochastic
// matrices. Pacific Journal of Mathematics, 21(2): 1967.

template <typename T>
class RowMaxColMaxMatrixScale: public MatrixScale<T>
{
public:
	inline RowMaxColMaxMatrixScale(const SolutionManager::ScaleOpt& scale);
	virtual ~RowMaxColMaxMatrixScale();

protected:
	virtual bool ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale);
	virtual std::ostream& vReport(std::ostream& os) const;

private:
	const doublereal dTol;
	doublereal maxNormRow, maxNormCol;
	std::vector<doublereal> normRow, normCol;
};

VectorHandler& MatrixScaleBase::ScaleRightHandSide(VectorHandler& bVH) const
{
	if (!rowScale.empty()) {
		ScaleVector(bVH, rowScale);
	}

	return bVH;
}

VectorHandler& MatrixScaleBase::ScaleSolution(VectorHandler& xVH) const
{
	if (!colScale.empty()) {
		ScaleVector(xVH, colScale);
	}

	return xVH;
}

MatrixHandler::Norm_t MatrixScaleBase::GetCondNumNorm()const
{
	switch (uFlags & SolutionManager::SCALEF_COND_NUM) {
	case SolutionManager::SCALEF_COND_NUM_1:
		return MatrixHandler::NORM_1;

	case SolutionManager::SCALEF_COND_NUM_INF:
		return MatrixHandler::NORM_INF;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

bool MatrixScaleBase::bReport() const
{
	return uFlags & (SolutionManager::SCALEF_WARN | SolutionManager::SCALEF_VERBOSE);
}

void MatrixScaleBase::Prepare(const MatrixHandler& mh, integer& nrows, integer& ncols)
{
	PrepareRows(mh, nrows);
	PrepareCols(mh, ncols);
}

template <typename T>
MatrixScale<T>::MatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScaleBase(scale)
{

}

template <typename T>
MatrixScale<T>::~MatrixScale()
{

}

template <typename T>
T& MatrixScale<T>::ScaleMatrix(T& mh) const
{
	if (uFlags & SolutionManager::SCALEF_COND_NUM) {
		dCondBefore = mh.ConditionNumber(GetCondNumNorm());
	}

	const bool bScaleRows = !rowScale.empty();
	const bool bScaleCols = !colScale.empty();

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		doublereal dCoef = i->dCoef;

		if (bScaleRows) {
			dCoef *= rowScale[i->iRow];
		}

		if (bScaleCols) {
			dCoef *= colScale[i->iCol];
		}

		mh(i->iRow + 1, i->iCol + 1) = dCoef;
	}

	if (uFlags & SolutionManager::SCALEF_COND_NUM) {
		dCondAfter = mh.ConditionNumber(GetCondNumNorm());
	}

	return mh;
}

template <typename T>
bool MatrixScale<T>::ComputeScaleFactors(const T& mh)
{
	bOK = ComputeScaleFactors(mh, rowScale, colScale);

	return bOK;
}

template <typename T>
MatrixScale<T>* MatrixScale<T>::Allocate(const SolutionManager::ScaleOpt& scale)
{
	MatrixScale<T>* pMatScale = 0;

	switch (scale.algorithm) {
	case SolutionManager::SCALEA_ROW_MAX:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   RowMaxMatrixScale<T>,
							   RowMaxMatrixScale<T>(scale));
		break;

	case SolutionManager::SCALEA_ROW_SUM:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   RowSumMatrixScale<T>,
							   RowSumMatrixScale<T>(scale));
		break;

	case SolutionManager::SCALEA_COL_MAX:
			SAFENEWWITHCONSTRUCTOR(pMatScale,
								   ColMaxMatrixScale<T>,
								   ColMaxMatrixScale<T>(scale));
			break;

	case SolutionManager::SCALEA_COL_SUM:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   ColSumMatrixScale<T>,
							   ColSumMatrixScale<T>(scale));
		break;

	case SolutionManager::SCALEA_LAPACK:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   LapackMatrixScale<T>,
							   LapackMatrixScale<T>(scale));
		break;

	case SolutionManager::SCALEA_ITERATIVE:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   IterativeMatrixScale<T>,
							   IterativeMatrixScale<T>(scale));
		break;

	case SolutionManager::SCALEA_ROW_MAX_COL_MAX:
		SAFENEWWITHCONSTRUCTOR(pMatScale,
							   RowMaxColMaxMatrixScale<T>,
							   RowMaxColMaxMatrixScale<T>(scale));
		break;
	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pMatScale;
}

template <typename T>
RowSumMatrixScale<T>::RowSumMatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScale<T>(scale)
{

}

template <typename T>
RowSumMatrixScale<T>::~RowSumMatrixScale()
{

}

template <typename T>
bool RowSumMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	integer nrows;

	MatrixScaleBase::PrepareRows(mh, nrows);

	if (normRow.empty()) {
		normRow.resize(nrows, 0.);
	} else {
		std::fill(normRow.begin(), normRow.end(), 0.);
	}

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		normRow[i->iRow] += std::abs(i->dCoef);
	}

	for (int i = 0; i < nrows; ++i) {
		rowScale[i] = 1. / normRow[i];
	}

	return true;
}

template <typename T>
std::ostream& RowSumMatrixScale<T>::vReport(std::ostream& os) const
{
	return os;
}

template <typename T>
RowMaxMatrixScale<T>::RowMaxMatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScale<T>(scale)
{

}

template <typename T>
RowMaxMatrixScale<T>::~RowMaxMatrixScale()
{

}

template <typename T>
bool RowMaxMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	integer nrows;

	MatrixScaleBase::PrepareRows(mh, nrows);

	if (normRow.empty()) {
		normRow.resize(nrows, 0.);
	} else {
		std::fill(normRow.begin(), normRow.end(), 0.);
	}

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		const doublereal d = std::abs(i->dCoef);

		if (d > normRow[i->iRow]) {
			normRow[i->iRow] = d;
		}
	}

	for (int i = 0; i < nrows; ++i) {
		rowScale[i] = 1. / normRow[i];
	}

	return true;
}

template <typename T>
std::ostream& RowMaxMatrixScale<T>::vReport(std::ostream& os) const
{
	return os;
}

template <typename T>
ColSumMatrixScale<T>::ColSumMatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScale<T>(scale)
{

}

template <typename T>
ColSumMatrixScale<T>::~ColSumMatrixScale()
{

}

template <typename T>
bool ColSumMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	integer ncols;

	MatrixScaleBase::PrepareCols(mh, ncols);

	if (normCol.empty()) {
		normCol.resize(ncols, 0.);
	} else {
		std::fill(normCol.begin(), normCol.end(), 0.);
	}

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		normCol[i->iCol] += std::abs(i->dCoef);
	}

	for (int i = 0; i < ncols; ++i) {
		colScale[i] = 1. / normCol[i];
	}

	return true;
}

template <typename T>
std::ostream& ColSumMatrixScale<T>::vReport(std::ostream& os) const
{
	return os;
}

template <typename T>
ColMaxMatrixScale<T>::ColMaxMatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScale<T>(scale)
{

}

template <typename T>
ColMaxMatrixScale<T>::~ColMaxMatrixScale()
{

}

template <typename T>
bool ColMaxMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	integer ncols;

	MatrixScaleBase::PrepareCols(mh, ncols);

	if (normCol.empty()) {
		normCol.resize(ncols, 0.);
	} else {
		std::fill(normCol.begin(), normCol.end(), 0.);
	}

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		const doublereal d = std::abs(i->dCoef);

		if (d > normCol[i->iCol]) {
			normCol[i->iCol] = d;
		}
	}

	for (int i = 0; i < ncols; ++i) {
		colScale[i] = 1. / normCol[i];
	}

	return true;
}

template <typename T>
std::ostream& ColMaxMatrixScale<T>::vReport(std::ostream& os) const
{
	return os;
}

template <typename T>
LapackMatrixScale<T>::LapackMatrixScale(const SolutionManager::ScaleOpt& scale)
	:MatrixScale<T>(scale)
{
#if defined(HAVE_DLAMCH) || defined(HAVE_DLAMCH_)
	// Use dlamch according to Netlib's dgeequ.f
	SMLNUM = __FC_DECL__(dlamch)("S");
#else
	// According to Netlib's dlamch for x86-64 machines
	SMLNUM = std::numeric_limits<doublereal>::min();
#endif
	BIGNUM = 1./SMLNUM;

	rowcnd = colcnd = amax = -1;
}

template <typename T>
LapackMatrixScale<T>::~LapackMatrixScale()
{

}

template <typename T>
bool LapackMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	integer nrows, ncols;
	MatrixScale<T>::Prepare(mh, nrows, ncols);

	doublereal rcmin;
	doublereal rcmax;

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		doublereal d = std::abs(i->dCoef);
		if (d > rowScale[i->iRow]) {
			rowScale[i->iRow] = d;
		}
	}

	rcmin = BIGNUM;
	rcmax = 0.;
	for (std::vector<doublereal>::iterator i = rowScale.begin(); i != rowScale.end(); ++i) {
		if (*i > rcmax) {
			rcmax = *i;
		}
		if (*i < rcmin) {
			rcmin = *i;
		}
	}

	amax = rcmax;

	if (rcmin == 0.) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"null min row value in dgeequ");
	}

	for (std::vector<doublereal>::iterator i = rowScale.begin(); i != rowScale.end(); ++i) {
		*i = 1./(std::min(std::max(*i, SMLNUM), BIGNUM));
	}

	rowcnd = std::max(rcmin, SMLNUM)/std::min(rcmax, BIGNUM);

	for (typename T::const_iterator i = mh.begin(); i != mh.end(); ++i) {
		doublereal d = std::abs(i->dCoef)*rowScale[i->iRow];
		if (d > colScale[i->iCol]) {
			colScale[i->iCol] = d;
		}
	}

	rcmin = BIGNUM;
	rcmax = 0.;
	for (std::vector<doublereal>::iterator i = colScale.begin(); i != colScale.end(); ++i) {
		if (*i > rcmax) {
			rcmax = *i;
		}
		if (*i < rcmin) {
			rcmin = *i;
		}
	}

	if (rcmin == 0.) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"null min column value in dgeequ");
	}

	for (std::vector<doublereal>::iterator i = colScale.begin(); i != colScale.end(); ++i) {
		*i = 1./(std::min(std::max(*i, SMLNUM), BIGNUM));
	}

	colcnd = std::max(rcmin, SMLNUM)/std::min(rcmax, BIGNUM);

	return true;
}

template <typename T>
std::ostream& LapackMatrixScale<T>::vReport(std::ostream& os) const
{
	if (amax < std::numeric_limits<doublereal>::epsilon()
		|| amax > 1./std::numeric_limits<doublereal>::epsilon())
	{
		os << "Warning: The matrix should be scaled\n";
	}

	if (colcnd >= 0.1) {
		os << "Warning: it is not worth scaling the columns\n";
	}

	if (rowcnd >= 0.1
		&& amax >= std::numeric_limits<doublereal>::epsilon()
		&& amax <= 1./std::numeric_limits<doublereal>::epsilon())
	{
		os << "Warning: it is not worth scaling the rows\n";
	}

	return os;
}

template <typename T>
IterativeMatrixScale<T>::IterativeMatrixScale(const SolutionManager::ScaleOpt& scale)
:MatrixScale<T>(scale),
 iMaxIter(scale.iMaxIter),
 dTol(scale.dTol),
 iIterTaken(-1),
 maxNormR(-1.),
 maxNormC(-1.)
{

}

template <typename T>
IterativeMatrixScale<T>::~IterativeMatrixScale()
{

}

template <typename T>
bool IterativeMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	const integer nrows = mh.iGetNumRows();
	const integer ncols = mh.iGetNumCols();

	if (nrows <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative row number");
	}

	if (ncols <= 0) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative column number");
	}

	if (rowScale.empty()) {
		rowScale.resize(nrows);
		normR.resize(nrows);
		DR.resize(nrows);
		std::fill(rowScale.begin(), rowScale.end(), 1.);
	} else if (rowScale.size() != static_cast<size_t>(nrows)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "row number mismatch");
	} else {
		// Use the scale values from the last Newton iteration
	}

	if (colScale.empty()) {
		colScale.resize(ncols);
		normC.resize(ncols);
		DC.resize(ncols);
		std::fill(colScale.begin(), colScale.end(), 1.);
	} else if (colScale.size() != static_cast<size_t>(ncols)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "column number mismatch");
	} else {
		// Use the scale values from the last Newton iteration
	}

	int i;
	bool bConverged = false, bFirstTime = true;

	while (true) {
		for (i = 1; i <= iMaxIter; ++i) {
			std::fill(DR.begin(), DR.end(), 0.);
			std::fill(DC.begin(), DC.end(), 0.);

			for (typename T::const_iterator im = mh.begin(); im != mh.end(); ++im) {
				doublereal d = im->dCoef;

				if (d == 0.) {
					continue;
				}

				d = std::abs(d * rowScale[im->iRow] * colScale[im->iCol]);

				if (d > DR[im->iRow]) {
					DR[im->iRow] = d;
				}

				if (d > DC[im->iCol]) {
					DC[im->iCol] = d;
				}
			}

			for (int j = 0; j < nrows; ++j) {
				rowScale[j] /= sqrt(DR[j]);
			}

			for (int j = 0; j < ncols; ++j) {
				colScale[j] /= sqrt(DC[j]);
			}

			std::fill(normR.begin(), normR.end(), 0.);
			std::fill(normC.begin(), normC.end(), 0.);

			for (typename T::const_iterator im = mh.begin(); im != mh.end(); ++im) {
				doublereal d = im->dCoef;

				if (d == 0.) {
					continue;
				}

				d = std::abs(d * rowScale[im->iRow] * colScale[im->iCol]);

				if (d > normR[im->iRow]) {
					normR[im->iRow] = d;
				}

				if (d > normC[im->iCol]) {
					normC[im->iCol] = d;
				}
			}

			ASSERT(normR.size() > 0);
			ASSERT(normC.size() > 0);

			maxNormR = 0.;

			for (std::vector<doublereal>::const_iterator ir = normR.begin();
				 ir != normR.end(); ++ir) {
				maxNormR = std::max(maxNormR, std::abs(1. - *ir));
			}

			maxNormC = 0.;

			for (std::vector<doublereal>::const_iterator ic = normC.begin();
				 ic != normC.end(); ++ic) {
				maxNormC = std::max(maxNormC, std::abs(1. - *ic));
			}

			if (maxNormR < dTol && maxNormC < dTol) {
				bConverged = true;
				break;
			}
		}

		if (bConverged) {
			break;	// Scale factors have been computed successfully!
		} else if (bFirstTime) {
			// No convergence: Reset the scale factors and try again!
			std::fill(rowScale.begin(), rowScale.end(), 1.);
			std::fill(colScale.begin(), colScale.end(), 1.);
			bFirstTime = false;
		} else {
			// Still no convergence: Bail out!
			break;
		}
	}

	iIterTaken = i;

	return bConverged;
}

template <typename T>
std::ostream& IterativeMatrixScale<T>::vReport(std::ostream& os) const
{
	if (!MatrixScaleBase::bOK) {
		os << "Warning: matrix scale did not converge\n";
	}

	os << "row scale: " << maxNormR << std::endl
	   << "col scale: " <<  maxNormC << std::endl
	   << "iter scale: " << iIterTaken << std::endl;

	return os;
}

template <typename T>
RowMaxColMaxMatrixScale<T>::RowMaxColMaxMatrixScale(const SolutionManager::ScaleOpt& scale)
:MatrixScale<T>(scale),
 dTol(scale.dTol),
 maxNormRow(-1.),
 maxNormCol(-1.)
{

}

template <typename T>
RowMaxColMaxMatrixScale<T>::~RowMaxColMaxMatrixScale()
{

}

template <typename T>
bool RowMaxColMaxMatrixScale<T>::ComputeScaleFactors(const T& mh, std::vector<doublereal>& rowScale, std::vector<doublereal>& colScale)
{
	const integer nrows = mh.iGetNumRows();
	const integer ncols = mh.iGetNumCols();

	if (nrows <= 0) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative row number");
	}

	if (ncols <= 0) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS,
			"invalid null or negative column number");
	}

	if (rowScale.empty()) {
		rowScale.resize(nrows, 1.);
		normRow.resize(nrows);
	} else if (rowScale.size() != static_cast<size_t>(nrows)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "row number mismatch");
	}

	if (colScale.empty()) {
		colScale.resize(ncols, 1.);
		normCol.resize(ncols);
	} else if (colScale.size() != static_cast<size_t>(ncols)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS, "column number mismatch");
	}

	std::fill(normRow.begin(), normRow.end(), 0.);

	for (typename T::const_iterator im = mh.begin(); im != mh.end(); ++im) {
		doublereal d = im->dCoef;

		if (d == 0.) {
			continue;
		}

		d = std::abs(d);

		if (d > normRow[im->iRow]) {
			normRow[im->iRow] = d;
		}
	}

	for (integer i = 0; i < nrows; ++i) {
		rowScale[i] = 1. / normRow[i];
	}

	std::fill(normCol.begin(), normCol.end(), 0.);

	for (typename T::const_iterator im = mh.begin(); im != mh.end(); ++im) {
		doublereal d = im->dCoef;

		if (d == 0.) {
			continue;
		}

		d = std::abs(rowScale[im->iRow] * d);

		if (d > normCol[im->iCol]) {
			normCol[im->iCol] = d;
		}
	}

	for (integer i = 0; i < ncols; ++i) {
		colScale[i] = 1. / normCol[i];
	}

	if (!this->bReport()) {
		// Test for convergence will not be checked
		// It would be a waste of time to compute the row and column norms
		return true;
	}

	std::fill(normRow.begin(), normRow.end(), 0.);
	std::fill(normCol.begin(), normCol.end(), 0.);

	for (typename T::const_iterator im = mh.begin(); im != mh.end(); ++im) {
		doublereal d = im->dCoef;

		if (d == 0.) {
			continue;
		}

		d = std::abs(colScale[im->iCol] * rowScale[im->iRow] * d);

		if (d > normRow[im->iRow]) {
			normRow[im->iRow] = d;
		}

		if (d > normCol[im->iCol]) {
			normCol[im->iCol] = d;
		}
	}

	maxNormRow = 0.;

	for (std::vector<doublereal>::const_iterator ir = normRow.begin();
		 ir != normRow.end(); ++ir) {
		maxNormRow = std::max(maxNormRow, std::abs(1. - *ir));
	}

	maxNormCol = 0.;

	for (std::vector<doublereal>::const_iterator ic = normCol.begin();
		 ic != normCol.end(); ++ic) {
		maxNormCol = std::max(maxNormCol, std::abs(1. - *ic));
	}

	return (maxNormRow < dTol && maxNormCol < dTol);
}

template <typename T>
std::ostream& RowMaxColMaxMatrixScale<T>::vReport(std::ostream& os) const
{
	if (!MatrixScaleBase::bOK) {
		os << "Warning: matrix scale did not converge\n";
	}

	os << "row scale: " << maxNormRow << std::endl
	   << "col scale: " << maxNormCol << std::endl;

	return os;
}

#endif // DGEEQU_H
