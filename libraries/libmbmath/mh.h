/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

/*****************************************************************************
 *                                                                           *
 *                            SOLUTION MANAGER                               *
 *                                                                           *
 *****************************************************************************/

/* Pierangelo Masarati */


#ifndef MH_H
#define MH_H

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

#include "ac/f2c.h"

/* per il debugging */
#include "myassert.h"
#include "mynewmem.h"
#include "except.h"

#include "vh.h"
#include "sp_gradient_base.h"

class SubMatrixHandler;
class VariableSubMatrixHandler;

/* MatrixHandler - begin */

class MatrixHandler {
public:
	enum Norm_t { NORM_1, NORM_INF };

	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrRebuildMatrix : public MBDynErrBase {
	public:
		ErrRebuildMatrix(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrMatrixIsSingular : public MBDynErrBase {
	public:
		ErrMatrixIsSingular(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

public:
	virtual ~MatrixHandler(void);

#ifdef DEBUG
	/* Usata per il debug */
	virtual void IsValid(void) const = 0;
#endif /* DEBUG */

	/* Ridimensiona la matrice */
	virtual void Resize(integer, integer) = 0;

	/* Ridimensiona ed inizializza. */
	virtual void ResizeReset(integer, integer);

	/* Restituisce un puntatore all'array di reali della matrice */
	virtual inline const doublereal* pdGetMat(void) const;
	virtual inline doublereal* pdGetMat(void);

	/* Restituisce un puntatore al vettore delle righe */
	virtual inline integer* piGetRows(void) const;

	/* Restituisce un puntatore al vettore delle colonne */
	virtual inline integer* piGetCols(void) const;

	/* Impacchetta la matrice; restituisce il numero di elementi 
	 * diversi da zero */
	virtual integer PacMat(void);

	/* Resetta la matrice ecc. */
	virtual void Reset(void) = 0;

	/* Inserisce un coefficiente */
	virtual void
	PutCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Incrementa un coefficiente - se non esiste lo crea */
	virtual void
	IncCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Decrementa un coefficiente - se non esiste lo crea */
	virtual void
	DecCoef(integer iRow, integer iCol, const doublereal& dCoef);

	/* Restituisce un coefficiente - zero se non e' definito */
	virtual const doublereal&
	dGetCoef(integer iRow, integer iCol) const;

	virtual const doublereal&
	operator () (integer iRow, integer iCol) const = 0;

	virtual doublereal&
	operator () (integer iRow, integer iCol) = 0;

	/* dimensioni */
	virtual integer iGetNumRows(void) const = 0;
	virtual integer iGetNumCols(void) const = 0;
	
	/* Overload di = */
	virtual MatrixHandler& operator = (const MatrixHandler& MH);

	/* Overload di += usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator += (const SubMatrixHandler& SubMH);

	/* Overload di -= usato per l'assemblaggio delle matrici */
	virtual MatrixHandler& operator -= (const SubMatrixHandler& SubMH);

	/* Overload di += usato per l'assemblaggio delle matrici
	 * questi li vuole ma non so bene perche'; forse per la doppia
	 * derivazione di VariableSubMatrixHandler */
	virtual MatrixHandler&
	operator += (const VariableSubMatrixHandler& SubMH);
	virtual MatrixHandler&
	operator -= (const VariableSubMatrixHandler& SubMH);

	/* */
	virtual MatrixHandler& ScalarMul(const doublereal& d);

        enum MatPrintFormat {
              MAT_PRINT_FULL,
              MAT_PRINT_TRIPLET,
              MAT_PRINT_SPCONVERT
        };
             
        virtual std::ostream& Print(std::ostream& os, MatPrintFormat eFormat = MAT_PRINT_FULL) const;

        typedef void EnumerateNzCallback(integer, integer, doublereal);
                
        virtual void EnumerateNz(const std::function<EnumerateNzCallback>& func) const;
                
        /* Matrix Matrix product */
protected:
	virtual MatrixHandler&
	MatMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatMul_base(void (MatrixHandler::*op)(integer iRow, integer iCol,
				const doublereal& dCoef),
			MatrixHandler& out, const MatrixHandler& in) const;

public:
	virtual MatrixHandler&
	MatMatMul(MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatMul(MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatMatIncMul(MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatIncMul(MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatMatDecMul(MatrixHandler& out, const MatrixHandler& in) const;
	virtual MatrixHandler&
	MatTMatDecMul(MatrixHandler& out, const MatrixHandler& in) const;

	/* Matrix Vector product */
protected:
	virtual VectorHandler&
	MatVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul_base(void (VectorHandler::*op)(integer iRow,
				const doublereal& dCoef),
			VectorHandler& out, const VectorHandler& in) const;

        template <typename T>
        static void IteratorScale(T& oMH,
				  const std::vector<doublereal>& oRowScale,
				  const std::vector<doublereal>& oColScale);
public:
	virtual VectorHandler&
	MatVecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatVecIncMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecIncMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatVecDecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual VectorHandler&
	MatTVecDecMul(VectorHandler& out, const VectorHandler& in) const;
	virtual doublereal ConditionNumber(enum Norm_t eNorm = NORM_1) const;
	virtual doublereal Norm(enum Norm_t eNorm = NORM_1) const;
        virtual void Scale(const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale);
        virtual MatrixHandler* Copy() const=0;
        virtual bool AddItem(integer iRow, const sp_grad::SpGradient& oItem);
};

/* Restituisce un puntatore all'array di reali della matrice */
inline const doublereal*
MatrixHandler::pdGetMat(void) const
{
	return NULL;
}

inline doublereal*
MatrixHandler::pdGetMat(void)
{
	return NULL;
}

/* Restituisce un puntatore al vettore delle righe */
inline integer*
MatrixHandler::piGetRows(void) const
{
	return NULL;
}

/* Restituisce un puntatore al vettore delle colonne */
inline integer*
MatrixHandler::piGetCols(void) const
{
	return NULL;
}

template <typename T>
void MatrixHandler::IteratorScale(T& oMH, const std::vector<doublereal>& oRowScale, const std::vector<doublereal>& oColScale)
{
     static_assert(std::is_base_of<MatrixHandler, T>::value, "argument oMH must be a matrix handler");
     
     const bool bScaleRows = !oRowScale.empty();
     const bool bScaleCols = !oColScale.empty();

     ASSERT(!bScaleRows || oRowScale.size() == static_cast<size_t>(oMH.iGetNumRows()));
     ASSERT(!bScaleCols || oColScale.size() == static_cast<size_t>(oMH.iGetNumCols()));

     for (const auto& oItem: oMH) {
	  doublereal dCoef = oItem.dCoef;

	  if (bScaleRows) {
	       dCoef *= oRowScale[oItem.iRow];
	  }

	  if (bScaleCols) {
	       dCoef *= oColScale[oItem.iCol];
	  }

	  oMH(oItem.iRow + 1, oItem.iCol + 1) = dCoef;
     }
}

extern std::ostream&
operator << (std::ostream& out, const MatrixHandler& MH);

/* MatrixHandler - end */

#endif /* MH_H */

