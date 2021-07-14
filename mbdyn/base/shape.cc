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

/* Classe di forme pluridimensionali */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "mbpar.h"
#include "dataman.h"
#include "shape_impl.h"

Shape::~Shape(void)
{
	NO_OP;
}

ShapeOwner::ShapeOwner(const Shape* pS)
: pShape(pS)
{
	ASSERT(pShape != NULL);
}

ShapeOwner::~ShapeOwner(void)
{
	ASSERT(pShape != NULL);
	if(pShape != NULL) {
		SAFEDELETE(pShape);
	}
}

doublereal
ShapeOwner::dGet(doublereal d) const
{
	return pShape->dGet(d);
}

doublereal
ShapeOwner::dGet(doublereal d1, doublereal d2) const
{
	return pShape->dGet(d1, d2);
}

const Shape *
ShapeOwner::pGetShape(void) const
{
	return pShape;
}

Shape1D::~Shape1D(void)
{
	NO_OP;
}

ConstShape1D::ConstShape1D(doublereal d)
: dConst(d)
{
	NO_OP;
}

ConstShape1D::~ConstShape1D(void)
{
	NO_OP;
}

doublereal
ConstShape1D::dGet(doublereal d, doublereal) const
{
	return dConst;
}

std::ostream&
ConstShape1D::Restart(std::ostream& out) const
{
	return out << "const, " << dConst;
}

PiecewiseConstShape1D::PiecewiseConstShape1D(int n,
	doublereal *x, doublereal *v)
: nPoints(n), pdX(x), pdV(v)
{
	ASSERT(nPoints > 0);
	ASSERT(pdX != NULL);
	ASSERT(pdV != NULL);
}

PiecewiseConstShape1D::~PiecewiseConstShape1D(void)
{
	SAFEDELETEARR(pdX);
	SAFEDELETEARR(pdV);
}

doublereal
PiecewiseConstShape1D::dGet(doublereal d, doublereal) const
{
	if (d <= pdX[0]) {
		return pdV[0];
	}

	for (int i = 1; i < nPoints; i++) {
		if (d < pdX[i]) {
			return pdV[i - 1];
		}
	}

	return pdV[nPoints - 1];
}

std::ostream&
PiecewiseConstShape1D::Restart(std::ostream& out) const
{
	out << "piecewise const, " << nPoints;

	for (int i = 0; i < nPoints; i++) {
		out << ", " << pdX[i] << ", " << pdV[i];
	}

	return out;
}

LinearShape1D::LinearShape1D(doublereal d0, doublereal d1)
: dShift(d0), dSlope(d1)
{
	NO_OP;
}

LinearShape1D::~LinearShape1D(void)
{
	NO_OP;
}

doublereal
LinearShape1D::dGet(doublereal d, doublereal) const
{
	return dShift + dSlope*d;
}

std::ostream&
LinearShape1D::Restart(std::ostream& out) const
{
	return out << "linear, " << dShift << ", " << dSlope;
}

PiecewiseLinearShape1D::PiecewiseLinearShape1D(int n,
	doublereal *x, doublereal *v)
: nPoints(n), pdX(x), pdV(v)
{
	ASSERT(nPoints > 0);
	ASSERT(pdX != NULL);
	ASSERT(pdV != NULL);
}

PiecewiseLinearShape1D::~PiecewiseLinearShape1D(void)
{
	SAFEDELETEARR(pdX);
	SAFEDELETEARR(pdV);
}

doublereal
PiecewiseLinearShape1D::dGet(doublereal d, doublereal) const
{
	if (d <= pdX[0]) {
		return pdV[0];
	}

	for (int i = 1; i < nPoints; i++) {
		if (d <= pdX[i]) {
			doublereal dl = pdX[i] - pdX[i-1];
			doublereal dw1 = pdX[i] - d;
			doublereal dw2 = d - pdX[i - 1];
			return (pdV[i]*dw2 + pdV[i - 1]*dw1)/dl;
		}
	}

	return pdV[nPoints - 1];
}

std::ostream&
PiecewiseLinearShape1D::Restart(std::ostream& out) const
{
	out << "piecewise linear, " << nPoints;

	for (int i = 0; i < nPoints; i++) {
		out << ", " << pdX[i] << ", " << pdV[i];
	}

	return out;
}

ParabolicShape1D::ParabolicShape1D(doublereal d0, doublereal d1, doublereal d2)
: da0(d0), da1(d1), da2(d2)
{
	NO_OP;
}

ParabolicShape1D::~ParabolicShape1D(void)
{
	NO_OP;
}

doublereal
ParabolicShape1D::dGet(doublereal d, doublereal) const
{
	return da0 + (da1 + da2*d)*d;
}

std::ostream&
ParabolicShape1D::Restart(std::ostream& out) const
{
	return out << "parabolic, " << da0 << ", "
		<< da1 << ", " << da2;
}

Shape2D::~Shape2D(void)
{
	NO_OP;
}

ConstShape2D::ConstShape2D(doublereal d)
: dConst(d)
{
	NO_OP;
}

ConstShape2D::~ConstShape2D(void)
{
	NO_OP;
}

doublereal
ConstShape2D::dGet(doublereal d, doublereal) const
{
	return dConst;
}

std::ostream&
ConstShape2D::Restart(std::ostream& out) const
{
	return out << "const, " << dConst;
}

BilinearShape2D::BilinearShape2D(doublereal d0, doublereal d1x,
	doublereal d1y, doublereal d1xy)
: da0(d0), da1x(d1x), da1y(d1y), da1xy(d1xy)
{
	NO_OP;
}

BilinearShape2D::~BilinearShape2D(void)
{
	NO_OP;
}

doublereal
BilinearShape2D::dGet(doublereal dx, doublereal dy) const
{
	return da0 + da1x*dx + da1y*dy + da1xy*dx*dy;
}

std::ostream&
BilinearShape2D::Restart(std::ostream& out) const
{
	return out << "bilinear, " << da0 << ", "
		<< da1x << ", " << da1y << ", " << da1xy;
}

/* Legge una shape1D;
 * NOTA: il proprietario del puntatore alla Shape la deve distruggere */

Shape*
ReadShape(MBDynParser& HP)
{
	DEBUGCOUTFNAME("ReadShape");

	const char* sKeyWords[] = {
		"const",
		"piecewise" "const",
		"linear",
		"piecewise" "linear",
		"parabolic",
		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		SHAPECONST = 0,
		PIECEWISECONST,
		LINEAR,
		PIECEWISELINEAR,
		PARABOLIC,
		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* lettura del tipo di drive */
	KeyWords CurrKeyWord;
	if ((CurrKeyWord = KeyWords(HP.IsKeyWord())) == UNKNOWN) {
		CurrKeyWord = SHAPECONST;
	}

#ifdef DEBUG
	if (CurrKeyWord >= 0) {
		std::cout << "shape type: " << sKeyWords[CurrKeyWord] << std::endl;
	}
#endif /* DEBUG */

	Shape* pS = NULL;

	switch (CurrKeyWord) {
	/* forma costante */
	case SHAPECONST: {
		/* lettura dei dati specifici */
		doublereal dConst = HP.GetReal();
		DEBUGLCOUT(MYDEBUG_INPUT, "Const value: " << dConst << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			ConstShape1D,
			ConstShape1D(dConst));
	} break;

	/* forma lineare */
	case LINEAR: {
		/* lettura dei dati specifici */
		doublereal da0;
		doublereal da1;
		if (HP.IsKeyWord("coefficients")) {
			da0 = HP.GetReal();
			da1 = HP.GetReal();
		} else {
			doublereal dm = HP.GetReal();
			doublereal dp = HP.GetReal();
			da0 = (dp + dm)/2.;
			da1 = (dp - dm)/2;
		}

		DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: "
			<< da0 << ", " << da1 << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			LinearShape1D,
			LinearShape1D(da0, da1));
	} break;

	/* forma lineare a tratti (costante al di fuori del dominio definito) */
	case PIECEWISECONST:
	case PIECEWISELINEAR: {
		const char *sType = 0;
		switch (CurrKeyWord) {
		case PIECEWISECONST:
			sType = "piecewise const";
			break;

		case PIECEWISELINEAR:
			sType = "piecewise linear";
			break;

		default:
			ASSERT(0);
			break;
		}

		int np = HP.GetInt();
		if (np <= 0) {
			silent_cerr("Illegal number of points " << np
       				<< " for " << sType << " shape at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal *px = NULL;
		doublereal *pv = NULL;

		SAFENEWARR(px, doublereal, np);
		SAFENEWARR(pv, doublereal, np);

		px[0] = HP.GetReal();
		if (px[0] < -1. || px[0] > 1.) {
			silent_cerr("Illegal value " << px[0]
				<< " for first point abscissa (must be -1. < x < 1.) "
				"in " << sType << " shape at line "
				<< HP.GetLineData() << std::endl);
			SAFEDELETEARR(px);
			SAFEDELETEARR(pv);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pv[0] = HP.GetReal();

		for (int i = 1; i < np; i++) {
			px[i] = HP.GetReal();
			if (px[i] <= px[i-1] || px[i] > 1.) {
				silent_cerr("Illegal value " << px[i]
					<< " for point " << i + 1 << " abscissa "
					"(must be " << px[i - 1] << " < x < 1.) "
					"in " << sType << " shape at line "
					<< HP.GetLineData() << std::endl);
				SAFEDELETEARR(px);
				SAFEDELETEARR(pv);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pv[i] = HP.GetReal();
		}

		/* allocazione e creazione */
		switch (CurrKeyWord) {
		case PIECEWISECONST:
			SAFENEWWITHCONSTRUCTOR(pS,
				PiecewiseConstShape1D,
				PiecewiseConstShape1D(np, px, pv));
			break;

		case PIECEWISELINEAR:
			SAFENEWWITHCONSTRUCTOR(pS,
				PiecewiseLinearShape1D,
				PiecewiseLinearShape1D(np, px, pv));
			break;

		default:
			ASSERT(0);
			break;
		}

	} break;

	/* forma lineare */
	case PARABOLIC: {
		/* lettura dei dati specifici */
		doublereal dm = HP.GetReal();
		doublereal da0 = HP.GetReal();
		doublereal dp = HP.GetReal();
		doublereal da1 = (dp - dm)/2.;
		doublereal da2 = (dp + dm)/2. - da0;
		DEBUGLCOUT(MYDEBUG_INPUT, "Coefficients: "
			<< da0 << ", " << da1 << ", " << da2 << std::endl);

		/* allocazione e creazione */
		SAFENEWWITHCONSTRUCTOR(pS,
			ParabolicShape1D,
			ParabolicShape1D(da0, da1, da2));
	} break;

	/* Non c'e' default in quanto all'inizio il default e' stato messo
	 * pari a SHAPECONST */
	default:
		ASSERTMSG(0, "You shouldn't have reached this point");
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(pS != NULL);
	return pS;
} /* ReadShape */

