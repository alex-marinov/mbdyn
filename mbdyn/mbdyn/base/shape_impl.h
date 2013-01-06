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

/* Classe di forme pluridimensionali */

#ifndef SHAPE_IMPL_H
#define SHAPE_IMPL_H

#include "shape.h"

/* Esempi sono:
 * - la forma costante
 * - la forma lineare
 * - la forma parabolica
 */

class PiecewiseConstShape1D : public Shape1D {
protected:
	int nPoints;
	doublereal *pdX;
	doublereal *pdV;

public:
	PiecewiseConstShape1D(int n, doublereal *x, doublereal *v);

	~PiecewiseConstShape1D(void);

	doublereal dGet(doublereal d, doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

class LinearShape1D : public Shape1D {
protected:
	doublereal dShift;
	doublereal dSlope;

public:
	LinearShape1D(doublereal d0, doublereal d1);

	~LinearShape1D(void);

	doublereal dGet(doublereal d, doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

class PiecewiseLinearShape1D : public Shape1D {
protected:
	int nPoints;
	doublereal *pdX;
	doublereal *pdV;

public:
	PiecewiseLinearShape1D(int n, doublereal *x, doublereal *v);

	~PiecewiseLinearShape1D(void);

	doublereal dGet(doublereal d, doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

class ParabolicShape1D : public Shape1D {
protected:
	doublereal da0;
	doublereal da1;
	doublereal da2;

public:
	ParabolicShape1D(doublereal d0, doublereal d1, doublereal d2);

	~ParabolicShape1D(void);

	doublereal dGet(doublereal d, doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

/* Esempi sono:
 * - la forma costante
 * - la forma bilineare
 */

class ConstShape2D : public Shape2D {
protected:
	doublereal dConst;

public:
	ConstShape2D(doublereal d);

	~ConstShape2D(void);

	doublereal dGet(doublereal /* d */ , doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

class BilinearShape2D : public Shape2D {
protected:
	doublereal da0;
	doublereal da1x;
	doublereal da1y;
	doublereal da1xy;

public:
	BilinearShape2D(doublereal d0, doublereal d1x,
		doublereal d1y, doublereal d1xy);

	~BilinearShape2D(void);

	doublereal dGet(doublereal dx, doublereal dy) const;

	std::ostream& Restart(std::ostream& out) const;
};

#endif // SHAPE_IMPL_H

