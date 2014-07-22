/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef SHAPE_H
#define SHAPE_H

#include "ac/f2c.h"

#include "myassert.h"
#include "mynewmem.h"
#include "output.h"


/* Classe base:
 * una class Shape puo' restituire un reale in funzione di un insieme di reali
 * per un reale solo si ha una forma monodimensionale, per due una forma
 * bidimensionale.
 */

class Shape {
public:
	virtual ~Shape(void);

	virtual doublereal dGet(doublereal d1, doublereal d2 = 0.) const = 0;
	virtual std::ostream& Restart(std::ostream& out) const = 0;
};

/* Possessore di puntatore a shape.
 * Questo oggetto consente un uso comodo di shapes senza alcuna conoscenza
 * sulla loro natura.
 * Viene costruito con un puntatore a Shape di natura qualsiasi.
 * Lui consente l'accesso in lettura
 * ed assicura la corretta distruzione della shape.
 */

class ShapeOwner {
protected:
	const Shape* pShape;

public:
	ShapeOwner(const Shape* pS);

	virtual ~ShapeOwner(void);

	virtual doublereal dGet(doublereal d) const;

	virtual doublereal dGet(doublereal d1, doublereal d2) const;

	virtual const Shape* pGetShape(void) const;
};

/* Classe base delle forme monodimensionali:
 * una forma monodimensionale restituisce un reale in funzione
 * di un altro reale.
 */

class Shape1D : public Shape {
public:
	virtual ~Shape1D(void);
};

// public because needed outside
class ConstShape1D : public Shape1D {
protected:
	doublereal dConst;

public:
	ConstShape1D(doublereal d);

	~ConstShape1D(void);

	doublereal dGet(doublereal /* d */ , doublereal = 0.) const;

	std::ostream& Restart(std::ostream& out) const;
};

/* Classe base delle forme bidimensionali:
 * una forma monodimensionale restituisce un reale in funzione
 * di due reali.
 */

class Shape2D : public Shape {
public:
	virtual ~Shape2D(void);
};

extern Shape* ReadShape(MBDynParser& HP);

#endif // SHAPE_H

