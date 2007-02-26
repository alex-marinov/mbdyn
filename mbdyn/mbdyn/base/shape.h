/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#include <ac/f2c.h>

#include <myassert.h>
#include <mynewmem.h>
#include <output.h>


/* Classe base:
 * una class Shape puo' restituire un reale in funzione di un insieme di reali
 * per un reale solo si ha una forma monodimensionale, per due una forma
 * bidimensionale.
 */

class Shape {
 public:
   virtual ~Shape(void) { 
      NO_OP; 
   };
   
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
   ShapeOwner(const Shape* pS) : pShape(pS) { 
      ASSERT(pShape != NULL);
   };
   
   virtual ~ShapeOwner(void) {
      ASSERT(pShape != NULL);
      if(pShape != NULL) {
	 SAFEDELETE(pShape);
      }
   };
   
   virtual doublereal dGet(doublereal d) const {
      return pShape->dGet(d);
   };
   
   virtual doublereal dGet(doublereal d1, doublereal d2) const {
      return pShape->dGet(d1, d2);
   };
   
   virtual const Shape* pGetShape(void) const { 
      return pShape;
   };
};



/* Classe base delle forme monodimensionali:
 * una forma monodimensionale restituisce un reale in funzione
 * di un altro reale.
 */

class Shape1D : public Shape {
 public:
   virtual ~Shape1D(void) { NO_OP; };
};

/* Esempi sono:
 * - la forma costante
 * - la forma lineare
 * - la forma parabolica
 */

class ConstShape1D : public Shape1D {
 protected:
   doublereal dConst;
 public:
   ConstShape1D(doublereal d) : dConst(d) { 
      NO_OP;
   };
   
   ~ConstShape1D(void) { 
      DEBUGCOUT("Const 1D shape deleted" << std::endl); 
   };
   
   doublereal dGet(doublereal /* d */ , doublereal = 0.) const { 
      return dConst; 
   };
   
   std::ostream& Restart(std::ostream& out) const { 
      return out << "const, " << dConst; 
   };
};

class LinearShape1D : public Shape1D {
 protected:
   doublereal dShift;
   doublereal dSlope;
 public:
   LinearShape1D(doublereal d0, doublereal d1) : dShift(d0), dSlope(d1) { 
      NO_OP; 
   };
   
   ~LinearShape1D(void) { 
      DEBUGCOUT("Linear 1D shape deleted" << std::endl); 
   };
   
   doublereal dGet(doublereal d, doublereal = 0.) const { 
      return dShift+dSlope*d; 
   };
   
   std::ostream& Restart(std::ostream& out) const { 
      return out << "linear, " << dShift << ", " << dSlope; 
   };
};

class PiecewiseLinearShape1D : public Shape1D {
 protected:
   int nPoints;
   doublereal *pdX;
   doublereal *pdV;
 public:
   PiecewiseLinearShape1D(int n, doublereal *x, doublereal *v)
     : nPoints(n), pdX(x), pdV(v) {
      ASSERT(nPoints > 0);
      ASSERT(pdX != NULL);
      ASSERT(pdV != NULL);
   };

   ~PiecewiseLinearShape1D(void) {
      SAFEDELETEARR(pdX);
      SAFEDELETEARR(pdV);
   };

   doublereal dGet(doublereal d, doublereal = 0.) const {
      if (d <= pdX[0]) {
	 return pdV[0];
      }
      for (int i = 1; i < nPoints; i++) {
	 if (d <= pdX[i]) {
	    doublereal dl = pdX[i]-pdX[i-1];
	    return (pdV[i]*(d-pdX[i-1])+pdV[i-1]*(pdX[i]-d))/dl;
	 }
      }

      return pdV[nPoints-1];
   };

   std::ostream& Restart(std::ostream& out) const {
      out << "piecewise linear, " << nPoints;
      for (int i = 0; i < nPoints; i++) {
	 out << ", " << pdX[i] << ", " << pdV[i];
      }
      return out;
   };
};

class ParabolicShape1D : public Shape1D {
 protected:
   doublereal da0;
   doublereal da1;
   doublereal da2;
 public:
   ParabolicShape1D(doublereal d0, doublereal d1, doublereal d2)
     : da0(d0), da1(d1), da2(d2) { 
	NO_OP; 
     };
   
   ~ParabolicShape1D(void) { 
      DEBUGCOUT("Parabolic 1D shape deleted" << std::endl); 
   };
   
   doublereal dGet(doublereal d, doublereal = 0.) const { 
      return da0+(da1+da2*d)*d; 
   };
   
   std::ostream& Restart(std::ostream& out) const { 
      return out << "parabolic, " << da0 << ", " << da1 << ", " << da2; 
   };
};

/* Classe base delle forme bidimensionali:
 * una forma monodimensionale restituisce un reale in funzione
 * di due reali.
 */

class Shape2D : public Shape {
 public:
   virtual ~Shape2D(void) { NO_OP; };
};

/* Esempi sono:
 * - la forma costante
 * - la forma bilineare
 */

class ConstShape2D : public Shape2D {
 protected:
   doublereal dConst;
 public:
   ConstShape2D(doublereal d) : dConst(d) { 
      NO_OP;
   };
   
   ~ConstShape2D(void) { 
      DEBUGCOUT("Const 2D shape deleted" << std::endl); 
   };
   
   doublereal dGet(doublereal /* d */ , doublereal = 0.) const { 
      return dConst; 
   };
   
   std::ostream& Restart(std::ostream& out) const { 
      return out << "const, " << dConst;
   };
};

class BilinearShape2D : public Shape2D {
 protected:
   doublereal da0;
   doublereal da1x;
   doublereal da1y;
   doublereal da1xy;
 public:
   BilinearShape2D(doublereal d0, doublereal d1x, 
		   doublereal d1y, doublereal d1xy)
     : da0(d0), da1x(d1x), da1y(d1y), da1xy(d1xy) { 
	NO_OP; 
     };
   
   ~BilinearShape2D(void) { 
      NO_OP;
   };
   
   doublereal dGet(doublereal dx, doublereal dy) const { 
      return da0+da1x*dx+da1y*dy+da1xy*dx*dy; 
   };
   
   std::ostream& Restart(std::ostream& out) const { 
      return out << "bilinear, " << da0 << ", "
	<< da1x << ", " << da1y << ", " << da1xy; 
   };
};

#endif /* SHAPE_H */

