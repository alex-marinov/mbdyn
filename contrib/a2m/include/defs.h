/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <map.h>
#include <string.h>


// Tipi derivati da tipi più semplice
typedef long double doublereal;

// Costanti matematiche trigonometriche
const double pi=3.141592653589793;
const double rad=360/(2*pi);
const doublereal dRaDegr = (180/pi);

// Tipo Label , o identificatore di Marker
typedef unsigned int Id;

// tipo booleano
enum Boolean { N, Y };

// tipo angolo
enum angle_type {
   DEGREE, RADIANS, GRAD
};

// tipo coordinata per iterazione su Jacobiano
 enum coord_type {
   _XC, _YC, _ZC, _PSIC, _THETAC, _PHIC
};

// Tipo vincolo
enum Joint {
   _CONVEL, _CYLINDRICAL, _FIXED, _HOOKE, _PLANAR, _RACKPIN, _REVOLUTE,
   _SCREW, _SPHERICAL, _TRANSLATIONAL,_UNIVERSAL
};

// tipo frizione relativa al vincolo
enum Friction {
   _ON, _OFF, _PRELOAD_ONLY
};

// tipo primitiva del vincolo
enum Joint_Primitive {
   _ATPOINT,_INLINE,_INPLANE,_ORIENTATION,_PARALLEL_AXES,_PERPENDICULAR
};

// tipo di smorzatore viscoelastico (SPRINGDAMPER)
enum Direction_Mode {
   _TRANSLATION, _ROTATION
};

// Struttura angolo : ANGOLO - tipo     Es: 109 DEGREE
struct Angle {
   Angle() : value(0),ref(RADIANS) {}
   Angle( float v, angle_type r ) : value(v), ref(r) {}
   ~Angle() {}
   float value;
   angle_type ref;
   friend ostream& operator << (ostream& ostr, const Angle& r)
     {
	cout << " " << r.value << " ";
	if (r.ref==DEGREE) cout << "DEGREES";
	else cout << "RADIANS";
	return ostr;
     }
   Angle& operator=(const Angle& r)
     {
	value=r.value;
	ref=r.ref;
     }
   void Set_Angle (double a, angle_type b)
     {
	value=a;
	ref=b;
     }
   void ConvToRadians ()
     {
	if (ref==DEGREE) {
	   value=value/rad; ref=RADIANS;
	}
     }
   void ConvToDegrees ()
     {
	if (ref==RADIANS) {
	   value=value*rad; ref=DEGREE;
	}
     }
};

/* Struttura stringa */
struct Strn {
   char *data;
   int len;
};

/* Definizione di una terna d'angoli, utile per angoli di Eulero*/
struct Euler {
   Euler () : Xangle(0,RADIANS),Yangle(0,RADIANS),Zangle(0,RADIANS) {}
   Euler (float a,angle_type at,float b,angle_type bt,float c,angle_type ct) : 
          Xangle(a,at),Yangle(b,bt),Zangle(c,ct) {}
   ~Euler () {}
   Angle Xangle,Yangle,Zangle;
   
   Euler& operator = (const Euler& r)
     {
	Xangle=r.Xangle;
	Yangle=r.Yangle;
	Zangle=r.Zangle;
     }
   friend ostream& operator << (ostream& ostr, const Euler& r)
     {
	cout << r.Xangle << ", " << r.Yangle << ", " << r.Zangle << endl;
	return ostr;
     }
   void Set_Euler (Angle a, Angle b, Angle c)
     {
	Xangle=a;
	Yangle=b;
	Zangle=c;
     }
   void Set (double A1,angle_type R1,
	     double A2,angle_type R2,
	     double A3,angle_type R3)
     {
	Xangle.value=A1; Xangle.ref=R1;
	Yangle.value=A2; Yangle.ref=R2;
	Zangle.value=A3; Zangle.ref=R3;
     }
   void ConvToRadians (void)
     {
	if (Xangle.ref==DEGREE) {
	   Xangle.ref=RADIANS;
	   Xangle.value=Xangle.value/57.3;
	}
	if (Yangle.ref==DEGREE) {
	   Yangle.ref=RADIANS;
	   Yangle.value=Yangle.value/57.3;
	}
	if (Zangle.ref==DEGREE) {
	   Zangle.ref=RADIANS;
	   Zangle.value=Zangle.value/57.3;
	}
     }
   void ConvToDegrees (void)
     {
	if (Xangle.ref==RADIANS) {
	   Xangle.ref=DEGREE;
	   Xangle.value=Xangle.value*57.3;
	}
	if (Yangle.ref==RADIANS) {
	   Yangle.ref=DEGREE;
	   Yangle.value=Yangle.value*57.3;
	}
	if (Zangle.ref==RADIANS) {
	   Zangle.ref=DEGREE;
	   Zangle.value=Zangle.value*57.3;
	}
     }
};

typedef map < char* , char* , less<char*> > formula_map;
typedef formula_map::value_type formula_entry;
typedef formula_map::iterator p_formula_entry;

#endif
