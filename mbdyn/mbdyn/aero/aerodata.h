/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifndef AERODATA_H
#define AERODATA_H

#include <myassert.h>
#include <myf2c.h>
extern "C" {
#include <aerodc81.h>
}
#include <withlab.h>

/* C81Data - begin */

class C81Data : public WithLabel, public c81_data {
public:
   	C81Data(unsigned int uLabel);
};

/* C81Data - end */


/* AeroData - begin */

class AeroData {
protected:
   	integer unsteadyflag;
   	doublereal VAM[6];
   	doublereal Omega;
   
public:
   	AeroData(integer u = 0);
   	virtual ~AeroData(void);
   
   	virtual ostream& Restart(ostream& out) const = 0;   
   	void SetAirData(const doublereal& rho, const doublereal& c);
   
   	void SetSectionData(const doublereal& chord,
			    const doublereal& forcepoint,
			    const doublereal& velocitypoint,
			    const doublereal& twist,
			    const doublereal& omega = 0.);
   
   	virtual int
	GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA) = 0;
};

/* AeroData - end */


/* STAHRAeroData - begin */

class STAHRAeroData : public AeroData {
protected:
   	integer profile;
   
public: 
   	STAHRAeroData(integer p, integer u);
   
   	ostream& Restart(ostream& out) const;   
   	int GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA);
};

/* STAHRAeroData - end */


/* C81AeroData - begin */

class C81AeroData : public AeroData {
protected:
   	integer profile;
   	const c81_data* data;
   
public: 
   	C81AeroData(integer p, integer u, const c81_data* d);

   	ostream& Restart(ostream& out) const;
   	int GetForces(doublereal* W, doublereal* TNG, doublereal* OUTA);
};

/* C81AeroData - end */

#endif /* AERODATA_H */

