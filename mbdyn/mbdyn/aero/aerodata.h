/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifndef AERODATA_H
#define AERODATA_H

#include <ac/f2c.h>

#include <myassert.h>
#include <withlab.h>
#include <drive.h>

extern "C" {
#include <aerodc81.h>
}

/* C81Data - begin */

class C81Data : public WithLabel, public c81_data {
public:
   	C81Data(unsigned int uLabel);
};

/* C81Data - end */


/* AeroMemory - begin */

class AeroMemory {
private:
	doublereal	*a;
	doublereal	*t;
	integer		iPoints;
	DriveCaller	*pTime;
	int		numUpdates;

protected:
	virtual int StorageSize(void) const = 0;

public:
	AeroMemory(DriveCaller *pt);
	virtual ~AeroMemory(void);

	void Predict(int i, doublereal alpha, doublereal &alf1, 
			doublereal &alf2);
	void Update(int i);
	void SetNumPoints(int i);
};

/* Memory - end */


/* AeroData - begin */

class AeroData : public AeroMemory {
public:
	enum UnsteadyModel {
		STEADY = 0,
		HARRIS = 1,
		BIELAWA = 2,

		LAST
	};
	
	enum { ALF1 = 8, ALF2 = 9 };

protected:
   	UnsteadyModel unsteadyflag;
   	doublereal VAM[6];
   	doublereal Omega;

	int StorageSize(void) const;
   
public:
   	AeroData(UnsteadyModel u = STEADY, DriveCaller *pt = NULL);
   	virtual ~AeroData(void);
   
   	virtual std::ostream& Restart(std::ostream& out) const = 0;   
	std::ostream& RestartUnsteady(std::ostream& out) const;
   	void SetAirData(const doublereal& rho, const doublereal& c);
   
   	virtual void SetSectionData(const doublereal& abscissa,
			    const doublereal& chord,
			    const doublereal& forcepoint,
			    const doublereal& velocitypoint,
			    const doublereal& twist,
			    const doublereal& omega = 0.);
   
   	virtual int
	GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA) = 0;
	inline AeroData::UnsteadyModel Unsteady(void) const;
};


inline AeroData::UnsteadyModel
AeroData::Unsteady(void) const
{
	return unsteadyflag;
}

/* AeroData - end */


/* STAHRAeroData - begin */

class STAHRAeroData : public AeroData {
protected:
   	integer 	profile;
   
public: 
   	STAHRAeroData(AeroData::UnsteadyModel u, integer p, 
			DriveCaller *ptime = NULL);
	virtual ~STAHRAeroData(void);
   
	std::ostream& Restart(std::ostream& out) const;   
   	int GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA);
};

/* STAHRAeroData - end */


/* C81AeroData - begin */

class C81AeroData : public AeroData {
protected:
   	integer profile;
   	const c81_data* data;
   
public: 
   	C81AeroData(AeroData::UnsteadyModel u, integer p, const c81_data* d,
			DriveCaller *ptime = NULL);

	std::ostream& Restart(std::ostream& out) const;
   	int GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA);
};

/* C81AeroData - end */


/* C81MultipleAeroData - begin */

class C81MultipleAeroData : public AeroData {
protected:
   	integer nprofiles;
	integer *profiles;
	doublereal *upper_bounds;
	const c81_data** data;
	integer curr_data;

public: 
   	C81MultipleAeroData(
			AeroData::UnsteadyModel u,
			integer np,
			integer *p,
			doublereal *ub,
			const c81_data** d,
			DriveCaller *ptime = NULL);
   	~C81MultipleAeroData(void);

	std::ostream& Restart(std::ostream& out) const;
   	void SetSectionData(const doublereal& abscissa,
			    const doublereal& chord,
			    const doublereal& forcepoint,
			    const doublereal& velocitypoint,
			    const doublereal& twist,
			    const doublereal& omega = 0.);
   
   	int GetForces(int i, doublereal* W, doublereal* TNG, doublereal* OUTA);
};

/* C81AeroData - end */

#endif /* AERODATA_H */

