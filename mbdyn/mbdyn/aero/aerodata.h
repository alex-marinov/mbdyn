/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#include "ac/f2c.h"

#include "myassert.h"
#include "withlab.h"
#include "drive.h"
#include "dofown.h"
#include "matvec6.h"
#include "matvec3n.h"
#include "shape.h"

#include "aerodc81.h"

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
	int		numUpdates;

protected:
	DriveCaller	*pTime;

	virtual int StorageSize(void) const = 0;

public:
	AeroMemory(DriveCaller *pt);
	virtual ~AeroMemory(void);

	void Predict(int i, doublereal alpha,
		doublereal &alf1, doublereal &alf2);
	void Update(int i);
	void SetNumPoints(int i);
	int GetNumPoints(void) const;
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

	enum {
		VX	= 0,
		VY	= 1,
		VZ	= 2,

		WX	= 3,
		WY	= 4,
		WZ	= 5,

		FX	= 0,
		FY	= 1,
		FZ	= 2,

		MX	= 3,
		MY	= 4,
		MZ	= 5
	};

protected:
	UnsteadyModel unsteadyflag;
	vam_t VAM;
	doublereal Omega;

	int StorageSize(void) const;

	int GetForcesJacForwardDiff_int(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
	int GetForcesJacCenteredDiff_int(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);

public:
	AeroData(int i_p, int i_dim,
		UnsteadyModel u = STEADY, DriveCaller *pt = 0);
	virtual ~AeroData(void);

	virtual std::ostream& Restart(std::ostream& out) const = 0;
	std::ostream& RestartUnsteady(std::ostream& out) const;
	virtual void SetAirData(const doublereal& rho, const doublereal& c);

	virtual void SetSectionData(const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega = 0.);

	virtual int
	GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	virtual int
	GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);

	// aerodynamic models with internal states
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual void
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr,
		integer iFirstIndex, integer iFirstSubIndex,
		int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	virtual void
	AssJac(FullSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr,
		integer iFirstIndex, integer iFirstSubIndex,
		const Mat3xN& vx, const Mat3xN& wx, Mat3xN& fq, Mat3xN& cq,
		int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
	virtual void
	AfterConvergence(int i,
		const VectorHandler& X, const VectorHandler& XP);

	AeroData::UnsteadyModel Unsteady(void) const;
};

/* AeroData - end */

extern void
ReadAeroData(DataManager* pDM, MBDynParser& HP, int dim,
	Shape** ppChord, Shape** ppForce,
	Shape** ppVelocity, Shape** ppTwist,
	Shape** ppTipLoss,
	integer* piNumber, DriveCaller** ppDC,
	AeroData** aerodata);

#endif /* AERODATA_H */

