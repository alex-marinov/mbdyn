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

#ifndef AERODATA_IMPL_H
#define AERODATA_IMPL_H

#include "aerodata.h"

/* STAHRAeroData - begin */

class STAHRAeroData : public AeroData {
protected:
	integer profile;

public:
	STAHRAeroData(
		int i_p, int i_dim,
		AeroData::UnsteadyModel u, integer p,
		DriveCaller *ptime = 0);
	virtual ~STAHRAeroData(void);

	std::ostream& Restart(std::ostream& out) const;
	int GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	int GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
};

/* STAHRAeroData - end */


/* C81AeroData - begin */

class C81AeroData : public AeroData {
protected:
	integer profile;
	const c81_data* data;

public:
	C81AeroData(
		int i_p, int i_dim,
		AeroData::UnsteadyModel u, integer p, const c81_data* d,
		DriveCaller *ptime = 0);
	virtual ~C81AeroData(void);

	virtual std::ostream& Restart(std::ostream& out) const;
	virtual int GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	virtual int GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
};

/* C81AeroData - end */


/* C81MultipleAeroData - begin */

struct C81AirfoilStation {
	integer profile;
	const c81_data *data;
	doublereal upper_bound;
};

class C81MultipleAeroData : public AeroData {
protected:
	std::vector<unsigned> profiles;
	std::vector<doublereal> upper_bounds;
	std::vector<const c81_data *> data;
	integer curr_data;

public:
	C81MultipleAeroData(
		int i_p, int i_dim,
		AeroData::UnsteadyModel u,
		std::vector<unsigned>& p,
		std::vector<doublereal>& ub,
		std::vector<const c81_data*>& d,
		DriveCaller *ptime = 0);
	~C81MultipleAeroData(void);

	std::ostream& Restart(std::ostream& out) const;
	void SetSectionData(const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega = 0.);

	int GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	int GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
};

/* C81MultipleAeroData - end */

/* C81InterpolatedAeroData - begin */

class C81InterpolatedAeroData : public AeroData {
protected:
	std::vector<unsigned> profiles;
	std::vector<doublereal> upper_bounds;
	std::vector<const c81_data *> data;

	std::vector<c81_data> i_data;

public:
	C81InterpolatedAeroData(
		int i_p, int i_dim,
		AeroData::UnsteadyModel u,
		std::vector<unsigned>& p,
		std::vector<doublereal>& ub,
		std::vector<const c81_data *>& d,
		doublereal dcltol,
		DriveCaller *ptime = 0);
	~C81InterpolatedAeroData(void);

	std::ostream& Restart(std::ostream& out) const;
	void SetSectionData(const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega = 0.);

	int GetForces(int i, const doublereal* W, doublereal* TNG, outa_t& OUTA);
	int GetForcesJac(int i, const doublereal* W, doublereal* TNG, Mat6x6& J, outa_t& OUTA);
};

/* C81InterpolatedAeroData - end */

/* TheodorsenAeroData - begin */

class TheodorsenAeroData : public AeroData {
protected:
	integer iParam;
	doublereal d14, d34;
	doublereal chord;
	doublereal a;
	doublereal A1, A2, b1, b2;
	doublereal *alpha_pivot, *dot_alpha_pivot, *dot_alpha, *ddot_alpha;
	doublereal *cfx_0, *cfy_0, *cmz_0;
	doublereal *clalpha;
	doublereal *prev_alpha_pivot, *prev_dot_alpha;
	doublereal prev_time;

	AeroData *pAeroData;

public:
	TheodorsenAeroData(
		int i_p, int i_dim,
		AeroData *pa, DriveCaller *ptime = 0);
	virtual ~TheodorsenAeroData(void);

	virtual std::ostream& Restart(std::ostream& out) const;
	virtual void SetAirData(const doublereal& rho, const doublereal& c);

	virtual void SetSectionData(const doublereal& abscissa,
		const doublereal& chord,
		const doublereal& forcepoint,
		const doublereal& velocitypoint,
		const doublereal& twist,
		const doublereal& omega = 0.);


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
	AfterConvergence( int i, const VectorHandler& X, const VectorHandler& XP );
};

/* TheodorsenAeroData - end */

#endif // AERODATA_IMPL_H

