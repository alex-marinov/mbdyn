/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/* Copyright (C) 2003 Marco Morandini*/

#ifndef FRICTION_H
#define FRICTION_H

#include "ScalarFunctions.h"
#include "simentity.h"
#include "JacSubMatrix.h"

/** Base class for friction models
 */
class BasicFriction : public SimulationEntity{
public:
/*
 * 	unsigned int iGetNumDof(void) const;
 * 	DofOrder::Order GetDofType(unsigned int i) const;
 * 	DofOrder::Order GetEqType (unsigned int i) const;
 * 	void SetValue(VectorHandler&X, VectorHandler&XP) const;
 * 	void BeforePredict(VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&) const;
 * 	void AfterPredict(VectorHandler&X, VectorHandler&XP);
 * 	void Update(const VectorHandler&XCurr, const VectorHandler&XPrimeCurr);
 * 	void AfterConvergence(VectorHandler&X, VectorHandler&XP);
 */
/** Return last computed friction coefficient
 */
	virtual doublereal fc(void) const = 0;
/** Compute self residual and friction coefficient
 */
	virtual void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP) = 0;
/** Compute self jacobian and friction coefficient derivatives
 */
	virtual void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableRowVector& dfc,
		const unsigned int startdof,
		const doublereal dCoef,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const = 0;
};

/** Base class for friction shape coefficient
 */
class BasicShapeCoefficient {
public:
/** Return last computed shape coefficient
 */
	virtual doublereal Sh_c(void) const = 0;
/** Compute the shape coefficient
 */
	virtual doublereal Sh_c(
		const doublereal f,
		const doublereal F,
		const doublereal v) = 0;
/** Compute derivatives of the shape coefficient
 */
	virtual void dSh_c(
		ExpandableRowVector& dShc,
		const doublereal f,
		const doublereal F,
		const doublereal v,
		const ExpandableRowVector& dfc,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const = 0;
};

/** A friction model based on 
 * "Dupont Pierre, Hayward, Vincent, Armstrong Brian
 * and Altpeter Friedhelm, Single state elasto-plastic friction models,
 * IEEE Transactions on Automatic Control, scheduled for June 2002"
 */
class ModLugreFriction : public BasicFriction {
private:
	const doublereal sigma0;
	const doublereal sigma1;
	const doublereal sigma2;
	const doublereal kappa;
	doublereal f;
	const DifferentiableScalarFunction & fss;
	doublereal alpha(const doublereal z,
		const doublereal v) const;
	doublereal alphad_v(const doublereal z,
		const doublereal v) const;
	doublereal alphad_z(const doublereal z,
		const doublereal v) const;
public:
	ModLugreFriction(
		const doublereal sigma0,
		const doublereal sigma1,
		const doublereal sigma2,
		const doublereal kappa,
		const BasicScalarFunction *const f);
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	doublereal fc(void) const;
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP);
	void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableRowVector& dfc,
		const unsigned int startdof,
		const doublereal dCoef,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const;
};

/** Simple, low load shape coefficient for revolute hinge (PlaneHingeJoint)
 */
class SimplePlaneHingeJointSh_c : public BasicShapeCoefficient {
private:
	doublereal r;
	doublereal shc;
public:
	SimplePlaneHingeJointSh_c(const doublereal rr);
	virtual doublereal Sh_c(void) const;
	doublereal Sh_c(
		const doublereal f,
		const doublereal F,
		const doublereal v);
	void dSh_c(
		ExpandableRowVector& dShc,
		const doublereal f,
		const doublereal F,
		const doublereal v,
		const ExpandableRowVector& dfc,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const;
};


#endif /* FRICTION_H */

