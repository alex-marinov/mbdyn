/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003
 *
 * Marco Morandini	<morandini@aero.polimi.it>
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
private:
	void SetValue(VectorHandler&X, VectorHandler&XP) const{};
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
/** Set Initial Values
 */
	virtual void SetValue(VectorHandler&X, 
		VectorHandler&XP, 
		const unsigned int solution_startdof) const;
/** Return last computed friction coefficient
 */
	virtual doublereal fc(void) const = 0;
/** Compute self residual and friction coefficient
 */
	virtual void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
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
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const = 0;
	virtual void AfterConvergence(
		const doublereal F,
		const doublereal v,
		const VectorHandler&X, 
		const VectorHandler&XP,
		const unsigned int solution_startdof) {};
};

/** Base class for friction shape coefficient
 */
class BasicShapeCoefficient {
public:
	virtual ~BasicShapeCoefficient(void){};
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
	const DifferentiableScalarFunction & fss;
	doublereal f;
	doublereal alpha(const doublereal z,
		const doublereal v) const;
	doublereal alphad_v(const doublereal z,
		const doublereal v) const;
	doublereal alphad_z(const doublereal z,
		const doublereal v) const;
	const doublereal fs(const doublereal&v) const;
	const doublereal fsd(const doublereal&v) const;
public:
	ModLugreFriction(
		const doublereal sigma0,
		const doublereal sigma1,
		const doublereal sigma2,
		const doublereal kappa,
		const BasicScalarFunction *const ff);
	void SetValue(VectorHandler&X, 
		VectorHandler&XP, 
		const unsigned int solution_startdof) const;
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	doublereal fc(void) const;
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP);
	void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableRowVector& dfc,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const;
};

/** A Coulomb model based on 
 * Morandini's ideas
 */
class DiscreteCoulombFriction : public BasicFriction {
private:
	enum tr_type{
		null,
		from_sticked_to_sliding,
		from_sticking_to_sliding,
		from_sliding_to_sticked,
		from_sliding_to_sticking};
	enum status_type{
		sticked,
		sticking,
		sliding};
	logical converged_sticked;
	status_type status;
	tr_type transition_type;	
	doublereal converged_v;
	logical first_iter;
	logical first_switch;
	doublereal previous_switch_v;
	doublereal current_velocity;
	
	const DifferentiableScalarFunction & fss;
	doublereal f;
public:
	DiscreteCoulombFriction(
		const BasicScalarFunction *const ff);
	void SetValue(VectorHandler&X, 
		VectorHandler&XP, 
		const unsigned int solution_startdof) const;
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	doublereal fc(void) const;
	void AfterConvergence(
		const doublereal F,
		const doublereal v,
		const VectorHandler&X, 
		const VectorHandler&XP,
		const unsigned int solution_startdof);
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP);
	void AssJac(
		FullSubMatrixHandler& WorkMat,
		ExpandableRowVector& dfc,
		const unsigned int startdof,
		const unsigned int solution_startdof,
		const doublereal dCoef,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const;
};


/** Simple shape coefficient: 1.
 */
class SimpleShapeCoefficient : public BasicShapeCoefficient {
public:
	virtual doublereal Sh_c(void) const;
	virtual doublereal Sh_c(
		const doublereal f,
		const doublereal F,
		const doublereal v);
	virtual void dSh_c(
		ExpandableRowVector& dShc,
		const doublereal f,
		const doublereal F,
		const doublereal v,
		const ExpandableRowVector& dfc,
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

//---------------------------------------

BasicFriction *const ParseFriction(MBDynParser& HP,
	DataManager * pDM);

BasicShapeCoefficient *const ParseShapeCoefficient(MBDynParser& HP);

#endif /* FRICTION_H */

