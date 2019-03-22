/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2017
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
#include "elem.h"
#include "JacSubMatrix.h"

/** Base class for friction models
 */
class BasicFriction : public SimulationEntity{
private:
	void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0) { NO_OP; };
public:
/*
 * 	unsigned int iGetNumDof(void) const;
 * 	DofOrder::Order GetDofType(unsigned int i) const;
 * 	DofOrder::Order GetEqType (unsigned int i) const;
 * 	void SetValue(VectorHandler&X, VectorHandler&XP);
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
	virtual void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP,
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
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
		const VectorHandler& XP) 
			throw(Elem::ChangedEquationStructure) = 0;
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
	void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP, 
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
	unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
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
		const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/;
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
	doublereal saved_sliding_velocity;
	doublereal saved_sliding_friction;
	doublereal sigma2;
	doublereal vel_ratio;
	doublereal current_friction_force;
	
	const DifferentiableScalarFunction & fss;
	doublereal f;
public:
	DiscreteCoulombFriction(
		const BasicScalarFunction *const ff,
		const doublereal s2,
		const doublereal vr);
	void SetValue(DataManager *pDM,
			VectorHandler&X, VectorHandler&XP, 
			SimulationEntity::Hints *ph = 0,
			const unsigned int solution_startdof = 0);
	unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
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
		const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/;
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
	doublereal shc;
public:
	SimplePlaneHingeJointSh_c(void);
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

/** Simple, low load shape coefficient for revolute hinge (PlaneHingeJoint)
 */
class ScrewJointSh_c : public BasicShapeCoefficient {
private:
	doublereal shc;
	
	const doublereal radius;
	const doublereal half_thread_angle;
	const doublereal pitch_angle;
	const doublereal sec_half_thread_angle;
	const doublereal tg_pitch;
	const doublereal tg_pitch2;
public:
	ScrewJointSh_c(
		const doublereal r,
		const doublereal hta
	): 
		radius(r), 
		half_thread_angle(hta), 
		pitch_angle(0.),
		sec_half_thread_angle(1. / std::cos(half_thread_angle)),
		tg_pitch(0.),
		tg_pitch2(0.) {
	};
	doublereal ComputePitchAngle(const doublereal pitch) {
		const_cast<doublereal&>(tg_pitch) = pitch / (2. * M_PI * radius);
		const_cast<doublereal&>(tg_pitch2) = tg_pitch * tg_pitch;
		const_cast<doublereal&>(pitch_angle) = std::atan(tg_pitch);
		const_cast<doublereal&>(half_thread_angle) = std::atan(
			std::cos(pitch_angle) * std::tan(half_thread_angle)
		);
		return radius / std::cos(pitch_angle);
	};
	virtual doublereal Sh_c(void) const {return shc;};
	doublereal Sh_c(
		const doublereal f,
		const doublereal F,
		const doublereal v) {
		
		shc  = radius * (f * sec_half_thread_angle * (1 + tg_pitch2)) /
			(f * sec_half_thread_angle * tg_pitch - 1.);
		return shc;
	};
	void dSh_c(
		ExpandableRowVector& dShc,
		const doublereal f,
		const doublereal F,
		const doublereal v,
		const ExpandableRowVector& dfc,
		const ExpandableRowVector& dF,
		const ExpandableRowVector& dv) const {
		doublereal dsh_fc = radius * (-sec_half_thread_angle * (1 + tg_pitch2) ) /
			std::pow(f * sec_half_thread_angle * tg_pitch - 1., 2.);
		dShc.ReDim(1);
		dShc.Set(dsh_fc, 1);
		dShc.Link(1, &dfc);
	};
};

//---------------------------------------

BasicFriction *const ParseFriction(MBDynParser& HP,
	DataManager * pDM);

BasicShapeCoefficient *const ParseShapeCoefficient(MBDynParser& HP);

#endif /* FRICTION_H */

