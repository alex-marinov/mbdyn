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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>

#include "mbpar.h"
#include "datamanforward.h"
#include "friction.h"
#include "submat.h"

int sign(const doublereal x) {
	if (x >= 0.) {
		return 1;
	} else if (x < 0.) {
		return -1;
	}
	return 0;
};

void
BasicFriction::SetValue(DataManager *pDM,
		VectorHandler&X, VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	NO_OP;
}

ModLugreFriction::ModLugreFriction(
		const doublereal s0,
		const doublereal s1,
		const doublereal s2,
		const doublereal k,
		const BasicScalarFunction *const ff) : 
sigma0(s0),
sigma1(s1),
sigma2(s2),
kappa(k),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*ff)),
f(0.)
{
	NO_OP;
}

void
ModLugreFriction::SetValue(DataManager *pDM,
		VectorHandler&X, 
		VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	X.PutCoef(solution_startdof+1,f/sigma0);
}

unsigned int ModLugreFriction::iGetNumDof(void) const {
	return 1;
};

std::ostream&
ModLugreFriction::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: ModLugreFriction state" << std::endl;
}

void
ModLugreFriction::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "ModLugreFriction state";
}

std::ostream&
ModLugreFriction::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: ModLugreFriction equation" << std::endl;
}

void
ModLugreFriction::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "ModLugreFriction equation";
}

DofOrder::Order ModLugreFriction::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetDofType");
	return DofOrder::DIFFERENTIAL;
};

DofOrder::Order ModLugreFriction::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

const doublereal ModLugreFriction::fs(const doublereal&v) const {
	return fss(std::abs(v))*sign(v);
};

const doublereal ModLugreFriction::fsd(const doublereal&v) const {
	return fss.ComputeDiff(std::abs(v),1)*sign(v);
};

doublereal ModLugreFriction::alpha(const doublereal z,
	const doublereal v) const {
	
	doublereal zss = fs(v)/sigma0;
	doublereal zba = kappa*zss;

	if (sign(v)==sign(z)) {
		if (std::abs(z) <= std::abs(zba)) {
		} else if ((std::abs(zba)<=std::abs(z)) &&
			(std::abs(z)<=std::abs(zss))) {
			return 0.5*std::sin(M_PI*sigma0*z/(fs(v)*(1.-kappa))
				-M_PI*(1.+kappa)/(2*(1.-kappa)))+0.5;
		} else {
			return 1.;
		}
	}
	return 0.;
};

doublereal ModLugreFriction::alphad_v(const doublereal z,
	const doublereal v) const {

	doublereal zss = fs(v)/sigma0;
	doublereal zba = kappa*zss;
	doublereal der;

	if (sign(v)==sign(z)) {
		if (std::fabs(z) <= std::fabs(zba)) {
			der = 0.;
		} else if ((std::fabs(zba) <= std::fabs(z)) &&
			(std::fabs(z) <= std::fabs(zss))) {
			der = -M_PI/2*std::cos(M_PI*sigma0*z/fs(v)/(1.-kappa)
				-M_PI*(1.+kappa)/2./(1.-kappa))
				*sigma0*z/std::pow(fs(v),2)/(1.-kappa)
				*fsd(v);
		} else {
			der = 0.;
		}
	} else {
		der = 0.;
	}
	return der;
};

doublereal ModLugreFriction::alphad_z(const doublereal z,
	const doublereal v) const {
	
	doublereal zss = fs(v)/sigma0;
	doublereal zba = kappa*zss;
	doublereal der;
	
	if (sign(v)==sign(z)) {
		if (std::fabs(z) <= std::fabs(zba)) {
			der = 0;
		} else if ((std::fabs(zba) <= std::fabs(z)) && 
			(std::fabs(z) <= std::fabs(zss))) {
			der = 0.5*M_PI*sigma0/fs(v)/(1.-kappa)
				*std::cos(M_PI*sigma0*z/fs(v)/(1.-kappa)
				-M_PI*(1.+kappa)/2./(1.-kappa));
		} else {
			der = 0;
		}
	} else {
		der = 0;
	}
	return der;
};

doublereal ModLugreFriction::fc(void) const {
	return f;
};

void ModLugreFriction::AssRes(
	SubVectorHandler& WorkVec,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal F,
	const doublereal v,
	const VectorHandler& X,
	const VectorHandler& XP) throw(Elem::ChangedEquationStructure) {
	doublereal z = X(solution_startdof+1);
	doublereal zp = XP(solution_startdof+1);
	f = sigma0*z + sigma1*zp + sigma2*v;
	WorkVec.IncCoef(startdof+1,zp-v+alpha(z,v)*v*z/fs(v)*sigma0);
};

void ModLugreFriction::AssJac(
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
	const ExpandableRowVector& dv) const {

	doublereal z = X(solution_startdof+1);
	//doublereal zp = XP(solution_startdof+1);
/*
 * 	attrito
 */
 	dfc.ReDim(2);
	dfc.Set(sigma0*dCoef+sigma1, 1, startdof+1);
	dfc.Set(sigma2, 2); dfc.Link(2, &dv);
	
/*
 * 	z
 */
 	doublereal alph = alpha(z,v);
	doublereal fsc = fs(v);
	WorkMat.IncCoef(startdof+1,startdof+1,-1.);
	WorkMat.IncCoef(startdof+1,startdof+1,
		-alphad_z(z,v)*v*z/fsc*sigma0*dCoef-
		alph*v/fsc*sigma0*dCoef);
	dv.Add(WorkMat,startdof+1,
		1.-
		alphad_v(z,v)*v*z/fsc*sigma0-
		alph*z/fsc*sigma0+
		alph*v*z/(fsc*fsc)*fsd(v)*sigma0);
//	std::cout << alphad_z(z,v) << std::endl;
/*
 * 	callback: dfc[] = df/d{F,v,(z+dCoef,zp)}
 */
};

//----------------------
DiscreteCoulombFriction::DiscreteCoulombFriction(
		const BasicScalarFunction *const ff,
		const doublereal s2,
		const doublereal vr) :
converged_sticked(true),
status(sticked),
transition_type(null),
converged_v(0),
first_iter(true),
first_switch(true),
previous_switch_v(0),
current_velocity(0),
sigma2(s2),
vel_ratio(vr),
current_friction_force(0),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*ff)),
f(0)
{
	NO_OP;
}

void
DiscreteCoulombFriction::SetValue(DataManager *pDM,
		VectorHandler&X, 
		VectorHandler&XP, 
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	X.PutCoef(solution_startdof+1,f);
}

unsigned int DiscreteCoulombFriction::iGetNumDof(void) const {
	return 1;
};

std::ostream&
DiscreteCoulombFriction::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: DiscreteCoulombFriction state" << std::endl;
}

void
DiscreteCoulombFriction::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "DiscreteCoulombFriction state";
}

std::ostream&
DiscreteCoulombFriction::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: DiscreteCoulombFriction equation" << std::endl;
}

void
DiscreteCoulombFriction::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "DiscreteCoulombFriction equation";
}

DofOrder::Order DiscreteCoulombFriction::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetDofType");
	return DofOrder::ALGEBRAIC;
};

DofOrder::Order DiscreteCoulombFriction::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

doublereal DiscreteCoulombFriction::fc(void) const {
	return current_friction_force;
};

void DiscreteCoulombFriction::AfterConvergence(
	const doublereal F,
	const doublereal v,
	const VectorHandler&X, 
	const VectorHandler&XP, 
	const unsigned int solution_startdof) {
	f = X(solution_startdof+1);
// 	std::cerr << " ** ";
// 	std::cerr << f << " " << v << " " << status << " " << transition_type << " ";
	converged_v = v;
	current_velocity = v;
	previous_switch_v = v;
	transition_type = null;
	first_iter = true;
	first_switch = true;
	if (status == sticking) {
//* 		std::cerr << "sticking" << std::endl;
		status = sticked;
	} else if (status == sliding) {
//* 		std::cerr << "sliding" << std::endl;
	} else {
//* 		std::cerr << "sticked" << std::endl;
	}
// 	std::cerr << status << " " << transition_type << std::endl;
//* 	std::cerr << "CONVERGENZA; v = " << v << "; f = " << f << std::endl;
};


void DiscreteCoulombFriction::AssRes(
	SubVectorHandler& WorkVec,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal F,
	const doublereal v,
	const VectorHandler& X,
	const VectorHandler& XP) throw(Elem::ChangedEquationStructure) {
// 	std::cerr << "Chimata residuo. Status:" << status << std::endl;
	f = X(solution_startdof+1);
//	if ((std::fabs(f)-fss(0) > 1.0E-6) && (first_iter == false)) {
//* 	std::cerr << "Attrito: " << f << " " << (std::fabs(f)-fss(0))/fss(0) << " - " << std::endl;
//*	std::cerr << "v: " << v << std::endl; 
	transition_type = null;
// 	std::cerr << "Attrito: " << f << " " << std::fabs(f)/fss(0) << " - " << std::endl;
// 	if ((std::fabs(std::fabs(f)-fss(0))/fss(0) > 1.0E-6)) {
	if (std::fabs(f)-fss(0) > 1.0E-6*fss(0)) {
		//unconditionally switch to sliding
		if (status == sticked) {
			transition_type = from_sticked_to_sliding;
//* 			std::cerr << "switch to sliding from sticked: " << transition_type << std::endl;
		} else if (status == sticking) {
			transition_type = from_sticking_to_sliding;
//* 			std::cerr << "switch to sliding from sticking: " << transition_type << std::endl;
		} else if (status == sliding) {
			//do nothing
// 			std::cerr << "DiscreteCoulombFriction::AssRes message:\n"
// 				<< "you shold not go here1! What's wrong?\n"
// 				<< "status: " << status << "\n"
// 				<< "transition: " << transition_type << "\n"
// 				<< "error: " << fabs(f)-fss(0)
// 				<< std::endl;
		} else {
			silent_cerr("DiscreteCoulombFriction::AssRes() "
					"logical error1" << std::endl);
		}
		status = sliding;
	}
 	//else 
	if (status == sliding) {
// 		std::cerr << "v*current_velocity: " << v*current_velocity << std::endl;
		if (v*current_velocity < 0.) {
//* 			std::cerr << "sono dentro; v: " << v << 
//*				" previous_switch_v: " << previous_switch_v << std::endl;
//* 			std::cerr << "current velocity: " << current_velocity << std::endl;
			if (((transition_type != from_sticked_to_sliding) ||
				(transition_type != from_sticking_to_sliding)) &&
				((std::fabs(v-current_velocity) < std::fabs(previous_switch_v)) ||
					(first_switch == true))) {
//* 				std::cerr << "Passo a sticking $$" << std::endl;;
				first_switch = false;
				status = sticking;
				transition_type = from_sliding_to_sticking;
				previous_switch_v = vel_ratio*(v-current_velocity);
				saved_sliding_velocity = v;
				saved_sliding_friction = f;
// 				std::cerr << "switch to sticking: " << transition_type << std::endl;
			} 
		}
 	}

	switch (status) {
	case sticking: {
		//switch to sticking: null velocity at the end of time step
//		std::cerr << "sono qui1" << std::endl;
		current_friction_force = f;
		WorkVec.IncCoef(startdof+1,v);
		break;
	}
	case sliding: {
		//still sliding
		//cur_sticking = false;
		switch (transition_type) {
		case from_sticked_to_sliding: {
//			std::cerr << "sono qui2" << std::endl;
			current_friction_force = sign(f)*fss(v)+sigma2*v;
			break;
		}
		case from_sticking_to_sliding: {
			current_friction_force = sign(saved_sliding_friction)*fss(v)+sigma2*v;
			break;
		}
		default: {
			if (std::fabs(v) > 0.) {
//				std::cerr << "sono qui3" << std::endl;
//				std::cerr << "v: " << v << std::endl;
				if (sign(v) == sign(current_velocity)) {
					current_friction_force = sign(v)*fss(v)+sigma2*v;
//*					std::cerr << "quixxxx" << std::endl;
				} else {
					current_friction_force = sign(f)*fss(v)+sigma2*v;
//*					std::cerr << "quiyyyy" << std::endl;
				}
			} else {
				//limit the force value while taking the sticking force direction
				current_friction_force = sign(f)*fss(v)+sigma2*v;
			}
		 	if (std::abs(v) < std::abs(current_velocity) && !first_iter) {
//*				std::cerr << "Aggiorno current velocity" << std::endl;
				current_velocity = v;
		 	}
			break;
		}
		}
		//save friction force value in the (algebric) state
//		std::cerr << "sono qui5" << std::endl;
		WorkVec.IncCoef(startdof+1,f-current_friction_force);
		break;
	}
	case sticked: {
		current_friction_force = f;
		WorkVec.IncCoef(startdof+1,v);
		break;
	}
	default: {
		silent_cerr("DiscreteCoulombFriction::AssRes() "
			"logical error" << std::endl);
	}
	}
	//update status
	first_iter = false;
	if (transition_type != null) {
		throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
//	current_velocity = v;
};

void DiscreteCoulombFriction::AssJac(
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
	const ExpandableRowVector& dv) const {
//* 	std::cerr << "Chimata jacobiano. Status:" << status << std::endl;
	switch (status) {
	case sticking: 
	case sticked: {
		//null velocity at the end of time step
		//WorkVec.IncCoef(startdof+1,v);
		dv.Sub(WorkMat,startdof+1);
		dfc.ReDim(1);
		dfc.Set(1.,1,startdof+1);
		break;
	}
	case sliding: {
		//still sliding
		////printf("here1\n");
		//cur_sticking = false;
		//save friction force value in the (algebric) state
		//WorkVec.IncCoef(startdof+1,f-current_friction_force);
		WorkMat.IncCoef(startdof+1,startdof+1,-1);
		dv.Add(WorkMat,startdof+1,
			sign(current_friction_force)*fss.ComputeDiff(v)+sigma2);
		dfc.ReDim(1);
		dfc.Set(sign(current_friction_force-sigma2*v)*fss.ComputeDiff(v)+sigma2,1); dfc.Link(1, &dv);
		break;
	}
	default: {
		silent_cerr("DiscreteCoulombFriction::AssJac() "
			"logical error" << std::endl);
	}
	}
};



//------------------------
doublereal SimpleShapeCoefficient::Sh_c(void) const {
	return 1.;
}
	
doublereal SimpleShapeCoefficient::Sh_c(
	const doublereal f,
	const doublereal F,
	const doublereal v) {
	return f;
};

void SimpleShapeCoefficient::dSh_c(
	ExpandableRowVector& dShc,
	const doublereal f,
	const doublereal F,
	const doublereal v,
	const ExpandableRowVector& dfc,
	const ExpandableRowVector& dF,
	const ExpandableRowVector& dv) const {
		dShc.ReDim(1);
		dShc.Set(1. ,1);
		dShc.Link(1,&dfc);
};

SimplePlaneHingeJointSh_c::SimplePlaneHingeJointSh_c()
 {};
	
doublereal SimplePlaneHingeJointSh_c::Sh_c(void) const {
	return shc;
}
	
doublereal SimplePlaneHingeJointSh_c::Sh_c(
	const doublereal f,
	const doublereal F,
	const doublereal v) {
	shc = f/std::sqrt(1.+f*f);
	return shc;
};

void SimplePlaneHingeJointSh_c::dSh_c(
	ExpandableRowVector& dShc,
	const doublereal f,
	const doublereal F,
	const doublereal v,
	const ExpandableRowVector& dfc,
	const ExpandableRowVector& dF,
	const ExpandableRowVector& dv) const {
		doublereal dsh_fc = 1./std::sqrt(1.+f*f)-0.5*std::pow(1.+f*f,-3./2.)*f;
// 		dShc.ReDim(2);
// 		dShc.Set(0.,1);
// 		dShc.Link(1,&dF);
// 		dShc.Set(dsh_fc,2);
// 		dShc.Link(2,&dfc);
		dShc.ReDim(1);
		dShc.Set(dsh_fc,1);
		dShc.Link(1,&dfc);
};


//---------------------------------------

BasicFriction *const ParseFriction(MBDynParser& HP,
	DataManager * pDM) 
{
   const char* sKeyWords[] = { 
      "modlugre",
      "discrete" "coulomb",
      NULL
   };
	enum KeyWords { 
	     MODLUGRE = 0,
	     DISCRETECOULOMB,
	     LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case MODLUGRE: {
		doublereal sigma0 = HP.GetReal();
		doublereal sigma1 = HP.GetReal();
		doublereal sigma2 = HP.GetReal();
		doublereal kappa = HP.GetReal();
		const BasicScalarFunction*const sf = 
			ParseScalarFunction(HP, pDM);
		return new ModLugreFriction(sigma0, sigma1, sigma2, kappa, sf);
		break;
	}
	case DISCRETECOULOMB: {
		const BasicScalarFunction*const sf = 
			ParseScalarFunction(HP, pDM);
		doublereal sigma2 = 0.;
		doublereal vel_ratio = 0.8;
		if (HP.IsKeyWord("sigma2")) {
			sigma2 = HP.GetReal();
		}
		if (HP.IsKeyWord("velocity" "ratio")) {
			vel_ratio = HP.GetReal();
		}
		return new DiscreteCoulombFriction(sf,sigma2, vel_ratio);
		break;
	}
	default: {
		silent_cerr("ParseFriction(): unrecognized friction type "
				"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};

BasicShapeCoefficient *const ParseShapeCoefficient(MBDynParser& HP) {
	const char* sKeyWords[] = {
		"simple",
		"simple" "plane" "hinge",
		"screw" "joint",
 		NULL
	};
	enum KeyWords { 
		SIMPLE = 0,
		SIMPLEPLANEHINGE,
		SCREWJOINT,
		LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case SIMPLE: {
		return new SimpleShapeCoefficient();
		break;
	}
	case SIMPLEPLANEHINGE: {
		return new SimplePlaneHingeJointSh_c();
		break;
	}
	case SCREWJOINT: {
		doublereal radius(0.), hta(0.);
		if (HP.IsKeyWord("radius")) {
			radius = HP.GetReal();
		} else {
			pedantic_cerr("ScrewJointShapeCoefficient: missing keyword \"radius\" at line "
				<< HP.GetLineData());
		}
		if (HP.IsKeyWord("half" "thread" "angle")) {
			hta = HP.GetReal();
		} else {
			pedantic_cerr("ScrewJointShapeCoefficient: missing keyword \"half thread angle\" at line "
				<< HP.GetLineData());
		}
		return new ScrewJointSh_c(radius, hta);
		break;
	}
	default: {
		silent_cerr("ParseShapeCoefficient(): "
			"unrecognized shape coefficient type "
			"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};
