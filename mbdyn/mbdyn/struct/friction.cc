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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include "friction.h"
#include "submat.h"

int sign(const doublereal x) {
	if (x > 0.) {
		return 1;
	} else if (x < 0.) {
		return -1;
	}
	return 0;
};

ModLugreFriction::ModLugreFriction(
		const doublereal s0,
		const doublereal s1,
		const doublereal s2,
		const doublereal k,
		const BasicScalarFunction *const f) : 
sigma0(s0),
sigma1(s1),
sigma2(s2),
kappa(k),
f(0.),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*f)) {};

unsigned int ModLugreFriction::iGetNumDof(void) const {
	return 2;
};

DofOrder::Order ModLugreFriction::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetDofType");
	if (i == 0) {
		return DofOrder::ALGEBRAIC;
	} else {
		return DofOrder::DIFFERENTIAL;
	}
};

DofOrder::Order ModLugreFriction::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

doublereal ModLugreFriction::alpha(const doublereal z,
	const doublereal v) const {
	
	doublereal zss = fss(v)/sigma0;
	doublereal zba = kappa*zss;

	if (sign(v)==sign(z)) {
		if (std::abs(z) <= std::abs(zba)) {
		} else if ((std::abs(zba)<=std::abs(z)) && 
			(std::abs(z)<=std::abs(zss))) {
			return 0.5*std::sin(M_PI*sigma0*z/(fss(v)*(1.-kappa))
				-M_PI*(1.+kappa)/(2*(1.-kappa)))+0.5;
		} else {
			return 1.;
		}
	} 
	return 0.;
};

doublereal ModLugreFriction::alphad_v(const doublereal z,
	const doublereal v) const {

	doublereal zss = fss(v)/sigma0;
	doublereal zba = kappa*zss;
	doublereal der;

	if (sign(v)==sign(z)) {
		if (std::fabs(z) <= std::fabs(zba)) {
			der = 0.;
		} else if ((std::fabs(zba) <= std::fabs(z)) &&
			(std::fabs(z) <= std::fabs(zss))) {
			der = -M_PI/2*std::cos(M_PI*sigma0*z/fss(v)/(1.-kappa)
				-M_PI*(1.+kappa)/2./(1.-kappa))
				*sigma0*z/std::pow(fss(v),2)/(1.-kappa)
				*fss.ComputeDiff(v,1);
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
	
	doublereal zss = fss(v)/sigma0;
	doublereal zba = kappa*zss;
	doublereal der;
	
	if (sign(v)==sign(z)) {
		if (std::fabs(z) <= std::fabs(zba)) {
			der = 0;
		} else if ((std::fabs(zba) <= std::fabs(z)) && 
			(std::fabs(z) <= std::fabs(zss))) {
			der = 0.5*M_PI*sigma0/fss(v)/(1.-kappa)
				*std::cos(M_PI*sigma0*z/fss(v)/(1.-kappa)
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
	const doublereal F,
	const doublereal v,
	const VectorHandler& X,
	const VectorHandler& XP) {
	doublereal z = X.dGetCoef(startdof+1);
	doublereal zp = XP.dGetCoef(startdof+1);	
	f = X.dGetCoef(startdof);
	WorkVec.fPutCoef(startdof+1,
		sigma0*z+
		sigma1*zp+
		sigma2*v-f);
	WorkVec.fPutCoef(startdof+2,-zp+v-alpha(z,v)*v*z/fss(v)*sigma0);
};

void ModLugreFriction::AssJac(
	FullSubMatrixHandler& WorkMat,
	ExpandableRowVector& dfc,
	const unsigned int startdof,
	const doublereal dCoef,
	const doublereal F,
	const doublereal v,
	const VectorHandler& X,
	const VectorHandler& XP,
	const ExpandableRowVector& dF,
	const ExpandableRowVector& dv) const {

	doublereal z = X.dGetCoef(startdof+2);
	doublereal zp = XP.dGetCoef(startdof+2);	
/*
 * 	prima equazione
 */
	WorkMat.fPutCoef(startdof+1,startdof+2,-sigma0*dCoef-sigma1);
	dv.Add(WorkMat,startdof+1,-sigma2);
	//f: algebrico
	WorkMat.fPutCoef(startdof+1,startdof+1,1.);
/*
 * 	seconda equazione
 */
 	doublereal alph = alpha(z,v);
	doublereal fs = fss(v);
	WorkMat.fPutCoef(startdof+2,startdof+2,1.);
	dv.Add(WorkMat,startdof+2,
		-1.+
		alphad_v(z,v)*v+
		z/fs*sigma0+
		alph*z/fs*sigma0-
		alph*v*z/(fs*fs)*fss.ComputeDiff(v,1)*sigma0);
/*
 * 	callback: dfc[] = df/d{F,v,(z+dCoef,zp)}
 */
 	dfc.ReDim(3);
	dfc.Set(0.,1);  dfc.Link(1,&dF);
	dfc.Set(sigma2,2); dfc.Link(2,&dv);
	dfc.Set(sigma0*dCoef+sigma1,3,startdof+2);
};


SimplePlaneHingeJointSh_c::SimplePlaneHingeJointSh_c(const doublereal rr): 
	r(rr) {};
	
doublereal SimplePlaneHingeJointSh_c::Sh_c(
	const doublereal f,
	const doublereal F,
	const doublereal v) {
	shc = r/std::sqrt(1.+f*f);
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
		doublereal dsh_fc = r*-0.5*std::pow(1.+f*f,-3./2.)*f;
		dShc.ReDim(2);
		dShc.Set(0.,1);
		dShc.Link(1,&dF);
		dShc.Set(dsh_fc,2);
		dShc.Link(1,&dfc);
};
