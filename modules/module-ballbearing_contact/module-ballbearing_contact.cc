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

/*
 AUTHOR: Reinhard Resch <r.resch@a1.net>
        Copyright (C) 2013(-2013) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <limits>
#include <iostream>
#include <iomanip>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <cstring>
#include <ctime>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include <userelem.h>

#include "module-ballbearing_contact.h"
#ifdef USE_AUTODIFF
#include <vector>
#include <gradient.h>
#include <matvec.h>
#include <matvecass.h>

using namespace grad;

class BallBearingContact: virtual public Elem, public UserDefinedElem
{
public:
	BallBearingContact(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~BallBearingContact(void);
	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		   doublereal dCoef,
		   const VectorHandler& XCurr,
		   const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	template <typename T>
	inline void
	AssRes(GradientAssVec<T>& WorkVec,
	       doublereal dCoef,
	       const GradientVectorHandler<T>& XCurr,
	       const GradientVectorHandler<T>& XPrimeCurr,
	       enum FunctionCall func);
	virtual void
	AfterConvergence(const VectorHandler& X,
					const VectorHandler& XP);
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetNumDof(void) const;
	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix, bool bInitial) const;
	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix, bool bInitial) const;
	virtual DofOrder::Order GetDofType(unsigned int) const;
	virtual DofOrder::Order GetEqType(unsigned int i) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

  private:
   	static const index_type iNumADVars = 14;
   	const DataManager* const pDM;
	const StructNode* pNode1;
    typedef Matrix<doublereal, 2, 2> Mat2x2;
    doublereal R;
    doublereal k;
    Mat2x2 Mk, Mk2, inv_Mk, inv_Mk_sigma0, Ms, Ms2, sigma0, sigma1;
	doublereal gamma, vs, beta;
    doublereal dStictionStateEquScale;
    doublereal dStictionStateDofScale;
    doublereal tPrev;
    doublereal tCurr;
    doublereal dtMax;
    doublereal dFMax;
    doublereal dFMin;
    bool bEnableFriction;
    bool bFirstRes;

    struct Washer {
    	Washer()
    	: pNode2(0),
    	  dPrev(0),
    	  dCurr(0),
    	  dd_dtPrev(0),
    	  dd_dtCurr(0),
    	  FnPrev(0),
    	  FnCurr(0)
    	{

    	}
    	const StructNode* pNode2;
		Matrix<doublereal, 3, 3> Rt2;
		Vector<doublereal, 3> o2;
		LocalDofMap dof;
		doublereal dPrev;
		doublereal dCurr;
		doublereal dd_dtPrev;
		doublereal dd_dtCurr;
		doublereal FnPrev;
		doublereal FnCurr;
		Vector<doublereal, 2> z;
		Vector<doublereal, 2> zP;
		Vector<doublereal, 2> uP;
        Vector<doublereal, 2> tau;
    };

    std::vector<Washer> washers;

    inline void CheckTimeStep(Washer& w, doublereal Fn, doublereal d, doublereal dn_dt);
    inline void CheckTimeStep(Washer& w, const Gradient<iNumADVars>&, const Gradient<iNumADVars>&, const Gradient<iNumADVars>&);

    static const int iNumPrivData = 11;

	static const struct PrivData {
		char szPattern[9];
	} rgPrivData[iNumPrivData];
};

const struct BallBearingContact::PrivData
BallBearingContact::rgPrivData[iNumPrivData] = {
	{"d[%u]"},
	{"dP[%u]"},
	{"F[%u]"},
	{"z1[%u]"},
	{"z2[%u]"},
	{"zP1[%u]"},
	{"zP2[%u]"},
	{"uP1[%u]"},
	{"uP2[%u]"},
        {"tau1[%u]"},
        {"tau2[%u]"}
};

BallBearingContact::BallBearingContact(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	pDM(pDM),
	pNode1(0),
	gamma(0.),
	vs(0.),
	beta(0.),
	dStictionStateEquScale(0.),
	dStictionStateDofScale(0.),
	tPrev(-std::numeric_limits<doublereal>::max()),
	tCurr(-std::numeric_limits<doublereal>::max()),
	dtMax(std::numeric_limits<doublereal>::max()),
	dFMax(std::numeric_limits<doublereal>::max()),
	dFMin(0.),
	bEnableFriction(false),
	bFirstRes(true)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"\n"
			"Module: 	ballbearing contact\n"
			"\n"
			"	ball bearing contact,\n"
			"		ball, (label) <node1>,\n"
			"		ball radius, (Scalar) <R>,\n"
			"		washers, (Integer) <N>,\n"
			"			(label) <node2>,\n"
			"			[offset, (Vec3) <o2>,]\n"
			"			[orientation, (OrientationMatrix) <Rt2>,\n"
			"			...\n"
			"		{elastic modulus, <E> |\n"
			"		 elastic modulus ball, (Scalar) <E1>,\n"
			"		 poisson ratio ball, (Scalar) <nu1>\n"
			"		 elastic modulus washer, (Scalar) <E2>,\n"
			"		 poisson ratio washer, (Scalar) <nu2>},\n"
			"		 [damping factor, (Scalar) <ks>,]\n"
			"		{coulomb friction coefficient, (Scalar) <mu> |\n"
			"		 coulomb friction coefficient x, (Scalar) <mux>,\n"
			"		 coulomb friction coefficient y, (Scalar) <muy>},\n"
			"		[{static friction coefficient, (Scalar) <mus> |\n"
			"		  static friction coefficient x, (Scalar) <musx>,\n"
			"		  static friction coefficient y, (Scalar) <musy>},]\n"
			"		[sliding velocity coefficient, (Scalar) <vs>,]\n"
			"		[sliding velocity exponent, (Scalar) <gamma>,]\n"
			"		{micro slip stiffness, (Scalar) <sigma0> |\n"
			"		 micro slip stiffness x, (Scalar) <sigma0x>,\n"
			"		 micro slip stiffness y, (Scalar) <sigma0y>}\n"
			"		[,{micro slip damping, (Scalar) <sigma1>, |\n"
			"		   micro slip damping x, (Scalar) <sigma1x>,\n"
			"		   micro slip damping y, (Scalar) <sigma1y>}]\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("ball")) {
		silent_cerr("ball bearing contact" << GetLabel()
				<< "): keyword \"ball\" expected at line "
				<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	if (!HP.IsKeyWord("ball" "radius") ) {
		silent_cerr("ball bearing contact(" << GetLabel()
					<< "): keyword \"ball radius\" expected at line "
					<< HP.GetLineData() << std::endl);

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	R = HP.GetReal();

	if (!HP.IsKeyWord("washers")) {
		silent_cerr("ball bearing contact" << GetLabel()
				<< "): keyword \"washers\" expected at line "
				<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const integer N = HP.GetInt();

	washers.resize(N);

	for (std::vector<Washer>::iterator i = washers.begin(); i != washers.end(); ++i) {
		i->pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

		if ( HP.IsKeyWord("offset") ) {
			i->o2 = HP.GetPosRel(ReferenceFrame(i->pNode2));
		} else {
			i->o2 = Zero3;
		}

		if (HP.IsKeyWord("orientation")) {
			i->Rt2 = HP.GetRotRel(ReferenceFrame(i->pNode2));
		} else {
			i->Rt2 = Eye3;
		}
	}

	doublereal E;

	if ( HP.IsKeyWord("elastic" "modulus") ) {
		E = HP.GetReal();
	} else {
		if ( !HP.IsKeyWord("elastic" "modulus" "ball") ) {
			silent_cerr("ball bearing contact(" << GetLabel()
						<< "): keyword \"elastic modulus\" or "
						   "\"elastic modulus ball\" expected at line "
						<< HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal E1 = HP.GetReal();

		if ( !HP.IsKeyWord("poisson" "ratio" "ball") ) {
			silent_cerr("ball bearing contact(" << GetLabel()
					<<  "): keyword \"poisson ratio ball\" expected at line "
					<< HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal nu1 = HP.GetReal();

		if ( !HP.IsKeyWord("elastic" "modulus" "washer") ) {
			silent_cerr("ball bearing contact(" << GetLabel()
					<< "): keyword \"elastic modulus washer\" expected at line "
					<< HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal E2 = HP.GetReal();

		if ( !HP.IsKeyWord("poisson" "ratio" "washer") ) {
			silent_cerr("ball bearing contact(" << GetLabel()
					<< "): keyword \"poisson ratio washer\" expected at line "
					<< HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal nu2 = HP.GetReal();

		E = 1 / ((1 - nu1 * nu1) / E1 + (1 - nu2 * nu2) / E2);
	}

	k = 4./3. * E * sqrt(R);

	if (HP.IsKeyWord("damping" "factor")) {
		const doublereal ks = HP.GetReal();
		beta = 3./2. * ks * k;
	} else {
		beta = 0.;
	}

	if (HP.IsKeyWord("coulomb" "friction" "coefficient")
			|| HP.IsKeyWord("coulomb" "friction" "coefficient" "x")) {
		bEnableFriction = true;

		const doublereal mukx = HP.GetReal();

		doublereal muky;

		if (HP.IsKeyWord("coulomb" "friction" "coefficient" "y")) {
			muky = HP.GetReal();
		} else {
			muky = mukx;
		}

		doublereal musx, musy;

		if (HP.IsKeyWord("static" "friction" "coefficient")
				|| HP.IsKeyWord("static" "friction" "coefficient" "x")) {
			musx = HP.GetReal();

			if (HP.IsKeyWord("static" "friction" "coefficient" "y")) {
				musy = HP.GetReal();
			} else {
				musy = musx;
			}
		} else {
			musx = mukx;
			musy = muky;
		}

		if (HP.IsKeyWord("sliding" "velocity" "coefficient")) {
			vs = HP.GetReal();
		} else {
			vs = 1.;
		}

		if (HP.IsKeyWord("sliding" "velocity" "exponent")) {
			gamma = HP.GetReal();
		} else {
			gamma = 1.;
		}

		if (!(HP.IsKeyWord("micro" "slip" "stiffness") || HP.IsKeyWord("micro" "slip" "stiffness" "x"))) {
			silent_cerr("ball bearing contact("
					<< GetLabel()
					<< "): keyword \"micro slip stiffness\" or \"micro slip stiffness x\" expected at line "
					<< HP.GetLineData() << std::endl);

			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		const doublereal sigma0x = HP.GetReal();

		doublereal sigma0y;

		if (HP.IsKeyWord("micro" "slip" "stiffness" "y")) {
			sigma0y = HP.GetReal();
		} else {
			sigma0y = sigma0x;
		}

		doublereal sigma1x, sigma1y;

		if (HP.IsKeyWord("micro" "slip" "damping") || HP.IsKeyWord("micro" "slip" "damping" "x")) {
			sigma1x = HP.GetReal();

			if (HP.IsKeyWord("micro" "slip" "damping" "y")) {
				sigma1y = HP.GetReal();
			} else {
				sigma1y = sigma1x;
			}
		} else {
			sigma1x = 0.;
			sigma1y = 0.;
		}

		Mk(1, 1) = mukx;
		Mk(2, 2) = muky;

		Mk2 = Mk * Mk;

		inv_Mk = Inv(Mk);

		Ms(1, 1) = musx;
		Ms(2, 2) = musy;

		Ms2 = Ms * Ms;

		sigma0(1, 1) = sigma0x;
		sigma0(2, 2) = sigma0y;

		sigma1(1, 1) = sigma1x;
		sigma1(2, 2) = sigma1y;

                sigma1 = Mat2x2(sigma1 * Inv(sigma0));

                inv_Mk_sigma0 = inv_Mk * sigma0;
                
		dStictionStateEquScale = HP.IsKeyWord("stiction" "state" "equation" "scale") ? HP.GetReal() : 1.;
		dStictionStateDofScale = HP.IsKeyWord("stiction" "state" "dof" "scale") ? HP.GetReal() : 1.;
	}

    if (HP.IsKeyWord("max" "force" "increment")) {
    	dFMax = HP.GetReal();

    	if (dFMax <= 0) {
    		silent_cerr("ball bearing contact(" << GetLabel()
    				<< "): max force increment must be greater than zero at line "
    				<< HP.GetLineData() << std::endl);
    		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    	}
    }

    if (HP.IsKeyWord("min" "force" "increment")) {
    	dFMin = HP.GetReal();

    	if (dFMin <= 0) {
    		silent_cerr("ball bearing contact(" << GetLabel()
    				<< "): min force increment must be greater than zero at line "
    				<< HP.GetLineData() << std::endl);
    		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    	}
    }

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

BallBearingContact::~BallBearingContact(void)
{

}

void
BallBearingContact::Output(OutputHandler& OH) const
{

}

void
BallBearingContact::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	integer iNumDof = 12;

	if (bEnableFriction) {
		iNumDof += 2;
	}

	*piNumRows = iNumDof * washers.size();
	*piNumCols = iNumDof * washers.size();
}

VariableSubMatrixHandler&
BallBearingContact::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	GradientAssVec<Gradient<iNumADVars> >::AssJac(this,
												  WorkMat.SetSparse(),
												  dCoef,
												  XCurr,
												  XPrimeCurr,
												  REGULAR_JAC,
												  0);
	return WorkMat;
}

SubVectorHandler&
BallBearingContact::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	GradientAssVec<doublereal>::AssRes(this,
									   WorkVec,
									   dCoef,
									   XCurr,
									   XPrimeCurr,
									   REGULAR_RES);

	return WorkVec;
}

template <typename T>
inline void BallBearingContact::AssRes(GradientAssVec<T>& WorkVec,
									  doublereal dCoef,
									  const GradientVectorHandler<T>& XCurr,
									  const GradientVectorHandler<T>& XPrimeCurr,
									  enum FunctionCall func)
{
	typedef Matrix<T, 3, 3> VMat3x3;
	typedef Vector<T, 3> VVec3;
	typedef Vector<T, 2> VVec2;

	VVec3 X1, X1P, omega1;
	VMat3x3 R2;
	VVec3 X2, X2P, omega2;
	VVec2 z, zP;
	VVec2 tau, Phi;
	T norm_Fn, kappa;

	integer offset = std::numeric_limits<integer>::min();
	const integer iFirstIndex = bEnableFriction ? iGetFirstIndex() : std::numeric_limits<integer>::min();
	const integer iFirstMomIndexNode1 = pNode1->iGetFirstMomentumIndex();
	dtMax = std::numeric_limits<doublereal>::max();

	for (std::vector<Washer>::iterator i = washers.begin(); i != washers.end(); ++i) {
		i->dof.Reset(func);
		pNode1->GetXCurr(X1, dCoef, func, &i->dof);
		pNode1->GetVCurr(X1P, dCoef, func, &i->dof);
		pNode1->GetWCurr(omega1, dCoef, func, &i->dof);

		i->pNode2->GetXCurr(X2, dCoef, func, &i->dof);
		i->pNode2->GetVCurr(X2P, dCoef, func, &i->dof);
		i->pNode2->GetRCurr(R2, dCoef, func, &i->dof);
		i->pNode2->GetWCurr(omega2, dCoef, func, &i->dof);

		if (bEnableFriction) {
			offset = 2 * (i - washers.begin());

			for (integer j = 1; j <= 2; ++j) {
				XCurr.dGetCoef(iFirstIndex + j + offset, z(j), dCoef, &i->dof);
				XPrimeCurr.dGetCoef(iFirstIndex + j + offset, zP(j), 1., &i->dof);
			}

			z *= dStictionStateDofScale;
			zP *= dStictionStateDofScale;

			for (int j = 1; j <= 2; ++j) {
				i->z(j) = dGetValue(z(j));
				i->zP(j) = dGetValue(zP(j));
			}
		}

		const VVec3 dX = X1 - X2;

		VVec3 v = Transpose(i->Rt2) * VVec3(Transpose(R2) * dX - i->o2);

		const T n = v(3);

		v(3) = 0.;

		const T dn_dt = Dot(R2 * i->Rt2.GetCol(3), X1P - X2P - Cross(omega2, dX));
		const T d = R - n;

		if (R > n) {
			norm_Fn = k * pow(d, 3./2.) - beta * sqrt(d) * dn_dt;
		} else {
			norm_Fn = 0.;
		}

		if (norm_Fn < 0.) {
			norm_Fn = 0.;
		}

		CheckTimeStep(*i, norm_Fn, d, dn_dt);

		if (bEnableFriction) {
			const VVec3 c1P = X1P - Cross(omega1, R2 * VVec3(i->Rt2.GetCol(3) * R));
			const VVec3 c2P = X2P + Cross(omega2, R2 * VVec3(i->o2 + i->Rt2 * v));
			
			VVec3 l = v;
			l(3) = R;
			const VVec3 Deltac = X1 - X2 - R2 * VVec3(i->o2 + i->Rt2 * l);
			
			const VVec2 uP = Transpose(SubMatrix<1, 3, 1, 2>(i->Rt2)) * VVec3(Transpose(R2) * VVec3(c1P - c2P - Cross(omega2, Deltac)));

			tau = Mk * VVec2(z + sigma1 * zP) * norm_Fn;

                        if (typeid(T) == typeid(doublereal)) {
			for (int j = 1; j <= 2; ++j) {
				i->uP(j) = dGetValue(uP(j));
                                i->tau(j) = dGetValue(tau(j));
                            }
			}

			const T Norm_uP = Norm(uP);

			if (Norm_uP == 0.) {
				kappa = 0.;
			} else {
				const T a0 = Norm(Mk2 * uP);
				const T a1 = a0 / Norm(Mk * uP);
				const T g = a1 + (Norm(Ms2 * uP) / Norm(Ms * uP) - a1) * exp(-pow(Norm_uP / vs, gamma));

				kappa = a0 / g;
			}

			Phi = (inv_Mk_sigma0 * VVec2(uP - inv_Mk * z * kappa) - zP) * dStictionStateEquScale;
		}

		const VVec3 R2_Rt2_e3_n = R2 * VVec3(i->Rt2.GetCol(3) * n);

		const VVec3 F1 = R2 * VVec3(i->Rt2 * VVec3(-tau(1), -tau(2), norm_Fn));
		const VVec3 M1 = -Cross(R2_Rt2_e3_n, F1);
		const VVec3 F2 = -F1;
		const VVec3 M2 = Cross(VVec3(dX - R2_Rt2_e3_n), F2);

		const integer iFirstMomIndexNode2 = i->pNode2->iGetFirstMomentumIndex();

		WorkVec.AddItem(iFirstMomIndexNode1 + 1, F1);
		WorkVec.AddItem(iFirstMomIndexNode1 + 4, M1);
		WorkVec.AddItem(iFirstMomIndexNode2 + 1, F2);
		WorkVec.AddItem(iFirstMomIndexNode2 + 4, M2);

		if (bEnableFriction) {
			for (int j = 1; j <= 2; ++j) {
				WorkVec.AddItem(iFirstIndex + j + offset, Phi(j));
			}
		}
	}
}

inline void BallBearingContact::CheckTimeStep(Washer& w, doublereal Fn, doublereal d, doublereal dn_dt)
{
	tCurr = pDM->dGetTime();
	w.dCurr = d;
	w.dd_dtCurr = -dn_dt;
	w.FnCurr = Fn;

	if (bFirstRes) {
		bFirstRes = false;
		return;
	}

	if (w.FnPrev < dFMin && w.FnCurr < dFMin) {
		// Bypass time step control
		// because the force is considered
		// to be too small!
		return;
	}

	doublereal a;

	if ((w.dCurr > 0 && w.dPrev < 0) || (w.dCurr < 0 && w.dPrev > 0)) {
		a = 0.;
	} else if (w.dCurr >= 0. && w.dPrev >= 0.) {
		const doublereal FnCurr = k * pow(w.dCurr, 3./2.);
		const doublereal FnPrev = k * pow(w.dPrev, 3./2.);
		const doublereal b = (FnPrev + copysign(dFMax, FnCurr - FnPrev)) / k;

		if (b >= 0.) {
			a = std::pow(b, 2. / 3.);
		} else {
			a = 0.;
		}
	} else {
		return;
	}

	const doublereal dd_dt2 = (w.dd_dtCurr - w.dd_dtPrev) / (tCurr - tPrev);
	const doublereal p = 2 * w.dd_dtPrev / dd_dt2;
	const doublereal q = 2 * (w.dPrev - a) / dd_dt2;
	const doublereal r = 0.25 * p * p - q;

	if (r >= 0) {
		const doublereal dt[2] = {
			-0.5 * p + sqrt(r),
			-0.5 * p - sqrt(r)
		};

		for (int i = 0; i < 2; ++i) {
			if (dt[i] > 0. && dt[i] < dtMax) {
				dtMax = dt[i];
			}
		}
	}
}

inline void BallBearingContact::CheckTimeStep(Washer& w, const Gradient<iNumADVars>&, const Gradient<iNumADVars>&, const Gradient<iNumADVars>&)
{
	// Do nothing
}

void BallBearingContact::AfterConvergence(const VectorHandler& X,
										  const VectorHandler& XP)
{
	tPrev = tCurr;

	for (std::vector<Washer>::iterator i = washers.begin(); i != washers.end(); ++i) {
		i->dPrev = i->dCurr;
		i->dd_dtPrev = i->dd_dtCurr;
		i->FnPrev = i->FnCurr;
	}
}

unsigned int
BallBearingContact::iGetNumPrivData(void) const
{
	return 1u + washers.size() * iNumPrivData;
}

unsigned int BallBearingContact::iGetPrivDataIdx(const char *s) const
{
	if (0 == strcmp(s, "max" "dt")) {
		return 1u;
	} else {
		for (int i = 0; i < iNumPrivData; ++i) {
			unsigned iWasher;

			if (1 != sscanf(s, rgPrivData[i].szPattern, &iWasher)) {
				continue;
			}

			if (iWasher <= 0 || iWasher > washers.size()) {
				return 0;
			}

			return (iWasher - 1) * iNumPrivData + i + 2u;
		}

		return 0;
	}
}

doublereal BallBearingContact::dGetPrivData(unsigned int i) const
{
	if (i == 1u) {
		return dtMax;
	}

	const div_t d = div(i - 2u, iNumPrivData);

	if (d.quot < 0 || size_t(d.quot) >= washers.size()) {
		silent_cerr("ballbearingcontact(" << GetLabel()
				<< "): invalid index " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	const Washer& w = washers[d.quot];

	switch (d.rem) {
	case 0:
		return w.dCurr;

	case 1:
		return w.dd_dtCurr;

	case 2:
		return w.FnCurr;

	case 3:
	case 4:
		return w.z(d.rem - 2);

	case 5:
	case 6:
		return w.zP(d.rem - 4);

	case 7:
	case 8:
		return w.uP(d.rem - 6);

        case 9:
        case 10:
            return w.tau(d.rem - 8);
            
	default:
		silent_cerr("ballbearingcontact(" << GetLabel()
				<< "): invalid index " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

int
BallBearingContact::iGetNumConnectedNodes(void) const
{
	return washers.size() + 1;
}

void
BallBearingContact::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.reserve(iGetNumConnectedNodes());
	connectedNodes.clear();
	connectedNodes.push_back(pNode1);

	for (std::vector<Washer>::const_iterator i = washers.begin(); i != washers.end(); ++i) {
		connectedNodes.push_back(i->pNode2);
	}
}

void
BallBearingContact::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{

}

std::ostream&
BallBearingContact::Restart(std::ostream& out) const
{
	return out;
}

unsigned int BallBearingContact::iGetNumDof(void) const
{
	return bEnableFriction ? 2 * washers.size() : 0;
}

std::ostream& BallBearingContact::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	if (bEnableFriction) {
		const integer iIndex = iGetFirstIndex();

		for (size_t i = 0; i < washers.size(); ++i) {
			for (int j = 1; j <= 2; ++j) {
				out << prefix << iIndex + j + 2 * i << ": stiction state z" << j << "[" << i << "]" << std::endl;
			}
		}
	}

	return out;
}

std::ostream& BallBearingContact::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	if (bEnableFriction) {
		const integer iIndex = iGetFirstIndex();

		for (size_t i = 0; i < washers.size(); ++i) {
			for (int j = 1; j <= 2; ++j) {
				out << prefix << iIndex + j + 2 * i << ": ode for stiction state z" << j << "[" << i << "]" << std::endl;
			}
		}
	}

	return out;
}

DofOrder::Order BallBearingContact::GetDofType(unsigned int) const
{
	return DofOrder::DIFFERENTIAL;
}

DofOrder::Order BallBearingContact::GetEqType(unsigned int i) const
{
	return DofOrder::DIFFERENTIAL;
}

unsigned int
BallBearingContact::iGetInitialNumDof(void) const
{
	return 0;
}

void
BallBearingContact::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
BallBearingContact::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
BallBearingContact::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}
#endif

bool ballbearing_contact_set(void)
{
#ifdef USE_AUTODIFF
	UserDefinedElemRead *rf = new UDERead<BallBearingContact>;

	if (!SetUDE("ball" "bearing" "contact", rf))
	{
		delete rf;
		return false;
	}

	return true;
#else
	return false;
#endif
}

#ifndef STATIC_MODULES

extern "C"
{

int module_init(const char *module_name, void *pdm, void *php)
{
	if (!ballbearing_contact_set())
	{
		silent_cerr("ballbearing_contact: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

}

#endif // ! STATIC_MODULE


