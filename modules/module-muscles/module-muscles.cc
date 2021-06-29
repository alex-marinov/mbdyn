/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * Author: Andrea Zanoni <andrea.zanoni@polimi.it>
 *         Pierangelo Masarati <masarati@aero.polimi.it>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "constltp_impl.h"

class MusclePennestriCL
: public ElasticConstitutiveLaw<doublereal, doublereal> {
protected:
	doublereal Li;
	doublereal L0;
	doublereal V0;
	doublereal F0;
	DriveOwner Activation;
	bool bActivationOverflow;
	bool bActivationOverflowWarn;
	doublereal a;
	doublereal aReq;	// requested activation: for output only
	virtual std::ostream& Restart_int(std::ostream& out) const {
		return out;
	};
#ifdef USE_NETCDF
	MBDynNcVar Var_dAct;
	MBDynNcVar Var_dActReq;
#endif // USE_NETCDF
public:
	MusclePennestriCL(const TplDriveCaller<doublereal> *pTplDC, doublereal dPreStress,
		doublereal Li, doublereal L0, doublereal V0, doublereal F0,
		const DriveCaller *pAct, bool bActivationOverflow, bool bActivationOverflowWarn)
	: ElasticConstitutiveLaw<doublereal, doublereal>(pTplDC, dPreStress),
	Li(Li), L0(L0), V0(V0), F0(F0),
	Activation(pAct), bActivationOverflow(bActivationOverflow), bActivationOverflowWarn(bActivationOverflowWarn)
#ifdef USE_NETCDFC
	, Var_dAct(0), Var_dActReq(0)
#endif // USE_NETCDFC
	{
		NO_OP;
	};

	virtual ~MusclePennestriCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriCL,
			MusclePennestriCL(pGetDriveCaller()->pCopy(),
				PreStress,
				Li, L0, V0, F0,
				Activation.pGetDriveCaller()->pCopy(),
				bActivationOverflow, 
				bActivationOverflowWarn));
		return pCL;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		out << "muscle pennestri"
			", initial length, " << Li
			<< ", reference length, " << L0
			<< ", reference velocity, " << V0
			<< ", reference force, " << F0
			<< ", activation, ", Activation.pGetDriveCaller()->Restart(out)
			<< ", activation check, " << bActivationOverflow
			<< ", warn, " << bActivationOverflowWarn;
			Restart_int(out)
			<< ", ", ElasticConstitutiveLaw<doublereal, doublereal>::Restart_int(out);
		return out;
	};

	virtual std::ostream& OutputAppend(std::ostream& out) const {
		return out << " " << a << " " << aReq;
	};

	virtual void NetCDFOutputAppend(OutputHandler& OH) const {
#ifdef USE_NETCDF
		OH.WriteNcVar(Var_dAct, a);
		OH.WriteNcVar(Var_dActReq, aReq);
#endif // USE_NETCDF
	}

	virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name) {
#ifdef USE_NETCDF
		ASSERT(OH.IsOpen(OutputHandler::NETCDF));
		if (OH.UseNetCDF(OutputHandler::LOADABLE)) 
		{
			Var_dAct = OH.CreateVar<doublereal>(name + ".a", 
				OutputHandler::Dimensions::Dimensionless, 
				"Muscular activation (effective value)");
			Var_dActReq = OH.CreateVar<doublereal>(name + ".aReq",  
				OutputHandler::Dimensions::Dimensionless,
				"Requested muscular activation");
		}
#endif // USE_NETCDF
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		aReq = Activation.dGet();
		a = aReq;
		if (a < 0.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation underflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 0.;
			}

		} else if (a > 1.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation overflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 1.;
			}
		}

		doublereal dxdEps = Li/L0;
		doublereal dvdEpsPrime = Li/V0;
		doublereal x = (1. + ConstitutiveLaw<doublereal, doublereal>::Epsilon)*dxdEps;
		doublereal v = ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime*dvdEpsPrime;

		doublereal f1 = std::exp(std::pow(x - 0.95, 2) - 40*std::pow(x - 0.95, 4));
		doublereal f2 = 1.6 - 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1., 4));
		doublereal f3 = 1.3*std::atan(0.1*std::pow(x - 0.22, 10));

		doublereal df1dx = f1*(2*(x - 0.95) - 4*40.*std::pow(x - 0.95, 3));
		doublereal df2dv = 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1, 4))*(2*0.1/std::pow(v - 1., 3) - 4*1.1/std::pow(v - 1., 5));
		doublereal df3dx = 1.3*std::pow(x - 0.22, 9)/(0.01*std::pow(x - 0.22, 20) + 1);

		ConstitutiveLaw<doublereal, doublereal>::F = PreStress + F0*(f1*f2*a + f3);
		ConstitutiveLaw<doublereal, doublereal>::FDE = F0*(df1dx*f2*a + df3dx)*dxdEps;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = F0*f1*df2dv*a*dvdEpsPrime;
	};
};

class MusclePennestriErgoCL
: public MusclePennestriCL {
public:
	MusclePennestriErgoCL(const TplDriveCaller<doublereal> *pTplDC, doublereal dPreStress,
		doublereal Li, doublereal L0, doublereal V0, doublereal F0,
		const DriveCaller *pAct, bool bActivationOverflow, bool bActivationOverflowWarn)
	: MusclePennestriCL(pTplDC, dPreStress, Li, L0, V0, F0, pAct, bActivationOverflow, bActivationOverflowWarn)
	{
		NO_OP;
	};

	virtual ~MusclePennestriErgoCL(void) {
		NO_OP;
	};

	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::ELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriErgoCL,
			MusclePennestriErgoCL(pGetDriveCaller()->pCopy(),
				PreStress,
				Li, L0, V0, F0,
				Activation.pGetDriveCaller()->pCopy(),
				bActivationOverflow, 
				bActivationOverflowWarn));
		return pCL;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		MusclePennestriCL::Update(Eps, 0.);
	};

protected:
	virtual std::ostream& Restart_int(std::ostream& out) const {
		out << ", ergonomy, yes";
		return out;
	};
};

class MusclePennestriReflexiveCL
: public MusclePennestriCL {
protected:
	DriveOwner Kp;
	DriveOwner Kd;
	DriveOwner ReferenceLength;
	
public:
	MusclePennestriReflexiveCL(const TplDriveCaller<doublereal> *pTplDC, doublereal dPreStress,
		doublereal Li, doublereal L0, doublereal V0, doublereal F0,
		const DriveCaller *pAct, bool bActivationOverflow, bool bActivationOverflowWarn,
		const DriveCaller *pKp, const DriveCaller *pKd, const DriveCaller *pReferenceLength)
	: MusclePennestriCL(pTplDC, dPreStress, Li, L0, V0, F0, pAct, bActivationOverflow, bActivationOverflowWarn),
	Kp(pKp), Kd(pKd), ReferenceLength(pReferenceLength)
	{
		NO_OP;
	};

	virtual ~MusclePennestriReflexiveCL(void) {
		NO_OP;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriReflexiveCL,
			MusclePennestriReflexiveCL(pGetDriveCaller()->pCopy(),
				PreStress,
				Li, L0, V0, F0,
				Activation.pGetDriveCaller()->pCopy(),
				bActivationOverflow,
				bActivationOverflowWarn,
				Kp.pGetDriveCaller()->pCopy(), Kd.pGetDriveCaller()->pCopy(),
				ReferenceLength.pGetDriveCaller()->pCopy()));
		return pCL;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal dxdEps = Li/L0;
		doublereal dvdEpsPrime = Li/V0;
		doublereal x = (1. + ConstitutiveLaw<doublereal, doublereal>::Epsilon)*dxdEps;
		doublereal v = ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime*dvdEpsPrime;

		doublereal dxRef = ReferenceLength.dGet()/L0;

		doublereal aRef = Activation.dGet();
		aReq = aRef + Kp.dGet()*(x - dxRef) + Kd.dGet()*v;
		a = aReq;

		if (aReq < 0.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation underflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 0.;
			}

		} else if (aReq > 1.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation overflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 1.;
			}
		}

		doublereal f1 = std::exp(std::pow(x - 0.95, 2) - 40*std::pow(x - 0.95, 4));
		doublereal f2 = 1.6 - 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1., 4));
		doublereal f3 = 1.3*std::atan(0.1*std::pow(x - 0.22, 10));

		doublereal df1dx = f1*(2*(x - 0.95) - 4*40.*std::pow(x - 0.95, 3));
		doublereal df2dv = 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1, 4))*(2*0.1/std::pow(v - 1., 3) - 4*1.1/std::pow(v - 1., 5));
		doublereal df3dx = 1.3*std::pow(x - 0.22, 9)/(0.01*std::pow(x - 0.22, 20) + 1);

		ConstitutiveLaw<doublereal, doublereal>::F = PreStress + F0*(f1*f2*a + f3);
		ConstitutiveLaw<doublereal, doublereal>::FDE = F0*((df1dx*aRef + f1*Kp.dGet())*f2 + df3dx)*dxdEps;
		ConstitutiveLaw<doublereal, doublereal>::FDEPrime = F0*f1*(df2dv*aRef + f2*Kd.dGet())*dvdEpsPrime;
	};
	 
protected:
	virtual std::ostream& Restart_int(std::ostream& out) const {
		out
			<< ", reflexive"
			<< ", proportional gain, ", Kp.pGetDriveCaller()->Restart(out)
			<< ", derivative gain, ", Kd.pGetDriveCaller()->Restart(out)
			<< ", reference length, ", ReferenceLength.pGetDriveCaller()->Restart(out);
		return out;
	};
};

class MusclePennestriReflexiveCLWithSRS
: public MusclePennestriCL {
protected:
	DriveOwner Kp;
	DriveOwner Kd;
	DriveOwner ReferenceLength;
	doublereal SRSGamma;
	doublereal SRSDelta;

public:
	MusclePennestriReflexiveCLWithSRS(const TplDriveCaller<doublereal> *pTplDC, doublereal dPreStress,
		doublereal Li, doublereal L0, doublereal V0, doublereal F0,
		const DriveCaller *pAct, bool bActivationOverflow, bool bActivationOverflowWarn,
		const DriveCaller *pKp, const DriveCaller *pKd, const DriveCaller *pReferenceLength,
		const doublereal SRSGamma, const doublereal SRSDelta)
	: MusclePennestriCL(pTplDC, dPreStress, Li, L0, V0, F0, pAct, bActivationOverflow, bActivationOverflowWarn),
	Kp(pKp), Kd(pKd), ReferenceLength(pReferenceLength), SRSGamma(SRSGamma), SRSDelta(SRSDelta)
	{
		NO_OP;
	};

	virtual ~MusclePennestriReflexiveCLWithSRS(void) {
		NO_OP;
	}

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		// pass parameters to copy constructor
		SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriReflexiveCLWithSRS,
			MusclePennestriReflexiveCLWithSRS(pGetDriveCaller()->pCopy(),
				PreStress,
				Li, L0, V0, F0,
				Activation.pGetDriveCaller()->pCopy(),
				bActivationOverflow,
				bActivationOverflowWarn,
				Kp.pGetDriveCaller()->pCopy(), Kd.pGetDriveCaller()->pCopy(),
				ReferenceLength.pGetDriveCaller()->pCopy(),
				SRSGamma, SRSDelta));
		return pCL;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime) {
		ConstitutiveLaw<doublereal, doublereal>::Epsilon = Eps - ElasticConstitutiveLaw<doublereal, doublereal>::Get();
		ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime = EpsPrime;

		doublereal dxdEps = Li/L0;
		doublereal dvdEpsPrime = Li/V0;
		doublereal x = (1. + ConstitutiveLaw<doublereal, doublereal>::Epsilon)*dxdEps;
		doublereal v = ConstitutiveLaw<doublereal, doublereal>::EpsilonPrime*dvdEpsPrime;

		doublereal dxRef = ReferenceLength.dGet()/L0;

		doublereal aRef = Activation.dGet();
		aReq = aRef + Kp.dGet()*(x - dxRef) + Kd.dGet()*v;
		a = aReq;
		
		if (aReq < 0.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation underflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 0.;
			}

		} else if (aReq > 1.) {
			if (bActivationOverflowWarn) {
				silent_cerr("MusclePennestriCL: activation overflow (aReq=" << aReq << ")" << std::endl);
			}
			if (bActivationOverflow) {
				a = 1.;
			}
		}

		doublereal f1 = std::exp(std::pow(x - 0.95, 2) - 40*std::pow(x - 0.95, 4));
		doublereal f2 = 1.6 - 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1., 4));
		doublereal f3 = 1.3*std::atan(0.1*std::pow(x - 0.22, 10));

		doublereal df1dx = f1*(2*(x - 0.95) - 4*40.*std::pow(x - 0.95, 3));
		doublereal df2dv = 1.6*std::exp(0.1/std::pow(v - 1., 2) - 1.1/std::pow(v - 1, 4))*(2*0.1/std::pow(v - 1., 3) - 4*1.1/std::pow(v - 1., 5));
		doublereal df3dx = 1.3*std::pow(x - 0.22, 9)/(0.01*std::pow(x - 0.22, 20) + 1);

		if ( (x - dxRef) > 0 ) {
			if ( (x - dxRef) < SRSDelta ) {
				ConstitutiveLaw<doublereal, doublereal>::F = PreStress + F0*(f1*(f2*a + SRSGamma*aRef*(x - dxRef)) + f3);
				ConstitutiveLaw<doublereal, doublereal>::FDE = 
					F0*(df1dx*aRef*(f2 + SRSGamma*(x - dxRef)) + f1*(SRSGamma*aRef + Kp.dGet()*f2) + df3dx)*dxdEps;
			} else {
				ConstitutiveLaw<doublereal, doublereal>::F = PreStress + F0*(f1*(f2*a + SRSGamma*aRef*SRSDelta) + f3);
				ConstitutiveLaw<doublereal, doublereal>::FDE = 
					F0*(df1dx*aRef*(f2 + SRSGamma*SRSDelta) + f1*(SRSGamma*aRef + Kp.dGet()*f2) + df3dx)*dxdEps;
			}
		} else {
			// No SRS contribution
			ConstitutiveLaw<doublereal, doublereal>::F = PreStress + F0*(f1*f2*a + f3);
			ConstitutiveLaw<doublereal, doublereal>::FDE = F0*((df1dx*aRef + f1*Kp.dGet())*f2 + df3dx)*dxdEps;
			ConstitutiveLaw<doublereal, doublereal>::FDEPrime = F0*f1*(df2dv*aRef + f2*Kd.dGet())*dvdEpsPrime;
		}
	};
};
		

/* specific functional object(s) */
struct MusclePennestriCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		if (HP.IsKeyWord("help")) {
			silent_cerr("MusclePennestriCL:\n"
				"        muscle Pennestri ,\n"
				"                [ initial length , <Li> , ]\n"
				"                reference length , <L0> ,\n"
				"                [ reference velocity , <V0> , ]\n"
				"                reference force , <F0> ,\n"
				"                activation , (DriveCaller) <activation>\n"
				"                [ , activation check , (bool)<activation_check> ]\n"
				"		 [ , warn , (bool)<warn_user> ]\n"
				"                [ , ergonomy , { yes | no } , ]\n"
				"                [ , reflexive , # only when ergonomy == no\n"
				"                        proportional gain , (DriveCaller) <kp> ,\n"
				"                        derivative gain , (DriveCaller) <kd> ,\n"
				"                        reference length, (DriveCaller) <lref> ]\n"
				"		  	[ , short range stiffness ]\n"
				" 		  	[ 	, gamma, (real) <gamma> ,]\n"
				"		   	[ 	, delta, (real) <delta> ,]\n"
				"                [ , prestress, <prestress> ]\n"
				"                [ , prestrain, (DriveCaller) <prestrain> ]\n"
				<< std::endl);

			if (!HP.IsArg()) {
				throw NoErr(MBDYN_EXCEPT_ARGS);
			}
		}

		bool bErgo(false);
		bool bGotErgo(false);
		if (HP.IsKeyWord("ergonomy")) {
			silent_cerr("MusclePennestriCL: deprecated, \"ergonomy\" "
					"at line " << HP.GetLineData()
					<< " should be at end of definition" << std::endl);
			bErgo = HP.GetYesNoOrBool(bErgo);
			bGotErgo = true;
		}

		doublereal Li = -1.;
		if (HP.IsKeyWord("initial" "length")) {
			Li = HP.GetReal();
			if (Li <= 0.) {
				silent_cerr("MusclePennestriCL: null or negative initial length "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		if (!HP.IsKeyWord("reference" "length")) {
			silent_cerr("MusclePennestriCL: \"reference length\" expected "
				"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal L0 = HP.GetReal();
		if (L0 <= 0.) {
			silent_cerr("MusclePennestriCL: null or negative reference length "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		// default (mm/s? m/s?)
		doublereal V0 = 2.5;
		if (HP.IsKeyWord("reference" "velocity")) {
			V0 = HP.GetReal();
			if (V0 <= 0.) {
				silent_cerr("MusclePennestriCL: null or negative reference velocity "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		if (!HP.IsKeyWord("reference" "force")) {
			silent_cerr("MusclePennestriCL: \"reference force\" expected "
				"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		doublereal F0 = HP.GetReal();
		if (F0 <= 0.) {
			silent_cerr("MusclePennestriCL: null or negative reference force "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (!HP.IsKeyWord("activation")) {
			silent_cerr("MusclePennestriCL: \"activation\" expected "
				"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		const DriveCaller *pAct = HP.GetDriveCaller();

		bool bActivationOverflow(false);
		if (HP.IsKeyWord("activation" "check")) {
			bActivationOverflow = HP.GetYesNoOrBool(bActivationOverflow);
		}

		bool bActivationOverflowWarn(false);
		if (HP.IsKeyWord("warn")) {
			bActivationOverflow = HP.GetYesNoOrBool(bActivationOverflowWarn);
		}

		// FIXME: "ergonomy" must be here
		if (!bGotErgo && HP.IsKeyWord("ergonomy")) {
			bErgo = HP.GetYesNoOrBool(bErgo);
		}

		bool bReflexive(false);
		DriveCaller* pKp(NULL);
		DriveCaller* pKd(NULL);
		const DriveCaller *pRefLen(0);
		// Short Range Stiffness
		bool bSRS(false);
		// Default values, from De Groote et Al, JBiomech 55 (2017)
		doublereal dSRSGamma = 280.; 
		doublereal dSRSDelta = 5.7e-3;
		if (HP.IsKeyWord("reflexive")) {
			if (bErgo) {
				silent_cerr("MusclePennestriCL: "
					"\"reflexive\" and \"ergonomy\" incompatible "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			bReflexive = true;

			if (!HP.IsKeyWord("proportional" "gain")) {
				silent_cerr("MusclePennestriCL: \"proportional gain\" expected "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pKp = HP.GetDriveCaller();

			if (!HP.IsKeyWord("derivative" "gain")) {
				silent_cerr("MusclePennestriCL: \"derivative gain\" expected "
					"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pKd = HP.GetDriveCaller();

			if (!HP.IsKeyWord("reference" "length")) {
				silent_cerr("MusclePennestriCL: \"reference length\" expected "
					"at line " << HP.GetLineData() << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			pRefLen = HP.GetDriveCaller();
		
			if (HP.IsKeyWord("short" "range" "stiffness")) {
				bSRS = true;
				
				if (HP.IsKeyWord("gamma")) {
					dSRSGamma = HP.GetReal();
				}

				if (HP.IsKeyWord("delta")) {
					dSRSDelta = HP.GetReal();
				}
			}
		}


		/* Prestress and prestrain */
		doublereal PreStress(0.);
		GetPreStress(HP, PreStress);
		TplDriveCaller<doublereal> *pTplDC = GetPreStrain<doublereal>(pDM, HP);

		if (Li == -1.) {
			Li = L0;
		}

		if (bErgo) {
			SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriErgoCL,
				MusclePennestriErgoCL(pTplDC, PreStress,
					Li, L0, V0, F0, pAct,
					bActivationOverflow, bActivationOverflowWarn));

		} else if (bReflexive) {
			if (bSRS) {
				SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriReflexiveCLWithSRS,
					MusclePennestriReflexiveCLWithSRS(pTplDC, PreStress,
						Li, L0, V0, F0, pAct,
						bActivationOverflow,
						bActivationOverflowWarn,
						pKp, pKd, pRefLen, dSRSGamma, dSRSDelta));
			} else {
				SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriReflexiveCL,
					MusclePennestriReflexiveCL(pTplDC, PreStress,
						Li, L0, V0, F0, pAct,
						bActivationOverflow,
						bActivationOverflowWarn,
						pKp, pKd, pRefLen));
			}
		} else {
			SAFENEWWITHCONSTRUCTOR(pCL, MusclePennestriCL,
				MusclePennestriCL(pTplDC, PreStress,
					Li, L0, V0, F0, pAct,
					bActivationOverflow, 
					bActivationOverflowWarn));
		}

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<doublereal, doublereal> *rf1D = new MusclePennestriCLR;
	if (!SetCL1D("muscle" "pennestri", rf1D)) {
		delete rf1D;

		silent_cerr("MusclePennestriCL: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

