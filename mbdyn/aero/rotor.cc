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

/* Rotore */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>

#include <sstream>

#ifdef USE_MPI
#include "mysleep.h"
const int mysleeptime = 300;

#ifdef MPI_PROFILING
extern "C" {
#include <mpe.h>
#include <stdio.h>
}
#endif /* MPI_PROFILING */
#include "mbcomm.h"
#endif /* USE_MPI */

#include "rotor.h"
#include "dataman.h"

static const doublereal dVTipTreshold = 1e-6;

/* Rotor - begin */

Rotor::Rotor(unsigned int uL, const DofOwner* pDO)
: Elem(uL, flag(0)),
InducedVelocityElem(uL, pDO, 0, 0, flag(0)),
pRotor(0), pGround(0),
dOmegaRef(0.), dRadius(0), dVTipRef(0.), dArea(0.),
dUMean(0.), dUMeanRef(0.), dUMeanPrev(0.),
iMaxIter(0),
iCurrIter(0),
dTolerance(0),
dEta(0),
bUMeanRefConverged(false),
Weight(), dWeight(0.),
dHoverCorrection(1.), dForwardFlightCorrection(1.),
RRotTranspose(::Zero3x3), RRot(::Eye3), RRot3(::Zero3),
VCraft(::Zero3),
dPsi0(0.), dSinAlphad(1.), dCosAlphad(0.),
dMu(0.), dLambda(1.), dChi(0.),
dVelocity(0.), dOmega(0.),
iNumSteps(0)
{
	NO_OP;
}

Rotor::Rotor(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, const StructNode* pG,
	ResForceSet **ppres,
	const doublereal& dR,
	unsigned int iMaxIt, const doublereal& dTol, const doublereal& dE,
	flag fOut)
: Elem(uL, fOut),
InducedVelocityElem(uL, pDO, 0, 0, flag(0)),
pRotor(0), pGround(0),
dOmegaRef(0.), dRadius(0), dVTipRef(0.), dArea(0.),
dUMean(0.), dUMeanRef(0.), dUMeanPrev(0.),
iMaxIter(0),
iCurrIter(0),
dTolerance(0),
dEta(0),
bUMeanRefConverged(false),
Weight(), dWeight(0.),
dHoverCorrection(1.), dForwardFlightCorrection(1.),
RRotTranspose(::Zero3x3), RRot(::Eye3), RRot3(::Zero3),
VCraft(::Zero3),
dPsi0(0.), dSinAlphad(1.), dCosAlphad(0.),
dMu(0.), dLambda(1.), dChi(0.),
dVelocity(0.), dOmega(0.),
iNumSteps(0)
{
	Init(pC, rrot, pR, pG, ppres, dR, iMaxIt, dTol, dE, fOut);
}

void
Rotor::Init(const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, const StructNode* pG,
	ResForceSet **ppres,
	const doublereal& dR,
	unsigned int iMaxIt, const doublereal& dTol, const doublereal& dE,
	flag fOut)
{
	InducedVelocityElem::pCraft = pC;
	InducedVelocityElem::ppRes = ppres;

	pRotor = pR;
	pGround = pG;
	dRadius = dR;
	iMaxIter = iMaxIt;
	dTolerance = dTol;
	dEta = dE;
	RRot = rrot;

	SetOutputFlag(fOut);

	Vec3 R3C(pCraft->GetRCurr()*RRot.GetVec(3));
	Vec3 R3R(pRotor->GetRCurr().GetVec(3));
	if (R3C.Dot(R3R) < 1. - std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("Rotor(" << GetLabel() << "): "
			"warning, possible misalignment "
			"of rotor node StructNode(" << pRotor->GetLabel() << ") "
			"axis Rr(3) = {" << R3R << "} "
			"and craft node StructNode(" << pCraft->GetLabel() << ") "
			"axis Rc*Rh(3) = {" << R3C << "} "
			<< "for Rotor(" << GetLabel() << ")"
			<< std::endl);
	}
}

Rotor::~Rotor(void)
{
	NO_OP;
}

void
Rotor::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	/* non mi ricordo a cosa serve! */
	iNumSteps++;

	/* updates the umean at the previous step */
	dUMeanPrev = dUMean;

	/*
	* FIXME: should go in AfterPredict ...
	* this way there's a one step delay in the weight value
	*/
	if (Weight.pGetDriveCaller() != 0) {
		dWeight = Weight.dGet();
		if (dWeight < 0.) {
			silent_cout("Rotor(" << GetLabel() << "): "
				"delay < 0.0; using 0.0" << std::endl);
			dWeight = 0.;

		} else if (dWeight > 1.) {
			silent_cout("Rotor(" << GetLabel() << "): "
				"delay > 1.0; using 1.0" << std::endl);
			dWeight = 1.;
		}
	}

	InducedVelocity::AfterConvergence(X, XP);
}

void
Rotor::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::ROTORS)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));
			std::ostringstream os;
			os << "elem.inducedvelocity." << GetLabel() << ".";
			std::string name = os.str();
			(void)OH.CreateVar(name, "Rotor");
			Var_f = OH.CreateVar<Vec3>(name + "f",
					OutputHandler::Dimensions::Force,
					"rotor force in x, y and z directions (lon, lat, thrust)");
			Var_m = OH.CreateVar<Vec3>(name + "m",
					OutputHandler::Dimensions::Moment,
					"rotor moment about x, y and z directions (pitch, roll, torque)");
			Var_dUMean = OH.CreateVar<doublereal>(name + "UMean",
					OutputHandler::Dimensions::Velocity,
					"mean inflow velocity");

			Var_dVelocity = OH.CreateVar<doublereal>(name + "VRef",
					OutputHandler::Dimensions::Velocity,
					"reference velocity (craft_node + airstream)");
			Var_dAlpha = OH.CreateVar<doublereal>(name + "Alpha",
					OutputHandler::Dimensions::rad,
					"rotor disk angle");
			Var_dMu = OH.CreateVar<doublereal>(name + "Mu",
					OutputHandler::Dimensions::Dimensionless,
					"advance parameter");
			Var_dLambda = OH.CreateVar<doublereal>(name + "Lambda",
					OutputHandler::Dimensions::Dimensionless,
					"inflow parameter");
			Var_dChi = OH.CreateVar<doublereal>(name + "Chi",
					OutputHandler::Dimensions::Dimensionless,
					"advance/inflow parameter");
			Var_dPsi0 = OH.CreateVar<doublereal>(name + "Psi0",
					OutputHandler::Dimensions::rad,
					"reference azimuthal direction");
			Var_bUMeanRefConverged = OH.CreateVar<integer>(name + "UMeanRefConverged",
					OutputHandler::Dimensions::Boolean,
					"boolean flag indicating reference induced velocity computation convergence");
			Var_iCurrIter = OH.CreateVar<integer>(name + "Iter",
					OutputHandler::Dimensions::Dimensionless,
					"number of iterations required for convergence");

		}
#endif // USE_NETCDF
	}
}

void
Rotor::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
#ifdef USE_MPI
		if (is_parallel && IndVelComm.Get_size() > 1) {
			if (IndVelComm.Get_rank() == 0) {
				Vec3 TmpF(pTmpVecR);
				Vec3 TmpM(pTmpVecR+3);
#ifdef USE_NETCDF
				if (OH.UseNetCDF(OutputHandler::ROTORS)) {
					OH.WriteNcVar(Var_f, RRotTranspose*TmpF);
					OH.WriteNcVar(Var_m, RRotTranspose*TmpM);
					OH.WriteNcVar(Var_dUMean, dUMean);
					OH.WriteNcVar(Var_dVelocity, dVelocity);
					OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
					OH.WriteNcVar(Var_dMu, dMu);
					OH.WriteNcVar(Var_dLambda, dLambda);
					OH.WriteNcVar(Var_dChi, dChi);
					OH.WriteNcVar(Var_dPsi0, dPsi0);
					OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
					OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
				}
#endif // USE_NETCDF
				if (OH.UseText(OutputHandler::ROTORS)) {
					OH.Rotors()
						<< std::setw(8) << GetLabel()	/* 1 */
						<< " " << RRotTranspose*TmpF /* 2-4 */
						<< " " << RRotTranspose*TmpM /* 5-7 */
						<< " " << dUMean 	/* 8 */
						<< " " << dVelocity	/* 9 */
						<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
						<< " " << dMu		/* 11 */
						<< " " << dLambda	/* 12 */
						<< " " << dChi		/* 13 */
						<< " " << dPsi0		/* 14 */
						<< " " << bUMeanRefConverged /* 15 */
						<< " " << iCurrIter	/* 16 */
						<< std::endl;

					for (int i = 0; ppRes && ppRes[i]; i++) {
						Vec3 TmpF(pTmpVecR+6+6*i);
						Vec3 TmpM(pTmpVecR+9+6*i);

						OH.Rotors()
							<< std::setw(8) << GetLabel()
							<< ":" << ppRes[i]->GetLabel()
							<< " " << TmpF
							<< " " << TmpM
							<< std::endl;
					}
				}
			}
		} else {
#ifdef USE_NETCDF
			if (OH.UseNetCDF(OutputHandler::ROTORS)) {
				OH.WriteNcVar(Var_f, RRotTranspose*Res.Force());
				OH.WriteNcVar(Var_m, RRotTranspose*Res.Moment());
				OH.WriteNcVar(Var_dUMean, dUMean);
				OH.WriteNcVar(Var_dVelocity, dVelocity);
				OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
				OH.WriteNcVar(Var_dMu, dMu);
				OH.WriteNcVar(Var_dLambda, dLambda);
				OH.WriteNcVar(Var_dChi, dChi);
				OH.WriteNcVar(Var_dPsi0, dPsi0);
				OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
				OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
			}
#endif // USE_NETCDF
			if (OH.UseText(OutputHandler::ROTORS)) {
				OH.Rotors()
					<< std::setw(8) << GetLabel()	/* 1 */
					<< " " << RRotTranspose*Res.Force()  /* 2-4 */
					<< " " << RRotTranspose*Res.Moment() /* 5-7 */
					<< " " << dUMean		/* 8 */
					<< " " << dVelocity		/* 9 */
					<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
					<< " " << dMu			/* 11 */
					<< " " << dLambda		/* 12 */
					<< " " << dChi			/* 13 */
					<< " " << dPsi0			/* 14 */
					<< " " << bUMeanRefConverged	/* 15 */
					<< " " << iCurrIter		/* 16 */
					<< std::endl;

				for (int i = 0; ppRes && ppRes[i]; i++) {
					OH.Rotors()
						<< std::setw(8) << GetLabel()
						<< ":" << ppRes[i]->GetLabel()
						<< " " << ppRes[i]->pRes->Force()
						<< " " << ppRes[i]->pRes->Moment()
						<< std::endl;
				}
			}
		}
#else /* !USE_MPI */
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::ROTORS)) {
			OH.WriteNcVar(Var_f, RRotTranspose*Res.Force());
			OH.WriteNcVar(Var_m, RRotTranspose*Res.Moment());
			OH.WriteNcVar(Var_dUMean, dUMean);
			OH.WriteNcVar(Var_dVelocity, dVelocity);
			OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
			OH.WriteNcVar(Var_dMu, dMu);
			OH.WriteNcVar(Var_dLambda, dLambda);
			OH.WriteNcVar(Var_dChi, dChi);
			OH.WriteNcVar(Var_dPsi0, dPsi0);
			OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
			OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::ROTORS)) {	
			OH.Rotors()
				<< std::setw(8) << GetLabel()	/* 1 */
				<< " " << RRotTranspose*Res.Force()	/* 2-4 */
				<< " " << RRotTranspose*Res.Moment()	/* 5-7 */
				<< " " << dUMean		/* 8 */
				<< " " << dVelocity		/* 9 */
				<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
				<< " " << dMu			/* 11 */
				<< " " << dLambda		/* 12 */
				<< " " << dChi			/* 13 */
				<< " " << dPsi0			/* 14 */
				<< " " << bUMeanRefConverged	/* 15 */
				<< " " << iCurrIter		/* 16 */
				<< std::endl;

			/* FIXME: check for parallel stuff ... */
			for (int i = 0; ppRes && ppRes[i]; i++) {
				OH.Rotors()
					<< std::setw(8) << GetLabel()
					<< ":" << ppRes[i]->GetLabel()
					<< " " << ppRes[i]->pRes->Force()
					<< " " << ppRes[i]->pRes->Moment()
					<< std::endl;
			}
		}
#endif /* !USE_MPI */
	}
}

/* Calcola la posizione azimuthale di un punto generico.
 * X e' il punto di cui e' chiesta la posizione azimuthale,
 *   nel sistema assoluto.
 * Gli viene sottratta la posizione del rotore (corpo attorno al quale avviene
 * la rotazione). Quindi il vettore risultante viene trasformato
 * nel riferimento del mozzo.
 * Si trascura l'eventuale componente fuori del piano, dalle altre due si
 * ricava l'angolo di rotazione relativa, a cui viene sommato l'angolo
 * che il corpo rotore forma con la direzione del vento relativo dPsi0,
 * calcolata in precedenza.
 */
doublereal
Rotor::dGetPsi(const Vec3& X) const
{
	doublereal dr, dp;
	GetPos(X, dr, dp);
	return dp;
}

/* Calcola la distanza di un punto dall'asse di rotazione in coordinate
 * adimensionali */
doublereal
Rotor::dGetPos(const Vec3& X) const
{
	doublereal dr, dp;
	GetPos(X, dr, dp);
	return dr;
}

void
Rotor::GetPos(const Vec3& X, doublereal& dr, doublereal& dp) const
{
	Vec3 XRel(RRotTranspose*(X - Res.Pole()));

	doublereal d1 = XRel(1);
	doublereal d2 = XRel(2);

	doublereal d = sqrt(d1*d1 + d2*d2);

	ASSERT(dRadius > std::numeric_limits<doublereal>::epsilon());
	dr = d/dRadius;

	dp = atan2(d2, d1) - dPsi0;
}

/* Calcola vari parametri geometrici
 * A partire dai corpi che identificano il velivolo ed il rotore
 */
void
Rotor::InitParam(bool bComputeMeanInducedVelocity)
{
	ASSERT(pCraft != 0);
	ASSERT(pRotor != 0);

	/* Trasposta della matrice di rotazione del rotore */
	RRotTranspose = pCraft->GetRCurr()*RRot;
	RRot3 = RRotTranspose.GetVec(3);
	RRotTranspose = RRotTranspose.Transpose();

	/* Posizione del rotore */
	Res.PutPole(pRotor->GetXCurr());

	/* Velocita' angolare del rotore */
	dOmega = (pRotor->GetWCurr() - pCraft->GetWCurr()).Norm();

	/* Velocita' di traslazione del velivolo */
	VCraft = -pRotor->GetVCurr();
	Vec3 VTmp(Zero3);
	if (fGetAirVelocity(VTmp, pRotor->GetXCurr())) {
		VCraft += VTmp;
	}

	/* Velocita' nel sistema del velivolo (del disco?) decomposta */
	VTmp = RRotTranspose*VCraft;
	doublereal dV1 = VTmp(1);
	doublereal dV2 = VTmp(2);
	doublereal dV3 = VTmp(3);
	doublereal dVV = dV1*dV1 + dV2*dV2;
	doublereal dV = sqrt(dVV);
	
	/* Angolo di azimuth 0 del rotore */
	dPsi0 = atan2(dV2, dV1);

	/* Angolo di influsso */
	dVelocity = sqrt(dV3*dV3 + dVV);
	if (dVelocity > dVTipTreshold*dVTipRef) {
		dSinAlphad = -dV3/dVelocity;
		dCosAlphad = dV/dVelocity;

	} else {
		dSinAlphad = 1.;
		dCosAlphad = 0.;
	}

	if (!bComputeMeanInducedVelocity) {
		return;
	}

	// bail out when density is negligible
	doublereal dRho = dGetAirDensity(GetXCurr());
	if (dRho <= std::numeric_limits<doublereal>::epsilon()) {
		return;
	}

	// Thrust in rotor reference frame
	doublereal dT = RRot3*Res.Force();

	/* Parametri di influsso (usano il valore di dUMean al passo precedente) */
	doublereal dVTip = 0.;
	dMu = 0.;
	dLambda = 0.;
	dVTip = dOmega*dRadius;

	if (dVTip > dVTipTreshold*dVTipRef) {
		dMu = (dVelocity*dCosAlphad)/dVTip;
	}

	/* NOTE: dUMeanRef starts at the value it had previously */
	if (std::abs(dUMeanRef) < std::numeric_limits<doublereal>::epsilon()
		&& std::abs(dT) > std::numeric_limits<doublereal>::epsilon())
	{
		// first guess: start with the value in hover
		dUMeanRef = copysign(std::sqrt(std::abs(dT)/(2*dRho*dArea)), dT);
	}

	/* Ground effect */
	doublereal dGE = 1.;
	if (pGround) {
		doublereal p = pGround->GetRCurr().GetVec(3)*(pRotor->GetXCurr() - pGround->GetXCurr());
		doublereal z = p/(dRadius/4.);

		if (z < .25) {
			if (z <= 0.) {
				silent_cerr("warning, illegal negative "
					"normalized altitude "
					"z=" << z << std::endl);
			}

			z = .25;
		}

		/*
		 * According to I.C. Cheeseman & N.E. Bennett,
		 * "The Effect of Ground on a Helicopter Rotor
		 * in Forward Flight", NASA TR-3021, 1955:
		 *
		 * U = Uref * ( 1 - ( R / ( 4 * z ) )^2 )
		 *
		 * We need to make R / ( 4 * z ) <= 1, so
		 * we must enforce z >= R / 4.
		 */
		dGE -= 1./(z*z);
	}

	bUMeanRefConverged = false;
	if (dVTip > dVTipTreshold*dVTipRef) {
		doublereal dCt = dT/(dRho*dArea*dVTip*dVTip);
		doublereal dLambda0 = dVelocity*dSinAlphad/dVTip;
		doublereal dLambdaInd = dUMeanRef/dVTip;

		bool bRetrying = false;
retry:;

		for (iCurrIter = 0; iCurrIter < iMaxIter; iCurrIter++) {
			dLambda = dLambda0 + dLambdaInd;

			doublereal dRef0 = dMu*dMu + dLambda*dLambda;
			doublereal dRef1 = 2.*sqrt(dRef0);
			doublereal dRef2 = dRef1*dRef0;
			doublereal dF = 0.;
			if (dRef1 > std::numeric_limits<doublereal>::epsilon()) {
				dF = dLambdaInd - dCt/dRef1;
				doublereal dFPrime = 1. + (dCt/dRef2)*dLambda;
				doublereal dDelta = dF/dFPrime;
				dLambdaInd -= dEta*dDelta;
			}

#if 0
			std::cerr
				<< " iter=" << iCurrIter
				<< " lambda=" << dLambda
				<< " dRef1=" << dRef1
				<< " dF=" << dF
				<< " Uref=" << dUMeanRef
				<< std::endl;
#endif

			if (std::abs(dF) <= dTolerance) {
				bUMeanRefConverged = true;
				break;
			}
		}

		// NOTE: the induced velocity must have the same sign
		// of the thrust coefficient
		if (dLambdaInd*dCt < 0.) {
			if (!bRetrying) {
				bRetrying = true;

				// first guess: revert the sign of induced velocity
				dLambdaInd = -dLambdaInd;
				goto retry;
			}

			silent_cerr("Rotor(" << GetLabel() << "): "
				"induced velocity and thrust signs differ "
				"(lambda_u=" << dLambdaInd << ", Ct=" << dCt << ")"
				<< std::endl);
		}

		dLambda = dLambda0 + dLambdaInd;

		// this is updated to serve as starting point for next iteration
		dUMeanRef = dLambdaInd*dVTip;

		/* if no convergence, simply accept the current value
		 * very forgiving choice, though */
#if 0
		if (iCurrIter == iMaxIter) {
			silent_cerr("unable to compute mean induced velocity "
				"for Rotor(" << GetLabel() << ")" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
#endif
	}

	//if (dMu == 0. && dLambda == 0.) {
	if (std::abs(dMu) < std::numeric_limits<doublereal>::epsilon() ) {
		dChi = 0.;

	} else {
		dChi = atan2(dMu, dLambda);
	}

	/*
	 * From Claudio Monteggia:
	 *
	 *                        Ct
	 * Um = -------------------------------------
	 *       sqrt( lambda^2 / KH^4 + mu^2 / KF^2)
	 */
	if (dVTip > dVTipTreshold*dVTipRef) {
		doublereal dMuTmp = dMu/dForwardFlightCorrection;
		doublereal dLambdaTmp = dLambda/(dHoverCorrection*dHoverCorrection);
		doublereal dRef = 2*sqrt(dMuTmp*dMuTmp + dLambdaTmp*dLambdaTmp);
		doublereal dCt = dT/(dRho*dArea*dVTip*dVTip);		
		if (dRef > std::abs(dCt)) {
			// NOTE: an inflow velocity larger than VTip
			// makes little sense
			dUMean = (1. - dWeight)*dGE*dVTip*dCt/dRef + dWeight*dUMeanPrev;

		} else {
			dUMean = dUMeanPrev;
		}

	} else {
		dUMean = dUMeanPrev;
	}
}

std::ostream&
Rotor::Restart(std::ostream& out) const
{
	return out << "  rotor: " << GetLabel() << ", "
		<< pCraft->GetLabel() << ", " << pRotor->GetLabel()
		<< ", induced velocity: ";
}

/* Rotor - end */


/* NoRotor - begin */

NoRotor::NoRotor(unsigned int uLabel,
	const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	NO_OP;
}

NoRotor::NoRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	ResForceSet **ppres,
	const doublereal& dR,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	Init(pCraft, rrot, pRotor, ppres, dR, fOut);
}

void
NoRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	ResForceSet **ppres,
	const doublereal& dR,
	flag fOut)
{
	Rotor::Init(pCraft, rrot, pRotor, 0, ppres, dR, 0, 0., 0., fOut);

#ifdef USE_MPI
	if (is_parallel && bToBeOutput()) {
		SAFENEWARR(pBlockLenght, int, 3);
		SAFENEWARR(pDispl, MPI::Aint, 3);
		for (int i = 0; i < 3; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
                        pDispl[i] = MPI::Get_address(const_cast<doublereal*>(&(Res.Pole().pGetVec()[i])));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
			MPI::Datatype(MPI::DOUBLE.Create_hindexed(3, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}

NoRotor::~NoRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}

/* assemblaggio residuo */
SubVectorHandler&
NoRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering NoRotor::AssRes()" << std::endl);

	flag out = fToBeOutput();

#ifdef USE_MPI
	if (out) {
		ExchangeLoads(out);
	}

	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		if (out) {
			/* Calcola parametri vari */
			Rotor::InitParam(false);

#ifdef USE_MPI
			ExchangeVelocity();
#endif /* USE_MPI */
		}
	}

	ResetForce();
	WorkVec.Resize(0);

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	return WorkVec;
}

std::ostream&
NoRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "no;" << std::endl;
}

// Somma alla trazione il contributo di forza di un elemento generico
void
NoRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
	// Non gli serve in quanto non calcola velocita' indotta.
	// Solo se deve fare l'output lo calcola
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif // USE_MPI

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
}

/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3
NoRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return Zero3;
}

/* NoRotor - end */


/* UniformRotor - begin */

UniformRotor::UniformRotor(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	NO_OP;
}

UniformRotor::UniformRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	Init(pCraft, rrot, pRotor, pGround, ppres, dOR, dR,
		pdW, iMaxIt, dTol, dE, dCH, dCFF, fOut);
}

void
UniformRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
{
	ASSERT(dOR > 0.);
	ASSERT(dR > 0.);
	ASSERT(pdW != 0);

	Rotor::Init(pCraft, rrot, pRotor, pGround, ppres, dR, iMaxIt, dTol, dE, fOut);

	dOmegaRef = dOR;
	dVTipRef = dOmegaRef*dRadius;
	dArea = M_PI*dRadius*dRadius;
	Weight.Set(pdW);
	dHoverCorrection = dCH;
	dForwardFlightCorrection = dCFF;

#ifdef USE_MPI
	if (is_parallel) {
		SAFENEWARR(pBlockLenght, int, 7);
		SAFENEWARR(pDispl, MPI::Aint, 7);
		for (int i = 0; i < 7; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
			pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));
		}
		pDispl[3] = MPI::Get_address(&dUMeanPrev);
		for (int i = 4; i <= 6; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(&(Res.Pole().pGetVec()[i-4])));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
				MPI::Datatype(MPI::DOUBLE.Create_hindexed(7, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}

UniformRotor::~UniformRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}

/* assemblaggio residuo */
SubVectorHandler&
UniformRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering UniformRotor::AssRes()" << std::endl);

#ifdef USE_MPI
	ExchangeLoads(fToBeOutput());
	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		/* Calcola parametri vari */
		Rotor::InitParam();

#ifdef DEBUG
		// Prova:
#if 0
		Vec3 XTmp(2.,2.,0.);
		doublereal dPsiTmp = dGetPsi(XTmp);
		doublereal dXTmp = dGetPos(XTmp);
		std::cout
			<< "X rotore:  " << pRotor->GetXCurr() << std::endl
			<< "V rotore:  " << VCraft << std::endl
			<< "X punto:   " << XTmp << std::endl
			<< "Omega:     " << dOmega << std::endl
			<< "Velocita': " << dVelocity << std::endl
			<< "Psi0:      " << dPsi0 << std::endl
			<< "Psi punto: " << dPsiTmp << std::endl
			<< "Raggio:    " << dRadius << std::endl
			<< "r punto:   " << dXTmp << std::endl
			<< "mu:        " << dMu << std::endl
			<< "lambda:    " << dLambda << std::endl
			<< "cos(ad):   " << dCosAlphad << std::endl
			<< "sin(ad):   " << dSinAlphad << std::endl
			<< "UMean:     " << dUMean << std::endl;
#endif
#endif /* DEBUG */
	}

#ifdef USE_MPI
	ExchangeVelocity();
#endif /* USE_MPI */

	ResetForce();

	/* Non tocca il residuo */
	WorkVec.Resize(0);

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	return WorkVec;
}

std::ostream&
UniformRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "uniform, " << dRadius << ", ",
		Weight.pGetDriveCaller()->Restart(out)
		<< ", correction, " << dHoverCorrection
		<< ", " << dForwardFlightCorrection << ';' << std::endl;
}

/* Somma alla trazione il contributo di forza di un elemento generico */
void
UniformRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	/* Solo se deve fare l'output calcola anche il momento */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else {
		Res.AddForce(F);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}

/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3
UniformRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	return RRot3*dUMeanPrev;
};

UniformRotor2::UniformRotor2(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
UniformRotor(uLabel, pDO)
{
	NO_OP;
}

UniformRotor2::UniformRotor2(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
: Elem(uLabel, fOut),
UniformRotor(uLabel, pDO, pCraft, rrot, pRotor, pGround, ppres, dOR, dR, pdW, iMaxIt, dTol, dE, dCH, dCFF, fOut)
{
	NO_OP;
}

UniformRotor2::~UniformRotor2(void)
{
	NO_OP;
}

bool
UniformRotor2::bSectionalForces(void) const
{
	return true;
}

/* Somma alla trazione il contributo di forza di un elemento generico */
void
UniformRotor2::AddSectionalForce(Elem::Type type,
		const Elem *pEl, unsigned uPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W)
{
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	/* Solo se deve fare l'output calcola anche il momento */
	if (bToBeOutput()) {
		Vec3 FTmp(F*dW);
		Vec3 MTmp(M*dW);
		Res.AddForces(FTmp, MTmp, X);
		InducedVelocity::AddForce(pEl, 0, FTmp, MTmp, X);

	} else {
		Res.AddForce(F*dW);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}

/* UniformRotor - end */


/* GlauertRotor - begin */

GlauertRotor::GlauertRotor(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
gtype(GlauertRotor::UNKNOWN)
{
	NO_OP;
}

GlauertRotor::GlauertRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	GlauertRotor::Type gtype,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
gtype(gtype)
{
	Init(pCraft, rrot, pRotor, pGround, ppres, dOR, dR,
		pdW, iMaxIt, dTol, dE, dCH, dCFF, fOut);
}

void
GlauertRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
{
	ASSERT(dOR > 0.);
	ASSERT(dR > 0.);
	ASSERT(pdW != 0);

	Rotor::Init(pCraft, rrot, pRotor, pGround, ppres, dR, iMaxIt, dTol, dE, fOut);

	dOmegaRef = dOR;
	dVTipRef = dOmegaRef*dRadius;
	dArea = M_PI*dRadius*dRadius;
	Weight.Set(pdW);
	dHoverCorrection = dCH;
	dForwardFlightCorrection = dCFF;

#ifdef USE_MPI
	if (is_parallel) {
		SAFENEWARR(pBlockLenght, int, 20);
		SAFENEWARR(pDispl, MPI::Aint, 20);
		for (int i = 0; i < 20; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
			pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));
		}
		pDispl[3] = MPI::Get_address(&dUMeanPrev);
		pDispl[4] = MPI::Get_address(&dLambda);
		pDispl[5] = MPI::Get_address(&dMu);
		pDispl[6] = MPI::Get_address(&dChi);
		pDispl[7] = MPI::Get_address(&dPsi0);
		for (int i = 8; i <= 10; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(&(Res.Pole().pGetVec()[i-8])));
		}
		for (int i = 11; i < 20; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(&(RRotTranspose.pGetMat()[i-11])));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
				MPI::Datatype(MPI::DOUBLE.Create_hindexed(20, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}


GlauertRotor::~GlauertRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}


/* assemblaggio residuo */
SubVectorHandler&
GlauertRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering GlauertRotor::AssRes()" << std::endl);

#ifdef USE_MPI
	ExchangeLoads(fToBeOutput());
	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		/* Calcola parametri vari */
		Rotor::InitParam();
	}

#ifdef USE_MPI
	ExchangeVelocity();
#endif /* USE_MPI */

	ResetForce();

	/* Non tocca il residuo */
	WorkVec.Resize(0);

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	return WorkVec;
}

std::ostream&
GlauertRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "Glauert, " << dRadius << ", ",
		Weight.pGetDriveCaller()->Restart(out)
		<< ", correction, " << dHoverCorrection
		<< ", " << dForwardFlightCorrection << ';' << std::endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void
GlauertRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	/* Solo se deve fare l'output calcola anche il momento */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else {
		Res.AddForce(F);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}


/*
 * Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale
 */
Vec3
GlauertRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	if (dUMeanPrev == 0.) {
		return Zero3;
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	if (std::abs(dLambda) < 1.e-9) {
		return RRot3*dUMeanPrev;
	}

	doublereal dr, dp;
	GetPos(X, dr, dp);

	// recall that tan(dChi) = dMu/dLambda

	doublereal k1, k2 = 0.;
	switch (gtype) {
	case GLAUERT:
		// NASA-TM-102219: attributed to Drees et al.
		k1 = 4./3.*(1. - 1.8*dMu*dMu)*tan(dChi/2.);
		break;

#if 0
	case WHEATLEY:
		// NASA-TM-102219
		k1 = 0.5;
		break;
#endif

	case COLEMAN_ET_AL:
		// NASA-TM-102219
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006
		// Wayne Johnson, "Rotorcraft Aeromechanics", 2013
		k1 = tan(dChi/2.);
		break;

	case DREES_1:
	case DREES_2:
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006
		// Wayne Johnson, "Rotorcraft Aeromechanics", 2013
		// k1 = 4./3.*(1 - cos(dChi) - 1.8*dMu*dMu)/sin(dChi); // risk of division by zero...
		k1 = 4./3.*(tan(dChi/2.) - 1.8*dLambda*dMu/cos(dChi));
		k2 = -2.*dMu;
		break;

	case PAYNE:
		// NASA-TM-102219
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006
		// k1 = 4./3.*(dMu/dLambda/(1.2 + dMu/dLambda));
		// reduce risk of division by zero
		k1 = 4./3.*dMu/(1.2*dLambda + dMu);
		break;

	case WHITE_AND_BLAKE:
		// NASA-TM-102219
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006
		// Wayne Johnson, "Rotorcraft Aeromechanics", 2013
		k1 = sqrt(2.)*sin(dChi);
		break;

	case PITT_AND_PETERS:
		// NASA-TM-102219
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006: "23" instead of "32"...
		k1 = 15.*M_PI/32.*tan(dChi/2.);
		break;

	case HOWLETT:
		// NASA-TM-102219
		// Gordon J. Leishman, "Principles of Helicopter Aerodynamics", 2nd Ed., 2006
		k1 = pow(sin(dChi), 2);
		break;

#if 0
	case DREES_2: {
		// Jianhua Zhang, "Active-Passive Hybrid Optimization of Rotor Blades With Trailing Edge Flaps", PhD Thesis, PSU, 2001
		// in the end, it is identical to Drees' as of Wayne Johnson and Gordon J. Leishman
		// FIXME: what if dMu ~ 0?
		doublereal dLdM = dLambda/dMu;
		// k1 = 4./3.*((1. - 1.8*dMu*dMu)*sqrt(1. + dLdM*dLdM - dLdM));
		k1 = 4./3.*((1. - 1.8*dMu*dMu)*sqrt(1. + dLdM*dLdM) - dLdM);
		k2 = -2.*dMu;
		} break;
#endif

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	doublereal dd = 1. + dr*(k1*cos(dp) + k2*sin(dp));

	return RRot3*(dd*dUMeanPrev);
};

/* GlauertRotor - end */


/* ManglerRotor - begin */

ManglerRotor::ManglerRotor(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	NO_OP;
}

ManglerRotor::ManglerRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO)
{
	Init(pCraft, rrot, pRotor, pGround, ppres, dOR, dR, pdW, iMaxIt, dTol, dE, dCH, dCFF, fOut);
}

void
ManglerRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	DriveCaller *pdW,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	flag fOut)
{
	ASSERT(dOR > 0.);
	ASSERT(dR > 0.);
	ASSERT(pdW != 0);

	Rotor::Init(pCraft, rrot, pRotor, pGround, ppres, dR, iMaxIt, dTol, dE, fOut);

	dOmegaRef = dOR;
	dVTipRef = dOmegaRef*dRadius;
	dArea = M_PI*dRadius*dRadius;
	Weight.Set(pdW);
	dHoverCorrection = dCH;
	dForwardFlightCorrection = dCFF;

#ifdef USE_MPI
	if (is_parallel) {
		SAFENEWARR(pBlockLenght, int, 18);
		SAFENEWARR(pDispl, MPI::Aint, 18);
		for (int i = 0; i < 18; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
			pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));
		}
		pDispl[3] = MPI::Get_address(&dUMeanPrev);
		pDispl[4] = MPI::Get_address(&dSinAlphad);
		pDispl[5] = MPI::Get_address(&dPsi0);
		for (int i = 6; i <= 8; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(&(Res.Pole().pGetVec()[i-6])));
		}
		for (int i = 9; i < 18; i++) {
			pDispl[i] = MPI::Get_address(&(RRotTranspose.pGetMat()[i-9]));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
				MPI::Datatype(MPI::DOUBLE.Create_hindexed(18, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}


ManglerRotor::~ManglerRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}


/* assemblaggio residuo */
SubVectorHandler&
ManglerRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering ManglerRotor::AssRes()" << std::endl);

#ifdef USE_MPI
	ExchangeLoads(fToBeOutput());
	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		/* Calcola parametri vari */
		Rotor::InitParam();

#ifdef DEBUG
		// Prova:
#if 0
		Vec3 XTmp(2.,2.,0.);
		doublereal dPsiTmp = dGetPsi(XTmp);
		doublereal dXTmp = dGetPos(XTmp);
		Vec3 IndV = GetInducedVelocity(XTmp);
		std::cout
			<< "X rotore:  " << pRotor->GetXCurr() << std::endl
			<< "V rotore:  " << VCraft << std::endl
			<< "X punto:   " << XTmp << std::endl
			<< "Omega:     " << dOmega << std::endl
			<< "Velocita': " << dVelocity << std::endl
			<< "Psi0:      " << dPsi0 << std::endl
			<< "Psi punto: " << dPsiTmp << std::endl
			<< "Raggio:    " << dRadius << std::endl
			<< "r punto:   " << dXTmp << std::endl
			<< "mu:        " << dMu << std::endl
			<< "lambda:    " << dLambda << std::endl
			<< "cos(ad):   " << dCosAlphad << std::endl
			<< "sin(ad):   " << dSinAlphad << std::endl
			<< "UMean:     " << dUMean << std::endl
			<< "iv punto:  " << IndV << std::endl;
#endif
#endif /* DEBUG */
	}

#ifdef USE_MPI
	ExchangeVelocity();
#endif /* USE_MPI */

	ResetForce();

	/* Non tocca il residuo */
	WorkVec.Resize(0);

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	return WorkVec;
}

std::ostream&
ManglerRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "Mangler, " << dRadius << ", ",
		Weight.pGetDriveCaller()->Restart(out)
		<< ", correction, " << dHoverCorrection
		<< ", " << dForwardFlightCorrection << ';' << std::endl;
}

/* Somma alla trazione il contributo di forza di un elemento generico */
void
ManglerRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	/* Solo se deve fare l'output calcola anche il momento */
	if (bToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else {
		Res.AddForce(F);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}


/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3
ManglerRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	if (dUMeanPrev == 0.) {
		return ::Zero3;
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	doublereal dr, dp;
	GetPos(X, dr, dp);

	doublereal dr2 = dr*dr;
	doublereal dm2 = 1. - dr2;
	doublereal dm = 0.;
	if (dm2 > 0.) {
		dm = sqrt(dm2);
	}
	doublereal da = 1. + std::abs(dSinAlphad);
	if (da > 1.e-9) {
		da = (1. - std::abs(dSinAlphad))/da;
	}
	if (da > 0.) {
		da = sqrt(da);
	}

	// c_0
	doublereal dd = 15./4.*dm*dr2;

	// c_1
	doublereal dc = -15./256.*M_PI*(9.*dr2 - 4.)*dr*da;
	dd -= 4.*dc*cos(dp);

	// c_3
	dc = 45./256.*M_PI*pow(da*dr, 3);
	dd -= 4.*dc*cos(3.*dp);

	// c_[5:2:infty] = 0

	// c_[2:2:infty] (truncated at 10)
	for (int i = 2; i <= 10; i += 2) {
		dc = pow(-1., i/2 - 1)*15./8.
			*((dm + i)/(i*i - 1.)*(3. - 9.*dr2 + i*i) + 3.*dm)/(i*i - 9.)
			*pow(((1. - dm)/(1. + dm))*da, i/2.);
		dd -= 4.*dc*cos(i*dp);
	}

	return RRot3*(dd*dUMeanPrev);
};

/* ManglerRotor - end */



/*
 * According to
 * Dale M. Pitt, David A. Peters, 
 * Theoretical Prediction of Dynamic-Inflow Derivatives,
 * Vertica, 5, 1981, pp. 21-34
 *
 * 8/(3*pi): uncorrected value
 * 128/(75*pi): corrected value
 */

#if 0
static const doublereal dM11 = 8./(3.*M_PI);
#else
static const doublereal dM11 = 128./(75.*M_PI); 
#endif

static const doublereal dM22 = -16./(45.*M_PI);
static const doublereal dM33 = -16./(45.*M_PI);

/* DynamicInflowRotor - begin */

DynamicInflowRotor::DynamicInflowRotor(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
dVConst(0), dVSine(0), dVCosine(0),
dL11(0.), dL13(0.), dL22(0.), dL31(0.), dL33(0.)
{
	NO_OP;
}

DynamicInflowRotor::DynamicInflowRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	const doublereal& dVConstTmp,
	const doublereal& dVSineTmp,
	const doublereal& dVCosineTmp,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
dVConst(0), dVSine(0), dVCosine(0),
dL11(0.), dL13(0.), dL22(0.), dL31(0.), dL33(0.)
{
	Init(pCraft, rrot, pRotor, pGround, ppres, dOR, dR,
		iMaxIt, dTol, dE, dCH, dCFF,
		dVConstTmp, dVSineTmp, dVCosineTmp, fOut);
}

void
DynamicInflowRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	const doublereal& dVConstTmp,
	const doublereal& dVSineTmp,
	const doublereal& dVCosineTmp,
	flag fOut)
{
	ASSERT(dOR > 0.);
	ASSERT(dR > 0.);

	Rotor::Init(pCraft, rrot, pRotor, pGround, ppres, dR, iMaxIt, dTol, dE, fOut);

	dVConst = dVConstTmp;
	dVSine = dVSineTmp;
	dVCosine = dVCosineTmp;

	dOmegaRef = dOR;
	dVTipRef = dOmegaRef*dRadius;
	dArea = M_PI*dRadius*dRadius;

	dHoverCorrection = dCH;
	dForwardFlightCorrection = dCFF;

	/* Significa che valuta la velocita' indotta media al passo corrente */
	dWeight = 0.;

#ifdef USE_MPI
	if (is_parallel) {
		SAFENEWARR(pBlockLenght, int, 20);
		SAFENEWARR(pDispl, MPI::Aint, 20);
		for (int i = 0; i < 20; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
			pDispl[i] = MPI::Get_address(RRot3.pGetVec()+i);
		}
		pDispl[3] = MPI::Get_address(&dVConst);
		pDispl[4] = MPI::Get_address(&dVSine);
		pDispl[5] = MPI::Get_address(&dVCosine);
		pDispl[6] = MPI::Get_address(&dOmega);
		pDispl[7] = MPI::Get_address(&dPsi0);
		for (int i = 8; i <= 10; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(Res.Pole().pGetVec()+i-8));
		}
		for (int i = 11; i < 20; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(RRotTranspose.pGetMat()+i-11));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
				MPI::Datatype(MPI::DOUBLE.Create_hindexed(20, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}


DynamicInflowRotor::~DynamicInflowRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}

void
DynamicInflowRotor::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
	     if (OH.UseNetCDF(OutputHandler::ROTORS)) {
		ASSERT(OH.IsOpen(OutputHandler::NETCDF));
		/* The first part of the output is the same for Rotor and
		 * DynamicInflowRotor
		 */
		Rotor::OutputPrepare(OH);

		std::ostringstream os;
		os << "elem.inducedvelocity." << GetLabel() << ".";
		std::string name = os.str();
		Var_dVConst = OH.CreateVar<doublereal>(name + "VConst",
				OutputHandler::Dimensions::Velocity,
				"constant inflow state");
		Var_dVSine = OH.CreateVar<doublereal>(name + "VSine",
				OutputHandler::Dimensions::Velocity,
				"sine inflow state (lateral)");
		Var_dVCosine= OH.CreateVar<doublereal>(name + "VCosine",
				OutputHandler::Dimensions::Velocity,
				"cosine inflow state (longitudinal)");
	     }
#endif // USE_NETCDF
	}
}

void
DynamicInflowRotor::Output(OutputHandler& OH) const
{
     	/*
	 * FIXME: posso usare dei temporanei per il calcolo della trazione
	 * totale per l'output, cosi' evito il giro dei cast
	 */
	if (bToBeOutput()) {
#ifdef USE_MPI
		if (is_parallel && IndVelComm.Get_size() > 1) {
			if (IndVelComm.Get_rank() == 0) {
				Vec3 TmpF(pTmpVecR), TmpM(pTmpVecR+3);
#ifdef USE_NETCDF
				if (OH.UseNetCDF(OutputHandler::ROTORS)) {
					OH.WriteNcVar(Var_f, RRotTranspose*TmpF);
					OH.WriteNcVar(Var_m, RRotTranspose*TmpM);
					OH.WriteNcVar(Var_dUMean, dUMean);
					OH.WriteNcVar(Var_dVelocity, dVelocity);
					OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
					OH.WriteNcVar(Var_dMu, dMu);
					OH.WriteNcVar(Var_dLambda, dLambda);
					OH.WriteNcVar(Var_dChi, dChi);
					OH.WriteNcVar(Var_dPsi0, dPsi0);
					OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
					OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
					OH.WriteNcVar(Var_dVConst, dVConst);
					OH.WriteNcVar(Var_dVSine, dVSine);
					OH.WriteNcVar(Var_dVCosine, dVCosine);
				}
#endif // USE_NETCDF
				if (OH.UseText(OutputHandler::ROTORS)) {
					OH.Rotors()
						<< std::setw(8) << GetLabel()	/* 1 */
						<< " " << RRotTranspose*TmpF /* 2-4 */
						<< " " << RRotTranspose*TmpM /* 5-7 */
						<< " " << dUMean	/* 8 */
						<< " " << dVelocity	/* 9 */
						<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
						<< " " << dMu		/* 11 */
						<< " " << dLambda	/* 12 */
						<< " " << dChi		/* 13 */
						<< " " << dPsi0		/* 14 */
						<< " " << bUMeanRefConverged /* 15 */
						<< " " << iCurrIter	/* 16 */
						<< " " << dVConst	/* 17 */
						<< " " << dVSine	/* 18 */
						<< " " << dVCosine	/* 19 */
						<< std::endl;

					for (int i = 0; ppRes && ppRes[i]; i++) {
						Vec3 TmpF(pTmpVecR+6+6*i);
						Vec3 TmpM(pTmpVecR+9+6*i);

						OH.Rotors()
							<< std::setw(8) << GetLabel()
							<< ":" << ppRes[i]->GetLabel()
							<< " " << TmpF
							<< " " << TmpM
							<< std::endl;
					}
				}
			}
		} else {
#ifdef USE_NETCDF
			if (OH.UseNetCDF(OutputHandler::ROTORS)) {
				OH.WriteNcVar(Var_f, RRotTranspose*Res.Force());
				OH.WriteNcVar(Var_m, RRotTranspose*Res.Moment());
				OH.WriteNcVar(Var_dUMean, dUMean);
				OH.WriteNcVar(Var_dVelocity, dVelocity);
				OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
				OH.WriteNcVar(Var_dMu, dMu);
				OH.WriteNcVar(Var_dLambda, dLambda);
				OH.WriteNcVar(Var_dChi, dChi);
				OH.WriteNcVar(Var_dPsi0, dPsi0);
				OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
				OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
				OH.WriteNcVar(Var_dVConst, dVConst);
				OH.WriteNcVar(Var_dVSine, dVSine);
				OH.WriteNcVar(Var_dVCosine, dVCosine);
			}
#endif // USE_NETCDF
			if (OH.UseText(OutputHandler::ROTORS)) {
				OH.Rotors()
					<< std::setw(8) << GetLabel()	/* 1 */
					<< " " << RRotTranspose*Res.Force()  /* 2-4 */
					<< " " << RRotTranspose*Res.Moment() /* 5-7 */
					<< " " << dUMean	/* 8 */
					<< " " << dVelocity	/* 9 */
					<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
					<< " " << dMu		/* 11 */
					<< " " << dLambda	/* 12 */
					<< " " << dChi		/* 13 */
					<< " " << dPsi0		/* 14 */
					<< " " << bUMeanRefConverged /* 15 */
					<< " " << iCurrIter	/* 16 */
					<< " " << dVConst	/* 17 */
					<< " " << dVSine	/* 18 */
					<< " " << dVCosine	/* 19 */
					<< std::endl;

				for (int i = 0; ppRes && ppRes[i]; i++) {
					OH.Rotors()
						<< std::setw(8) << GetLabel()
						<< ":" << ppRes[i]->GetLabel()
						<< " " << ppRes[i]->pRes->Force()
						<< " " << ppRes[i]->pRes->Moment()
						<< std::endl;
				}
			}
		}
#else /* !USE_MPI */
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::ROTORS)) {
			OH.WriteNcVar(Var_f, RRotTranspose*Res.Force());
			OH.WriteNcVar(Var_m, RRotTranspose*Res.Moment());
			OH.WriteNcVar(Var_dUMean, dUMean);
			OH.WriteNcVar(Var_dVelocity, dVelocity);
			OH.WriteNcVar(Var_dAlpha, atan2(dSinAlphad, dCosAlphad));
			OH.WriteNcVar(Var_dMu, dMu);
			OH.WriteNcVar(Var_dLambda, dLambda);
			OH.WriteNcVar(Var_dChi, dChi);
			OH.WriteNcVar(Var_dPsi0, dPsi0);
			OH.WriteNcVar(Var_bUMeanRefConverged, (int)bUMeanRefConverged);
			OH.WriteNcVar(Var_iCurrIter, (int)iCurrIter);
			OH.WriteNcVar(Var_dVConst, dVConst);
			OH.WriteNcVar(Var_dVSine, dVSine);
			OH.WriteNcVar(Var_dVCosine, dVCosine);
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::ROTORS)) {
			OH.Rotors()
				<< std::setw(8) << GetLabel()	/* 1 */
				<< " " << RRotTranspose*Res.Force()	/* 2-4 */
				<< " " << RRotTranspose*Res.Moment()	/* 5-7 */
				<< " " << dUMean	/* 8 */
				<< " " << dVelocity	/* 9 */
				<< " " << atan2(dSinAlphad,dCosAlphad)	/* 10 */
				<< " " << dMu		/* 11 */
				<< " " << dLambda	/* 12 */
				<< " " << dChi		/* 13 */
				<< " " << dPsi0		/* 14 */
				<< " " << bUMeanRefConverged /* 15 */
				<< " " << iCurrIter	/* 16 */
				<< " " << dVConst	/* 17 */
				<< " " << dVSine	/* 18 */
				<< " " << dVCosine	/* 19 */
				<< std::endl;

			for (int i = 0; ppRes && ppRes[i]; i++) {
				OH.Rotors()
					<< std::setw(8) << GetLabel()
					<< ":" << ppRes[i]->GetLabel()
					<< " " << ppRes[i]->pRes->Force()
					<< " " << ppRes[i]->pRes->Moment()
					<< std::endl;
			}
		}
#endif /* !USE_MPI */
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
DynamicInflowRotor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering DynamicInflowRotor::AssJac()" << std::endl);

	WorkMat.SetNullMatrix();

#ifdef USE_MPI
	if (is_parallel && IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		SparseSubMatrixHandler& WM = WorkMat.SetSparse();
		integer iFirstIndex = iGetFirstIndex();

		WM.ResizeReset(5, 0);

		WM.PutItem(1, iFirstIndex + 1, iFirstIndex + 1, dM11 + dCoef*dL11);
		WM.PutItem(2, iFirstIndex + 3, iFirstIndex + 1, dCoef*dL31);
		WM.PutItem(3, iFirstIndex + 2, iFirstIndex + 2, dM22 + dCoef*dL22);
		WM.PutItem(4, iFirstIndex + 1, iFirstIndex + 3, dCoef*dL13);
		WM.PutItem(5, iFirstIndex + 3, iFirstIndex + 3, dM33 + dCoef*dL33);
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
DynamicInflowRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
     	DEBUGCOUT("Entering DynamicInflowRotor::AssRes()" << std::endl);

#ifdef USE_MPI
     	ExchangeLoads(flag(1));

     	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
	   	/* Calcola parametri vari */
	   	Rotor::InitParam();

       	   	WorkVec.Resize(3);

	   	integer iFirstIndex = iGetFirstIndex();

	   	WorkVec.PutRowIndex(1, iFirstIndex + 1);
	   	WorkVec.PutRowIndex(2, iFirstIndex + 2);
	   	WorkVec.PutRowIndex(3, iFirstIndex + 3);

	   	dVConst = XCurr(iFirstIndex + 1);
	   	dVSine = XCurr(iFirstIndex + 2);
	   	dVCosine = XCurr(iFirstIndex + 3);

	   	doublereal dVConstPrime = XPrimeCurr(iFirstIndex + 1);
	   	doublereal dVSinePrime = XPrimeCurr(iFirstIndex + 2);
	   	doublereal dVCosinePrime = XPrimeCurr(iFirstIndex + 3);

		doublereal dCT = 0.;
		doublereal dCl = 0.;
		doublereal dCm = 0.;

	   	/*
		 * Attenzione: moltiplico tutte le equazioni per dOmega
		 * (ovvero, i coefficienti CT, CL e CM sono divisi
		 * per dOmega anziche' dOmega^2)
		 */

		dL11 = 0.;
		dL13 = 0.;
		dL22 = 0.;
		dL31 = 0.;
		dL33 = 0.;

	   	doublereal dDim = dGetAirDensity(GetXCurr())*dArea*dOmega*(dRadius*dRadius);
	   	if (dDim > std::numeric_limits<doublereal>::epsilon()) {
			/*
			 * From Claudio Monteggia:
			 *
			 *                        Ct
			 * Um = -------------------------------------
			 *       sqrt( lambda^2 / KH^4 + mu^2 / KF^2)
			 */
		 	doublereal dLambdaTmp
				= dLambda/(dHoverCorrection*dHoverCorrection);
		 	doublereal dMuTmp = dMu/dForwardFlightCorrection;

		 	doublereal dVT
				= sqrt(dLambdaTmp*dLambdaTmp + dMuTmp*dMuTmp);
		 	doublereal dVm = 0.;
		 	if (dVT > dVTipTreshold*dVTipRef) {
		       		dVm = (dMuTmp*dMuTmp + dLambdaTmp*(dLambdaTmp + dVConst))/dVT;
		 	}

			/*
			 * dUMean is just for output;
			 */
		   	dUMean = dVConst*dOmega*dRadius;

		   	/* Trazione nel sistema rotore */
		   	doublereal dT = RRot3*Res.Force();

		   	/* Momento nel sistema rotore-vento */
		   	doublereal dCosP = cos(dPsi0);
		   	doublereal dSinP = sin(dPsi0);
		   	Mat3x3 RTmp( dCosP, -dSinP, 0.,
		      			dSinP, dCosP, 0.,
		      			0., 0., 1.);

		   	Vec3 M(RTmp*(RRotTranspose*Res.Moment()));

		 	/* Thrust, roll and pitch coefficients */
		 	dCT = dT/dDim;
		 	dDim *= dRadius;
			// Note: Neda Taymourtash noted a possible sign error in the Cm; however, a sign error in the Cl is more likely
		 	dCl = - M(1)/dDim;
		 	// dCm = - M(2)/dDim;
		 	dCm = M(2)/dDim;

			if (dVT > std::numeric_limits<doublereal>::epsilon()
				&& dVm > std::numeric_limits<doublereal>::epsilon())
			{

				/* Matrix coefficients */
				/* FIXME: divide by 0? */
			 	doublereal dl11 = .5/dVT;
				/* FIXME: divide by 0? */
			 	doublereal d = 15./64.*M_PI*tan(dChi/2.);
				/* FIXME: divide by 0? */
			 	doublereal dl13 = d/dVm;
				/* FIXME: divide by 0? */
			 	doublereal dl31 = d/dVT;

			 	doublereal dCosChi2 = cos(dChi/2.);
			 	d = 2.*dCosChi2*dCosChi2;
				/* FIXME: divide by 0? */
			 	doublereal dl22 = -4./(d*dVm);
				/* FIXME: divide by 0? */
				doublereal dl33 = -4.*(d - 1)/(d*dVm);
	
				d = dl11*dl33 - dl31*dl13;
				/* FIXME: divide by 0? */
				dL11 = dOmega*dl33/d;
				dL31 = -dOmega*dl31/d;
				dL13 = -dOmega*dl13/d;
				dL33 = dOmega*dl11/d;
				dL22 = dOmega/dl22;
		   	}
		}
	
#ifdef DEBUG
	   	/* Prova: */
		static int i = -1;
		int iv[] = { 0, 1, 0, -1, 0 };
		if (++i == 4) {
			i = 0;
		}
	   	Vec3 XTmp(pRotor->GetXCurr()+pCraft->GetRCurr()*Vec3(dRadius*iv[i],dRadius*iv[i+1],0.));
	   	doublereal dPsiTmp, dXTmp;
	   	GetPos(XTmp, dXTmp, dPsiTmp);
	   	Vec3 IndV = GetInducedVelocity(GetElemType(), GetLabel(), 0, XTmp);
	   	std::cout
		 	<< "X rotore:  " << pRotor->GetXCurr() << std::endl
		 	<< "V rotore:  " << VCraft << std::endl
		 	<< "X punto:   " << XTmp << std::endl
		 	<< "Omega:     " << dOmega << std::endl
		 	<< "Velocita': " << dVelocity << std::endl
		 	<< "Psi0:      " << dPsi0 << " ("
				<< dPsi0*180./M_PI << " deg)" << std::endl
		 	<< "Psi punto: " << dPsiTmp << " ("
				<< dPsiTmp*180./M_PI << " deg)" << std::endl
		 	<< "Raggio:    " << dRadius << std::endl
		 	<< "r punto:   " << dXTmp << std::endl
		 	<< "mu:        " << dMu << std::endl
		 	<< "lambda:    " << dLambda << std::endl
		 	<< "cos(ad):   " << dCosAlphad << std::endl
		 	<< "sin(ad):   " << dSinAlphad << std::endl
		 	<< "UMean:     " << dUMean << std::endl
		 	<< "iv punto:  " << IndV << std::endl;
#endif /* DEBUG */

		WorkVec.PutCoef(1, dCT - dM11*dVConstPrime
				- dL11*dVConst - dL13*dVCosine);
	 	WorkVec.PutCoef(2, dCl - dM22*dVSinePrime
				- dL22*dVSine);
	 	WorkVec.PutCoef(3, dCm - dM33*dVCosinePrime
				- dL31*dVConst - dL33*dVCosine);

#ifdef USE_MPI
     	} else {
	   	WorkVec.Resize(0);
#endif /* USE_MPI */
     	}

#ifdef USE_MPI
	ExchangeVelocity();
#endif /* USE_MPI */

	/* Ora la trazione non serve piu' */
	ResetForce();

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

     	return WorkVec;
}

/* Relativo ai ...WithDofs */
void
DynamicInflowRotor::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}


/* Relativo ai ...WithDofs */
void
DynamicInflowRotor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		XP.PutCoef(iFirstIndex + iCnt, 0.);
	}

	X.PutCoef(iFirstIndex + 1, dVConst);
	X.PutCoef(iFirstIndex + 2, dVSine);
	X.PutCoef(iFirstIndex + 3, dVCosine);
}


/* Restart */
std::ostream&
DynamicInflowRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "dynamic inflow, " << dRadius
		<< ", correction, " << dHoverCorrection
		<< ", " << dForwardFlightCorrection << ';' << std::endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void
DynamicInflowRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
	/*
	 * Gli serve la trazione ed il momento rispetto al rotore,
	 * che si calcola da se'
	 */
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	Res.AddForces(F, M, X);
	if (bToBeOutput()) {
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}


/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3
DynamicInflowRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	doublereal dr, dpsi;
	GetPos(X, dr, dpsi);

	return RRot3*((dVConst + dr*(dVCosine*cos(dpsi) + dVSine*sin(dpsi)))*dRadius*dOmega);
};

/* DynamicInflowRotor - end */


/* PetersHeRotor - begin */

PetersHeRotor::PetersHeRotor(unsigned int uLabel, const DofOwner* pDO)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
dVConst(0), dVSine(0), dVCosine(0),
dL11(0.), dL13(0.), dL22(0.), dL31(0.), dL33(0.)
{
	NO_OP;
}

PetersHeRotor::PetersHeRotor(unsigned int uLabel,
	const DofOwner* pDO,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	const doublereal& dVConstTmp,
	const doublereal& dVSineTmp,
	const doublereal& dVCosineTmp,
	flag fOut)
: Elem(uLabel, flag(0)),
Rotor(uLabel, pDO),
dVConst(0), dVSine(0), dVCosine(0),
dL11(0.), dL13(0.), dL22(0.), dL31(0.), dL33(0.)
{
	Init(pCraft, rrot, pRotor, pGround, ppres, dOR, dR,
		iMaxIt, dTol, dE, dCH, dCFF,
		dVConstTmp, dVSineTmp, dVCosineTmp, fOut);
}

void
PetersHeRotor::Init(const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor,
	const StructNode* pGround,
	ResForceSet **ppres,
	const doublereal& dOR,
	const doublereal& dR,
	unsigned int iMaxIt,
	const doublereal& dTol,
	const doublereal& dE,
	const doublereal& dCH,
	const doublereal& dCFF,
	const doublereal& dVConstTmp,
	const doublereal& dVSineTmp,
	const doublereal& dVCosineTmp,
	flag fOut)
{
	ASSERT(dOR > 0.);
	ASSERT(dR > 0.);

	Rotor::Init(pCraft, rrot, pRotor, pGround, ppres, dR, iMaxIt, dTol, dE, fOut);

	dVConst = dVConstTmp;
	dVSine = dVSineTmp;
	dVCosine = dVCosineTmp;

	dOmegaRef = dOR;
	dVTipRef = dOmegaRef*dRadius;
	dArea = M_PI*dRadius*dRadius;

	dHoverCorrection = dCH;
	dForwardFlightCorrection = dCFF;

	/* Significa che valuta la velocita' indotta media al passo corrente */
	dWeight = 0.;

#ifdef USE_MPI
	if (is_parallel) {
		SAFENEWARR(pBlockLenght, int, 20);
		SAFENEWARR(pDispl, MPI::Aint, 20);
		for (int i = 0; i < 20; i++) {
			pBlockLenght[i] = 1;
		}
		for (int i = 0; i < 3; i++) {
			pDispl[i] = MPI::Get_address(RRot3.pGetVec()+i);
		}
		pDispl[3] = MPI::Get_address(&dVConst);
		pDispl[4] = MPI::Get_address(&dVSine);
		pDispl[5] = MPI::Get_address(&dVCosine);
		pDispl[6] = MPI::Get_address(&dOmega);
		pDispl[7] = MPI::Get_address(&dPsi0);
		for (int i = 8; i <= 10; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(Res.Pole().pGetVec()+i-8));
		}
		for (int i = 11; i < 20; i++) {
			pDispl[i] = MPI::Get_address(const_cast<doublereal*>(RRotTranspose.pGetMat()+i-11));
		}
		SAFENEWWITHCONSTRUCTOR(pIndVelDataType, MPI::Datatype,
				MPI::Datatype(MPI::DOUBLE.Create_hindexed(20, pBlockLenght, pDispl)));
		pIndVelDataType->Commit();
	}
#endif /* USE_MPI */
}


PetersHeRotor::~PetersHeRotor(void)
{
#ifdef USE_MPI
	SAFEDELETEARR(pBlockLenght);
	SAFEDELETEARR(pDispl);
	SAFEDELETE(pIndVelDataType);
#endif /* USE_MPI */
}


void
PetersHeRotor::Output(OutputHandler& OH) const
{
     	/*
	 * FIXME: posso usare dei temporanei per il calcolo della trazione
	 * totale per l'output, cosi' evito il giro dei cast
	 */
	if (bToBeOutput()) {
#ifdef USE_MPI
		if (is_parallel && IndVelComm.Get_size() > 1) {
			if (IndVelComm.Get_rank() == 0) {
				Vec3 TmpF(pTmpVecR), TmpM(pTmpVecR+3);

				OH.Rotors()
					<< std::setw(8) << GetLabel()	/* 1 */
					<< " " << RRotTranspose*TmpF /* 2-4 */
					<< " " << RRotTranspose*TmpM /* 5-7 */
					<< " " << dUMean	/* 8 */
					<< " " << dVelocity	/* 9 */
					<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
					<< " " << dMu		/* 11 */
					<< " " << dLambda	/* 12 */
					<< " " << dChi		/* 13 */
					<< " " << dPsi0		/* 14 */
					<< " " << bUMeanRefConverged /* 15 */
					<< " " << iCurrIter	/* 16 */
					<< " " << dVConst	/* 17 */
					<< " " << dVSine	/* 18 */
					<< " " << dVCosine	/* 19 */
					<< std::endl;

				for (int i = 0; ppRes && ppRes[i]; i++) {
					Vec3 TmpF(pTmpVecR+6+6*i);
					Vec3 TmpM(pTmpVecR+9+6*i);

					OH.Rotors()
						<< std::setw(8) << GetLabel()
						<< ":" << ppRes[i]->GetLabel()
						<< " " << TmpF
						<< " " << TmpM
						<< std::endl;
				}
			}
		} else {
			OH.Rotors()
				<< std::setw(8) << GetLabel()	/* 1 */
				<< " " << RRotTranspose*Res.Force()  /* 2-4 */
				<< " " << RRotTranspose*Res.Moment() /* 5-7 */
				<< " " << dUMean	/* 8 */
				<< " " << dVelocity	/* 9 */
				<< " " << atan2(dSinAlphad, dCosAlphad)	/* 10 */
				<< " " << dMu		/* 11 */
				<< " " << dLambda	/* 12 */
				<< " " << dChi		/* 13 */
				<< " " << dPsi0		/* 14 */
				<< " " << bUMeanRefConverged /* 15 */
				<< " " << iCurrIter	/* 16 */
				<< " " << dVConst	/* 17 */
				<< " " << dVSine	/* 18 */
				<< " " << dVCosine	/* 19 */
				<< std::endl;

			for (int i = 0; ppRes && ppRes[i]; i++) {
				OH.Rotors()
					<< std::setw(8) << GetLabel()
	    				<< ":" << ppRes[i]->GetLabel()
					<< " " << ppRes[i]->pRes->Force()
					<< " " << ppRes[i]->pRes->Moment()
					<< std::endl;
			}
		}

#else /* !USE_MPI */
		OH.Rotors()
			<< std::setw(8) << GetLabel()	/* 1 */
			<< " " << RRotTranspose*Res.Force()	/* 2-4 */
			<< " " << RRotTranspose*Res.Moment()	/* 5-7 */
			<< " " << dUMean	/* 8 */
			<< " " << dVelocity	/* 9 */
			<< " " << atan2(dSinAlphad,dCosAlphad)	/* 10 */
			<< " " << dMu		/* 11 */
			<< " " << dLambda	/* 12 */
			<< " " << dChi		/* 13 */
			<< " " << dPsi0		/* 14 */
			<< " " << bUMeanRefConverged /* 15 */
			<< " " << iCurrIter	/* 16 */
			<< " " << dVConst	/* 17 */
			<< " " << dVSine	/* 18 */
			<< " " << dVCosine	/* 19 */
			<< std::endl;

		for (int i = 0; ppRes && ppRes[i]; i++) {
			OH.Rotors()
				<< std::setw(8) << GetLabel()
    				<< ":" << ppRes[i]->GetLabel()
				<< " " << ppRes[i]->pRes->Force()
				<< " " << ppRes[i]->pRes->Moment()
				<< std::endl;
		}
#endif /* !USE_MPI */
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
PetersHeRotor::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering PetersHeRotor::AssJac()" << std::endl);

	WorkMat.SetNullMatrix();

#ifdef USE_MPI
	if (is_parallel && IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
		SparseSubMatrixHandler& WM = WorkMat.SetSparse();
		integer iFirstIndex = iGetFirstIndex();

		WM.ResizeReset(5, 0);

		WM.PutItem(1, iFirstIndex + 1, iFirstIndex + 1, dM11 + dCoef*dL11);
		WM.PutItem(2, iFirstIndex + 3, iFirstIndex + 1, dCoef*dL31);
		WM.PutItem(3, iFirstIndex + 2, iFirstIndex + 2, dM22 + dCoef*dL22);
		WM.PutItem(4, iFirstIndex + 1, iFirstIndex + 3, dCoef*dL13);
		WM.PutItem(5, iFirstIndex + 3, iFirstIndex + 3, dM33 + dCoef*dL33);
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
PetersHeRotor::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
     	DEBUGCOUT("Entering PetersHeRotor::AssRes()" << std::endl);

#ifdef USE_MPI
     	ExchangeLoads(flag(1));

     	if (!is_parallel || IndVelComm.Get_rank() == 0)
#endif /* USE_MPI */
	{
	   	/* Calcola parametri vari */
	   	Rotor::InitParam();

       	   	WorkVec.Resize(3);

	   	integer iFirstIndex = iGetFirstIndex();

	   	WorkVec.PutRowIndex(1, iFirstIndex + 1);
	   	WorkVec.PutRowIndex(2, iFirstIndex + 2);
	   	WorkVec.PutRowIndex(3, iFirstIndex + 3);

	   	dVConst = XCurr(iFirstIndex + 1);
	   	dVSine = XCurr(iFirstIndex + 2);
	   	dVCosine = XCurr(iFirstIndex + 3);

	   	doublereal dVConstPrime = XPrimeCurr(iFirstIndex + 1);
	   	doublereal dVSinePrime = XPrimeCurr(iFirstIndex + 2);
	   	doublereal dVCosinePrime = XPrimeCurr(iFirstIndex + 3);

		doublereal dCT = 0.;
		doublereal dCl = 0.;
		doublereal dCm = 0.;

	   	/*
		 * Attenzione: moltiplico tutte le equazioni per dOmega
		 * (ovvero, i coefficienti CT, CL e CM sono divisi
		 * per dOmega anziche' dOmega^2)
		 */

		dL11 = 0.;
		dL13 = 0.;
		dL22 = 0.;
		dL31 = 0.;
		dL33 = 0.;

	   	doublereal dDim = dGetAirDensity(GetXCurr())*dArea*dOmega*(dRadius*dRadius);
	   	if (dDim > std::numeric_limits<doublereal>::epsilon()) {
			/*
			 * From Claudio Monteggia:
			 *
			 *                        Ct
			 * Um = -------------------------------------
			 *       sqrt( lambda^2 / KH^4 + mu^2 / KF^2)
			 */
		 	doublereal dLambdaTmp
				= dLambda/(dHoverCorrection*dHoverCorrection);
		 	doublereal dMuTmp = dMu/dForwardFlightCorrection;

		 	doublereal dVT
				= sqrt(dLambdaTmp*dLambdaTmp + dMuTmp*dMuTmp);
		 	doublereal dVm = 0.;
		 	if (dVT > dVTipTreshold*dVTipRef) {
		       		dVm = (dMuTmp*dMuTmp + dLambdaTmp*(dLambdaTmp + dVConst))/dVT;
		 	}

			/*
			 * dUMean is just for output;
			 */
		   	dUMean = dVConst*dOmega*dRadius;

		   	/* Trazione nel sistema rotore */
		   	doublereal dT = RRot3*Res.Force();

		   	/* Momento nel sistema rotore-vento */
		   	doublereal dCosP = cos(dPsi0);
		   	doublereal dSinP = sin(dPsi0);
		   	Mat3x3 RTmp( dCosP, -dSinP, 0.,
		      			dSinP, dCosP, 0.,
		      			0., 0., 1.);

		   	Vec3 M(RTmp*(RRotTranspose*Res.Moment()));

		 	/* Thrust, roll and pitch coefficients */
		 	dCT = dT/dDim;
		 	dDim *= dRadius;
		 	dCl = - M(1)/dDim;
		 	dCm = - M(2)/dDim;

			if (dVT > std::numeric_limits<doublereal>::epsilon()
				&& dVm > std::numeric_limits<doublereal>::epsilon())
			{

				/* Matrix coefficients */
				/* FIXME: divide by 0? */
			 	doublereal dl11 = .5/dVT;
				/* FIXME: divide by 0? */
			 	doublereal d = 15./64.*M_PI*tan(dChi/2.);
				/* FIXME: divide by 0? */
			 	doublereal dl13 = d/dVm;
				/* FIXME: divide by 0? */
			 	doublereal dl31 = d/dVT;

			 	doublereal dCosChi2 = cos(dChi/2.);
			 	d = 2.*dCosChi2*dCosChi2;
				/* FIXME: divide by 0? */
			 	doublereal dl22 = -4./(d*dVm);
				/* FIXME: divide by 0? */
				doublereal dl33 = -4.*(d - 1)/(d*dVm);
	
				d = dl11*dl33 - dl31*dl13;
				/* FIXME: divide by 0? */
				dL11 = dOmega*dl33/d;
				dL31 = -dOmega*dl31/d;
				dL13 = -dOmega*dl13/d;
				dL33 = dOmega*dl11/d;
				dL22 = dOmega/dl22;
		   	}
		}
	
#ifdef DEBUG
	   	/* Prova: */
		static int i = -1;
		int iv[] = { 0, 1, 0, -1, 0 };
		if (++i == 4) {
			i = 0;
		}
	   	Vec3 XTmp(pRotor->GetXCurr()+pCraft->GetRCurr()*Vec3(dRadius*iv[i],dRadius*iv[i+1],0.));
	   	doublereal dPsiTmp, dXTmp;
	   	GetPos(XTmp, dXTmp, dPsiTmp);
	   	Vec3 IndV = GetInducedVelocity(GetElemType(), GetLabel(), 0, XTmp);
	   	std::cout
		 	<< "X rotore:  " << pRotor->GetXCurr() << std::endl
		 	<< "V rotore:  " << VCraft << std::endl
		 	<< "X punto:   " << XTmp << std::endl
		 	<< "Omega:     " << dOmega << std::endl
		 	<< "Velocita': " << dVelocity << std::endl
		 	<< "Psi0:      " << dPsi0 << " ("
				<< dPsi0*180./M_PI << " deg)" << std::endl
		 	<< "Psi punto: " << dPsiTmp << " ("
				<< dPsiTmp*180./M_PI << " deg)" << std::endl
		 	<< "Raggio:    " << dRadius << std::endl
		 	<< "r punto:   " << dXTmp << std::endl
		 	<< "mu:        " << dMu << std::endl
		 	<< "lambda:    " << dLambda << std::endl
		 	<< "cos(ad):   " << dCosAlphad << std::endl
		 	<< "sin(ad):   " << dSinAlphad << std::endl
		 	<< "UMean:     " << dUMean << std::endl
		 	<< "iv punto:  " << IndV << std::endl;
#endif /* DEBUG */

		WorkVec.PutCoef(1, dCT - dM11*dVConstPrime
				- dL11*dVConst - dL13*dVCosine);
	 	WorkVec.PutCoef(2, dCl - dM22*dVSinePrime
				- dL22*dVSine);
	 	WorkVec.PutCoef(3, dCm - dM33*dVCosinePrime
				- dL31*dVConst - dL33*dVCosine);

#ifdef USE_MPI
     	} else {
	   	WorkVec.Resize(0);
#endif /* USE_MPI */
     	}

#ifdef USE_MPI
	ExchangeVelocity();
#endif /* USE_MPI */

	/* Ora la trazione non serve piu' */
	ResetForce();

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Done();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

     	return WorkVec;
}

/* Relativo ai ...WithDofs */
void
PetersHeRotor::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}


/* Relativo ai ...WithDofs */
void
PetersHeRotor::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	integer iFirstIndex = iGetFirstIndex();

	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		XP.PutCoef(iFirstIndex + iCnt, 0.);
	}

	X.PutCoef(iFirstIndex + 1, dVConst);
	X.PutCoef(iFirstIndex + 2, dVSine);
	X.PutCoef(iFirstIndex + 3, dVCosine);
}


/* Restart */
std::ostream&
PetersHeRotor::Restart(std::ostream& out) const
{
	return Rotor::Restart(out) << "dynamic inflow, " << dRadius
		<< ", correction, " << dHoverCorrection
		<< ", " << dForwardFlightCorrection << ';' << std::endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void
PetersHeRotor::AddForce(const Elem *pEl, const StructNode *pNode,
	const Vec3& F, const Vec3& M, const Vec3& X)
{
	/*
	 * Gli serve la trazione ed il momento rispetto al rotore,
	 * che si calcola da se'
	 */
#ifdef USE_MPI
	if (ReqV != MPI::REQUEST_NULL) {
		while (!ReqV.Test()) {
			MYSLEEP(mysleeptime);
		}
	}
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_lock(&forces_mutex);
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	Res.AddForces(F, M, X);
	if (bToBeOutput()) {
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	pthread_mutex_unlock(&forces_mutex);
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
}


/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3
PetersHeRotor::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */

	doublereal dr, dp;
	GetPos(X, dr, dp);

	return RRot3*((dVConst + dr*(dVCosine*cos(dp) + dVSine*sin(dp)))*dRadius*dOmega);
};

/* PetersHeRotor - end */


/* Legge un rotore */
Elem *
ReadRotor(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel)
{
     	DEBUGCOUT("Entering ReadRotor()" << std::endl);

	/* demote to pedantic; syntax changed a long ago... */
     	pedantic_cout("WARNING: the syntax changed; use a comma ',' "
     			"instead of a colon ':' after the keyword "
     			"\"induced velocity\"" << std::endl);

     	const char* sKeyWords[] = {
	  	"induced" "velocity",
			"no",
			"uniform",
				"uniform" "sectional",
			"glauert",
			"mangler",
			"dynamic" "inflow",

		NULL
     	};

     	/* enum delle parole chiave */
     	enum KeyWords {
	  	UNKNOWN = -1,
		INDUCEDVELOCITY = 0,
			NO,
			UNIFORM,
				UNIFORM_SECTIONAL,
			GLAUERT,
			MANGLER,
			DYNAMICINFLOW,
			PETERS_HE,

		LASTKEYWORD
     	};

     	/* tabella delle parole chiave */
     	KeyTable K(HP, sKeyWords);

     	/* aircraft node */
	const StructNode* pCraft = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	/* rotor orientation with respect to aircraft */
     	Mat3x3 rrot(Eye3);
     	if (HP.IsKeyWord("orientation")) {
     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);

     	} else if (HP.IsKeyWord("hinge")) {
		silent_cerr("InducedVelocity(" << uLabel << "): deprecated keyword \"hinge\"; use \"orientation\" instead at line " << HP.GetLineData() << std::endl);

     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);
     	}

     	/* rotor node */
     	const StructNode* pRotor = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	KeyWords InducedType = NO;
     	if (HP.IsArg() && HP.IsKeyWord("induced" "velocity")) {
      		InducedType = KeyWords(HP.GetWord());
     	}

     	Elem* pEl = 0;
     	ResForceSet **ppres = 0;

     	switch (InducedType) {
	case NO: {
		DEBUGCOUT("No induced velocity is considered" << std::endl);

	 	doublereal dR = 0.;
	 	if (HP.IsKeyWord("radius")) {
	      		dR = HP.GetReal();
			if (dR <= 0) {
				silent_cerr("Rotor(" << uLabel << "): "
					"invalid radius " << dR
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
	 	}

	 	ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);

	 	SAFENEWWITHCONSTRUCTOR(pEl,
  				NoRotor,
  				NoRotor(uLabel, pDO, pCraft, rrot, pRotor,
  					ppres, dR, fOut));
		} break;

    	case UNIFORM:
    	case UNIFORM_SECTIONAL:
    	case GLAUERT:
    	case MANGLER:
    	case DYNAMICINFLOW: {
		GlauertRotor::Type gtype(GlauertRotor::GLAUERT);
		if (InducedType == GLAUERT && HP.IsKeyWord("type")) {
			if (HP.IsKeyWord("glauert")) {
				gtype = GlauertRotor::GLAUERT;

			} else if (HP.IsKeyWord("coleman")) {
				gtype = GlauertRotor::COLEMAN_ET_AL;

			} else if (HP.IsKeyWord("drees")) {
				gtype = GlauertRotor::DREES_1;

			} else if (HP.IsKeyWord("payne")) {
				gtype = GlauertRotor::PAYNE;

			} else if (HP.IsKeyWord("white" "and" "blake")) {
				gtype = GlauertRotor::WHITE_AND_BLAKE;

			} else if (HP.IsKeyWord("pitt" "and" "peters")) {
				gtype = GlauertRotor::PITT_AND_PETERS;

			} else if (HP.IsKeyWord("howlett")) {
				gtype = GlauertRotor::HOWLETT;

			} else if (HP.IsKeyWord("drees" "2")) {
				silent_cerr("warning, \"drees 2\" deprecated at line " << HP.GetLineData() << "; use \"drees\" instead" << std::endl);
				gtype = GlauertRotor::DREES_2;

			} else {
				silent_cerr("Rotor(" << uLabel << "): "
					"unknown variant of Glauert's induced velocity at line "
					<< HP.GetLineData() << std::endl);
	      			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}

		doublereal dOR = HP.GetReal();
	 	DEBUGCOUT("Reference rotation speed: " << dOR << std::endl);
	 	if (dOR <= 0.) {
	      		silent_cerr("Rotor(" << uLabel << "): "
				"invalid reference speed " << dOR
				<< " at line " << HP.GetLineData()
				<< std::endl);
	      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	 	}

	 	doublereal dR = HP.GetReal();
	 	DEBUGCOUT("Radius: " << dR << std::endl);
	 	if (dR <= 0.) {
	      		silent_cerr("Rotor(" << uLabel << "): "
				"invalid radius " << dR
				<< " at line " << HP.GetLineData()
				<< std::endl);
	      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	 	}

		// optional parameters
		bool bGotInitialValues(false);
	 	doublereal dVConst = 0.;
	 	doublereal dVSine = 0.;
	 	doublereal dVCosine = 0.;
	 	DriveCaller *pdW = 0;
		const StructNode *pGround = 0;
		unsigned iMaxIter = unsigned(-1);
		doublereal dTolerance = std::numeric_limits<double>::max();
		doublereal dEta = -1.;
	 	doublereal dCH = -1.;
	 	doublereal dCFF = -1.;

		while (HP.IsArg()) {
			if (HP.IsKeyWord("ground")) {
				/*
				 * ground node for ground effect modeling
				 */
				if (pGround != 0) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"ground\" node "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* ground node */
     				pGround = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

			} else if (HP.IsKeyWord("initial" "value")) {
	 			if (InducedType != DYNAMICINFLOW) {
					silent_cerr("Rotor(" << uLabel << "): "
						"invalid parameter \"initial value\" "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bGotInitialValues) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing \"initial value\" another time "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

		   		dVConst = HP.GetReal();
		   		dVSine = HP.GetReal();
		   		dVCosine = HP.GetReal();

				bGotInitialValues = true;

	      		} else if (HP.IsKeyWord("weight") || HP.IsKeyWord("delay")) {
	 			if (InducedType == DYNAMICINFLOW) {
					silent_cerr("Rotor(" << uLabel << "): "
						"invalid parameter \"delay\" "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (pdW != 0) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"delay\" driver "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

	      			/*
				 * Legge il coefficiente di peso della velocita'
				 * indotta ("weight" e' deprecato, si preferisce
				 * "delay")
				 *
				 * nota:
				 *
				 * U = U_n * ( 1 - dW ) + U_n-1 * dW
				 *
				 * quindi dW rappresenta il peso che si da'
				 * al valore al passo precedente; in questo modo
				 * si introduce un ritardo euristico (attenzione:
				 * il ritardo vero dipende dal passo temporale)
				 * che aiuta ad evitare problemi di convergenza.
				 * Se si desidera un ritardo "fisico", conviene
				 * provare il "Dynamic Inflow".
				 */
		   		pdW = HP.GetDriveCaller();

	      		} else if (HP.IsKeyWord("max" "iterations")) {
				if (iMaxIter != unsigned(-1)) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"max iterations\" value "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* max iterations when computing reference inflow velocity;
				 * after iMaxIter iterations, the current value is accepted
				 * regardless of convergence; thus, 1 reproduces original
				 * behavior */
				int i = HP.GetInt();
				if (i <= 0) {
					silent_cerr("illegal max iterations "
						<< i << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				iMaxIter = i;

	 		} else if (HP.IsKeyWord("tolerance")) {
				if (dTolerance != std::numeric_limits<double>::max()) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"tolerance\" value "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* tolerance when computing reference inflow velocity;
				 * when the difference in inflow velocity between two
				 * iterations is less than tolerance in module, the
				 * cycle breaks */

				dTolerance = HP.GetReal();
				if (dTolerance <= 0.) {
					silent_cerr("illegal tolerance "
						<< dTolerance << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("eta")) {
				if (dEta != -1.) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"eta\" relaxation factor "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* increment factor when computing reference inflow velocity;
				 * only a fraction dEta of the difference between two iterations
				 * is applied */

				dEta = HP.GetReal();
				if (dEta <= 0.) {
					silent_cerr("illegal eta "
						<< dEta << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("correction")) {
				if (dCH != -1.) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"correction\" factor "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

		 		/* Legge la correzione della velocita' indotta */
		 		dCH = HP.GetReal();
		 		DEBUGCOUT("Hover correction: " << dCH << std::endl);
		 		if (dCH <= 0.) {
		 			silent_cerr("Rotor(" << uLabel << "): "
						"illegal null or negative hover inflow correction "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		 		}

	 			dCFF = HP.GetReal();
				DEBUGCOUT("Forward-flight correction: " << dCFF << std::endl);
	 			if (dCFF <= 0.) {
	 				silent_cerr("Rotor(" << uLabel << "): "
						"illegal null or negative forward-flight inflow correction "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	 			}

			} else {
				break;
		 	}
		}

	 	ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);

		// check consistency and initialize defaults
		if (InducedType != DYNAMICINFLOW && pdW == 0) {
		   	SAFENEW(pdW, NullDriveCaller);
		}

		if (iMaxIter == unsigned(-1)) {
			iMaxIter = 1;

		} else {
			if (dTolerance == std::numeric_limits<double>::max()) {
				silent_cerr("Rotor(" << uLabel << "): "
					"warning, \"max iterations\" is meaningless with default tolerance"
					" at line " << HP.GetLineData()
					<< std::endl);
			}
		}

		if (dEta == -1.) {
			dEta = 1.;
		}

		if (dCH == -1.) {
	 		dCH = 1.;
	 		dCFF = 1.;
		}

		// create element
      		switch (InducedType) {
		case UNIFORM:
	  		DEBUGCOUT("Uniform induced velocity" << std::endl);
			SAFENEWWITHCONSTRUCTOR(pEl,
   					UniformRotor,
   					UniformRotor(uLabel, pDO, pCraft, rrot,
   						pRotor, pGround,
   						ppres, dOR, dR, pdW,
						iMaxIter, dTolerance, dEta,
						dCH, dCFF,
   						fOut));
	  		break;

		case UNIFORM_SECTIONAL:
	  		DEBUGCOUT("Uniform induced velocity" << std::endl);
			SAFENEWWITHCONSTRUCTOR(pEl,
   					UniformRotor2,
   					UniformRotor2(uLabel, pDO, pCraft, rrot,
   						pRotor, pGround,
   						ppres, dOR, dR, pdW,
						iMaxIter, dTolerance, dEta,
						dCH, dCFF,
   						fOut));
	  		break;

     		case GLAUERT:
	  		DEBUGCOUT("Glauert induced velocity" << std::endl);
	  		SAFENEWWITHCONSTRUCTOR(pEl,
   					GlauertRotor,
   					GlauertRotor(uLabel, pDO, pCraft, rrot,
   						pRotor, pGround,
   						ppres, dOR, dR, pdW,
						iMaxIter, dTolerance, dEta,
						dCH, dCFF, gtype,
   						fOut));
	  		break;

     		case MANGLER:
	  		DEBUGCOUT("Mangler induced velocity" << std::endl);

	  		SAFENEWWITHCONSTRUCTOR(pEl,
   					ManglerRotor,
   					ManglerRotor(uLabel, pDO, pCraft, rrot,
   						pRotor, pGround,
   						ppres, dOR, dR, pdW,
						iMaxIter, dTolerance, dEta,
						dCH, dCFF,
   						fOut));
	  		break;

     		case DYNAMICINFLOW:
	  		DEBUGCOUT("Dynamic inflow" << std::endl);

	  		SAFENEWWITHCONSTRUCTOR(pEl,
       					DynamicInflowRotor,
       					DynamicInflowRotor(uLabel, pDO,
						pCraft, rrot, pRotor,
						pGround, ppres,
						dOR, dR,
						iMaxIter, dTolerance, dEta,
						dCH, dCFF,
						dVConst, dVSine, dVCosine,
						fOut));
			break;

     		default:
			ASSERTMSG(0, "You shouldn't have reached this point");
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
      		}
	 	break;
	}

    	case PETERS_HE: {
		doublereal dOR = HP.GetReal();
	 	DEBUGCOUT("Reference rotation speed: " << dOR << std::endl);
	 	if (dOR <= 0.) {
	      		silent_cerr("Rotor(" << uLabel << "): "
				"invalid reference speed " << dOR
				<< " at line " << HP.GetLineData()
				<< std::endl);
	      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	 	}

	 	doublereal dR = HP.GetReal();
	 	DEBUGCOUT("Radius: " << dR << std::endl);
	 	if (dR <= 0.) {
	      		silent_cerr("Rotor(" << uLabel << "): "
				"invalid radius " << dR
				<< " at line " << HP.GetLineData()
				<< std::endl);
	      		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	 	}

		// optional parameters
		bool bGotInitialValues(false);
	 	doublereal dVConst = 0.;
	 	doublereal dVSine = 0.;
	 	doublereal dVCosine = 0.;
		const StructNode *pGround = 0;
		unsigned iMaxIter = unsigned(-1);
		doublereal dTolerance = std::numeric_limits<double>::max();
		doublereal dEta = -1.;
	 	doublereal dCH = -1.;
	 	doublereal dCFF = -1.;

		while (HP.IsArg()) {
			if (HP.IsKeyWord("ground")) {
				/*
				 * ground node for ground effect modeling
				 */
				if (pGround != 0) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"ground\" node "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* ground node */
     				pGround = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

			} else if (HP.IsKeyWord("initial" "value")) {
				if (bGotInitialValues) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing \"initial value\" another time "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

		   		dVConst = HP.GetReal();
		   		dVSine = HP.GetReal();
		   		dVCosine = HP.GetReal();

				bGotInitialValues = true;

	      		} else if (HP.IsKeyWord("max" "iterations")) {
				if (iMaxIter != unsigned(-1)) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"max iterations\" value "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* max iterations when computing reference inflow velocity;
				 * after iMaxIter iterations, the current value is accepted
				 * regardless of convergence; thus, 1 reproduces original
				 * behavior */
				int i = HP.GetInt();
				if (i <= 0) {
					silent_cerr("illegal max iterations "
						<< i << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				iMaxIter = i;

	 		} else if (HP.IsKeyWord("tolerance")) {
				if (dTolerance != std::numeric_limits<double>::max()) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"tolerance\" value "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* tolerance when computing reference inflow velocity;
				 * when the difference in inflow velocity between two
				 * iterations is less than tolerance in module, the
				 * cycle breaks */

				dTolerance = HP.GetReal();
				if (dTolerance <= 0.) {
					silent_cerr("illegal tolerance "
						<< dTolerance << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("eta")) {
				if (dEta != -1.) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"eta\" relaxation factor "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				/* increment factor when computing reference inflow velocity;
				 * only a fraction dEta of the difference between two iterations
				 * is applied */

				dEta = HP.GetReal();
				if (dEta <= 0.) {
					silent_cerr("illegal eta "
						<< dEta << " for Rotor(" << uLabel << ")");
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

			} else if (HP.IsKeyWord("correction")) {
				if (dCH != -1.) {
					silent_cerr("Rotor(" << uLabel << "): "
						"providing another \"correction\" factor "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

		 		/* Legge la correzione della velocita' indotta */
		 		dCH = HP.GetReal();
		 		DEBUGCOUT("Hover correction: " << dCH << std::endl);
		 		if (dCH <= 0.) {
		 			silent_cerr("Rotor(" << uLabel << "): "
						"illegal null or negative hover inflow correction "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		 		}

	 			dCFF = HP.GetReal();
				DEBUGCOUT("Forward-flight correction: " << dCFF << std::endl);
	 			if (dCFF <= 0.) {
	 				silent_cerr("Rotor(" << uLabel << "): "
						"illegal null or negative forward-flight inflow correction "
						"at line " << HP.GetLineData()
						<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	 			}

			} else {
				break;
		 	}
		}

	 	ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);

		if (iMaxIter == unsigned(-1)) {
			iMaxIter = 1;

		} else {
			if (dTolerance == std::numeric_limits<double>::max()) {
				silent_cerr("Rotor(" << uLabel << "): "
					"warning, \"max iterations\" is meaningless with default tolerance"
					<< std::endl);
			}
		}

		if (dEta == -1.) {
			dEta = 1.;
		}

		if (dCH == -1.) {
	 		dCH = 1.;
	 		dCFF = 1.;
		}

		// create element
      		DEBUGCOUT("Dynamic inflow" << std::endl);

  		SAFENEWWITHCONSTRUCTOR(pEl,
			PetersHeRotor,
			PetersHeRotor(uLabel, pDO,
				pCraft, rrot, pRotor,
				pGround, ppres,
				dOR, dR,
				iMaxIter, dTolerance, dEta,
				dCH, dCFF,
				dVConst, dVSine, dVCosine,
				fOut));
	 	break;
	}

	default:
		silent_cerr("Rotor(" << uLabel << "): "
			"unknown induced velocity type at line "
	       		<< HP.GetLineData() << std::endl);
	 	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("Rotor(" << uLabel << "): "
			"semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(pEl != 0);
	return pEl;
} /* End of DataManager::ReadRotor() */

