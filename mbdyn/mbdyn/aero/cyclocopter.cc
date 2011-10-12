/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* Elementi di rotore */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>

#include "cyclocopter.h"
#include "dataman.h"

#define EXTENDED_PLOT

/* CyclocopterInflow - begin */

CyclocopterInflow::CyclocopterInflow(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, 
	flag fOut)
: Elem(uL, fOut),
InducedVelocityElem(uL, pDO, pC, ppres, fOut),
pRotor(pR),
RRot(rrot)
{
	NO_OP;	
}

CyclocopterInflow::~CyclocopterInflow(void)
{
	NO_OP;
}

InducedVelocity::Type
CyclocopterInflow::GetInducedVelocityType(void) const
{
	return InducedVelocity::CYCLOCOPTER;
}

void
CyclocopterInflow::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;	
}

void
CyclocopterInflow::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {

                OH.Rotors()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << "0."                	 /* 9 */
                        << " " << "0."                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << "0."                	 /* 14 */
                        << " " << "0."   		 /* 15 */
                        << " " << "0."           	 /* 16 */
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


}

std::ostream&
CyclocopterInflow::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

void
CyclocopterInflow::SetInitialValue(VectorHandler& X)
{
	NO_OP;
}

void
CyclocopterInflow::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

void 
CyclocopterInflow::SetFilterCoefficients( const doublereal dOmegaFilter, 
	const doublereal dDeltaT ) 
{
 /* Butterworth discrete low-pass filter coefficients */
	if( dDeltaT > 0 && dOmegaFilter > 0 ) {
		doublereal dTmp = 4. + 2*sqrt(2)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter;
		a1 = (-8.+2*dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;
		a2 = (4. -2*sqrt(2)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;

		dTmp = dOmegaFilter*dOmegaFilter*dDeltaT*dDeltaT/dTmp;
		b0 = dTmp;
		b1 = 2*dTmp;
		b2 = dTmp;
	} else {
		a1 = 0.;
		a2 = 0.;
		b0 = 1.;
		b1 = 0.;
		b2 = 0.;
	}
}

/* CyclocopterInflow - end */


/* CyclocopterNoInflow - begin */

CyclocopterNoInflow::CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, 
	flag fOut)
: Elem(uL, fOut),
CyclocopterInflow(uL, pDO, pC, rrot, pR, ppres, fOut)
{
	dUindMean = 0.;
	NO_OP;	
}

CyclocopterNoInflow::~CyclocopterNoInflow(void)
{
	NO_OP;
}

SubVectorHandler&
CyclocopterNoInflow::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	if (fToBeOutput() ){
		RRotorTranspose = pCraft->GetRCurr()*RRot;
		RRotorTranspose = RRotorTranspose.Transpose();
	}
	
	ResetForce();
	WorkVec.Resize(0);	

	return WorkVec;
}

void
CyclocopterNoInflow::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}
}

Vec3
CyclocopterNoInflow::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	
	return Zero3;
}

/* CyclocopterNoInflow - end */

/* CyclocopterUniform1D - begin */

CyclocopterUniform1D::CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, 
	const unsigned int& iFlagAve, const doublereal& dR,
	const doublereal& dL, const doublereal& dOmegaFilter,
	const doublereal& dDeltaT, DriveCaller *pdW, 
	flag fOut)
: Elem(uL, fOut),
CyclocopterInflow(uL, pDO, pC, rrot, pR, ppres, fOut)
{
	ASSERT(dR > 0.);	
	ASSERT(dL > 0.);	
	ASSERT(pdW != 0);	

	iFlagAverage = iFlagAve;
	dRadius = dR;
	dSpan = dL;
	dArea = 2*dRadius*dSpan;

	Weight.Set(pdW);
	dWeight = 0.;
	
	dUindMean = 0.;
	dUindMeanPrev = 0.;

	bFlagIsFirstBlade = 1;
	
	dAzimuth = 0.;
	dAzimuthPrev = 0.;
	
	dTz = 0.;
	dTzMean = 0.;

	F = Zero3;
	FMean = Zero3;
	FMeanOut = Zero3;

	iStepCounter = 0;

	SetFilterCoefficients( dOmegaFilter, dDeltaT );
	
	/* ingresso del filtro */
	Uk = 0.;
	Uk_1 = 0.;
	Uk_2 = 0.;
	/* uscita del filtro */
	Yk = 0.;
	Yk_1 = 0.;
	Yk_2 = 0.;

}

CyclocopterUniform1D::~CyclocopterUniform1D(void)
{
	NO_OP;
}

void
CyclocopterUniform1D::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {

                OH.Rotors()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << FMeanOut                	 /* 14 */
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


}

void
CyclocopterUniform1D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = 1;
	/* calcolo la forza media sul giro generata dal rotore */
	dTzMean += dTz;
	FMean = FMean + F;
	iStepCounter++;
	//if( (dAzimuth>0. && dAzimuthPrev<0.) || (dAzimuth<0. && dAzimuthPrev>0.) ){
	if( (dAzimuth>0. && dAzimuthPrev<0.) ){
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (iFlagAverage == 1) {
			dTzMean = dTzMean/iStepCounter;
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = copysign(std::sqrt(std::abs(dTzMean)/(2*dRho*dArea)), dTzMean);
			dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
			dTzMean = 0.;
		}
		FMean = Zero3;
		iStepCounter = 0;
	}
	dAzimuthPrev = dAzimuth;
	
	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;
		
	dUindMeanPrev = dUindMean;

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

SubVectorHandler&
CyclocopterUniform1D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity (Moble version)*/
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRot3 = RRotor.GetVec(3);
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	dTz= RRot3*Res.Force();
	if (iFlagAverage == 0){
		/* filtro le forze */
		Uk = dTz;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		dTz = Yk;	
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = copysign(std::sqrt(std::abs(dTz)/(2*dRho*dArea)), dTz);

		dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
	} 
		 
	
	ResetForce();
	WorkVec.Resize(0);	

	return WorkVec;
}

void
CyclocopterUniform1D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	
	/* colcolo la posizione azimutale della prima pala */
	if (bFlagIsFirstBlade == 1 ) {
		Vec3 XRel(RRotorTranspose*(X-pRotor->GetXCurr()));
		doublereal d1 = XRel.dGet(2);
		doublereal d2 = XRel.dGet(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = 0;
	}


	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform1D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return RRot3*dUindMean;
}

/* CyclocopterUniform1D - end */

/* CyclocopterUniform2D - begin */

CyclocopterUniform2D::CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, 
	const unsigned int& iFlagAve, const doublereal& dR,
	const doublereal& dL, const doublereal& dOmegaFilter,
	const doublereal& dDeltaT, DriveCaller *pdW, 
	flag fOut)
: Elem(uL, fOut),
CyclocopterInflow(uL, pDO, pC, rrot, pR, ppres, fOut)
{
	ASSERT(dR > 0.);	
	ASSERT(dL > 0.);	
	ASSERT(pdW != 0);	

	iFlagAverage = iFlagAve;

	dRadius = dR;
	dSpan = dL;
	dArea = 2*dRadius*dSpan;

	Weight.Set(pdW);
	dWeight = 0.;
	
	dUind = Zero3;
	dUindPrev = Zero3;
	dUindMean = 0.;

	bFlagIsFirstBlade = 1;
	
	dAzimuth = 0.;
	dAzimuthPrev = 0.;
	
	F = Zero3;
	FMean = Zero3;
	FMeanOut = Zero3;

	iStepCounter = 0;
	
	SetFilterCoefficients( dOmegaFilter, dDeltaT );

	/* ingresso del filtro */
	Uk = Zero3;
	Uk_1 = Zero3;
	Uk_2 = Zero3;
	/* uscita del filtro */
	Yk = Zero3;
	Yk_1 = Zero3;
	Yk_2 = Zero3;

}

CyclocopterUniform2D::~CyclocopterUniform2D(void)
{
	NO_OP;
}

void
CyclocopterUniform2D::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {

                OH.Rotors()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << dUind                	 /* 11-13 */
                        << " " << FMeanOut                	 /* 14-16 */
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

}

void
CyclocopterUniform2D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = 1;
#if 0
	if (iFlagAverage == 1) {
		/* calcolo la forza media sul giro generata dal rotore */
		FMean = FMean + F;
		iStepCounter++;
		//if( (dAzimuth>0. && dAzimuthPrev<0.) || (dAzimuth<0. && dAzimuthPrev>0.) ){
		if( (dAzimuth>0. && dAzimuthPrev<0.) ){
			FMean = FMean/iStepCounter;
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT= sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema 
	 		* rotore */
			dUind = 0.;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

			FMean = 0.;
			iStepCounter = 0;
		}
	}
#endif
	/* calcolo la forza media sul giro generata dal rotore */
	FMean = FMean + F;
	iStepCounter++;
	if( (dAzimuth>0. && dAzimuthPrev<0.) ){
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (iFlagAverage == 1) {
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT= sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema 
	 		* rotore */
			dUind = Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
		}

		FMean = Zero3;
		iStepCounter = 0;
	}


	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

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

SubVectorHandler&
CyclocopterUniform2D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	if (iFlagAverage == 0){
		/* filtro le forze */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;	
		/* Forza nel piano normale all'asse di rotazione */
		doublereal dT= sqrt(F(2)*F(2) + F(3)*F(3));
		/* Velocità indotta: calcolata in base alla dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = sqrt(dT/(2*dRho*dArea));
		/* Componenti della velocità indotta nel sistema 
	 	* rotore */
		dUind = Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);		

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
	}
	

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterUniform2D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	
	/* colcolo la posizione azimutale della prima pala */
	//if (bFlagIsFirstBlade == 1 && iFlagAverage == 1) {
	if (bFlagIsFirstBlade == 1) {
		Vec3 XRel(RRotorTranspose*(X-pRotor->GetXCurr()));
		doublereal d1 = XRel.dGet(2);
		doublereal d2 = XRel.dGet(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = 0;
	}

	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform2D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	//printf("%f %f %f\n",dUind(1),dUind(2),dUind(3));
	return RRotor*dUind;
}

/* CyclocopterUniform2D - end */

/* CyclocopterPolimi - begin */

CyclocopterPolimi::CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, 
	const unsigned int& iFlagAve, const doublereal& dR,
	const doublereal& dL, const doublereal& dOmegaFilter,
	const doublereal& dDeltaT, DriveCaller *pdW, 
	flag fOut)
: Elem(uL, fOut),
CyclocopterInflow(uL, pDO, pC, rrot, pR, ppres, fOut)
{
	ASSERT(dR > 0.);	
	ASSERT(dL > 0.);	
	ASSERT(pdW != 0);	

	iFlagAverage = iFlagAve;

	dRadius = dR;
	dSpan = dL;
	dArea = 2*dRadius*dSpan;

	Weight.Set(pdW);
	dWeight = 0.;
	
	dUind = Zero3;
	dUindPrev = Zero3;
	dUindMean = 0.;

	dXi = 0.;

	bFlagIsFirstBlade = 1;
	
	dAzimuth = 0.;
	dAzimuthPrev = 0.;
	
	F = Zero3;
	FMean = Zero3;
	FMeanOut = Zero3;

	iStepCounter = 0;
	
	SetFilterCoefficients( dOmegaFilter, dDeltaT );
	
	/* ingresso del filtro */
	Uk = Zero3;
	Uk_1 = Zero3;
	Uk_2 = Zero3;
	/* uscita del filtro */
	Yk = Zero3;
	Yk_1 = Zero3;
	Yk_2 = Zero3;
	
}

CyclocopterPolimi::~CyclocopterPolimi(void)
{
	NO_OP;
}

void
CyclocopterPolimi::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = 1;
	/* calcolo la forza media sul giro generata dal rotore */
	FMean = FMean + F;;
	iStepCounter++;
	//if( (dAzimuth>0. && dAzimuthPrev<0.) || (dAzimuth<0. && dAzimuthPrev>0.) ){
	if( (dAzimuth>0. && dAzimuthPrev<0.) ){
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (iFlagAverage == 1) {
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT= sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema 
	 		* rotore */
			dUind = Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
			/* angolo di cui è ruotata la trazione */
			dXi = atan2(FMean(3), FMean(2)) - M_PI/2.;
		}

		FMean = Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

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
CyclocopterPolimi::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {

                OH.Rotors()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean              /* 8 */
                        << " " << dUind                	 /* 9 -11*/
                        << " " << dXi                	 /* 12 */
                        << " " << dAzimuth               /* 13 */
                        << " " << FMeanOut                	 /* 14-16 */
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


}

SubVectorHandler&
CyclocopterPolimi::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	if (iFlagAverage == 0){
		/* filtro le forze */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;	
		/* Forza nel piano normale all'asse di rotazione */
		doublereal dT= sqrt(F(2)*F(2) + F(3)*F(3));
		/* Velocità indotta: calcolata in base alla dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = sqrt(dT/(2*dRho*dArea));
		/* Componenti della velocità indotta nel sistema 
	 	* rotore */
		dUind = Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
		/* angolo di cui è ruotata la trazione */
		dXi = atan2(F(3), F(2)) - M_PI/2.;
	}
	
	ResetForce();
	WorkVec.Resize(0);	

	return WorkVec;
}

void
CyclocopterPolimi::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	
	/* colcolo la posizione azimutale della prima pala */
	if (bFlagIsFirstBlade == 1) {
		Vec3 XRel(RRotorTranspose*(X-pRotor->GetXCurr()));
		doublereal d1 = XRel.dGet(2);
		doublereal d2 = XRel.dGet(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = 0;
	}

	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F,M,X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterPolimi::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	Vec3 XRel(RRotorTranspose*(X-pRotor->GetXCurr()));

	doublereal d1 = XRel.dGet(2);
	doublereal d2 = XRel.dGet(3);

	/* dPsi0 non serve a nulla perchè uso l'angolo
	 * relativo: (dp-dXi)!!! */
	doublereal dpp = atan2(d2, d1);

	doublereal r = sqrt(d1*d1+d2*d2)*cos(dpp-dXi);
	
	return RRotor*((dUind*(M_PI/2.))*cos((M_PI/2.)*(r/dRadius)));

}

/* CyclocopterPolimi - end */

/* CyclocopterKARI - begin */

CyclocopterKARI::CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
	const StructNode* pC, const Mat3x3& rrot,
	const StructNode* pR, ResForceSet **ppres, flag fOut)
: Elem(uL, fOut),
InducedVelocityElem(uL, pDO, pC, ppres, fOut),
pRotor(pR),
RRot(rrot)
{	
	NO_OP;
}

CyclocopterKARI::~CyclocopterKARI(void)
{
	NO_OP;
}

InducedVelocity::Type
CyclocopterKARI::GetInducedVelocityType(void) const
{
	return InducedVelocity::CYCLOCOPTER;
}

void
CyclocopterKARI::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterKARI::Output(OutputHandler& OH) const
{
	NO_OP;
}

std::ostream&
CyclocopterKARI::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

void
CyclocopterKARI::SetInitialValue(VectorHandler& X)
{
	NO_OP;
}

SubVectorHandler&
CyclocopterKARI::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{	
	return WorkVec;
}

#if 0
void
CyclocopterKARI::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	NO_OP;
}
#endif

Vec3
CyclocopterKARI::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return Zero3;
}

void
CyclocopterKARI::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

/* CyclocopterKARI - end */

Elem*
ReadCyclocopter(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO, 
	unsigned int uLabel,
	const StructNode* pCraft,
	const Mat3x3& rrot,
	const StructNode* pRotor)
{
	Elem *pEl = 0;

	const char* sKeyWords[] = {
		"type",
		"no",
		"uniform1D",
		"uniform2D",
		"polimi",
		"KARI",
		NULL
	};

	enum KeyWords {
		UNKNOWN = -1,
		type = 0,
		NO,
		uniform1D,
		uniform2D,
		polimi,
		KARI,

		LASTKEYWORD
	};

	KeyTable K(HP, sKeyWords);

	KeyWords CyclocopterInducedType = NO;
	if (HP.IsArg() && HP.IsKeyWord("type")) {
        	CyclocopterInducedType = KeyWords(HP.GetWord());
	}

	switch( CyclocopterInducedType ) {
	case NO: {
		ResForceSet **ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);
		pEl = new CyclocopterNoInflow(uLabel, pDO,
			pCraft, rrot, pRotor,
  			ppres, fOut);
 
		break;
	}
	case uniform1D:	
	case uniform2D:	
	case polimi:	{
		unsigned int iFlagAve = HP.GetInt();
		if ( (iFlagAve != 0) && (iFlagAve != 1) ) {
			silent_cerr("Illegal input "
				"for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal dR = HP.GetReal();
		if (dR <= 0.) {
			silent_cerr("Illegal null or negative radius"
				"for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		doublereal dL = HP.GetReal();
		if (dL <= 0.) {
			silent_cerr("Illegal null or negative blade"
				"length for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	
		DriveCaller *pdW = 0;
		if (HP.IsKeyWord("delay")) {
			pdW = HP.GetDriveCaller();
		} else {
			SAFENEW( pdW, NullDriveCaller);
		}

		doublereal dOmegaFilter = 0.;
		if (HP.IsKeyWord("omegacut")) {
			dOmegaFilter = HP.GetReal();
			if (dOmegaFilter <= 0){
				silent_cerr("Illegal null or negative filter"
					"cut frequency for rotor" << uLabel
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		} else {
			dOmegaFilter = 0.;
		}

		doublereal dDeltaT = 0.;
		if (HP.IsKeyWord("timestep")) {
			dDeltaT = HP.GetReal();
			if (dDeltaT <= 0){
				silent_cerr("Illegal null or negative time"
					"step for rotor" << uLabel
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		} else {
			dDeltaT = 0.;
		}

     		ResForceSet **ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);

		switch( CyclocopterInducedType ) {
		case uniform1D: {
			pEl = new CyclocopterUniform1D(uLabel, pDO,
				pCraft, rrot, pRotor,
  				ppres, iFlagAve, dR, dL,
				dOmegaFilter, dDeltaT, pdW, fOut);
			break;
		}
		case uniform2D: {
			pEl = new CyclocopterUniform2D(uLabel, pDO,
				pCraft, rrot, pRotor,
  				ppres, iFlagAve, dR, dL,
				dOmegaFilter, dDeltaT, pdW, fOut);
			break;
		}
		case polimi: {
			pEl = new CyclocopterPolimi(uLabel, pDO,
				pCraft, rrot, pRotor,
  				ppres, iFlagAve, dR, dL,
				dOmegaFilter, dDeltaT, pdW, fOut);
			break;
		}
		}
		break;
		
	}
	case KARI: {
     		ResForceSet **ppres = ReadResSets(pDM, HP);

	 	flag fOut = pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY);
		pEl = new CyclocopterKARI(uLabel, pDO,
			pCraft, rrot, pRotor,
  			ppres, fOut);
		break;
	}
	default:
		silent_cerr("Rotor(" << uLabel << "): "
			"unknown cyclocopter inflow model "
			"at line " << HP.GetLineData()
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pEl;
}

