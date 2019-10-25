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
 */


/* Bezier Rod */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>
#include <limits>

#include "dataman.h"
#include "gauss.h" 

#include "rodbezj.h"

/* Rod Bezier - begin */

/* Costructor */
RodBezier::RodBezier(unsigned int uL, const DofOwner* pDO,
			const ConstitutiveLaw1D* pCL,
			const StructNode* pN1, const StructNode* pN2,
			const Vec3& fOTmp, const Vec3& fATmp, 
			const Vec3& fBTmp, const Vec3& fITmp,
			doublereal dLength, bool bFromNodes, 
			const unsigned int iIntOrder, const unsigned int iIntSegments,
			flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw1DOwner(pCL),
pNode1(pN1),
pNode2(pN2),
fO(fOTmp),
fA(fATmp),
fB(fBTmp),
fI(fITmp),
dL0(dLength),
l1(Zero3),
l2(Zero3),
dElle(0.),
dEpsilon(0.),
dEpsilonPrime(0.),
iIntOrd(iIntOrder),
iIntSeg(iIntSegments)
#ifdef USE_NETCDFC
,
Var_F2(0),
Var_l(0),
Var_l2(0),
Var_l1(0),
Var_v(0)
#endif // USE_NETCDFC
{
	/* Check initial data consistency */
	ASSERT(pN1 != NULL);
	ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN2 != NULL);
	ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);

	/* Integration sites and weights (for curve length) 
	 * Gauss-Legendre quadrature used. */

	gdi = new GaussDataIterator(iIntOrder);
	
	const Mat3x3& R1(pN1->GetRRef());
	const Mat3x3& R2(pN2->GetRRef());
	
	l1 = R1*(fATmp - fOTmp);
	l2 = R2*(fBTmp - fITmp);

	l1 = l1/l1.Norm();
	l2 = l2/l2.Norm();

	if (bFromNodes) {
		
		const Vec3 b1 = pN1->GetXCurr() + R1*fO; 
		const Vec3 b2 = pN1->GetXCurr() + R1*fA; 
		const Vec3 b3 = pN2->GetXCurr() + R2*fB; 
		const Vec3 b4 = pN2->GetXCurr() + R2*fI; 
	
		Vec3 p;
		doublereal q, up;
		doublereal s = 1.0/(2.0*iIntSegments);

		for (unsigned int jj = 0; jj < iIntSegments; jj++) {
			q = s*(2*jj + 1);

			PntWght gp = gdi->GetFirst();
			do {
				up = s*gp.dGetPnt() + q;
				p = b1*Ap1(up) + b2*Ap2(up) + b3*Ap3(up) + b4*Ap4(up);
				dLength = dLength + s*gp.dGetWght()*p.Norm();

			} while (gdi->fGetNext(gp)); 
		}

		dL0 = dLength;

	}

	ASSERT(dLength > std::numeric_limits<doublereal>::epsilon());

}

/* Destructor */
RodBezier::~RodBezier(void)
{
	NO_OP;
}

/* Restart file contribute */
std::ostream&
RodBezier::Restart(std::ostream& out) const
{
	Joint::Restart(out); 
	out << ", rod bezier, "
		<< pNode1->GetLabel() << ", "
		<< "position, reference, node, "; 
		fO.Write(out, ", ") << ", "
		<< "position, reference, node, "; fA.Write(out, ", ") << "," << std::endl
		<< pNode2->GetLabel() << ", "
		<< "position, reference, node, "; fO.Write(out, ", ") << ", "
		<< "position, reference, node, "; fA.Write(out, ", ") << "," << std::endl
		<< dL0 << "," << std::endl
		<< "integration order, " << iIntOrd << ", "
		<< "integration segments, " << iIntSeg << ", ";
	return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
RodBezier::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw1DOwner::AfterConvergence(dEpsilon, dEpsilonPrime);
}

/* Jacobian matrix contribute */
VariableSubMatrixHandler&
RodBezier::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering RodBezier::AssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Resizes and resets the work matrix */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Gets the nodes indexes */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Sets the matrix indexes */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	const Mat3x3& R1(pNode1->GetRRef());
	const Mat3x3& R2(pNode2->GetRRef());
	const Vec3 fOTmp(R1*fO);
	const Vec3 fATmp(R1*fA);
	const Vec3 fBTmp(R2*fB);
	const Vec3 fITmp(R2*fI);
	
	const Vec3& Omega1(pNode1->GetWRef());
	const Vec3& Omega2(pNode2->GetWRef());

	/* Force and slopes */
	doublereal dF = GetF();
	doublereal dFDE = GetFDE();
	doublereal dFDEPrime = GetFDEPrime();

	const Vec3 b1 = pNode1->GetXCurr() + R1*fO; 
	const Vec3 b2 = pNode1->GetXCurr() + R1*fA; 
	const Vec3 b3 = pNode2->GetXCurr() + R2*fB; 
	const Vec3 b4 = pNode2->GetXCurr() + R2*fI; 
	
	Vec3 p;
	doublereal q, up, pNorm, dKC, dAp1, dAp2, dAp3, dAp4;
	doublereal s = 1.0/(2.0*iIntSeg);

	Vec3 kcx1(Zero3);
	Vec3 kcg1(Zero3);
	Vec3 kcx2(Zero3);
	Vec3 kcg2(Zero3);

	Vec3 TmpV;
	
	for (unsigned int jj = 0; jj < iIntSeg; jj++) {
		q = s*(2*jj + 1);
		
		PntWght gp = gdi->GetFirst();
		do { 
			up = s*gp.dGetPnt() + q;
			
			dAp1 = Ap1(up);
			dAp2 = Ap2(up);
			dAp3 = Ap3(up);
			dAp4 = Ap4(up);

			p = b1*dAp1 + b2*dAp2 +	b3*dAp3 + b4*dAp4;
			pNorm = p.Norm();
			dKC = (s*gp.dGetWght())/(2*dL0*pNorm);
			
			kcx1 = kcx1 + p*(dFDE*dKC*dCoef*(dAp1 + dAp2));
			
			kcg1 = kcg1 + p.Cross(fOTmp*dAp1 + fATmp*dAp2)*dFDE*dKC*dCoef;
			
			kcx2 = kcx2 + p*dFDE*dKC*dCoef*(dAp3 + dAp4);
			
			kcg2 = kcg2 + p.Cross(fBTmp*dAp3 +fITmp*dAp4)*dFDE*dKC*dCoef;

			if (dFDEPrime != 0.) {
				kcx1 = kcx1 + p*dFDEPrime*dKC*(dAp1 + dAp2);
				
				TmpV = p*2. - p.Cross(fOTmp*dAp1 + fATmp*dAp2);
				kcg1 = kcg1 - Omega1.Cross(TmpV)*dFDEPrime*dKC*dCoef;
				
				kcx2 = kcx2 + p*dFDEPrime*dKC*(dAp3 + dAp4);

				TmpV = p*2. - p.Cross(fBTmp*dAp3 + fITmp*dAp4);
				kcg2 = kcg2 + Omega2.Cross(TmpV)*dFDEPrime*dKC*dCoef;
			}

		} while (gdi->fGetNext(gp));
	}

	
	/* Forces on node 1 from dx1 */
	Mat3x3 F1x1 = l1.Tens(kcx1);
	WM.Add(1, 1, F1x1);

	/* Moments on node 1 from dx1 */
	Mat3x3 M1x1 = Mat3x3(MatCross, fOTmp)*F1x1;
	WM.Sub(3 + 1, 1, M1x1);

	/* Forces on node 2 from dx1 */ 
	Mat3x3 F2x1 = l2.Tens(kcx1);
	WM.Add(6 + 1, 1, F2x1);

	/* Moments on node 2 from dx1 */
	Mat3x3 M2x1 = Mat3x3(MatCross, fITmp)*F2x1;
	WM.Sub(9 + 1, 1, M2x1);

	/* Forces on node 1 from dg1 */
	Mat3x3 F1g1 = l1.Tens(kcg1) - (Mat3x3(MatCross, l1) - Mat3x3(MatCross, l1)*l1.Tens())*dCoef*dF;
	WM.Add(1, 3 + 1, F1g1);

	/* Moments on node 1 from dg1 */
	Mat3x3 M1g1 = Mat3x3(MatCross, fOTmp)*F1g1 - Mat3x3(MatCrossCross, l1, fOTmp)*dCoef*dF;
	WM.Sub(3 + 1, 3 + 1, M1g1);

	/* Forces on node 2 from dg1 */
	Mat3x3 F2g1 = l2.Tens(kcg1);
	WM.Add(6 + 1, 3 + 1, F2g1);

	/* Moments on node 2 from dg1 */
	Mat3x3 M2g1 = Mat3x3(MatCross, fITmp)*F2g1;
	WM.Sub(9 + 1, 3 + 1, M2g1);

	/* Forces on node 1 from dx2 */
	Mat3x3 F1x2 = l1.Tens(kcx2);
	WM.Add(1, 6 + 1, F1x2);

	/* Moments on node 1 from dx2 */
	Mat3x3 M1x2 = Mat3x3(MatCross, fOTmp)*F1x2;
	WM.Sub(3 + 1, 6 + 1, M1x2);

	/* Forces on node 2 from dx2 */
	Mat3x3 F2x2 = l2.Tens(kcx2);
	WM.Add(6 + 1, 6 + 1, F2x2);

	/* Moments on node 2 from dx2 */
	Mat3x3 M2x2 = Mat3x3(MatCross, fITmp)*F2x2;
	WM.Sub(9 + 1, 6 + 1, M2x2);

	/* Forces on node 1 from dg2 */
	Mat3x3 F1g2 = l1.Tens(kcg2);
	WM.Add(1, 9 + 1, F1g2);

	/* Moments on node 1 from dg2 */
	Mat3x3 M1g2 = Mat3x3(MatCross, fOTmp)*F1g2;
	WM.Sub(3 + 1, 9 + 1, M1g2);

	/* Forces on node 2 from dg2 */
	Mat3x3 F2g2 = l2.Tens(kcg2) - (Mat3x3(MatCross, l2) - Mat3x3(MatCross, l2)*l2.Tens())*dF*dCoef;
	WM.Add(6 + 1, 9 + 1, F2g2);

	/* Moments on node 2 from dg2 */
	Mat3x3 M2g2 = Mat3x3(MatCross, fITmp)*F2g2 - Mat3x3(MatCrossCross, l2, fITmp)*dCoef*dF;
	WM.Sub(9 + 1, 9 + 1, M2g2);

	return WorkMat;
}

/* Residual contribute */
SubVectorHandler&
RodBezier::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("RodBezier::AssRes()" << std::endl);

	/* Resizes and resets the work matrix */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Gets the nodes indexes */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Sets the indexes */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

void
RodBezier::AssVec(SubVectorHandler& WorkVec)
{
	DEBUGCOUT("RodBezier::AssVec()" << std::endl);

	/* Data */
	const Mat3x3& R1(pNode1->GetRCurr());
	const Mat3x3& R2(pNode2->GetRCurr());
	Vec3 fOTmp(R1*fO);
	Vec3 fATmp(R1*fA);
	Vec3 fBTmp(R2*fB);
	Vec3 fITmp(R2*fI);

	const Vec3& x1(pNode1->GetXCurr());
	const Vec3& x2(pNode2->GetXCurr());

	const Vec3& v1(pNode1->GetVCurr());
	const Vec3& v2(pNode2->GetVCurr());
	const Vec3& Omega1(pNode1->GetWCurr());
	const Vec3& Omega2(pNode2->GetWCurr());

	dElle = 0.;
	doublereal dEllePrime = 0.;
	Vec3 p, pprime;
	doublereal q, up, dAp1, dAp2, dAp3, dAp4;
	doublereal s = 1.0/(2.0*iIntSeg);

	for (unsigned int jj = 0; jj < iIntSeg; jj++) {
		q = s*(2*jj + 1); 
		
		PntWght gp = gdi->GetFirst();
		do {
			up = s*gp.dGetPnt() + q;

			dAp1 = Ap1(up);
			dAp2 = Ap2(up);
			dAp3 = Ap3(up);
			dAp4 = Ap4(up);

			p = (x1 + fOTmp)*dAp1 + 
			    (x1 + fATmp)*dAp2 + 
			    (x2 + fBTmp)*dAp3 + 
			    (x2 + fITmp)*dAp4;

			dElle = dElle + s*gp.dGetWght()*p.Norm();
			pprime = (v1 + Omega1.Cross(fOTmp))*dAp1 + 
				 (v1 + Omega1.Cross(fATmp))*dAp2 +
				 (v2 + Omega2.Cross(fBTmp))*dAp3 + 
				 (v2 + Omega2.Cross(fITmp))*dAp4;
			dEllePrime = dEllePrime + .5*s*gp.dGetWght()/p.Norm()*p.Dot(pprime);
		} while (gdi->fGetNext(gp));
	}

	/* Check that the distance is not null */
	if (dElle <= std::numeric_limits<doublereal>::epsilon()) {
		silent_cerr("RodBezier(" << GetLabel() << "): "
			"null length of element." << std::endl);
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	/* Strain */
	dEpsilon = dElle/dL0 - 1.;

	/* Strain velocity */
	dEpsilonPrime  = dEllePrime/dL0;

	/* Amplitude of force from CL */
	bool ChangeJac(false);
	try {
		ConstitutiveLaw1DOwner::Update(dEpsilon, dEpsilonPrime);

	} catch (Elem::ChangedEquationStructure& err) {
		ChangeJac = true;
	}

	doublereal dF = GetF();

	/* Force vectors */
	Vec3 F1(l1*dF);
	Vec3 F2(l2*dF);

	WorkVec.Add(1, F1);
	WorkVec.Add(3 + 1, fOTmp.Cross(F1));
	WorkVec.Add(6 + 1, F2);
	WorkVec.Add(9 + 1, fITmp.Cross(F2));

	if (ChangeJac) {
		throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
	}
}

void
RodBezier::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("rod bezier", OH, name);
			
			Var_F2 = OH.CreateVar<Vec3>(name + "F2", "N",
				"force on Node 2 (x, y, z)");

			Var_l = OH.CreateVar<doublereal>(name + "l", "m",
				"length of the element");

			Var_l1 = OH.CreateVar<Vec3>(name + "l1", "-",
				"node 1 reference unit vector (x, y, z)");
			
			Var_l2 = OH.CreateVar<Vec3>(name + "l2", "-",
				"node 2 reference unit vector (x, y, z)");

			Var_v = OH.CreateVar<doublereal>(name + "v", "m/s",
				"length rate of change");
		}
#endif // USE_NETCDF
	}
}

void
RodBezier::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		ASSERT(dElle > std::numeric_limits<doublereal>::epsilon());
		doublereal dF = GetF();
		Vec3 FTmp = Vec3(dF, 0., 0.);

		if (OH.UseText(OutputHandler::JOINTS)) {
			std::ostream& out = OH.Joints();

			Joint::Output(out, "RodBezier", GetLabel(),
					FTmp, Zero3, l1*dF, Zero3)
				<< " " << l2*dF << " " << dElle << " " << l1 << " " 
				<< " " << l2 << " " << " " << dEpsilonPrime*dL0,
				ConstitutiveLaw1DOwner::OutputAppend(out) << std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, FTmp, Zero3, l1*dF, Zero3);
			OH.WriteNcVar(Var_F2, l2*dF);
			OH.WriteNcVar(Var_l, dElle);
			OH.WriteNcVar(Var_l1, l1);
			OH.WriteNcVar(Var_l2, l2);
			OH.WriteNcVar(Var_v, dEpsilonPrime*dL0);
		}
#endif // USE_NETCDF

	}
}

/* Jacobian contribute during initial assembly */
VariableSubMatrixHandler&
RodBezier::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("Entering RodBezier::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode1FirstVelIndex + iCnt);

		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(18 + iCnt, iNode2FirstVelIndex + iCnt);
	}

	const Mat3x3& R1(pNode1->GetRRef());
	const Mat3x3& R2(pNode2->GetRRef());
	const Vec3 fOTmp(R1*fO);
	const Vec3 fATmp(R1*fA);
	const Vec3 fBTmp(R2*fB);
	const Vec3 fITmp(R2*fI);

	const Vec3& x1(pNode1->GetXCurr());
	const Vec3& x2(pNode2->GetXCurr());

	//const Vec3& v1(pNode1->GetVCurr());
	//const Vec3& v2(pNode2->GetVCurr());
	const Vec3& Omega1(pNode1->GetWRef());
	const Vec3& Omega2(pNode2->GetWRef());

	/* Forza e slopes */
	doublereal dF = GetF();
	doublereal dFDE = GetFDE();
	doublereal dFDEPrime = GetFDEPrime();
	
	const Vec3 b1 = x1 + R1*fO; 
	const Vec3 b2 = x1 + R1*fA; 
	const Vec3 b3 = x2 + R2*fB; 
	const Vec3 b4 = x2 + R2*fI; 

	Vec3 p;
	doublereal q, up, pNorm, dKC, dAp1, dAp2, dAp3, dAp4;
	doublereal s = 1.0/(2.0*iIntSeg);

	Vec3 kx1(Zero3);
	Vec3 kg1(Zero3);
	Vec3 kx2(Zero3);
	Vec3 kg2(Zero3);

	Vec3 cg1(Zero3);
	Vec3 cxp1(Zero3);
	Vec3 com1(Zero3);
	Vec3 cg2(Zero3);
	Vec3 cxp2(Zero3);
	Vec3 com2(Zero3);

	Vec3 TmpV;

	for (unsigned int jj = 0; jj < iIntSeg; jj++) {
		q = s*(2*jj + 1);

		PntWght gp = gdi->GetFirst();
		do {
			up = s*gp.dGetPnt() + q;
			
			dAp1 = Ap1(up);
			dAp2 = Ap2(up);
			dAp3 = Ap3(up);
			dAp4 = Ap4(up);

			p = b1*dAp1 + b2*dAp2 + b3*dAp3 + b4*dAp4;
			pNorm = p.Norm();
			dKC = (s*gp.dGetWght())/(2*dL0*pNorm);
			
			kx1 = kx1 + p*dFDE*dKC*(dAp1 + dAp2);
			
			TmpV = fOTmp*dAp1 + fATmp*dAp2;
			kg1 = kg1 - TmpV.Cross(p)*dFDE*dKC;
			
			kx2 = kx2 + p*dFDE*dKC*(dAp3 + dAp4);
			
			TmpV = fBTmp*dAp3 + fITmp*dAp4;
			kg2 = kg2 - TmpV.Cross(p)*dFDE*dKC;

			if (dFDEPrime != 0.) {
				TmpV = Omega1.Cross(p);
				cg1 = cg1 - TmpV.Cross(fATmp*dAp2 - fOTmp*dAp1)*dFDEPrime*dKC;
				
				cxp1 = cxp1 + p*dFDEPrime*dKC*(dAp1 + dAp2);

				com1 = com1 - p.Cross(fATmp*dAp2 - fOTmp*dAp1)*dFDEPrime*dKC;

				TmpV = Omega2.Cross(p);
				cg2 = cg2 + TmpV.Cross(fITmp*dAp4 - fBTmp*dAp3)*dFDEPrime*dKC;
				
				cxp2 = cxp2 + p*dFDEPrime*dKC*(dAp3 + dAp4);

				com2 = com2 - p.Cross(fITmp*dAp4 - fBTmp*dAp3)*dFDEPrime*dKC;

			}

		} while (gdi->fGetNext(gp));
	}

	/* Force on node 1 from dx1 */
	Mat3x3 F1x1 = l1.Tens(kx1);
	WM.Add(1, 1, F1x1);

	/* Force on node 1 from dg1 */
	Vec3 Tmp1 = kg1;
	if (dFDEPrime != 0.) {
		Tmp1 += cg1;
	}
	Mat3x3 F1g1 = l1.Tens(Tmp1) - (Mat3x3(MatCross, l1) - Mat3x3(MatCross, l1)*l1.Tens())*dF;
	WM.Add(1, 3 + 1, F1g1);

	/* Force on node 1 from dxp1 */
	Mat3x3 F1xp1;
	if (dFDEPrime != 0.) {
		F1xp1 = l1.Tens(cxp1);
		WM.Add(1, 6 + 1, F1xp1);
	}

	/* Force on node 1 from dom1 */
	Mat3x3 F1om1;
	if (dFDEPrime != 0.) {
		F1om1 = l1.Tens(com1);
		WM.Add(1, 9 + 1, F1om1);
	}

	/* Force on node 1 from dx2 */
	Mat3x3 F1x2 = l1.Tens(kx2);
	WM.Add(1, 12 + 1, F1x2);

	/* Force on node 1 from dg2 */
	Vec3 Tmp2 = kg2;
	if (dFDEPrime != 0.) {
		Tmp2 += cg2;
	}
	Mat3x3 F1g2 = l1.Tens(Tmp2);
	WM.Add(1, 15 + 1, F1g2);

	/* Force on node 1 from dxp2 */
	Mat3x3 F1xp2;
	if (dFDEPrime != 0.) {
		F1xp2 = l1.Tens(cxp2);
		WM.Add(1, 18 + 1, F1xp2);
	}

	/* Force on node 1 from dom2 */
	Mat3x3 F1om2;
	if (dFDEPrime != 0.) {
		F1om2 = l1.Tens(com2);
		WM.Add(1, 21 + 1, F1om2);
	}

	/* Moment on node 1 from dx1 */
	Mat3x3 M1x1 = Mat3x3(MatCross, fOTmp)*F1x1;
	WM.Sub(3 + 1, 1, M1x1);

	/* Moment on node 1 from dg1 */
	Mat3x3 M1g1 = Mat3x3(MatCross, fOTmp)*F1g1 + Mat3x3(MatCrossCross, l1, fOTmp);
	WM.Sub(3 + 1, 3 + 1, M1g1);

	/* Moment on node 1 from dxp1 */
	if (dFDEPrime != 0.) {
		Mat3x3 M1xp1 = Mat3x3(MatCross, fOTmp)*F1xp1;
		WM.Sub(3 + 1, 6 + 1, M1xp1);
	}

	/* Moment on node 1 from dom1 */
	if (dFDEPrime != 0.) {	
		Mat3x3 M1om1 = Mat3x3(MatCross, fOTmp)*F1om1;
		WM.Sub(3 + 1, 9 + 1, M1om1);
	}

	/* Moment on node 1 from dx2 */
	Mat3x3 M1x2 = Mat3x3(MatCross, fOTmp)*F1x2;
	WM.Sub(3 + 1, 12 + 1, M1x2);

	/* Moment on node 1 from dg2 */
	Mat3x3 M1g2 = Mat3x3(MatCross, fOTmp)*F1g2;
	WM.Sub(3 + 1, 15 + 1, M1g2);

	/* Moment on node 1 from dxp2 */
	if (dFDEPrime != 0.) {
		Mat3x3 M1xp2 = Mat3x3(MatCross, fOTmp)*F1xp2;
		WM.Sub(3 + 1, 18 + 1, M1xp2);
	}

	/* Moment on node 1 from dom2 */
	if (dFDEPrime != 0.) {
		Mat3x3 M1om2 = Mat3x3(MatCross, fOTmp)*F1om2;
		WM.Sub(3 + 1, 21 + 1, M1om2);
	}

	/* Force on node 2 from dx1 */
	Mat3x3 F2x1 = l2.Tens(kx1);
	WM.Add(6 + 1, 1, F1x1);

	/* Force on node 2 from dg1 */
	Tmp1 = kg1;
	Mat3x3 F2g1 = l2.Tens(Tmp1);
	WM.Add(6 + 1, 3 + 1, F2g1);

	/* Force on node 2 from dxp1 */
	Mat3x3 F2xp1;
	if (dFDEPrime != 0.) {
		F2xp1 = l2.Tens(cxp1);
		WM.Add(6 + 1, 6 + 1, F2xp1);
	}

	/* Force on node 2 from dom1 */
	Mat3x3 F2om1;
	if (dFDEPrime != 0.) {
		F2om1 = l2.Tens(com1);
		WM.Add(6 + 1, 9 + 1, F2om1);
	}

	/* Force on node 2 from dx2 */
	Mat3x3 F2x2 = l2.Tens(kx2);
	WM.Add(6 + 1, 12 + 1, F2x2);

	/* Force on node 2 from dg2 */

	Mat3x3 F2g2 = l2.Tens(Tmp2) - (Mat3x3(MatCross, l2) - Mat3x3(MatCross, l2)*l2.Tens())*dF;
	WM.Add(6 + 1, 15 + 1, F2g2);

	/* Force on node 2 from dxp2 */
	Mat3x3 F2xp2;
	if (dFDEPrime != 0.) {
		F2xp2 = l2.Tens(cxp2);
		WM.Add(6 + 1, 18 + 1, F2xp2);
	}

	/* Force on node 2 from dom2 */
	Mat3x3 F2om2;
	if (dFDEPrime != 0.) {
		F2om2 = l2.Tens(com2);
		WM.Add(6 + 1, 21 + 1, F2om2);
	}

	/* Moment on node 2 from dx2 */
	Mat3x3 M2x1 = Mat3x3(MatCross, fITmp)*F2x1;
	WM.Sub(9 + 1, 1, M2x1);

	/* Moment on node 2 from dg1 */
	Mat3x3 M2g1 = Mat3x3(MatCross, fITmp)*F2g1;
	WM.Sub(9 + 1, 3 + 1, M1g1);

	/* Moment on node 2 from dxp1 */
	if (dFDEPrime != 0.) {
		Mat3x3 M2xp1 = Mat3x3(MatCross, fITmp)*F2xp1;
		WM.Sub(9 + 1, 6 + 1, M2xp1);
	}

	/* Moment on node 2 from dom1 */
	if (dFDEPrime != 0.) {	
		Mat3x3 M2om1 = Mat3x3(MatCross, fITmp)*F2om1;
		WM.Sub(9 + 1, 9 + 1, M2om1);
	}

	/* Moment on node 2 from dx2 */
	Mat3x3 M2x2 = Mat3x3(MatCross, fITmp)*F2x2;
	WM.Sub(9 + 1, 12 + 1, M2x2);

	/* Moment on node 2 from dg2 */
	Mat3x3 M2g2 = Mat3x3(MatCross, fITmp)*F2g2 + Mat3x3(MatCrossCross, l2, fITmp);
	WM.Sub(9 + 1, 15 + 1, M2g2);

	/* Moment on node 2 from dxp2 */
	if (dFDEPrime != 0.) {
		Mat3x3 M2xp2 = Mat3x3(MatCross, fITmp)*F2xp2;
		WM.Sub(9 + 1, 18 + 1, M2xp2);
	}

	/* Moment on node 2 from dom2 */
	if (dFDEPrime != 0.) {
		Mat3x3 M2om2 = Mat3x3(MatCross, fITmp)*F2om2;
		WM.Sub(9 + 1, 21 + 1, M2om2);
	}

	return WorkMat;
}

/* Residual contribute during initial assembly */
SubVectorHandler&
RodBezier::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUT("RodBezier::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Inverse dynamics capable element */
bool
RodBezier::bInverseDynamics(void) const
{
	return true;
}

/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
RodBezier::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	ASSERT(bIsErgonomy());
	return AssJac(WorkMat, 1., XCurr, XCurr);
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
RodBezier::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr, 
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	DEBUGCOUT("Entering RodBezier::AssRes()" << std::endl);

	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS
		|| (iOrder == InverseDynamics::POSITION && bIsErgonomy()));
	
	return AssRes(WorkVec, 1., XCurr, XPrimeCurr);
}

/* Inverse Dyamics update */
void
RodBezier::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	NO_OP;
}

/* Inverse Dynamics after convergence */
void
RodBezier::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
}

unsigned int
RodBezier::iGetNumPrivData(void) const
{
	return 3 + ConstitutiveLaw1DOwner::iGetNumPrivData();
}

unsigned int
RodBezier::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "F") == 0) {
		return 1;
	}

	if (strcmp(s, "L") == 0) {
		return 2;
	}

	if (strcmp(s, "LPrime") == 0) {
		return 3;
	}

	size_t l = STRLENOF("constitutiveLaw.");
	if (strncmp(s, "constitutiveLaw.", l) == 0) {
		return 3 + ConstitutiveLaw1DOwner::iGetPrivDataIdx(&s[l]);
	}

	/* error; handle later */
	return 0;
}

doublereal
RodBezier::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	switch (i) {
	case 1:
		return GetF();

	case 2:
		return dElle;

	case 3:
		return dL0*dEpsilonPrime;
	}

	i -= 3;
	ASSERT(i <= ConstitutiveLaw1DOwner::iGetNumPrivData());

	return ConstitutiveLaw1DOwner::dGetPrivData(i);
}

/* RodBezier - end */

