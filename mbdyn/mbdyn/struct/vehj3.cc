/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

/* Cerniera deformabile */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "dataman.h"
#include "vehj3.h"

#include "matvecexp.h"
#include "Rot.hh"

/* DeformableJoint - begin */

/* Costruttore non banale */
DeformableJoint::DeformableJoint(unsigned int uL,
		const DofOwner* pDO, 
		const ConstitutiveLaw6D* pCL,
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& f1Tmp,
		const Vec3& f2Tmp,
		const Mat3x3& R1,
		const Mat3x3& R2,
		flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::DEFORMABLEJOINT, pDO, fOut), 
ConstitutiveLaw6DOwner(pCL),
pNode1(pN1), pNode2(pN2), f1(f1Tmp), f2(f2Tmp), R1h(R1), R2h(R2)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}

   
/* Distruttore */
DeformableJoint::~DeformableJoint(void)
{
	NO_OP;
}

   
/* Contributo al file di restart */
std::ostream&
DeformableJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", deformable joint, " 
     << pNode1->GetLabel() << ", reference, node, ",
     f1.Write(out, ", ") << ", hinge, reference, node, 1, ", 
     (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, ",
     f2.Write(out, ", ") << ", hinge, reference, node, 1, ",
     (R2h.GetVec(1)).Write(out, ", ")
       << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ", ";
   return pGetConstLaw()->Restart(out) << ';' << std::endl;
}


void
DeformableJoint::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		Mat3x3 R1(pNode1->GetRCurr());
      		Vec3 F(GetF().GetVec1());
		Vec3 M(GetF().GetVec2());
		Joint::Output(OH.Joints(), "DeformableJoint", GetLabel(),
				F, M, R1*F, R1*(R1h*M)) << std::endl;	
   }     
}

unsigned int
DeformableJoint::iGetNumPrivData(void) const
{
	return 18 + ConstitutiveLaw6DOwner::iGetNumPrivData();
}

unsigned int
DeformableJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	
	unsigned idx = 0;
	
	switch (s[0]) {
	case 'd':
		break;

	case 'r':
		idx += 3;
		break;

	case 'v':
		idx += 6;
		break;

	case 'w':
		idx += 9;
		break;

	case 'F':
		idx += 12;
		break;

	case 'M':
		idx += 15;
		break;

	default:
	{
		size_t l = sizeof("constitutiveLaw.") - 1;
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			return 18 + ConstitutiveLaw6DOwner::iGetPrivDataIdx(s + l);
		}
		return 0;
	}
	}
	
	switch (s[1]) {
	case 'x':
		idx += 1;
		break;
	case 'y':
		idx += 2;
		break;
	case 'z':
		idx += 3;
		break;
	default:
		return 0;
	}
	
	if (s[2] != '\0') {
		return 0;
	}
	
	return idx;
}

doublereal
DeformableJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());
	
	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		Vec3 f1Tmp(pNode1->GetRCurr()*f1);
		Vec3 f2Tmp(pNode2->GetRCurr()*f2);
		Mat3x3 R1T((pNode1->GetRCurr()*R1h).Transpose());
		Vec3 d(R1T*(pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr() - f1Tmp));

		return d(i);
	}
	
	case 4:
	case 5:
	case 6:
	{
		Mat3x3 R1T((pNode1->GetRCurr()*R1h).Transpose());
		Mat3x3 R2(pNode2->GetRCurr()*R2h);

		Vec3 v(RotManip::VecRot(R1T*R2));

		return v(i - 3);
	}
	
	case 7:
	case 8:
	case 9:
	{
		Vec3 f2Tmp(pNode2->GetRCurr()*f2);
		Mat3x3 R1T(pNode1->GetRCurr().Transpose());
		Vec3 v(R1T*(pNode2->GetVCurr() - pNode1->GetVCurr()
					+ (pNode2->GetXCurr() - pNode1->GetXCurr()).Cross(pNode1->GetWCurr())
					- f2Tmp.Cross(pNode2->GetWCurr() - pNode1->GetWCurr())));

		return v(i - 6);
	}
	
	case 10:
	case 11:
	case 12:
	{
		Mat3x3 R1T((pNode1->GetRCurr()*R1h).Transpose());
		Vec3 W1(pNode1->GetWCurr());
		Vec3 W2(pNode2->GetWCurr());
		Vec3 w = R1T*(W2 - W1);

		return w(i - 9);
	}

	case 13:
	case 14:
	case 15:
	case 16:
	case 17:
	case 18:
		return GetF()(i - 12);
	
	default:
		return ConstitutiveLaw6DOwner::dGetPrivData(i - 18);
	}
}

/* DeformableJoint - end */


/* ElasticJoint - begin */

ElasticJoint::ElasticJoint(unsigned int uL, 
		const DofOwner* pDO, 
		const ConstitutiveLaw6D* pCL,
		const StructNode* pN1, 
		const StructNode* pN2,
		const Vec3& f1Tmp,
		const Vec3& f2Tmp,
		const Mat3x3& R1,
		const Mat3x3& R2,
		flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableJoint(uL, pDO, pCL, pN1, pN2, f1Tmp, f2Tmp, R1, R2, fOut),
ThetaRef(0.), k(0.)
{
	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	Mat3x3 R(pNode1->GetRRef()*R1h);
	FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE(), R);
}


ElasticJoint::~ElasticJoint(void)
{
	NO_OP;
}

void
ElasticJoint::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(k);
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ElasticJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();
   
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex+iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssMat(WM, dCoef);

	return WorkMat;
}


void
ElasticJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)   
{      
	Vec3 f1Tmp(pNode1->GetRRef()*f1);
	Vec3 f2Tmp(pNode2->GetRRef()*f2);
	Vec3 d1(pNode2->GetXCurr()+f2Tmp-pNode1->GetXCurr());

	/* D11 */
	Mat3x3 DTmp = FDE.GetMat11()*dCoef;

	WM.Add(1, 1, DTmp);
	WM.Sub(1, 6 + 1, DTmp);
	WM.Sub(6 + 1, 1, DTmp);
	WM.Add(6 + 1, 6 + 1, DTmp);

	/* D11 * (f1 x) */
	Mat3x3 FTmp = DTmp*Mat3x3(f1Tmp);

	WM.Sub(1, 4, FTmp);
	WM.Add(6 + 1, 4, FTmp);
	
	/* D11 * (f2 x) */
	FTmp = DTmp*Mat3x3(f2Tmp);

	WM.Add(1, 6 + 4, FTmp);
	WM.Sub(6 + 1, 6 + 4, FTmp);

	/* (f1 x) * D11 */
	FTmp = Mat3x3(f1Tmp)*DTmp;

	WM.Add(4, 1, FTmp);
	WM.Sub(4, 6 + 1, FTmp);
	
	/* (f2 x) * D11 */
	FTmp = Mat3x3(f2Tmp)*DTmp;

	WM.Add(6 + 4, 1, FTmp);
	WM.Sub(6 + 4, 6 + 1, FTmp);
	

	/* D12 - (F1 x) */
	DTmp = (FDE.GetMat12() - Mat3x3(F.GetVec1()))*dCoef;

	WM.Add(1, 4, DTmp);
	WM.Sub(1, 6 + 4, DTmp);
	WM.Sub(6 + 1, 4, DTmp);
	WM.Add(6 + 1, 6 + 4, DTmp);

	/* (f1 x) * D11 */
	FTmp = Mat3x3(f1Tmp)*DTmp;

	WM.Add(4, 4, FTmp);
	WM.Sub(4, 6 + 4, FTmp);
	
	/* (f2 x) * D11 */
	FTmp = Mat3x3(f2Tmp)*DTmp;

	WM.Add(6 + 4, 4, FTmp);
	WM.Sub(6 + 4, 6 + 4, FTmp);


	/* D21 */
	DTmp = FDE.GetMat21()*dCoef;

	WM.Add(4, 1, DTmp);
	WM.Sub(4, 6 + 1, DTmp);
	WM.Sub(6 + 4, 1, DTmp);
	WM.Add(6 + 4, 6 + 1, DTmp);

	/* D21 * (f1 x) */
	FTmp = DTmp*Mat3x3(f1Tmp);

	WM.Sub(4, 4, FTmp);
	WM.Add(6 + 4, 4, FTmp);
	
	/* D21 * (f2 x) */
	FTmp = DTmp*Mat3x3(f2Tmp);

	WM.Add(4, 6 + 4, FTmp);
	WM.Sub(6 + 4, 6 + 4, FTmp);


	/* D22 - (F2 x) */
	DTmp = (FDE.GetMat22() - Mat3x3(F.GetVec2()))*dCoef;

	WM.Add(4, 4, DTmp);
	WM.Sub(4, 6 + 4, DTmp);
	WM.Sub(6 + 4, 4, DTmp);
	WM.Add(6 + 4, 6 + 4, DTmp);


	/* ((F1 x f1) x) */
	FTmp = Mat3x3(F.GetVec1().Cross(f1Tmp*dCoef));

	WM.Add(4, 4, FTmp);
	WM.Sub(6 + 4, 4, FTmp);

	/* ((F1 x f2) x) */
	FTmp = Mat3x3(F.GetVec1().Cross(f2Tmp*dCoef));

	WM.Sub(4, 6 + 4, FTmp);
	WM.Add(6 + 4, 6 + 4, FTmp);
}

/* assemblaggio residuo */
SubVectorHandler& 
ElasticJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex+iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

   
void
ElasticJoint::AssVec(SubVectorHandler& WorkVec)
{   
	Mat3x3 R1(pNode1->GetRCurr()*R1h);
	Vec3 f1Tmp(pNode1->GetRCurr()*f1);
	Vec3 f2Tmp(pNode2->GetRCurr()*f2);

	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Mat3x3 R1T(R1.Transpose());
		Mat3x3 R2(pNode2->GetRCurr()*R2h);
		Vec3 d1(pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr());   

		ThetaCurr = RotManip::VecRot(R1T*R2);

		k = Vec6(R1T*(d1 - f1Tmp), ThetaCurr);
		
		ConstitutiveLaw6DOwner::Update(k);
	}

	F = MultRV(GetF(), R1);
   
	WorkVec.Add(1, F.GetVec1());
	WorkVec.Add(4, f1Tmp.Cross(F.GetVec1()) + F.GetVec2());
	WorkVec.Sub(6 + 1, F.GetVec1());   
	WorkVec.Sub(6 + 4, f2Tmp.Cross(F.GetVec1()) + F.GetVec2());
}

void
ElasticJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Mat3x3 R1(pNode1->GetRRef()*R1h);
	Mat3x3 R2(pNode2->GetRRef()*R2h);
	Mat3x3 R1T(R1.Transpose());

	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
	ThetaCurr = ThetaRef = RotManip::VecRot(R1T*R2);

	/* Calcola l'inversa di Gamma di ThetaRef */
	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);

	Vec3 f1Tmp(pNode1->GetRRef()*f1);
	Vec3 f2Tmp(pNode2->GetRRef()*f2);
	Vec3 d1(pNode2->GetXCurr() + f2Tmp - pNode1->GetXCurr());   

	k = Vec6(R1T*(d1 - f1Tmp), ThetaRef);

	/* Aggiorna il legame costitutivo */
	ConstitutiveLaw6DOwner::Update(k);

	/* Chiede la matrice tangente di riferimento e la porta
	 * nel sistema globale */
	/* FIXME: horrible */
        FDE = MultRMRt(ConstitutiveLaw6DOwner::GetFDE()*Mat6x6(Eye3, Zero3x3, Zero3x3, GammaRefm1), R1);

	bFirstRes = true;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ElasticJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssMat(WM, 1.);

	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ElasticJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex+iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* ElasticJoint - end */

#if 0

/* ViscousJoint - begin */

ViscousJoint::ViscousJoint(unsigned int uL, 
		const DofOwner* pDO, 
		const ConstitutiveLaw6D* pCL,
		const StructNode* pN1, 
		const StructNode* pN2,
		const Vec3& f1Tmp,
		const Vec3& f2Tmp,
		const Mat3x3& R1,
		const Mat3x3& R2,
		flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableJoint(uL, pDO, pCL, pN1, pN2, f1Tmp, f2Tmp, R1, R2, fOut)
{
   NO_OP;
   
   /* Temporary */
   silent_cerr("DeformableHingeJoint(" << GetLabel() << "): "
	   "warning, this element is not implemented yet" << std::endl);
   throw ErrNotImplementedYet();
}


ViscousJoint::~ViscousJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscousJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1.Transpose());
   Mat3x3 Tmp(FDEPrime-FDEPrime*Mat3x3(Omega2*dCoef));
   
   WM.Add(4, 4, Tmp);
   WM.Sub(1, 4, Tmp);
   
   Tmp += Mat3x3(F*dCoef);
   WM.Add(1, 1, Tmp);   
   WM.Sub(4, 1, Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscousJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }  
 
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   kPrime = R1.Transpose()*(Omega2-Omega1);
   IncrementalUpdate(Zero3, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);   
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscousJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }  

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1.Transpose());      
   Mat3x3 Tmp(Mat3x3(F)-FDEPrime*Mat3x3(Omega2-Omega1));
   
   WM.Add(1, 1, Tmp);
   WM.Sub(4, 1, Tmp);
   
   WM.Add(1, 4, FDEPrime);
   WM.Add(4, 10, FDEPrime);
   
   WM.Sub(1, 10, FDEPrime);   
   WM.Sub(4, 4, FDEPrime);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscousJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   k = R1.Transpose()*(g2*(4./(4.+g2.Dot()))
			  -g1*(4./(4.+g1.Dot())));
   kPrime = R1.Transpose()*(Omega2-Omega1);
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);   
   
   return WorkVec;
}

/* ViscousJoint - end */


/* ViscoElasticJoint - begin */

ViscoElasticJoint::ViscoElasticJoint(unsigned int uL, 
		const DofOwner* pDO, 
		const ConstitutiveLaw6D* pCL,
		const StructNode* pN1, 
		const StructNode* pN2,
		const Vec3& f1Tmp,
		const Vec3& f2Tmp,
		const Mat3x3& R1,
		const Mat3x3& R2, 
		flag fOut)
: Elem(uL, Elem::JOINT, fOut), 
DeformableJoint(uL, pDO, pCL, pN1, pN2, f1Tmp, f2Tmp, R1, R2, fOut)
{
   /* Temporary */
   silent_cerr("DeformableHingeJoint(" << GetLabel() << "): "
	   "warning, this element is not implemented yet" << std::endl);
   throw ErrNotImplementedYet();
}


ViscoElasticJoint::~ViscoElasticJoint(void)
{
   NO_OP;
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
ViscoElasticJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);	
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  

   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDE(R1*GetFDE()*(R1t*dCoef));
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1t);
   
   Mat3x3 Tmp(FDEPrime-FDEPrime*Mat3x3(Omega2*dCoef)+FDE);
         
   WM.Add(4, 4, Tmp);
   WM.Sub(1, 4, Tmp);
   
   Tmp += Mat3x3(F*dCoef);
   WM.Add(1, 1, Tmp);   
   WM.Sub(4, 1, Tmp);

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& 
ViscoElasticJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   /* Recupera gli indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
   }  
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   k = R1t*(g2*(4./(4.+g2.Dot())) -g1*(4./(4.+g1.Dot())));
   kPrime = R1t*(Omega2-Omega1);
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);   
   
   return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ViscoElasticJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& /* XCurr */ )
{
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      
      WM.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }  
   
   Mat3x3 R1(pNode1->GetRRef()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());

   Vec3 F(R1*GetF());
   Mat3x3 FDE(R1*GetFDE()*R1t);
   Mat3x3 FDEPrime(R1*GetFDEPrime()*R1t);
   
   Mat3x3 Tmp(Mat3x3(F)-FDEPrime*Mat3x3(Omega2-Omega1)+FDE);
   
   WM.Add(1, 1, Tmp);
   WM.Sub(4, 1, Tmp);
   
   WM.Add(4, 7, FDE);
   WM.Sub(1, 7, FDE);
   
   WM.Add(1, 4, FDEPrime);
   WM.Add(4, 10, FDEPrime);
   
   WM.Sub(1, 10, FDEPrime);   
   WM.Sub(4, 4, FDEPrime);

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ViscoElasticJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Recupera gli indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   
   /* Setta gli indici della matrice */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
   }  
   
   Mat3x3 R1(pNode1->GetRCurr()*R1h);
   Mat3x3 R1t(R1.Transpose());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   Vec3 g1(pNode1->GetgCurr());
   Vec3 g2(pNode2->GetgCurr());
   
   /* Aggiornamento: k += R1^T(G(g2)*g2-G(g1)*g1) */
   k = R1t*(g2*(4./(4.+g2.Dot())) - g1*(4./(4.+g1.Dot())));
   kPrime = R1t*(Omega2-Omega1);
   IncrementalUpdate(k, kPrime);
   
   Vec3 F(R1*GetF());
   
   WorkVec.Add(1, F);
   WorkVec.Sub(4, F);   
   
   return WorkVec;
}

/* ViscoElasticJoint - end */
#endif
