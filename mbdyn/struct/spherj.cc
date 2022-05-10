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

/* Giunti sferici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "spherj.h"
#include "Rot.hh"


/* SphericalHingeJoint - begin */

/* Costruttore non banale */
SphericalHingeJoint::SphericalHingeJoint(unsigned int uL, const DofOwner* pDO,
					 const StructNode* pN1, 
					 const StructNode* pN2,
					 const Vec3& dTmp1, const Mat3x3& RTmp1h,
					 const Vec3& dTmp2, const Mat3x3& RTmp2h,
					 const OrientationDescription& od,
					 flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut),
pNode1(pN1), pNode2(pN2), 
d1(dTmp1), R1h(RTmp1h),
d2(dTmp2), R2h(RTmp2h), 
F(Zero3),
od(od)
{
   NO_OP;
}


/* Distruttore banale */
SphericalHingeJoint::~SphericalHingeJoint(void)
{
   NO_OP;
};


/* Contributo al file di restart */
std::ostream& SphericalHingeJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", spherical hinge, "
     << pNode1->GetLabel() << ", reference, node, ",
     d1.Write(out, ", ")  << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "       
     << pNode2->GetLabel() << ", reference, node, ",
     d2.Write(out, ", ") << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
SphericalHingeJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering SphericalHingeJoint::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(54, 1);
   
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   Vec3 dTmp1(pNode1->GetRRef()*d1);
   Vec3 dTmp2(pNode2->GetRRef()*d2);
   
   
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere consistente con la posizione dei due nodi:
    *      x2 + d2 = x1 + d1
    * 
    * con: d2 = R2 * d2_0
    *      d1 = R1 * d1_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     -d1 /\ F    per il nodo 1,
    *      d2 /\ F    per il nodo 2
    *
    * 
    *         x1   g1        x2     g2        F
    * Q1 |  0      0         0      0         I    | | x1 |   | -F           |
    * G1 |  0      cF/\d1/\  0      0         d1/\ | | g1 |   | -d1/\F       |
    * Q2 |  0      0         0     -cF/\d2/\ -I    | | x2 | = |  F           |
    * G2 |  0      0         0      0        -d2/\ | | g2 |   |  d2/\F       |
    * F  | -c*I    c*d1/\    c*I   -c*d2/\    0    | | F  |   |  x1+d1-x2-d2 |
    * 
    * con d1 = R1*d01, d2 = R2*d02, c = dCoef
    */

   /* Moltiplico la forza per il coefficiente del metodo.
    * Nota: F, la reazione vincolare, e' stata aggiornata da AssRes */
   Vec3 FTmp = F*dCoef;
   
   /* termini di reazione sul nodo 1 */
   WM.PutDiag(1, iNode1FirstMomIndex, iFirstReactionIndex, 1.);
   WM.PutCross(4, iNode1FirstMomIndex+3, iFirstReactionIndex, dTmp1);
   
   WM.PutMat3x3(10, iNode1FirstMomIndex+3,
		 iNode1FirstPosIndex+3, Mat3x3(MatCrossCross, FTmp, dTmp1));

   /* termini di reazione sul nodo 2 */
   WM.PutDiag(19, iNode2FirstMomIndex, iFirstReactionIndex, -1.);
   WM.PutCross(22, iNode2FirstMomIndex+3, iFirstReactionIndex, -dTmp2);

   WM.PutMat3x3(28, iNode2FirstMomIndex+3,
		 iNode2FirstPosIndex+3, Mat3x3(MatCrossCross, FTmp, -dTmp2));
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   WM.PutDiag(37, iFirstReactionIndex, iNode1FirstPosIndex, -1.);
   WM.PutCross(40, iFirstReactionIndex, iNode1FirstPosIndex+3, dTmp1);
      
   /* termini di vincolo dovuti al nodo 1 */
   WM.PutDiag(46, iFirstReactionIndex, iNode2FirstPosIndex, 1.);
   WM.PutCross(49, iFirstReactionIndex, iNode2FirstPosIndex+3, -dTmp2);

   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& SphericalHingeJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr, 
					      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering SphericalHingeJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);   
   }
   
   F = Vec3(XCurr, iFirstReactionIndex+1);
   
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   
   Vec3 dTmp1(pNode1->GetRCurr()*d1);
   Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
   WorkVec.Sub(1, F);
   WorkVec.Sub(4, dTmp1.Cross(F));
   WorkVec.Add(7, F);
   WorkVec.Add(10, dTmp2.Cross(F));
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
   WorkVec.Add(13, (x1+dTmp1-x2-dTmp2)/dCoef);

   return WorkVec;
}

			    
DofOrder::Order SphericalHingeJoint::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i >=0 and i < iGetNumDof(), 
		"INDEX ERROR in SphericalHingeJoint::GetEqType");
	return DofOrder::ALGEBRAIC;
}

void
SphericalHingeJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Spherical hinge", OH, name);

			Var_Phi = OH.CreateRotationVar(name, "", od, 
					"relative orientation, in joint reference frame");

		}
#endif // USE_NETCDF
	}
}

/* Output (da mettere a punto) */
void SphericalHingeJoint::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) {
		Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);
		Mat3x3 RTmp(R1Tmp.MulTM(pNode2->GetRCurr()*R2h));
		Vec3 E;
		switch (od) {
		case EULER_123:
			E = MatR2EulerAngles123(RTmp)*dRaDegr;
			break;

		case EULER_313:
			E = MatR2EulerAngles313(RTmp)*dRaDegr;
			break;

		case EULER_321:
			E = MatR2EulerAngles321(RTmp)*dRaDegr;
			break;

		case ORIENTATION_VECTOR:
			E = RotManip::VecRot(RTmp);
			break;

		case ORIENTATION_MATRIX:
			break;

		default:
			/* impossible */
			break;
		}
      
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, (R1Tmp.MulTV(F)), Zero3, F, Zero3);
			switch (od) {
			case EULER_123:
			case EULER_313:
			case EULER_321:
			case ORIENTATION_VECTOR:
				OH.WriteNcVar(Var_Phi, E);
				break;

			case ORIENTATION_MATRIX:
				OH.WriteNcVar(Var_Phi, RTmp);
				break;

			default:
				/* impossible */
				break;
			}
		}
#endif // USE_NETCDF
		if (OH.UseText(OutputHandler::JOINTS)) {
				Joint::Output(OH.Joints(), "SphericalHinge", GetLabel(),
				R1Tmp.MulTV(F), Zero3, F, Zero3)
				<< " ";

				switch (od) {
				case EULER_123:
				case EULER_313:
				case EULER_321:
				case ORIENTATION_VECTOR:
					OH.Joints() << E;
					break;

				case ORIENTATION_MATRIX:
					OH.Joints() << RTmp;
					break;

				default:
					/* impossible */
					break;
				}

				OH.Joints() << std::endl;
		}
   }   
}

void
SphericalHingeJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::OffsetHint<1> *>(pjh)) {
				const Mat3x3& R1(pNode1->GetRCurr());
				Vec3 dTmp2(pNode2->GetRCurr()*d2);
   
				d1 = R1.MulTV(pNode2->GetXCurr() + dTmp2 - pNode1->GetXCurr());

			} else if (dynamic_cast<Joint::OffsetHint<2> *>(pjh)) {
				const Mat3x3& R2(pNode2->GetRCurr());
				Vec3 dTmp1(pNode1->GetRCurr()*d1);
   
				d2 = R2.MulTV(pNode1->GetXCurr() + dTmp1 - pNode2->GetXCurr());

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}
}

Hint *
SphericalHingeJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "offset{" /*}*/, STRLENOF("offset{" /*}*/)) == 0) {
		s += STRLENOF("offset{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::OffsetHint<1>;

		case '2':
			return new Joint::OffsetHint<2>;
		}
	}

	return 0;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
SphericalHingeJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering SphericalHingeJoint::InitialAssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza, nodo 1 */
      WM.PutCoef(iCnt, 24+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 1 */
      WM.PutCoef(6+iCnt, 27+iCnt, 1.);
      
      /* Contributo di forza all'equazione della forza, nodo 2 */
      WM.PutCoef(12+iCnt, 24+iCnt, -1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza, nodo 2 */
      WM.PutCoef(18+iCnt, 27+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 1 */
      WM.PutCoef(24+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo, nodo 1 */
      WM.PutCoef(27+iCnt, 6+iCnt, -1.);
      
      /* Equazione di vincolo, nodo 2 */
      WM.PutCoef(24+iCnt, 12+iCnt, 1.);
      
      /* Derivata dell'equazione di vincolo, nodo 2 */
      WM.PutCoef(27+iCnt, 18+iCnt, 1.);
   }
   
   /* Recupera i dati */
   const Mat3x3& R1(pNode1->GetRRef());
   const Mat3x3& R2(pNode2->GetRRef());
   const Vec3& Omega1(pNode1->GetWRef());
   const Vec3& Omega2(pNode2->GetWRef());
   /* F e' stata aggiornata da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);

   /* Matrici F/\d1/\, -F/\d2/\ */
   Mat3x3 FWedged1Wedge(MatCrossCross, F, d1Tmp);
   Mat3x3 FWedged2Wedge(MatCrossCross, F, -d2Tmp);
   
   /* Matrici (omega1/\d1)/\, -(omega2/\d2)/\ */
   Mat3x3 O1Wedged1Wedge(MatCross, Omega1.Cross(d1Tmp));
   Mat3x3 O2Wedged2Wedge(MatCross, d2Tmp.Cross(Omega2));
   
   /* Equazione di momento, nodo 1 */
   WM.Add(4, 4, FWedged1Wedge);
   WM.Add(4, 25, Mat3x3(MatCross, d1Tmp));
   
   /* Equazione di momento, nodo 2 */
   WM.Add(16, 16, FWedged2Wedge);
   WM.Sub(16, 25, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di momento, nodo 1 */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega1))*Mat3x3(MatCross, d1Tmp));
   WM.Add(10, 10, FWedged1Wedge);
   WM.Add(10, 25, O1Wedged1Wedge);
   WM.Add(10, 28, Mat3x3(MatCross, d1Tmp));
   
   /* Derivata dell'equazione di momento, nodo 2 */
   WM.Sub(22, 16, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega2))*Mat3x3(MatCross, d2Tmp));
   WM.Add(22, 22, FWedged2Wedge);
   WM.Add(22, 25, O2Wedged2Wedge);
   WM.Sub(22, 28, Mat3x3(MatCross, d2Tmp));
      
   /* Equazione di vincolo */
   WM.Add(25, 4, Mat3x3(MatCross, d1Tmp));
   WM.Sub(25, 16, Mat3x3(MatCross, d2Tmp));
   
   /* Derivata dell'equazione di vincolo */
   WM.Add(28, 4, O1Wedged1Wedge);
   WM.Add(28, 10, Mat3x3(MatCross, d1Tmp));
   WM.Add(28, 16, O2Wedged2Wedge);
   WM.Sub(28, 22, Mat3x3(MatCross, d2Tmp));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
SphericalHingeJoint::InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering SphericalHingeJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
   
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iNode2FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(24+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   const Vec3& x1(pNode1->GetXCurr());
   const Vec3& x2(pNode2->GetXCurr());
   const Vec3& v1(pNode1->GetVCurr());
   const Vec3& v2(pNode2->GetVCurr());
   const Mat3x3& R1(pNode1->GetRCurr());
   const Mat3x3& R2(pNode2->GetRCurr());
   const Vec3& Omega1(pNode1->GetWCurr());
   const Vec3& Omega2(pNode2->GetWCurr());
   F = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);

   /* Vettori omega1/\d1, -omega2/\d2 */
   Vec3 O1Wedged1(Omega1.Cross(d1Tmp));
   Vec3 O2Wedged2(Omega2.Cross(d2Tmp));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(d1Tmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(d1Tmp)-O1Wedged1.Cross(F));
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(13, F);
   WorkVec.Add(16, d2Tmp.Cross(F)); 
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(19, FPrime);
   WorkVec.Add(22, d2Tmp.Cross(FPrime)+O2Wedged2.Cross(F));
   
   /* Equazione di vincolo */
   WorkVec.Add(25, x1+d1Tmp-x2-d2Tmp);
   
   /* Deivata dell'equazione di vincolo */
   WorkVec.Add(28, v1+O1Wedged1-v2-O2Wedged2);
      
   return WorkVec;
}

const MBUnits::Dimensions
SphericalHingeJoint::GetEquationDimension(integer index) const {
	// DOF == 3
	MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = MBUnits::Dimensions::Length;
			break;
		case 2:
			dimension = MBUnits::Dimensions::Length;
			break;
      case 3:
			dimension = MBUnits::Dimensions::Length;
			break;
	}

	return dimension;
}

std::ostream&
SphericalHingeJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"relative position constraints" << std::endl;

	return out;
}

/* SphericalHingeJoint - end */


/* PinJoint - begin */
/* Costruttore non banale */
PinJoint::PinJoint(unsigned int uL, const DofOwner* pDO,	       
		   const StructNode* pN,
		   const Vec3& X0Tmp, const Vec3& dTmp, flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), pNode(pN), X0(X0Tmp), d(dTmp), 
F(Zero3)
{
   NO_OP;
}


/* Distruttore banale */
PinJoint::~PinJoint(void)
{
   NO_OP;
};


/* Contributo al file di restart */
std::ostream& PinJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", pin, "
     << pNode->GetLabel() << ", reference, node, ",
     d.Write(out, ", ") << ", reference, global, ",
     X0.Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PinJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		 doublereal dCoef,
		 const VectorHandler& /* XCurr */ ,
		 const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PinJoint::AssJac()" << std::endl);
      
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(27, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   const Mat3x3& R(pNode->GetRRef());
   Vec3 dTmp(R*d);
      
   /* 
    * L'equazione di vincolo afferma che il punto in cui si trova la
    * cerniera deve essere fissato:
    *      x + d = x0
    * 
    * con: d = R * d_0
    * 
    * La forza e' data dalla reazione vincolare F, nel sistema globale
    * La coppia dovuta all'eccentricita' e' data rispettivamente da:
    *     d /\ F
    *
    * 
    *       x      g         F
    * Q1 |  0      0         I   | | x |   | -F          |
    * G1 |  0      cF/\d1/\  d/\ | | g |   | -d/\F       |
    * F  |  I      d/\       0   | | F |   |  (x+d-x0)/c |
    * 
    * con d = R*d_0, c = dCoef
    */

   /* termini di reazione sul nodo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(iCnt, iFirstMomentumIndex+iCnt, 
		  iFirstReactionIndex+iCnt, 1.);
   }   
   WM.PutCross(4, iFirstMomentumIndex+3,
		iFirstReactionIndex, dTmp);
      
   /* Nota: F, la reazione vincolare, e' stata aggiornata da AssRes */
   
   /* Termini diagonali del tipo: c*F/\d/\Delta_g 
    * nota: la forza e' gia' moltiplicata per dCoef */      
   WM.PutMat3x3(10, iFirstMomentumIndex+3,
		 iFirstPositionIndex+3, Mat3x3(MatCrossCross, F*dCoef, dTmp));

   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* termini di vincolo dovuti al nodo 1 */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutItem(18+iCnt, iFirstReactionIndex+iCnt, 
		  iFirstPositionIndex+iCnt, -1.);
   }
   WM.PutCross(22, iFirstReactionIndex,
		iFirstPositionIndex+3, dTmp);
         
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& PinJoint::AssRes(SubVectorHandler& WorkVec,
					      doublereal dCoef,
					      const VectorHandler& XCurr,
					      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PinJoint::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
     
   
   /* Indici del vincolo */
   for(int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }

   F = Vec3(XCurr, iFirstReactionIndex+1);
   
   const Vec3& x(pNode->GetXCurr());
   const Mat3x3& R(pNode->GetRCurr());
   
   Vec3 dTmp(R*d);
   
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);
   WorkVec.Add(7, (x+dTmp-X0)/dCoef);

   return WorkVec;
}

DofOrder::Order PinJoint::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i >=0 and i < iGetNumDof(), 
		"INDEX ERROR in PinJoint::GetEqType");
	return DofOrder::ALGEBRAIC;
}


void
PinJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Spherical pin", OH, name);

			Var_Phi = OH.CreateVar<Vec3>(name + "E",
				MBUnits::Dimensions::deg,
				"node orientation (E123)");

		}
#endif // USE_NETCDF
	}
}


/* Output (da mettere a punto) */
void PinJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "Pin", GetLabel(), F, Zero3, F, Zero3) 
				<< " " << MatR2EulerAngles(pNode->GetRCurr())*dRaDegr << std::endl;
			OH.Joints() << std::endl;
		}
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, Zero3, F, Zero3);
			OH.WriteNcVar(Var_Phi, MatR2EulerAngles(pNode->GetRCurr())*dRaDegr);
		}
#endif // USE_NETCDF
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PinJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering PinJoint::InitialAssJac()" << std::endl);

   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);

   /* Equazioni: vedi joints.dvi */
    
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;

   /* Setto gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Matrici identita' */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      /* Contributo di forza all'equazione della forza */
      WM.PutCoef(iCnt, 12+iCnt, 1.);
      
      /* Contrib. di der. di forza all'eq. della der. della forza */
      WM.PutCoef(6+iCnt, 15+iCnt, 1.);
      
      /* Equazione di vincolo */
      WM.PutCoef(12+iCnt, iCnt, -1.);
      
      /* Derivata dell'equazione di vincolo */
      WM.PutCoef(15+iCnt, 6+iCnt, -1.);
   }
   
   /* Recupera i dati */
   const Mat3x3& R(pNode->GetRRef());
   const Vec3& Omega(pNode->GetWRef());
   /* F e' stata aggiornata da InitialAssRes */
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Matrici F/\d/\ */
   Mat3x3 FWedgedWedge(MatCrossCross, F, dTmp);
   
   /* Matrici (omega/\d)/\ */
   Mat3x3 OWedgedWedge(MatCross, Omega.Cross(dTmp));
   
   /* Equazione di momento */
   WM.Add(4, 4, FWedgedWedge);
   WM.Add(4, 13, Mat3x3(MatCross, dTmp));
   
   /* Derivata dell'equazione di momento */
   WM.Add(10, 4, (Mat3x3(MatCross, FPrime) + Mat3x3(MatCrossCross, F, Omega))*Mat3x3(MatCross, dTmp));
   WM.Add(10, 10, FWedgedWedge);
   WM.Add(10, 13, OWedgedWedge);
   WM.Add(10, 16, Mat3x3(MatCross, dTmp));
   
   /* Equazione di vincolo */
   WM.Add(13, 4, Mat3x3(MatCross, dTmp));
   
   /* Derivata dell'equazione di vincolo */
   WM.Add(16, 4, OWedgedWedge);
   WM.Add(16, 10, Mat3x3(MatCross, dTmp));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PinJoint::InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PinJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
   
   /* Indici */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati */
   const Vec3& x(pNode->GetXCurr());
   const Vec3& v(pNode->GetVCurr());
   const Mat3x3& R(pNode->GetRCurr());
   const Vec3& Omega(pNode->GetWCurr());
   F = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 FPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanza nel sistema globale */
   Vec3 dTmp(R*d);

   /* Vettori omega/\d */
   Vec3 OWedged(Omega.Cross(dTmp));
   
   /* Equazioni di equilibrio */
   WorkVec.Sub(1, F);
   WorkVec.Add(4, F.Cross(dTmp)); /* Sfrutto il fatto che F/\d = -d/\F */
   
   /* Derivate delle equazioni di equilibrio */
   WorkVec.Sub(7, FPrime);
   WorkVec.Add(10, FPrime.Cross(dTmp)-OWedged.Cross(F));
   
   /* Equazione di vincolo */
   WorkVec.Add(13, x+dTmp-X0);
   
   /* Derivata dell'equazione di vincolo */
   WorkVec.Add(16, v+OWedged);
      
   return WorkVec;
}

const MBUnits::Dimensions
PinJoint::GetEquationDimension(integer index) const {
	// DOF == 3
   MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = MBUnits::Dimensions::Length;
			break;
		case 2:
			dimension = MBUnits::Dimensions::Length;
			break;
      case 3:
			dimension = MBUnits::Dimensions::Length;
			break;
	}

	return dimension;
}

std::ostream&
PinJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"position constraints" << std::endl;

	return out;
}
/* PinJoint - end */
