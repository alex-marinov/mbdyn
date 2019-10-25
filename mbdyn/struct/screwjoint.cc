/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2009
 *
 * Marco Morandini	<morandini@aero.polimi.it>
 * Mauro Manetti	<manetti@aero.polimi.it>
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

/* Cerniera pilotata */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "screwjoint.h"
#include "Rot.hh"
#include "hint_impl.h"

// #define _GNU_SOURCE 1
// #include <fenv.h>
// static void __attribute__ ((constructor))
// trapfpe ()
// {
//   /* Enable some exceptions.  At startup all exceptions are masked.  */
// 
//   feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
// }

/* ScrewJoint - begin */

void unwrap(doublereal & theta,
		const doublereal & thetacurr,
		doublereal & thetaprev,
		integer & ntheta)
{
	if (std::abs(thetacurr-thetaprev) > M_PI) {
		if (thetacurr > 0.) {
			ntheta -= 1;
		} else {
			ntheta += 1;
		}
	}
	theta = ntheta * 2 * M_PI + thetacurr;
	thetaprev = thetacurr;
}



/* Costruttore non banale */
ScrewJoint::ScrewJoint(unsigned int uL,
				const DofOwner* pDO,
				const StructNode* pN1,
				const StructNode* pN2,
				const Vec3& f1Tmp, 
				const Vec3& f2Tmp,
				const Mat3x3& R1t,
//				const Mat3x3& R2t,
				const doublereal& p, 
				flag fOut,
				BasicShapeCoefficient *const sh,
				BasicFriction *const f
			)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
nTheta(0),
dThetaPrev(0.), dThetaCurr(0.), dTheta(0.), dD(0.), dTheta0(0.), dD0(0.),
dLambda(0.),
dPitch(p),
Theta(Zero3),
pNode1(pN1), pNode2(pN2), 
R1h(R1t), //R2h(R2t),
f1(f1Tmp), f2(f2Tmp),
GammaInv(Eye3),
F1(Zero3), C1(Zero3), C2(Zero3),
Sh_c(sh),
fc(f)
#ifdef USE_NETCDFC
,
Var_dTheta(0),
Var_Theta(0),
Var_vrel(0),
Var_fc(0),
Var_MFR(0)
#endif // USE_NETCDFC
{
	ASSERT(pNode1 != NULL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);

	Mat3x3 R1(pNode1->GetRCurr());
	Mat3x3 R2(pNode2->GetRCurr());
	Mat3x3 R1R1h(R1*R1h);
//	Mat3x3 R2R2h(R2*R2h);
	
	R1f1 = R1 * f1;
	R2f2 = R2 * f2;
	
	e1hz = R1R1h.GetVec(3);
	const_cast<doublereal&>(dTheta0) = dThetaPrev = dThetaCurr = dTheta =
		RotManip::VecRot(R1.MulTM(R2)).Dot(e1hz);
	D = pNode2->GetXCurr() + R2f2 - pNode1->GetXCurr() - R1f1;
	const_cast<doublereal&>(dD0) = dD = 
		D.Dot(e1hz);
	if (fc != 0) {
		ScrewJointSh_c *const ts = dynamic_cast<ScrewJointSh_c *const>(Sh_c);
		if (ts == 0) {
			std::cerr << "Screw joint " << uL << " fatal input error: with friction"
				" you can define only a \"screw joint\" shape function" <<
				std::endl;
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		} else {
			cos_pitch_angle_r = ts->ComputePitchAngle(dPitch);
		}
	}
}


/* Distruttore */
ScrewJoint::~ScrewJoint(void)
{
	NO_OP;
}


// /* Contributo al file di restart */
// std::ostream& //TODO
// ScrewJoint::Restart(std::ostream& out) const
// {
// 	Joint::Restart(out) << ", screw hinge, "
// 		<< pNode1->GetLabel() << ", reference, node, 1, ",
// 		(R1h.GetVec(3)).Write(out, ", ")
// 		<< ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
// 		<< pNode2->GetLabel() << ", reference, node, 1, ",
// 		(R2h.GetVec(3)).Write(out, ", ")
// 		<< ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << 
// 		';' << std::endl;
// 	return out;
// }

void
ScrewJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Screw Joint", OH, name);

			Var_dTheta = OH.CreateVar<doublereal>(name + "dTheta", "rad",
				"screw angle magnitude [deg]");
			Var_Theta = OH.CreateVar<Vec3>(name + "Theta", "-",
				"screw axis (x, y, z)");
			if (fc) {
				Var_vrel = OH.CreateVar<doublereal>(name + "vRel", "m/s",
					"contact point sliding velocity");
				Var_fc = OH.CreateVar<doublereal>(name + "fc", "-",
					"friction coefficient");
				Var_MFR = OH.CreateVar<doublereal>(name + "MFR", "Nm",
					"friction moment");

			}
		}
#endif // USE_NETCDF
	}
}

void //TODO
ScrewJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 R1(pNode1->GetRCurr()*R1h);
		Vec3 FTilde(F1 * dLambda);
		Vec3 MTilde(C1 * dLambda + e1hz*M3diff);
		Vec3 F(R1*F1 * dLambda);
		Vec3 M(R1*(C1 * dLambda + e1hz*M3diff));
		
		if (OH.UseText(OutputHandler::JOINTS)) {
			std::ostream &of = Joint::Output(OH.Joints(), "ScrewJoint", GetLabel(),
					FTilde, MTilde, F, M)
				<< " " << dTheta*dRaDegr << " " << dD
				<< " " << Theta << " " << D;
			if (fc) {
				of << " " << vrel << " " << fc->fc() << " " << e1hz*M3diff;
			}	
			of << std::endl;
		}
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, FTilde, MTilde, F, M);
			OH.WriteNcVar(Var_dTheta, dD);
			OH.WriteNcVar(Var_Theta, D);

			if (fc) {
				OH.WriteNcVar(Var_vrel, vrel);
				OH.WriteNcVar(Var_fc, fc->fc());
				OH.WriteNcVar(Var_MFR, e1hz*M3diff);
			}
		}
#endif // USE_NETCDF
	}
}


// void //TODO
// ScrewJoint::SetValue(DataManager *pDM,
// 		VectorHandler& X, VectorHandler& XP,
// 		SimulationEntity::Hints *ph)
// {
// 	if (ph) {
// 		for (unsigned i = 0; i < ph->size(); i++) {
// 			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);
// 
// 			if (pjh) {
// 				if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
// 					(Mat3x3&)R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2h;
// 
// 				} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
// 					(Mat3x3&)R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1h;
// 
// 				} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
// 					/* TODO */
// 				}
// 				continue;
// 			}
// 		}
// 	}
// }

void
ScrewJoint::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{

	R1 = pNode1->GetRCurr();
	R2 = pNode2->GetRCurr();
	R1R1h = R1*R1h;
//	R2R2h = R2*R2h;
	
	R1f1 = R1 * f1;
	R2f2 = R2 * f2;
	
	Theta = RotManip::VecRot(R1.Transpose()*R2);
	e1hz = R1R1h.GetVec(3);

	if (fc) {
		Vec3 Omega1(pNode1->GetWCurr());
		Vec3 Omega2(pNode2->GetWCurr());

		Vec3 OmegaTheta = R1.MulTV(Omega2 - Omega1);
		doublereal vrel = e1hz.Dot(OmegaTheta) * cos_pitch_angle_r;

		//compute
		doublereal modF, v;
			 //relative velocity
		v = vrel;
			 //reaction norm
		dLambda = X(iGetFirstIndex()+1);
		modF = std::abs(dLambda);
		fc->AfterConvergence(modF,v,X,XP,iGetFirstIndex()+1);
	}
}

// Hint * //TODO
// ScrewJoint::ParseHint(DataManager *pDM, const char *s) const
// {
// 	if (strncasecmp(s, "hinge{" /*}*/ , STRLENOF("hinge{" /*}*/ )) == 0)
// 	{
// 		s += STRLENOF("hinge{" /*}*/ );
// 
// 		if (strcmp(&s[1], /*{*/ "}") != 0) {
// 			return 0;
// 		}
// 
// 		switch (s[0]) {
// 		case '1':
// 			return new Joint::HingeHint<1>;
// 
// 		case '2':
// 			return new Joint::HingeHint<2>;
// 		}
// 	}
// 
// 	return SimulationEntity::ParseHint(pDM, s);
// }

std::ostream&
ScrewJoint::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 1 << ": "
			"reaction force [f]" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 1 << ": "
				"reaction force derivatives [fP]" << std::endl;
	}
	iIndex += 1;

	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction dof(s)" << std::endl
				<< "        ", fc->DescribeDof(out, prefix, bInitial);
		}
	}

	return out;
}

static const char *dof[] = { "reaction force f", "reaction force derivative fP" };
static const char *eq[] = { "screw constraint g", "screw constraint derivative w" };

void //TODO
ScrewJoint::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	ASSERTMSGBREAK((i == 1) || (i == -1), "INDEX ERROR in ScrewJoint::DescribeDof");
	if (i == -1) {
		if (bInitial) {
			iend = 1;

		} else {
			iend = 2;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "ScrewJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << dof[i];

			desc[i] = os.str();
		}

	} else {
		os << ": " << dof[i];
		desc[0] = os.str();
	}
}

std::ostream&
ScrewJoint::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 1 << ": "
			"screw constraint [g]" << std::endl;

	if (bInitial) {
		iIndex += 1;
		out
			<< prefix << iIndex + 1 << "->" << iIndex + 1 << ": "
				"screw constraint derivative [gP]" << std::endl;
	}

	iIndex += 1;
	if (fc) {
		integer iFCDofs = fc->iGetNumDof();
		if (iFCDofs > 0) {
			out << prefix << iIndex + 1;
			if (iFCDofs > 1) {
				out << "->" << iIndex + iFCDofs;
			}
			out << ": friction equation(s)" << std::endl
				<< "        ", fc->DescribeEq(out, prefix, bInitial);
		}
	}

	return out;

}

void //TODO
ScrewJoint::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	int iend = 1;
	ASSERTMSGBREAK((i == 1) || (i == -1), "INDEX ERROR in ScrewJoint::DescribeEq");
	if (i == -1) {
		if (bInitial) {
			iend = 2;

		} else {
			iend = 1;
		}
	}
	desc.resize(iend);

	std::ostringstream os;
	os << "ScrewJoint(" << GetLabel() << ")";

	if (i == -1) {
		std::string name = os.str();
		for (i = 0; i < iend; i++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": " << eq[i];

			desc[i] = os.str();
		}

	} else {
		os << ": " << eq[i];
		desc[0] = os.str();
	}
}

/* Dati privati (aggiungere magari le reazioni vincolari) */
unsigned int //TODO
ScrewJoint::iGetNumPrivData(void) const
{
	return 8;
};

unsigned int //TODO
ScrewJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);
	ASSERT(s[0] != '\0');

	unsigned int idx = 0;

	switch (s[0]) {
	case 'F':
		break;
	case 'M':
		idx += 3;
		break;
	case 'd':
		idx = 7;
		break;
	case 'r':
		idx = 8;
		break;

	default:
		return 0;
	}

	if (s[1] == '\0' || s[2] != '\0') {
		if ( idx > 6 ) {
			return idx;
		} else {
			return 0;
		}
	}

	switch (s[1]) {
	case 'x':
		return idx + 1;
	case 'y':
		return idx + 2;
	case 'z':
		return idx + 3;
	}

	return 0;
}

doublereal //TODO
ScrewJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 8);

	switch (i) {
	case 1:
	case 2:
	case 3:
		return F1(i) * dLambda;

	case 4:
	case 5:
	case 6:
		return C1(i - 3) * dLambda;
	case 7:
		return dD;
	case 8:
		return dTheta;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ScrewJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering ScrewJoint::AssJac()" << std::endl);

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
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
	}

	AssMat(WM, dCoef, XCurr, XPrimeCurr);

	return WorkMat;
}


void//TODO
ScrewJoint::AfterPredict(VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
// 	/* Recupera i dati */
// 	R1Ref = pNode1->GetRRef()*R1h;
// 	Mat3x3 R1T = R1Ref.Transpose();
// 	Mat3x3 RD(R1T*pNode2->GetRRef()*R2h);
// 
// 	/* Calcola la deformazione corrente nel sistema locale (nodo a) */
// 	ThetaCurr = ThetaRef = RotManip::VecRot(RD);
// 
// 	/* Calcola l'inversa di Gamma di ThetaRef */
// 	Mat3x3 GammaRefm1 = RotManip::DRot_I(ThetaRef);
// 
// 	/* Contributo alla linearizzazione ... */
// 	RRef = GammaRefm1*R1T;
// 
// 	/* Flag di aggiornamento dopo la predizione */
// 	bFirstRes = true;
}


void//TODO
ScrewJoint::AssMat(FullSubMatrixHandler& WM, doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	
	doublereal dLambdadCoef = dLambda * dCoef;
	Mat3x3 R1GammaInvT = R1.MulMT(GammaInv);
	Mat3x3 HorridusMatrix = R1GammaInvT * (
					RotManip::Elle(-Theta, GammaInv.Transpose() * e1hz) * GammaInv.MulMT(R1)
				);
	Mat3x3 e1hzxdLambdadCoef(MatCross, e1hz * dLambdadCoef);
	//F1 - g1
	WM.Sub(1, 4, e1hzxdLambdadCoef );
	//F1 - lambda
	WM.Add(1, 13, e1hz);
	//F2 - g1
	WM.Add(7, 4, e1hzxdLambdadCoef );
	//F2 - lambda
	WM.Sub(7, 13, e1hz);
	
	//C1 - x1
	WM.Add(4, 1, e1hzxdLambdadCoef );
	//C1 - x2
	WM.Sub(4, 7, e1hzxdLambdadCoef );
	//C1 - g1
	WM.Add(4, 4, (
			(Theta.Cross(Mat3x3(MatCross, e1hz)) -
					e1hz.Cross(R1GammaInvT.Transpose()) +
					Mat3x3(MatCross, R1GammaInvT*e1hz) +
					HorridusMatrix +
					R1GammaInvT * Mat3x3(MatCross, e1hz)
			) * (dPitch / ( 2. * M_PI )) -
			Mat3x3(MatCrossCross, e1hz, R1f1) - 
			Mat3x3(MatCrossCross, D, e1hz) +
			Mat3x3(MatCrossCross, e1hz, R1f1) -
			Mat3x3(MatCrossCross, R1f1, e1hz)
		) * dLambdadCoef
	);
	//C1 - g2
	WM.Add(4, 10, (
			(
				e1hz.Cross(R1GammaInvT.Transpose()) -
				HorridusMatrix
			) * (dPitch / ( 2. * M_PI )) +
			Mat3x3(MatCrossCross, e1hz, R2f2)
		) * dLambdadCoef
	);
	//C1 - lambda
	WM.Add(4, 13, C1);
	//C2 - g1
	WM.Add(10, 4, (
			(
				- Mat3x3(MatCross, R1GammaInvT*e1hz) - 
				HorridusMatrix
			) * (dPitch / ( 2. * M_PI )) +
			Mat3x3(MatCrossCross, R2f2, e1hz)
		) * dLambdadCoef
	);
	//C2 - g2
	WM.Add(10, 10,(
			(
				HorridusMatrix
			) * (dPitch / ( 2. * M_PI )) -
			Mat3x3(MatCrossCross, e1hz, R2f2)
		) * dLambdadCoef
	); 
	//C2 - lambda
	WM.Add(10, 13, C2);
	//eq - x1
	WM.AddT(13, 1, F1);
	//eq - g1
	WM.AddT(13, 4, C1);
	//eq - x2
	WM.SubT(13, 7, F1);
	//eq - g2
	WM.AddT(13, 10, C2);


	if (fc) {
		integer iFirstReactionIndex = iGetFirstIndex();
		//retrive
			 //friction coef
		doublereal f = fc->fc();
			 //shape function
		doublereal shc = Sh_c->Sh_c();
			 //omega and omega rif
		Vec3 Omega1(pNode1->GetWCurr());
		Vec3 Omega2(pNode2->GetWCurr());

		Vec3 OmegaTheta = R1.MulTV(Omega2 - Omega1);
// 		doublereal vrel = e1hz.Dot(OmegaTheta) * cos_pitch_angle_r;


		Vec3 Omega1r(pNode1->GetWRef());
		Vec3 Omega2r(pNode2->GetWRef());	
		//compute
		doublereal modF, v;
			 //relative velocity
		doublereal F_sign = 1.;
		if (dLambda < 0.) {
			F_sign = -1;
		}
		v = vrel * F_sign;
			 //reaction norm
		modF = std::abs(dLambda);
			 //reaction moment
		//doublereal M3 = shc*modF*r;
		
		ExpandableRowVector dfc;
		ExpandableRowVector dF, dF0;
		ExpandableRowVector dv;
			 //variation of reaction force
		dF0.ReDim(1);
// 		if ((modF[0] == 0.) or (F.Norm() < preF)) {
		if ((modF == 0.)) {
			 dF0.Set(0., 1, 12+1);
		} else {
			 dF0.Set(dLambda > 0 ? 1. : -1., 1, 12+1);
		}
		dF.ReDim(1);
		dF.Set(0., 1, 12+1);
			 //variation of relative velocity
		dv.ReDim(6);
		
		Mat3x3 tempmat1(MatCross, R1.MulTV(Omega2 - Omega1));
		tempmat1 += R1.MulTM(Mat3x3(MatCross, Omega2 - Omega1));
		Vec3 tempvec(tempmat1.MulTV(e1hz));
		dv.Set(tempvec * cos_pitch_angle_r * (dCoef * F_sign) - 
			(e1hz + Omega1r.Cross(e1hz * dCoef)) * (cos_pitch_angle_r * F_sign), 1, 4);
		dv.Set((e1hz + Omega2r.Cross(e1hz * dCoef)) * (cos_pitch_angle_r * F_sign), 4, 6+4);

		ExpandableMatrix de1hz;
		de1hz.ReDim(3,1);
			de1hz.SetBlockDim(1, 3);
			de1hz.SetBlockIdx(1, 4);
		de1hz.Set(Mat3x3(MatCross, e1hz * (-dCoef)), 1, 1);

		//assemble friction states
		fc->AssJac(WM,dfc,12+1,iFirstReactionIndex+1,dCoef,modF,v,
				XCurr,XPrimeCurr,dF0,dv);

		ExpandableMatrix dM3diff;
		ExpandableRowVector dShc;
		//compute 
			 //variation of shape function
		Sh_c->dSh_c(dShc,f,modF,v,dfc,dF0,dv);
			 //variation of moment component
		dM3diff.ReDim(3, 3);
		dM3diff.SetBlockDim(1, 1);
		dM3diff.SetBlockDim(2, 1);
		dM3diff.SetBlockDim(3, 3);
		dM3diff.SetBlockIdx(3, 4);
		
		dM3diff.Set(e1hz*shc, 1, 1); dM3diff.Link(1, &dF);
		dM3diff.Set(e1hz*modF, 1, 2); dM3diff.Link(2, &dShc);
		dM3diff.Set(Mat3x3(MatCross, e1hz * (-dCoef * shc * modF)), 1, 3); 
		//dM3diff.Set(Mat3x3(1.) * shc * modF[0], 1, 3); 
// 		std::cerr << "Ci siamo\n";
		
		//dM3diff.Link(3, &de1hz);
		//assemble first node
			 //variation of moment component
// 		std::cerr << "Ci risiamo\n";
		dM3diff.Add(WM, 4, 1.);
		//assemble second node
			 //variation of moment component
		dM3diff.Sub(WM, 6+4, 1.);
	}

}


/* assemblaggio residuo */
SubVectorHandler&
ScrewJoint::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering ScrewJoint::AssRes()" << std::endl);
	
	
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
	integer iFirstReactionIndex = iGetFirstIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}
	for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
		WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
	}

	dLambda = XCurr(iFirstReactionIndex+1);

	AssVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;
}


void
ScrewJoint::AssVec(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr)
{
	integer iFirstReactionIndex = iGetFirstIndex();
	R1 = pNode1->GetRCurr();
	R2 = pNode2->GetRCurr();
	R1R1h = R1*R1h;
//	R2R2h = R2*R2h;
	
	R1f1 = R1 * f1;
	R2f2 = R2 * f2;
	
	Theta = RotManip::VecRot(R1.Transpose()*R2);
	e1hz = R1R1h.GetVec(3);
	dThetaCurr = Theta.Dot(e1hz);
	unwrap(dTheta, dThetaCurr, dThetaPrev, nTheta);
	GammaInv = RotManip::DRot_I(Theta);
	
	D = pNode2->GetXCurr() + R2f2 - pNode1->GetXCurr() - R1f1;
	dD = D.Dot(e1hz);

	F1 = e1hz;
	//Vec3 F2 = -F1;
	C1 = ( - ( Theta.Cross(e1hz) + R1.MulMT(GammaInv)*e1hz ) * (dPitch / ( 2. * M_PI )) + 
		D.Cross(e1hz) + R1f1.Cross(e1hz) );
	C2 = ( R1.MulMT(GammaInv) * e1hz * (dPitch / ( 2. * M_PI )) -
		R2f2.Cross(e1hz) );
		
	doublereal eq = dPitch / (2. * M_PI) * (dTheta - dTheta0) - (dD - dD0);
// 	std::cerr << "XXXXXXXXXXX " << nTheta << " " << dTheta << " " << dThetaCurr << " " << dThetaPrev << 
// 		" " << dTheta0 << " " << dD << " " << dD0 <<  
// 		" " << eq  << " " << dCoef << 
// 		" " << dPitch / (2. * M_PI) * (dTheta - dTheta0) << 
// 		" " << (dD - dD0) << std::endl;

// 	if (bFirstRes) {
// 		/* La rotazione e' gia' stata aggiornata da AfterPredict */
// 		bFirstRes = false;
// 
// 	} else {
// 		Mat3x3 R2(pNode2->GetRCurr()*R2h);
// 		ThetaCurr = RotManip::VecRot(R1.Transpose()*R2);
// 	}

	WorkVec.Sub(1, F1 * dLambda);
	WorkVec.Sub(3 + 1, C1 * dLambda);
	WorkVec.Add(6 + 1, F1 * dLambda);
	WorkVec.Sub(6 + 3 + 1, C2 * dLambda);
	ASSERT(dCoef != 0.);
	WorkVec.DecCoef(13, eq / dCoef);

	if (fc) {
		bool ChangeJac(false);
		
		Vec3 Omega1(pNode1->GetWCurr());
		Vec3 Omega2(pNode2->GetWCurr());
		Vec3 OmegaTheta = R1.MulTV(Omega2 - Omega1);
		vrel = e1hz.Dot(OmegaTheta) * cos_pitch_angle_r;
		
		
		doublereal modF, v;
		doublereal F_sign = 1.;
		if (dLambda < 0.) {
			F_sign = -1;
		}
		v = vrel * F_sign;
		modF = std::abs(dLambda);
		try {
			fc->AssRes(WorkVec,12+1,iFirstReactionIndex+1,modF,v,XCurr,XPrimeCurr);
		}
		catch (Elem::ChangedEquationStructure& err) {
			ChangeJac = true;
		}
		doublereal f = fc->fc();
// 		std::cerr << f << " " << vrel << " " << dLambda << std::endl;
		doublereal shc = Sh_c->Sh_c(f,modF,v);
		
		M3diff = shc*dLambda;
		WorkVec.Sub(4, e1hz*M3diff);
		WorkVec.Add(10, e1hz*M3diff);
//!!!!!!!!!!!!!!
//		M += e3a*M3;
		if (ChangeJac) {
			throw Elem::ChangedEquationStructure(MBDYN_EXCEPT_ARGS);
		}
	}

}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ScrewJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ScrewJoint::InitialAssJac()" << std::endl);

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

// 	/* Recupera gli indici */
// 	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
// 	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
// 	integer iFirstReactionIndex = iGetFirstIndex();
// 	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
// 	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
// 	integer iReactionPrimeIndex = iFirstReactionIndex + 3;
// 
// 	/* Setta gli indici della matrice */
// 	for (int iCnt = 1; iCnt <= 3; iCnt++) {
// 		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
// 		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
// 		WM.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
// 		WM.PutColIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
// 		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
// 		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
// 		WM.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
// 		WM.PutColIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
// 		WM.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
// 		WM.PutColIndex(12 + iCnt, iFirstReactionIndex + iCnt);
// 		WM.PutRowIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
// 		WM.PutColIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
// 	}
// 
// 	Mat3x3 Ra(pNode1->GetRRef()*R1h);
// 	Mat3x3 RaT(Ra.Transpose());
// 	Vec3 Wa(pNode1->GetWRef());
// 	Vec3 Wb(pNode2->GetWRef());
// 
// 	Mat3x3 MTmp(M);
// 	Mat3x3 MPrimeTmp(Ra*Vec3(XCurr, iReactionPrimeIndex + 1));
// 
// 	WM.Add(1, 1, MTmp);
// 	WM.Add(3 + 1, 3 + 1, MTmp);
// 	WM.Sub(6 + 1, 1, MTmp);
// 	WM.Sub(9 + 1, 3 + 1, MTmp);
// 
// 	MTmp = Mat3x3(Wa)*MTmp + MPrimeTmp;
// 	WM.Add(3 + 1, 1, MTmp);
// 	WM.Sub(9 + 1, 1, MTmp);
// 
// 	WM.Add(6 + 1, 12 + 1, Ra);
// 	WM.Add(9 + 1, 15 + 1, Ra);
// 	WM.Sub(1, 12 + 1, Ra);
// 	WM.Sub(3 + 1, 15 + 1, Ra);
// 
// 	MTmp = Mat3x3(Wa)*Ra;
// 	WM.Add(9 + 1, 12 + 1, MTmp);
// 	WM.Sub(3 + 1, 12 + 1, MTmp);
// 
// 	WM.Add(12 + 1, 6 + 1, RaT);
// 	WM.Sub(12 + 1, 1, RaT);
// 	WM.Sub(15 + 1, 3 + 1, RaT);
// 	WM.Add(15 + 1, 9 + 1, RaT);
// 	WM.Add(15 + 1, 1, RaT*Mat3x3(Wb - Wa));
// 
	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ScrewJoint::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering ScrewJoint::InitialAssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

// 	/* Recupera gli indici */
// 	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex() + 3;
// 	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex() + 3;
// 	integer iFirstReactionIndex = iGetFirstIndex();
// 	integer iNode1FirstVelIndex = iNode1FirstPosIndex + 6;
// 	integer iNode2FirstVelIndex = iNode2FirstPosIndex + 6;
// 	integer iReactionPrimeIndex = iFirstReactionIndex + 3;
// 
// 	/* Setta gli indici del vettore */
// 	for (int iCnt = 1; iCnt <= 3; iCnt++) {
// 		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
// 		WorkVec.PutRowIndex(3 + iCnt, iNode1FirstVelIndex + iCnt);
// 		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
// 		WorkVec.PutRowIndex(9 + iCnt, iNode2FirstVelIndex + iCnt);
// 		WorkVec.PutRowIndex(12 + iCnt, iFirstReactionIndex + iCnt);
// 		WorkVec.PutRowIndex(15 + iCnt, iReactionPrimeIndex + iCnt);
// 	}
// 
// 	Mat3x3 R1(pNode1->GetRCurr()*R1h);
// 	Mat3x3 R1T(R1.Transpose());
// 	Vec3 Wa(pNode1->GetWCurr());
// 	Vec3 Wb(pNode2->GetWCurr());
// 
// 	M = Vec3(XCurr, iFirstReactionIndex+1);
// 	Vec3 MPrime = Vec3(XCurr, iReactionPrimeIndex+1);
// 
// 	Vec3 MTmp(R1*M);
// 	Vec3 MPrimeTmp(Wa.Cross(MTmp) + R1*MPrime);
// 
// 	Mat3x3 R2(pNode2->GetRCurr()*R2h);
// 	ThetaCurr = RotManip::VecRot(R1T*R2);
// 
// 	Vec3 ThetaPrime = R1T*(Wb-Wa);
// 
// 	WorkVec.Add(1, MTmp);
// 	WorkVec.Add(3 + 1, MPrimeTmp);
// 	WorkVec.Sub(6 + 1, MTmp);
// 	WorkVec.Sub(9 + 1, MPrimeTmp);
// 	WorkVec.Add(12 + 1, Get() - ThetaCurr);
// 	if (bIsDifferentiable()) {
// 		ThetaPrime -= GetP();
// 	}
// 	WorkVec.Sub(15 + 1, ThetaPrime);

	return WorkVec;
}

/* ScrewJoint - end */
