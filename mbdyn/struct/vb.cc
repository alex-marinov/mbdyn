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

/* Cerniera deformabile */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "vb.h"

#include "matvecexp.h"
#include "Rot.hh"

/* Costruttore non banale */
ViscousBody::ViscousBody(unsigned int uL,
	const DofOwner* pDO,
	const ConstitutiveLaw6D* pCL,
	const StructNode* pN,
	const Vec3& tilde_f,
	const Mat3x3& tilde_Rh,
	const OrientationDescription& od,
	flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
ConstitutiveLaw6DOwner(pCL),
pNode(pN),
tilde_f(tilde_f),
tilde_Rh(tilde_Rh),
od(od),
tilde_kPrime(Zero6),
bFirstRes(false)
{
	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);

	Rh = pNode->GetRRef()*tilde_Rh;

	/*
	 * Chiede la matrice tangente di riferimento
	 * e la porta nel sistema globale
	 */
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), Rh);
}


/* Distruttore */
ViscousBody::~ViscousBody(void)
{
	NO_OP;
}


/* Contributo al file di restart */
std::ostream&
ViscousBody::Restart(std::ostream& out) const
{
	Joint::Restart(out) << ", viscous body, "
		<< pNode->GetLabel() << ", position, reference, node, ",
	tilde_f.Write(out, ", ") << ", orientation, reference, node, 1, ",
	(tilde_Rh.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (tilde_Rh.GetVec(2)).Write(out, ", ") << ", ";
	return pGetConstLaw()->Restart(out) << ';' << std::endl;
}

void
ViscousBody::OutputPrepare(OutputHandler& OH)
{	
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("viscous body", OH, name);

			Var_v = OH.CreateVar<Vec3>(name + "V",
				MBUnits::Dimensions::Velocity,
				"local relative linear velocity (x, y, z)");
			
			Var_omega = OH.CreateVar<Vec3>(name + "Omega",
				MBUnits::Dimensions::AngularVelocity,
				"local relative angular velocity (x, y, z)");
		}
#endif // USE_NETCDF
	}
}


void
ViscousBody::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		Mat3x3 Rh(pNode->GetRCurr()*tilde_Rh);
		Vec3 F(GetF().GetVec1());
		Vec3 M(GetF().GetVec2());

		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "ViscousBody", GetLabel(),
					F, M, Rh*F, Rh*M);

			OH.Joints() << " " << tilde_kPrime << " " << std::endl;
		}

#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, F, M, Rh*F, Rh*M);
			OH.WriteNcVar(Var_v, tilde_kPrime.GetVec1());
			OH.WriteNcVar(Var_omega, tilde_kPrime.GetVec2());
		}
#endif // USE_NETCDF

	}
}

void
ViscousBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	if (ph) {
		/* pass to constitutive law */
		ConstitutiveLaw6DOwner::SetValue(pDM, X, XP, ph);
	}
}

Hint *
ViscousBody::ParseHint(DataManager *pDM, const char *s) const
{
	return ConstitutiveLaw6DOwner::ParseHint(pDM, s);
}

unsigned int
ViscousBody::iGetNumPrivData(void) const
{
	return 15 + ConstitutiveLaw6DOwner::iGetNumPrivData();
}

unsigned int
ViscousBody::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	unsigned idx = 0;

	switch (s[0]) {
	case 'r':
		break;

	case 'v':
		idx += 3;
		break;

	case 'w':
		idx += 6;
		break;

	case 'F':
		idx += 9;
		break;

	case 'M':
		idx += 12;
		break;

	default:
	{
		size_t l = STRLENOF("constitutiveLaw.");
		if (strncmp(s, "constitutiveLaw.", l) == 0) {
			idx = ConstitutiveLaw6DOwner::iGetPrivDataIdx(&s[l]);
			if (idx > 0) {
				return 15 + idx;
			}
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
ViscousBody::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0);

	ASSERT(i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	{
		Mat3x3 RhT((pNode->GetRCurr()*tilde_Rh).Transpose());

		Vec3 tilde_Theta(RotManip::VecRot(RhT));

		return tilde_Theta(i);
	}

	case 4:
	case 5:
	case 6:
	{
		Vec3 f(pNode->GetRCurr()*tilde_f);
		Mat3x3 RhT(pNode->GetRCurr().Transpose());
		Vec3 tilde_dPrime(RhT*(pNode->GetVCurr() - f.Cross(pNode->GetWCurr())));

		return tilde_dPrime(i - 3);
	}

	case 7:
	case 8:
	case 9:
	{
		Mat3x3 RhT((pNode->GetRCurr()*tilde_Rh).Transpose());
		Vec3 tilde_Omega = RhT*(pNode->GetWCurr());

		return tilde_Omega(i - 6);
	}

	case 10:
	case 11:
	case 12:
	case 13:
	case 14:
	case 15:
		return GetF()(i - 9);

	default:
		return ConstitutiveLaw6DOwner::dGetPrivData(i - 15);
	}
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
ViscousBody::AssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	AssMats(WM, WM, dCoef);

	return WorkMat;
}

/* Jacobian matrix assembly - all but Elastic */
void
ViscousBody::AssMats(VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	FullSubMatrixHandler& WMA = WorkMatA.SetFull();
	FullSubMatrixHandler& WMB = WorkMatB.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WMA.ResizeReset(iNumRows, iNumCols);
	WMB.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WMA.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WMA.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);

		WMB.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WMB.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	AssMats(WMA, WMB, 1.);
}

/* assemblaggio residuo */
SubVectorHandler&
ViscousBody::AssRes(SubVectorHandler& WorkVec,
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
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* inverse dynamics capable element */
bool
ViscousBody::bInverseDynamics(void) const
{
	return true;
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
ViscousBody::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */,
	const VectorHandler& /* XPrimeCurr */, 
	const VectorHandler& /* XPrimePrimeCurr */, 
	InverseDynamics::Order iOrder)
{	
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS);

	bFirstRes = false;

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNodeFirstMomIndex = pNode->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
ViscousBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

void
ViscousBody::AssMats(FullSubMatrixHandler& WMA,
	FullSubMatrixHandler& WMB, doublereal dCoef)
{
	// common
	WMA.Sub(1, 4, Mat3x3(MatCross, F.GetVec1()*dCoef));
	WMA.Sub(4, 4, Mat3x3(MatCross, F.GetVec2()*dCoef));

	// viscous
	const Mat3x3& F_dPrime = FDEPrime.GetMat11();
	const Mat3x3& F_thetaPrime = FDEPrime.GetMat12();
	const Mat3x3& M_dPrime = FDEPrime.GetMat21();
	const Mat3x3& M_thetaPrime = FDEPrime.GetMat22();

	WMB.Add(1, 1, F_dPrime);
	Mat3x3 MTmp2(M_dPrime + f.Cross(F_dPrime));
	WMB.Add(4, 1, MTmp2);
	Mat3x3 MTmp(F_thetaPrime - F_dPrime*Mat3x3(MatCross, f));
	WMB.Add(1, 4, MTmp);
	WMB.Add(4, 4, M_thetaPrime - M_dPrime*Mat3x3(MatCross, f) + f.Cross(MTmp));

	MTmp = Mat3x3(MatCross, f.Cross(pNode->GetWCurr()*dCoef));
	WMA.Add(1, 4, F_dPrime*MTmp);
	WMA.Add(4, 4, MTmp2*MTmp);
}

void
ViscousBody::AssVec(SubVectorHandler& WorkVec)
{
	if (bFirstRes) {
		bFirstRes = false;

	} else {
		Rh = pNode->GetRCurr()*tilde_Rh;
		f = pNode->GetRCurr()*tilde_f;

		Vec3 tilde_v = Rh.MulTV(pNode->GetVCurr() + pNode->GetWCurr().Cross(f));
		Vec3 tilde_omega = Rh.MulTV(pNode->GetWCurr());

		tilde_kPrime = Vec6(tilde_v, tilde_omega);

		ConstitutiveLaw6DOwner::Update(Zero6, tilde_kPrime);
	}

	F = MultRV(ConstitutiveLaw6DOwner::GetF(), Rh);

	WorkVec.Sub(1, F.GetVec1());
	WorkVec.Sub(4, f.Cross(F.GetVec1()) + F.GetVec2());
}

void
ViscousBody::AfterPredict(VectorHandler& /* X */ ,
	VectorHandler& /* XP */ )
{
	/* Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE */

	/* Recupera i dati */
	Rh = pNode->GetRRef()*tilde_Rh;
	f = pNode->GetRRef()*tilde_f;
	Vec3 tilde_v = Rh.MulTV(pNode->GetVCurr() + pNode->GetWCurr().Cross(f));
	Vec3 tilde_omega = Rh.MulTV(pNode->GetWCurr());

	tilde_kPrime = Vec6(tilde_v, tilde_omega);

	ConstitutiveLaw6DOwner::Update(Zero6, tilde_kPrime);

	/* FIXME: we need to be able to regenerate FDE
	 * if the constitutive law throws ChangedEquationStructure */
	FDEPrime = MultRMRt(ConstitutiveLaw6DOwner::GetFDEPrime(), Rh);

	bFirstRes = true;
}

void
ViscousBody::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	ConstitutiveLaw6DOwner::AfterConvergence(Zero6, tilde_kPrime);
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler&
ViscousBody::InitialAssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& /* XCurr */ )
{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	// FIXME: wrong
	AssMats(WM, WM, 1.);

	return WorkMat;
}

const MBUnits::Dimensions
ViscousBody::GetEquationDimension(integer index) const {
	// DOF == 0
	return MBUnits::Dimensions::UnknownDimension;
}
/* ViscousBody - end */

