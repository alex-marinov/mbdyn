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

/* Elementi elettrici */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "strnode.h"
#include "elecnode.h"
#include "accelerometer.h"

/* Accelerometer - begin */

/* Costruttore */
Accelerometer::Accelerometer(unsigned int uL, const DofOwner* pDO,
	const StructNode* pS,
	const ScalarDifferentialNode* pA,
	const Vec3& TmpDir,
	doublereal dO, doublereal dT,
	doublereal dC, doublereal dK,
	flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pStrNode(pS), pAbsNode(pA),
Dir(TmpDir), dOmega(dO), dTau(dT), dCsi(dC), dKappa(dK)
{
	NO_OP;
}

/* Distruttore banale */
Accelerometer::~Accelerometer(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
Accelerometer::Restart(std::ostream& out) const
{
	Electric::Restart(out) << ", accelerometer, "
		<< pStrNode->GetLabel() << ", "
		<< pAbsNode->GetLabel() << ", "
		"reference, node, ", Dir.Write(out, ", ") << ", "
		"node, "
		<< dOmega << ", "
		<< dTau << ", "
		<< dCsi << ", "
		<< dKappa << ';'
		<< std::endl;
	return out;
}

/* Costruisce il contributo allo jacobiano */
VariableSubMatrixHandler&
Accelerometer::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering Accelerometer::AssJac()" << std::endl);

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(15, 0);

	/* Indici delle equazioni */
	integer iFirstPositionIndex = pStrNode->iGetFirstPositionIndex();
	integer iAbstractIndex = pAbsNode->iGetFirstIndex();
	integer iFirstIndex = iGetFirstIndex();

	/*
	 * | c      0 0                0 -a^T  c*xP^T*a/\ || Delta_vP  |
	 * | 0      1 c*O^2/T          0  0    0          || Delta_y1P |
	 * | 0      0 1+c*(2*C*O+1/T) -c  0    0          || Delta_y2P |
	 * |-K*O^2 -c c*O*(1+2*C/T)    1  0    0          || Delta_zP  | = res
	 *                                                 | Delta_xP  |
	 *                                                 | Delta_gP  |
	 *
	 * con: c  = dCoef
	 *      a  = RNode*Dir
	 *      xp = Velocita' del nodo
	 *      O  = dOmega
	 *      T  = dTau
	 *      C  = dCsi
	 *      K  = dKappa
	 *
	 * v e' la misura della veocita' del punto,
	 * y1 e y2 sono stati dell'acceleromtro; v, y1, y2 appartengono al DofOwner
	 * dell'elemento;
	 * z e' la variabile del nodo astratto;
	 * x e g sono posizione e parametri di rotazione del nodo strutturale
	 *
	 *
	 * Funzione di trasferimento dell'accelerometro nel dominio di Laplace:
	 *
	 *  e0                        T * s
	 * ----(s) = K * ---------------------------------
	 * acc.          (1 + T*s)*(1 + 2*C/O*s + s^2/O^2)
	 *
	 */

	/* Dinamica dell'accelerometro */
	WM.PutItem(1, iFirstIndex + 1, iFirstIndex + 1, dCoef);
	WM.PutItem(2, iAbstractIndex + 1, iFirstIndex + 1,
		-dKappa*dOmega*dOmega);
	WM.PutItem(3, iFirstIndex + 2, iFirstIndex + 2, 1.);
	WM.PutItem(4, iAbstractIndex + 1, iFirstIndex + 2, -dCoef);
	WM.PutItem(5, iFirstIndex + 2, iFirstIndex + 3,
		dCoef*dOmega*dOmega/dTau);
	WM.PutItem(6, iFirstIndex + 3, iFirstIndex + 3,
		1. + dCoef*(2.*dCsi*dOmega + 1./dTau));
	WM.PutItem(7, iAbstractIndex + 1, iFirstIndex + 3,
	dCoef*dOmega*(dOmega + 2.*dCsi/dTau));
	WM.PutItem(8, iFirstIndex + 3, iAbstractIndex + 1, -dCoef);
	WM.PutItem(9, iAbstractIndex + 1, iAbstractIndex + 1, 1.);

	/* Misura dell'accelerazione */
	Vec3 TmpDir((pStrNode->GetRRef())*Dir);
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutItem(9 + iCnt, iFirstIndex + 1,
			iFirstPositionIndex + iCnt,
			- TmpDir.dGet(iCnt));
	}

	Vec3 XP(pStrNode->GetVCurr());
	TmpDir = -TmpDir.Cross(XP);
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutItem(12 + iCnt, iFirstIndex + 1,
			iFirstPositionIndex + 3 + iCnt,
			dCoef*TmpDir.dGet(iCnt));
	}

	return WorkMat;
}

/* Costruisce il contributo al residuo */
SubVectorHandler&
Accelerometer::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering Accelerometer::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iAbstractIndex = pAbsNode->iGetFirstIndex();
	integer iFirstIndex = iGetFirstIndex();

	/*
	 *        | Delta_vP  |   | -v+a^T*xP                      |
	 *        | Delta_y1P |   | -y1P-O^2/T*y2                  |
	 *  jac * | Delta_y2P | = | -y2P+z-(2*C*O+1/T)*y2          |
	 *        | Delta_zP  |   | -zP+y1-O*(O+2*C/T)*y2+K*O^2*vP |
	 *        | Delta_xP  |
	 *        | Delta_gP  |
	 *
	 * per il significato dei termini vedi AssJac
	 */

	WorkVec.PutRowIndex(1, iFirstIndex + 1);
	WorkVec.PutRowIndex(2, iFirstIndex + 2);
	WorkVec.PutRowIndex(3, iFirstIndex + 3);
	WorkVec.PutRowIndex(4, iAbstractIndex + 1);

	Vec3 XP(pStrNode->GetVCurr());
	Mat3x3 R(pStrNode->GetRCurr());
	doublereal v = XCurr(iFirstIndex + 1);
	doublereal vp = XPrimeCurr(iFirstIndex + 1);
	doublereal y1 = XCurr(iFirstIndex + 2);
	doublereal y1p = XPrimeCurr(iFirstIndex + 2);
	doublereal y2 = XCurr(iFirstIndex + 3);
	doublereal y2p = XPrimeCurr(iFirstIndex + 3);
	doublereal z = XCurr(iAbstractIndex + 1);
	doublereal zp = XPrimeCurr(iAbstractIndex + 1);

	WorkVec.PutCoef(1, (R*Dir).Dot(XP) - v);
	WorkVec.PutCoef(2, -y1p - dOmega*dOmega/dTau*y2);
	WorkVec.PutCoef(3, -y2p + z - (2.*dCsi*dOmega + 1./dTau)*y2);
	WorkVec.PutCoef(4, -zp + y1 - dOmega*(dOmega + 2.*dCsi/dTau)*y2
		+ dKappa*dOmega*dOmega*vp);

	return WorkVec;
}

unsigned int
Accelerometer::iGetNumDof(void) const
{
	return 3;
}

DofOrder::Order
Accelerometer::GetDofType(unsigned int i) const
{
	ASSERT(i >= 0 && i < 3);
	return DofOrder::DIFFERENTIAL;
}

void
Accelerometer::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 4;
	*piNumCols = 10;
}

void
Accelerometer::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}

void
Accelerometer::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	doublereal v = (pStrNode->GetRCurr()*Dir).Dot(pStrNode->GetVCurr());
	X.PutCoef(iGetFirstIndex() + 1, v);
}

/* Accelerometer - end */


/* TranslAccel - begin */

/* Costruttore */
TranslAccel::TranslAccel(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pS,
	const ScalarDifferentialNode* pA,
	const Vec3& TmpDir,
	const Vec3& Tmpf,
	flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pStrNode(pS), pAbsNode(pA),
Dir(TmpDir), f(Tmpf)
{
	ASSERT(pStrNode != NULL);
	ASSERT(pStrNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pAbsNode != NULL);
	ASSERT(pAbsNode->GetNodeType() == Node::ABSTRACT);
}

/* Distruttore banale */
TranslAccel::~TranslAccel(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
TranslAccel::Restart(std::ostream& out) const
{
	Electric::Restart(out) << ", accelerometer, translational, "
		<< pStrNode->GetLabel() << ", "
		<< pAbsNode->GetLabel() << ", "
		"reference, node, ", Dir.Write(out, ", ") << ", "
		"reference, node, ", f.Write(out, ", ") << ';'
		<< std::endl;
	return out;
}

/* Costruisce il contributo allo jacobiano */
VariableSubMatrixHandler&
TranslAccel::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering TranslAccel::AssJac()" << std::endl);

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(9, 0);

	/* Indici delle equazioni */
	integer iFirstColIndex = pStrNode->iGetFirstColIndex();
	integer iAbstractIndex = pAbsNode->iGetFirstIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WM.PutItem(1, iAbstractIndex, iAbstractIndex, dCoef);
	WM.PutItem(2, iAbstractIndex, iFirstIndex, -1.);
	WM.PutItem(3, iFirstIndex, iFirstIndex, dCoef);

	Vec3 tmpf = pStrNode->GetRCurr()*f;
	Vec3 tmpd = pStrNode->GetRCurr()*Dir;
	Vec3 tmp = tmpf.Cross(tmpd);
	WM.PutItem(4, iFirstIndex, iFirstColIndex + 1, -tmpd.dGet(1));
	WM.PutItem(5, iFirstIndex, iFirstColIndex + 2, -tmpd.dGet(2));
	WM.PutItem(6, iFirstIndex, iFirstColIndex + 3, -tmpd.dGet(3));

	tmp = tmpd.Cross((pStrNode->GetVCurr()*(-dCoef) + tmpf));
	WM.PutItem(7, iFirstIndex, iFirstColIndex + 4, tmp.dGet(1));
	WM.PutItem(8, iFirstIndex, iFirstColIndex + 5, tmp.dGet(2));
	WM.PutItem(9, iFirstIndex, iFirstColIndex + 6, tmp.dGet(3));

	return WorkMat;
}

/* Costruisce il contributo al residuo */
SubVectorHandler&
TranslAccel::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering TranslAccel::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	WorkVec.Resize(2);

	integer iAbstractIndex = pAbsNode->iGetFirstIndex() + 1;
	integer iFirstIndex = iGetFirstIndex()+1;

	WorkVec.PutRowIndex(1, iAbstractIndex);
	WorkVec.PutRowIndex(2, iFirstIndex);

	Vec3 tmpf = pStrNode->GetRCurr()*f;
	Vec3 tmpd = pStrNode->GetRCurr()*Dir;

	doublereal v = XCurr(iFirstIndex);
	doublereal vp = XPrimeCurr(iFirstIndex);
	doublereal a = pAbsNode->dGetX();

	WorkVec.PutCoef(1, vp - a);
	WorkVec.PutCoef(2, tmpd.Dot((pStrNode->GetVCurr() + pStrNode->GetWCurr().Cross(tmpf))) - v);

	return WorkVec;
}

unsigned int
TranslAccel::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
TranslAccel::GetDofType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::DIFFERENTIAL;
}

void
TranslAccel::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 8;
}

void
TranslAccel::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}

void
TranslAccel::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	doublereal v =
		(pStrNode->GetRCurr()*Dir).Dot(pStrNode->GetVCurr()
			+ pStrNode->GetWCurr().Cross(pStrNode->GetRCurr()*f));
	X.PutCoef(iGetFirstIndex() + 1, v);
	XP.PutCoef(iGetFirstIndex() + 1, 0.);
}

/* TranslAccel - end */


/* RotAccel - begin */

/* Costruttore */
RotAccel::RotAccel(unsigned int uL,
	const DofOwner* pDO,
	const StructNode* pS,
	const ScalarDifferentialNode* pA,
	const Vec3& TmpDir,
	flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pStrNode(pS), pAbsNode(pA),
Dir(TmpDir)
{
	ASSERT(pStrNode != NULL);
	ASSERT(pStrNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pAbsNode != NULL);
	ASSERT(pAbsNode->GetNodeType() == Node::ABSTRACT);
}

/* Distruttore banale */
RotAccel::~RotAccel(void)
{
	NO_OP;
}

/* Contributo al file di restart */
std::ostream&
RotAccel::Restart(std::ostream& out) const
{
	Electric::Restart(out) << ", accelerometer, rotational, "
		<< pStrNode->GetLabel() << ", "
		<< pAbsNode->GetLabel() << ", "
		"reference, node, ", Dir.Write(out, ", ") << ';'
		<< std::endl;
	return out;
}

/* Costruisce il contributo allo jacobiano */
VariableSubMatrixHandler&
RotAccel::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering RotAccel::AssJac()" << std::endl);

	/* Casting di WorkMat */
	SparseSubMatrixHandler& WM = WorkMat.SetSparse();
	WM.ResizeReset(6, 0);

	/* Indici delle equazioni */
	integer iFirstColIndex = pStrNode->iGetFirstColIndex();
	integer iAbstractIndex = pAbsNode->iGetFirstIndex()+1;
	integer iFirstIndex = iGetFirstIndex()+1;

	WM.PutItem(1, iAbstractIndex, iAbstractIndex, dCoef);
	WM.PutItem(2, iAbstractIndex, iFirstIndex, -1.);
	WM.PutItem(3, iFirstIndex, iFirstIndex, dCoef);

	Vec3 tmp = pStrNode->GetRCurr()*Dir;
	WM.PutItem(4, iFirstIndex, iFirstColIndex+4, tmp.dGet(1));
	WM.PutItem(5, iFirstIndex, iFirstColIndex+5, tmp.dGet(2));
	WM.PutItem(6, iFirstIndex, iFirstColIndex+6, tmp.dGet(3));

	return WorkMat;
}

/* Costruisce il contributo al residuo */
SubVectorHandler&
RotAccel::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering RotAccel::AssRes()" << std::endl);

	/* Dimensiona e resetta la matrice di lavoro */
	WorkVec.Resize(2);

	integer iAbstractIndex = pAbsNode->iGetFirstIndex() + 1;
	integer iFirstIndex = iGetFirstIndex() + 1;

	WorkVec.PutRowIndex(1, iAbstractIndex);
	WorkVec.PutRowIndex(2, iFirstIndex);

	Vec3 tmp = pStrNode->GetRCurr()*Dir;

	doublereal v = XCurr(iFirstIndex);
	doublereal vp = XPrimeCurr(iFirstIndex);
	doublereal a = pAbsNode->dGetX();

	WorkVec.PutCoef(1, vp - a);
	WorkVec.PutCoef(2, tmp.Dot(pStrNode->GetWCurr()) - v);

	return WorkVec;
}

unsigned int
RotAccel::iGetNumDof(void) const
{
	return 1;
}

DofOrder::Order
RotAccel::GetDofType(unsigned int i) const
{
	ASSERT(i == 0);
	return DofOrder::DIFFERENTIAL;
}

void
RotAccel::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 2;
	*piNumCols = 5;
}

void
RotAccel::SetInitialValue(VectorHandler& /* X */ )
{
	NO_OP;
}

void
RotAccel::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* XP */ ,
	SimulationEntity::Hints *ph)
{
	doublereal v = (pStrNode->GetRCurr()*Dir).Dot(pStrNode->GetWCurr());
	X.PutCoef(iGetFirstIndex() + 1, v);
}

/* RotAccel - end */

