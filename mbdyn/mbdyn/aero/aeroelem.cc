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

/* Elementi aerodinamici */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "aeroelem.h"
#include "aerodata.h"
#include "dataman.h"
#include "shapefnc.h"
#include "drive_.h"

#include "aerodc81.h"
#include "c81data.h"
#include "Rot.hh"

AerodynamicOutput::AerodynamicOutput(flag f, int iNP)
{
	SetOutputFlag(f, iNP);
}

AerodynamicOutput::~AerodynamicOutput(void)
{
	NO_OP;
}

void
AerodynamicOutput::SetOutputFlag(flag f, int iNP)
{
	m_eOutput = f;

	if (IsOutput() && IsPGAUSS()) {
		pOutput.resize(iNP);

	} else {
		pOutput.resize(0);
	}
}

void
AerodynamicOutput::ResetIterator(void)
{
	if (IsOutput() && IsPGAUSS()) {
		ASSERT(!pOutput.empty());
		pTmpOutput = pOutput.begin();
	}
}

void
AerodynamicOutput::SetData(const Vec3& v, const doublereal* pd)
{
	ASSERT(IsPGAUSS());
	ASSERT(!pOutput.empty());
	ASSERT(pTmpOutput >= pOutput.begin());
	ASSERT(pTmpOutput < pOutput.end());

	pTmpOutput->alpha = 180./M_PI*atan2(-v.dGet(2), v.dGet(1));
	pTmpOutput->f = Vec3(pd[1], pd[0], pd[5]);

	// move iterator forward
	pTmpOutput++;
}

AerodynamicOutput::eOutput
AerodynamicOutput::GetOutput(void) const
{
	return eOutput(m_eOutput & AEROD_OUT_MASK);
}

bool
AerodynamicOutput::IsOutput(void) const
{
	return (m_eOutput & 0x1);
}

bool
AerodynamicOutput::IsSTD(void) const
{
	return GetOutput() == AEROD_OUT_STD;
}

bool
AerodynamicOutput::IsPGAUSS(void) const
{
	return GetOutput() == AEROD_OUT_PGAUSS;
}

bool
AerodynamicOutput::IsNODE(void) const
{
	return GetOutput() == AEROD_OUT_NODE;
}

/* AerodynamicBody - begin */

AerodynamicBody::AerodynamicBody(unsigned int uLabel,
	const DofOwner *pDO,
	const StructNode* pN, InducedVelocity* pR,
	const Vec3& fTmp, doublereal dS,
	const Mat3x3& RaTmp,
	const Shape* pC, const Shape* pF,
	const Shape* pV, const Shape* pT,
	integer iN, AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut),
DriveOwner(pDC),
AerodynamicOutput(fOut, iN),
aerodata(a),
pNode(pN),
pIndVel(pR),
fPassiveInducedVelocity(0),
f(fTmp),
dHalfSpan(dS/2.),
Ra(RaTmp),
Ra3(RaTmp.GetVec(3)),
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
GDI(iN),
pdOuta(NULL),
pvdOuta(NULL),
F(0.),
M(0.),
bJacobian(bUseJacobian)
{
	DEBUGCOUTFNAME("AerodynamicBody::AerodynamicBody");

	ASSERT(pNode != NULL);
	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
	ASSERT(aerodata != NULL);

#ifdef DEBUG
	if (pIndVel != NULL) {
		ASSERT(pIndVel->GetElemType() == Elem::INDUCEDVELOCITY);
	}
#endif /* DEBUG */

	SAFENEWARR(pdOuta, doublereal, 20*iN);
	SAFENEWARR(pvdOuta, doublereal*, iN);
	for (integer i = 20*iN; i-- > 0; ) {
		pdOuta[i] = 0.;
	}
	for (integer i = iN; i-- > 0; ) {
		pvdOuta[i] = pdOuta+20*i;
	}
}

AerodynamicBody::~AerodynamicBody(void)
{
	DEBUGCOUTFNAME("AerodynamicBody::~AerodynamicBody");

	SAFEDELETEARR(pvdOuta);
	SAFEDELETEARR(pdOuta);

	SAFEDELETE(aerodata);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
void
AerodynamicBody::SetOutputFlag(flag f)
{
	DEBUGCOUTFNAME("AerodynamicBody::SetOutputFlag");
	ToBeOutput::SetOutputFlag(f);
	AerodynamicOutput::SetOutputFlag(f, GDI.iGetNum());
}

/* Dimensioni del workspace */
void
AerodynamicBody::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = *piNumCols = 6 + iGetNumDof();
}

/* inherited from SimulationEntity */
unsigned int
AerodynamicBody::iGetNumDof(void) const
{
	return aerodata->iGetNumDof()*GDI.iGetNum();
}

std::ostream&
AerodynamicBody::DescribeDof(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(bInitial == false);

	integer iFirstIndex = iGetFirstIndex();
	integer iNumDof = aerodata->iGetNumDof();

	for (integer i = 0; i < GDI.iGetNum(); i++) {
		out 
			<< prefix << iFirstIndex + 1 << "->" << iFirstIndex + iNumDof
			<< ": internal states at point #" << i << " of " << GDI.iGetNum() << std::endl;

		iFirstIndex += iNumDof;
	}

	return out;
}

void
AerodynamicBody::DescribeDof(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	if (i < -1) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(bInitial == false);

	integer iNumDof = aerodata->iGetNumDof();

	std::ostringstream os;
	os << "AerodynamicBody(" << GetLabel() << ")";

	if (i == -1) {
		integer iTotDof = iNumDof*GDI.iGetNum();
		desc.resize(iTotDof);

		std::string name(os.str());

		for (integer s = 0; s < iTotDof; s++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": internal state #" << s % iNumDof << " of " << iNumDof
				<< " at point #" << s/iNumDof << " of " << GDI.iGetNum();
			desc[s] = os.str();
		}

		return;
	}

	desc.resize(1);
	os << ": internal state #" << i % iNumDof << " of " << iNumDof
		<< " at point #" << i/iNumDof << " of " << GDI.iGetNum();
	desc[0] = os.str();
}

std::ostream&
AerodynamicBody::DescribeEq(std::ostream& out,
	const char *prefix, bool bInitial) const
{
	ASSERT(bInitial == false);

	integer iFirstIndex = iGetFirstIndex();
	integer iNumDof = aerodata->iGetNumDof();

	for (integer i = 0; i < GDI.iGetNum(); i++) {
		out 
			<< prefix << iFirstIndex + 1 << "->" << iFirstIndex + iNumDof
			<< ": dynamics equations at point #" << i << " of " << GDI.iGetNum() << std::endl;

		iFirstIndex += iNumDof;
	}

	return out;
}

void
AerodynamicBody::DescribeEq(std::vector<std::string>& desc,
	bool bInitial, int i) const
{
	if (i < -1) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ASSERT(bInitial == false);

	integer iNumDof = aerodata->iGetNumDof();

	std::ostringstream os;
	os << "AerodynamicBody(" << GetLabel() << ")";

	if (i == -1) {
		integer iTotDof = iNumDof*GDI.iGetNum();
		desc.resize(iTotDof);

		std::string name(os.str());

		for (integer s = 0; s < iTotDof; s++) {
			os.str(name);
			os.seekp(0, std::ios_base::end);
			os << ": internal equation #" << s % iNumDof << " of " << iNumDof
				<< " at point #" << s/iNumDof << " of " << GDI.iGetNum();
			desc[s] = os.str();
		}

		return;
	}

	desc.resize(1);
	os << ": internal equation #" << i % iNumDof << " of " << iNumDof
		<< " at point #" << i/iNumDof << " of " << GDI.iGetNum();
	desc[0] = os.str();
}

DofOrder::Order
AerodynamicBody::GetDofType(unsigned int i) const
{
	ASSERT(aerodata->iGetNumDof() > 0);

	return aerodata->GetDofType(i % aerodata->iGetNumDof());
}

void
AerodynamicBody::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints* h)
{
	integer iNumDof = aerodata->iGetNumDof();
	if (iNumDof) {
		vx.Resize(6);
		wx.Resize(6);
		fq.Resize(iNumDof);
		cq.Resize(iNumDof);

		vx.Reset();
		wx.Reset();
		fq.Reset();
		cq.Reset();
	}
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBody::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBody::Restart");

	out << "  aerodynamic body: " << GetLabel() << ", "
		<< pNode->GetLabel();
	if (pIndVel != NULL) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f.Write(out, ", ")
		<< ", " << dHalfSpan*2.
		<< ", reference, node, 1, ", (Ra.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< ", " << GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

VariableSubMatrixHandler&
AerodynamicBody::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUT("Entering AerodynamicBody::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici delle varie incognite */
	integer iNodeFirstMomIndex = pNode->iGetFirstMomentumIndex();
	integer iNodeFirstPosIndex = pNode->iGetFirstPositionIndex();

	/* Setta gli indici delle equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNodeFirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNodeFirstPosIndex + iCnt);
	}

	/* Equations start here... */
	doublereal dW[6];

	/* Dati del nodo */
	const Vec3& Xn(pNode->GetXCurr());
	const Mat3x3& Rn(pNode->GetRCurr());
	const Vec3& Vn(pNode->GetVCurr());
	const Vec3& Wn(pNode->GetWCurr());

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR(Rn*Ra);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso la posizione del nodo perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	doublereal** pvd = pvdOuta;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 6;

		integer iOffset = iFirstEq - 6;

		for (int iCnt = 6 + 1; iCnt <= iNumRows; iCnt++) {
			WM.PutRowIndex(iCnt, iOffset + iCnt);
			WM.PutColIndex(iCnt, iOffset + iCnt);
		}
	}

	/* Ciclo sui punti di Gauss */
	PntWght PW = GDI.GetFirst();
	int iPnt = 0;
	do {
		doublereal dCsi = PW.dGetPnt();
		Vec3 Xr(Rn*(f + Ra3*(dHalfSpan*dCsi)));
		Vec3 Xnr = Xn + Xr;
		Vec3 Vr(Vn + Wn.Cross(Xr));

		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xnr)) {
			Vr -= VTmp;
		}

		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pIndVel != NULL) {
			Vr += pIndVel->GetInducedVelocity(Xnr);
		}

		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(dCsi) + dGet();
		aerodata->SetSectionData(dCsi,
			Chord.dGet(dCsi),
			ForcePoint.dGet(dCsi),
			VelocityPoint.dGet(dCsi),
			dTw,
			dOmega);

		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		Mat3x3 RRloc;
		if (dTw != 0.) {
			doublereal dCosT = cos(dTw);
			doublereal dSinT = sin(dTw);
			/* Assumo lo svergolamento positivo a cabrare */
			Mat3x3 RTw( dCosT, dSinT, 0.,
				-dSinT, dCosT, 0.,
				0.,    0.,    1.);
			RRloc = RR*RTw;

		} else {
			RRloc = RR;
		}

		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		VTmp = RRloc.MulTV(Vr);
		VTmp.PutTo(dW);

		Vec3 WTmp = RRloc.MulTV(Wn);
		WTmp.PutTo(&dW[3]);

		/* Funzione di calcolo delle forze aerodinamiche */
		doublereal Fa0[6];
		Mat6x6 JFa;

		/* Jacobian Assembly... */
		doublereal cc = dHalfSpan*PW.dGetWght();

		if (iNumDof) {
			// prepare (v/dot{x} + dCoef*v/x) and so
			Mat3x3 RRlocT(RRloc.Transpose());

			vx.PutMat3x3(1, RRlocT);
			vx.PutMat3x3(4, RRloc.MulTM(Mat3x3((Vr + Xr.Cross(Wn))*dCoef - Xr)));
			wx.PutMat3x3(4, RRlocT);

			// equations from iFirstEq on are dealt with by aerodata
			aerodata->AssJac(WM, dCoef, XCurr, XPrimeCurr,
                		iFirstEq, iFirstSubEq,
				vx, wx, fq, cq, iPnt, dW, Fa0, JFa, *pvd);

			// deal with (f/dot{q} + dCoef*f/q) and so
			integer iOffset = 6 + iPnt*iNumDof;
			for (integer iCol = 1; iCol <= fq.iGetNumCols(); iCol++) {
				Vec3 fqTmp((RRloc*fq.GetVec(iCol))*cc);
				Vec3 cqTmp((RRloc*cq.GetVec(iCol))*cc);

				WM.Sub(1, iOffset + iCol, fqTmp);
				WM.Sub(4, iOffset + iCol, Xr.Cross(fqTmp) + cqTmp);
			}

			// first equation
			iFirstEq += iNumDof;
			iFirstSubEq += iNumDof;

		} else {
			aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, *pvd);
		}

		/* OUTA */
		pvd++;

		// rotate force, couple and Jacobian matrix in absolute frame
		Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

		Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
		Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + Xr.Cross(fTmp));

		// compute submatrices (see tecman.pdf)
		Mat3x3 Mat21(Xr.Cross(JFaR.GetMat11()) + JFaR.GetMat21());

		Mat3x3 Mat12(JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(Xr));

		Mat3x3 Mat22(Xr.Cross(JFaR.GetMat12()) + JFaR.GetMat22() - Mat21*Mat3x3(Xr));

		Mat3x3 MatV((Vr + Xr.Cross(Wn))*dCoef);

		Mat12 += JFaR.GetMat11()*MatV - Mat3x3(fTmp);

		Mat22 += Mat21*MatV - Mat3x3(cTmp);

		// add (actually, sub) contributions, scaled by weight
		WM.Sub(1, 1, JFaR.GetMat11());
		WM.Sub(1, 4, Mat12);
		WM.Sub(4, 1, Mat21);
		WM.Sub(4, 4, Mat22);

		iPnt++;

	} while (GDI.fGetNext(PW));

	return WorkMat;
}


SubVectorHandler&
AerodynamicBody::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::AssRes");

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstIndex = pNode->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	AssVec(WorkVec, dCoef, XCurr, XPrimeCurr);

	return WorkVec;
}


SubVectorHandler&
AerodynamicBody::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::InitialAssRes");

	if (aerodata->iGetNumDof() > 0) {
		WorkVec.ResizeReset(0);
		return WorkVec;
	}

	integer iNumRows;
	integer iNumCols;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	integer iFirstIndex = pNode->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstIndex + iCnt);
	}

	AssVec(WorkVec, 1., XCurr, XCurr);

	return WorkVec;
}


/* assemblaggio residuo */
void
AerodynamicBody::AssVec(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("AerodynamicBody::AssVec");

	doublereal dTng[6];
	doublereal dW[6];

	/* Dati del nodo */
	const Vec3& Xn(pNode->GetXCurr());
	const Mat3x3& Rn(pNode->GetRCurr());
	const Vec3& Vn(pNode->GetVCurr());
	const Vec3& Wn(pNode->GetWCurr());

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR(Rn*Ra);

	/*
	 * Se l'elemento e' collegato ad un rotore,
 	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/* Resetta i dati */
	F.Reset();
	M.Reset();

	/*
	 * Dati "permanenti" (uso la posizione del nodo perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	doublereal** pvd = pvdOuta;

	ResetIterator();

	integer iNumDof = aerodata->iGetNumDof();
	integer iFirstEq = -1;
	integer iFirstSubEq = -1;
	if (iNumDof > 0) {
		iFirstEq = iGetFirstIndex();
		iFirstSubEq = 6;

		integer iOffset = iFirstEq - 6;
		integer iNumRows = WorkVec.iGetSize();
		for (int iCnt = 6 + 1; iCnt <= iNumRows; iCnt++) {
			WorkVec.PutRowIndex(iCnt, iOffset + iCnt);
		}
	}

	/* Ciclo sui punti di Gauss */
	PntWght PW = GDI.GetFirst();
	int iPnt = 0;
	do {
		doublereal dCsi = PW.dGetPnt();
		Vec3 Xr(Rn*(f + Ra3*(dHalfSpan*dCsi)));
		Vec3 Xnr = Xn + Xr;
		Vec3 Vr(Vn + Wn.Cross(Xr));

		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xnr)) {
	 		Vr -= VTmp;
		}

		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pIndVel != NULL) {
	 		Vr += pIndVel->GetInducedVelocity(Xnr);
		}

		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(dCsi) + dGet();
		aerodata->SetSectionData(dCsi,
			Chord.dGet(dCsi),
			ForcePoint.dGet(dCsi),
			VelocityPoint.dGet(dCsi),
			dTw,
			dOmega);

		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		Mat3x3 RRloc;
		if (dTw != 0.) {
			doublereal dCosT = cos(dTw);
			doublereal dSinT = sin(dTw);
			/* Assumo lo svergolamento positivo a cabrare */
			Mat3x3 RTw( dCosT, dSinT, 0.,
				   -dSinT, dCosT, 0.,
				    0.,    0.,    1.);
			RRloc = RR*RTw;

		} else {
			RRloc = RR;
		}

		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		VTmp = RRloc.MulTV(Vr);
		VTmp.PutTo(dW);

		Vec3 WTmp = RRloc.MulTV(Wn);
		WTmp.PutTo(&dW[3]);

		/* Funzione di calcolo delle forze aerodinamiche */
		if (iNumDof) {
			aerodata->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr,
                		iFirstEq, iFirstSubEq, iPnt, dW, dTng, *pvd);

			// first equation
			iFirstEq += iNumDof;
			iFirstSubEq += iNumDof;

		} else {
			aerodata->GetForces(iPnt, dW, dTng, *pvd);
		}

		/* OUTA */
		pvd++;

		// specific for Gauss points force output
		if (fToBeOutput() && IsPGAUSS()) {
	 		SetData(VTmp, dTng);
		}

		/* Dimensionalizza le forze */
		doublereal dWght = PW.dGetWght();
		Vec3 FTmp(RRloc*(Vec3(dTng)*(dHalfSpan*dWght)));
		F += FTmp;
		M += RRloc*(Vec3(&dTng[3])*(dHalfSpan*dWght));
		M += Xr.Cross(FTmp);

		iPnt++;

	} while (GDI.fGetNext(PW));

	/* Se e' definito il rotore, aggiungere il contributo alla trazione */
	if (pIndVel != NULL && !fPassiveInducedVelocity) {
		pIndVel->AddForce(GetLabel(), F, M, Xn);
	}

	/* Sommare il termine al residuo */
	WorkVec.Add(1, F);
	WorkVec.Add(4, M);
}

void
AerodynamicBody::AfterConvergence(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	/* Memoria in caso di forze instazionarie */
	switch (aerodata->Unsteady()) {
	case AeroData::STEADY:
		break;

	case AeroData::HARRIS:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case AeroData::BIELAWA:
		for (integer i = 0; i < GDI.iGetNum(); i++) {
	 		aerodata->Update(i);
		}
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBody::Output(OutputHandler& OH) const
{
	/* Output delle forze aerodinamiche F, M su apposito file */
	if (fToBeOutput()) {
		std::ostream& out = OH.Aerodynamic()
			<< std::setw(8) << GetLabel();

		switch (GetOutput()) {
		case AEROD_OUT_NODE:
			out << " " << std::setw(8) << pNode->GetLabel()
				<< " ", F.Write(out, " ") << " ", M.Write(out, " ");
			break;

		case AEROD_OUT_PGAUSS:
			ASSERT(!pOutput.empty());
			for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
				i != pOutput.end(); i++)
			{
	 			out << " " << i->alpha
					<< " " << i->f;
			}
			break;

		case AEROD_OUT_STD:
			for (int i = 0; i < GDI.iGetNum(); i++) {
	 			for (int j = 1; j <= 6; j++) {
	 				out << " " << pvdOuta[i][j];
	 			}
				out << " " << pvdOuta[i][OUTA_ALF1]
					<< " " << pvdOuta[i][OUTA_ALF2];
			}
			break;

		default:
			ASSERT(0);
			break;
		}

		out << std::endl;
	}
}

/* AerodynamicBody - end */


/* Legge dati aerodinamici */

static AeroData::UnsteadyModel
ReadUnsteadyFlag(MBDynParser& HP)
{
	if (HP.IsArg()) {
		AeroData::UnsteadyModel eInst = AeroData::STEADY;
		if (HP.IsKeyWord("unsteady")) {
			/*
			 * swallow "unsteady" keyword
			 */
			if (HP.IsKeyWord("steady")) {
				eInst = AeroData::STEADY;
			} else if (HP.IsKeyWord("harris")) {
				eInst = AeroData::HARRIS;
			} else if (HP.IsKeyWord("bielawa")) {
				eInst = AeroData::BIELAWA;
			} else {
				/* demote to pedantic, because the integer
				 * form allows to change unsteady model
				 * parametrically (while waiting for string
				 * vars) */
				pedantic_cerr("deprecated unsteady model "
					"given by integer number;"
					" use \"steady\", \"Harris\" or \"Bielawa\" "
					"instead, at line " << HP.GetLineData()
					<< std::endl);

				int i = HP.GetInt();
				if (i < AeroData::STEADY || i >= AeroData::LAST) {
					silent_cerr("illegal unsteady flag "
							"numeric value " << i
							<< " at line "
							<< HP.GetLineData()
							<< std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
				eInst = AeroData::UnsteadyModel(i);
			}
		}

		switch (eInst) {
		case AeroData::STEADY:
		case AeroData::BIELAWA:
			break;

		case AeroData::HARRIS:
			silent_cerr("\"Harris\" unsteady aerodynamics "
					"are not available at line "
					<< HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		default:
	 		silent_cerr("illegal unsteady flag at line "
					<< HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		/*
		 * unsteady flag
		 */
		return eInst;
	}

	/*
	 * default: no unsteady ...
	 */
	return AeroData::STEADY;
}

static void
ReadC81MultipleAeroData(DataManager* pDM, MBDynParser& HP, AeroData** aerodata)
{
	integer nProfiles = HP.GetInt();
	if (nProfiles <= 0) {
		silent_cerr("Need at least one profile at line "
				<< HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	integer *profiles = NULL;
	SAFENEWARR(profiles, integer, nProfiles+1);
	doublereal *upper_bounds = NULL;
	SAFENEWARR(upper_bounds, doublereal, nProfiles);
	const c81_data** data = NULL;
	SAFENEWARR(data, const c81_data*, nProfiles+1);

	for (int i = 0; i < nProfiles; i++) {
		profiles[i] = HP.GetInt();
		upper_bounds[i] = HP.GetReal();
		if (upper_bounds[i] <= -1.) {
			silent_cerr("upper bound " << i+1 << " = "
					<< upper_bounds[i]
					<< " too small at line "
					<< HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (upper_bounds[i] > 1.) {
			silent_cerr("upper bound " << i+1 << " = "
					<< upper_bounds[i]
					<< " too large at line "
					<< HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (i > 0 && upper_bounds[i] <= upper_bounds[i-1]) {
			silent_cerr("upper bound " << i+1 << " = "
					<< upper_bounds[i]
					<< " not in increasing order "
					"at line " << HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		data[i] = HP.GetC81Data(profiles[i]);
		if (data[i] == NULL) {
			silent_cerr("Unable to find profile "
					<< profiles[i] << " at line "
					<< HP.GetLineData()
					<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
 		DEBUGLCOUT(MYDEBUG_INPUT, "profile data " << i+1
				<< " is from file c81 " << profiles[i]
				<< std::endl);
	}

	if (upper_bounds[nProfiles-1] != 1.) {
		silent_cerr("warning: the last upper bound should be 1.0 "
				"at line " << HP.GetLineData()
				<< std::endl);
	}

	profiles[nProfiles] = -1;
	data[nProfiles] = NULL;

	AeroData::UnsteadyModel
		eInst = ReadUnsteadyFlag(HP);
	DriveCaller *ptime = NULL;
	if (eInst != AeroData::STEADY) {
		SAFENEWWITHCONSTRUCTOR(ptime,
				TimeDriveCaller,
				TimeDriveCaller(pDM->pGetDrvHdl()));
	}
	SAFENEWWITHCONSTRUCTOR(*aerodata,
			C81MultipleAeroData,
			C81MultipleAeroData(eInst,
				nProfiles, profiles,
				upper_bounds, data,
				ptime));
}

static void
ReadC81InterpolatedAeroData(DataManager* pDM, MBDynParser& HP, AeroData** aerodata)
{
	silent_cerr("C81InterpolatedAeroData not implemented yet" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

static void
ReadAeroData(DataManager* pDM, MBDynParser& HP,
	Shape** ppChord, Shape** ppForce,
	Shape** ppVelocity, Shape** ppTwist,
	integer* piNumber, DriveCaller** ppDC,
	AeroData** aerodata)
{
	DEBUGCOUTFNAME("ReadAeroData");

	/* Keywords */
	const char* sKeyWords[] = {
		"naca0012",
		"rae9671",
		"c81",

		"umd",		// FIXME

		NULL
	};

	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,
		NACA0012 = 0,
		RAE9671,
		C81,

		UMD,

		LASTKEYWORD
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	*ppChord = ReadShape(HP);
	*ppForce = ReadShape(HP);
	*ppVelocity = ReadShape(HP);
	*ppTwist = ReadShape(HP);

	*piNumber = HP.GetInt();
	if (*piNumber <= 0) {
		silent_cerr("need at least 1 Gauss integration point at line "
				<< HP.GetLineData()  << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	DEBUGLCOUT(MYDEBUG_INPUT, "Gauss points number: "
			<< *piNumber << std::endl);

	if (HP.IsKeyWord("control")) {
		/* Driver di un'eventuale controllo */
		*ppDC = HP.GetDriveCaller();

	} else {
		SAFENEW(*ppDC, NullDriveCaller);
	}

	if (HP.IsArg()) {
		switch (HP.IsKeyWord()) {
 		default:
	  		silent_cerr("unknown profile type \"" << HP.GetString()
					<< "\" at line " << HP.GetLineData()
					<< "; using default (NACA0012)"
					<< std::endl);

 		case NACA0012: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				   "profile is NACA0012" << std::endl);

			AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
			DriveCaller *ptime = NULL;
			if (eInst != AeroData::STEADY) {
				SAFENEWWITHCONSTRUCTOR(ptime, TimeDriveCaller,
						TimeDriveCaller(pDM->pGetDrvHdl()));
			}
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(eInst, 1, ptime));
	  		break;
 		}

 		case RAE9671: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				"profile is RAE9671" << std::endl);

			AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
			DriveCaller *ptime = NULL;
			if (eInst != AeroData::STEADY) {
				SAFENEWWITHCONSTRUCTOR(ptime, TimeDriveCaller,
						TimeDriveCaller(pDM->pGetDrvHdl()));
			}
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(eInst, 2, ptime));
	  		break;
 		}

	 	/*
		 * uso tabelle standard, che vengono cercate in base all'indice
		 * (sistemare)
		 */
 		case C81:
			if (HP.IsKeyWord("multiple")) {
				ReadC81MultipleAeroData(pDM, HP, aerodata);

			} else if (HP.IsKeyWord("interpolated")) {
				ReadC81InterpolatedAeroData(pDM, HP, aerodata);

			} else {
	  			integer iProfile = HP.GetInt();
		  		const c81_data* data = HP.GetC81Data(iProfile);

		  		DEBUGLCOUT(MYDEBUG_INPUT,
					"profile data is from file c81 "
					<< iProfile << std::endl);
				AeroData::UnsteadyModel
					eInst = ReadUnsteadyFlag(HP);
				DriveCaller *ptime = NULL;
				if (eInst != AeroData::STEADY) {
					SAFENEWWITHCONSTRUCTOR(ptime,
							TimeDriveCaller,
							TimeDriveCaller(pDM->pGetDrvHdl()));
				}
	  			SAFENEWWITHCONSTRUCTOR(*aerodata,
					C81AeroData,
					C81AeroData(eInst, iProfile,
						data, ptime));
			}
			break;

		case UMD:
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				UMDAeroData,
				UMDAeroData(0));
			break;
			
		}

	} else {
		/* FIXME: better abort! */
		silent_cerr("missing profile type at line "
				<< HP.GetLineData()
				<< "; using default (NACA0012)"
				<< std::endl);

		AeroData::UnsteadyModel eInst = ReadUnsteadyFlag(HP);
		SAFENEWWITHCONSTRUCTOR(*aerodata,
			STAHRAeroData, STAHRAeroData(eInst, 1));
	}
}

static InducedVelocity*
ReadInducedVelocity(DataManager *pDM, MBDynParser& HP,
	unsigned uLabel, const char *sElemType)
{
	bool bReadIV(false);
	InducedVelocity *pIndVel = 0;
	if (HP.IsKeyWord("rotor")) {
		silent_cerr(sElemType << "(" << uLabel << "): "
			"\"rotor\" keyword is deprecated; "
			"use \"induced velocity\" instead "
			"at line " << HP.GetLineData()
			<< std::endl);

		bReadIV = true;

	} else if (HP.IsKeyWord("induced" "velocity")) {
		bReadIV = true;
	}

	if (bReadIV) {
		unsigned int uIV = (unsigned int)HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Linked to InducedVelocity(" << uIV << ")" << std::endl);

		/*
		 * verifica di esistenza del rotore
		 * NOTA: ovviamente il rotore deve essere definito
		 * prima dell'elemento aerodinamico
		 */
		Elem* p = pDM->pFindElem(Elem::INDUCEDVELOCITY, uIV);
		if (p == NULL) {
			silent_cerr(sElemType << "(" << uLabel << "): "
				"InducedVelocity(" << uIV << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		pIndVel = dynamic_cast<InducedVelocity *>(p);
		ASSERT(pIndVel != 0);
	}

	return pIndVel;
}

Elem *
ReadAerodynamicBody(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBody");

	/* Nodo */
	StructNode* pNode = dynamic_cast<StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));

	InducedVelocity* pIndVel = ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBody");

	ReferenceFrame RF(pNode);
	Vec3 f(HP.GetPosRel(RF));

	DEBUGLCOUT(MYDEBUG_INPUT, "Offset: " << f << std::endl);

	Mat3x3 Ra(HP.GetRotRel(RF));

	doublereal dSpan = HP.GetReal();
	DEBUGLCOUT(MYDEBUG_INPUT, "Span: " << dSpan << std::endl);

	Shape* pChord = NULL;
	Shape* pForce = NULL;
	Shape* pVelocity = NULL;
	Shape* pTwist = NULL;

	integer iNumber = 0;
	DriveCaller* pDC = NULL;
	AeroData* aerodata = NULL;

	ReadAeroData(pDM, HP,
		&pChord, &pForce, &pVelocity, &pTwist,
		     &iNumber, &pDC, &aerodata);

	aerodata->SetNumPoints(iNumber);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = true;
	}

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBody(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = NULL;
	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBody,
		AerodynamicBody(uLabel, pDO, pNode, pIndVel, f, dSpan, Ra,
				pChord, pForce, pVelocity, pTwist,
				iNumber, aerodata, pDC, bUseJacobian, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Vec3 Ra3 = Ra.GetVec(3);
	doublereal dCm1 = pChord->dGet(-1.);
	doublereal dPm1 = pForce->dGet(-1.);
	doublereal dTm1 = pTwist->dGet(-1.);
	Vec3 Ram1 = Ra*Vec3(std::cos(dTm1), std::sin(dTm1), 0.);

	doublereal dCp1 = pChord->dGet(1.);
	doublereal dPp1 = pForce->dGet(1.);
	doublereal dTp1 = pTwist->dGet(1.);
	Vec3 Rap1 = Ra*Vec3(std::cos(dTp1), std::sin(dTp1), 0.);

	std::ostream& out = pDM->GetLogFile();
	out << "aero0: " << uLabel
		<< " " << pNode->GetLabel()
		<< " ", (f - Ra3*(dSpan/2.) + Ram1*(dPm1 - dCm1*3./4.)).Write(out, " ")
		<< " ", (f - Ra3*(dSpan/2.) + Ram1*(dPm1 + dCm1/4.)).Write(out, " ")
		<< " ", (f + Ra3*(dSpan/2.) + Rap1*(dPp1 - dCp1*3./4.)).Write(out, " ")
		<< " ", (f + Ra3*(dSpan/2.) + Rap1*(dPp1 + dCp1/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of DataManager::ReadAerodynamicBody() */


/* AerodynamicBeam - begin */

AerodynamicBeam::AerodynamicBeam(unsigned int uLabel,
	const DofOwner *pDO,
	const Beam* pB, InducedVelocity* pR,
	const Vec3& fTmp1,
	const Vec3& fTmp2,
	const Vec3& fTmp3,
	const Mat3x3& Ra1Tmp,
	const Mat3x3& Ra2Tmp,
	const Mat3x3& Ra3Tmp,
	const Shape* pC, const Shape* pF,
	const Shape* pV, const Shape* pT,
	integer iN, AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	flag fOut)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut),
DriveOwner(pDC),
AerodynamicOutput(fOut, 3*iN),
aerodata(a),
pBeam(pB),
pIndVel(pR),
fPassiveInducedVelocity(0),
f1(fTmp1),
f2(fTmp2),
f3(fTmp3),
Ra1(Ra1Tmp),
Ra2(Ra2Tmp),
Ra3(Ra3Tmp),
Ra1_3(Ra1Tmp.GetVec(3)),
Ra2_3(Ra2Tmp.GetVec(3)),
Ra3_3(Ra3Tmp.GetVec(3)),
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
GDI(iN),
pdOuta(NULL),
pvdOuta(NULL),
bJacobian(bUseJacobian)
{
	DEBUGCOUTFNAME("AerodynamicBeam::AerodynamicBeam");

	ASSERT(pBeam != NULL);
	ASSERT(pBeam->GetElemType() == Elem::BEAM);

	pNode1 = pBeam->pGetNode(1);
	pNode2 = pBeam->pGetNode(2);
	pNode3 = pBeam->pGetNode(3);

	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode3 != NULL);
	ASSERT(pNode3->GetNodeType() == Node::STRUCTURAL);

#ifdef DEBUG
	if (pIndVel != NULL) {
		ASSERT(pIndVel->GetElemType() == Elem::INDUCEDVELOCITY);
	}
#endif /* DEBUG */

	SAFENEWARR(pdOuta, doublereal, 20*3*iN);
	SAFENEWARR(pvdOuta, doublereal*, 3*iN);
	for (integer i = 20*3*iN; i-- > 0; ) {
		pdOuta[i] = 0.;
	}
	for (integer i = 3*iN; i-- > 0; ) {
		pvdOuta[i] = pdOuta+20*i;
	}
}

AerodynamicBeam::~AerodynamicBeam(void)
{
	DEBUGCOUTFNAME("AerodynamicBeam::~AerodynamicBeam");

	SAFEDELETEARR(pvdOuta);
	SAFEDELETEARR(pdOuta);

	SAFEDELETE(aerodata);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
void
AerodynamicBeam::SetOutputFlag(flag f)
{
	DEBUGCOUTFNAME("AerodynamicBeam::SetOutputFlag");

	ToBeOutput::SetOutputFlag(f);
	AerodynamicOutput::SetOutputFlag(f, 3*GDI.iGetNum());
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBeam::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Restart");
	out << "  aerodynamic beam: " << GetLabel()
		<< ", " << pBeam->GetLabel();
	if (pIndVel != NULL) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f1.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra1.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra1.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f2.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra2.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra2.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f3.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra3.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra3.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

static const doublereal d13 = 1./sqrt(3.);
static const doublereal pdsi3[] = { -1., -d13, d13 };
static const doublereal pdsf3[] = { -d13, d13, 1. };

/* Jacobian assembly */
VariableSubMatrixHandler&
AerodynamicBeam::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal  dCoef  ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AerodynamicBeam::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstMomentumIndex();
	integer iNode3FirstPosIndex = pNode3->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
		WM.PutColIndex(12 + iCnt, iNode3FirstPosIndex + iCnt);
	}

	doublereal dW[6];

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[3];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Xn[NODE3] = pNode3->GetXCurr();
	const Mat3x3& Rn3(pNode3->GetRCurr());
	const Vec3& Vn3(pNode3->GetVCurr());
	const Vec3& Wn3(pNode3->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);
	Vec3 f3Tmp(Rn3*f3);

	Vec3 fTmp[3];
	fTmp[NODE1] = f1Tmp;
	fTmp[NODE2] = f2Tmp;
	fTmp[NODE3] = f3Tmp;

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);
	Vec3 X3Tmp(Xn[NODE3] + f3Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));
	Vec3 V3Tmp(Vn3 + Wn3.Cross(f3Tmp));

	Vec3 Omega1Crossf1(Wn1.Cross(f1Tmp));
	Vec3 Omega2Crossf2(Wn2.Cross(f2Tmp));
	Vec3 Omega3Crossf3(Wn3.Cross(f3Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);
	Mat3x3 RR3(Rn3*Ra3);

	/*
	 * Parametri di rotazione dai nodi 1 e 3 al nodo 2 (nell'ipotesi
	 * che tale trasformazione non dia luogo ad una singolarita')
	 */

	Vec3 g1(ER_Rot::Param, RR2.MulTM(RR1));
	Vec3 g3(ER_Rot::Param, RR2.MulTM(RR3));

	Mat3x3 GammaInv1(ER_Rot::MatGm1, g1);
	Mat3x3 GammaInv3(ER_Rot::MatGm1, g3);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn[NODE2], rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	/* OUTA */
	doublereal** pvd = pvdOuta;
	int iPnt = 0;

	ResetIterator();

	for (int iNode = 0; iNode < LASTNODE; iNode++) {
		doublereal dsi = pdsi3[iNode];
		doublereal dsf = pdsf3[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		Mat3x3 WM_F[6], WM_M[6];
		WM_F[DELTAx1].Reset();
		WM_F[DELTAg1].Reset();
		WM_F[DELTAx2].Reset();
		WM_F[DELTAg2].Reset();
		WM_F[DELTAx3].Reset();
		WM_F[DELTAg3].Reset();
		WM_M[DELTAx1].Reset();
		WM_M[DELTAg1].Reset();
		WM_M[DELTAx2].Reset();
		WM_M[DELTAg2].Reset();
		WM_M[DELTAx3].Reset();
		WM_M[DELTAg3].Reset();

		//unsigned int iDelta_x1, iDelta_g1, iDelta_x2, iDelta_g2, iDelta_x3, iDelta_g3;
		//iDelta_x1 = 0;, iDelta_g1, iDelta_x2, iDelta_g2, iDelta_x3, iDelta_g3;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi3N(ds,
				Xn[NODE1], Xn[NODE2], Xn[NODE3]);

			doublereal dN1 = ShapeFunc3N(ds, 1);
			doublereal dN2 = ShapeFunc3N(ds, 2);
			doublereal dN3 = ShapeFunc3N(ds, 3);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2 + X3Tmp*dN3);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2 + V3Tmp*dN3);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2 + Wn3*dN3);
			Vec3 gr(g1*dN1 + g3*dN3);
			Mat3x3 Gamma(ER_Rot::MatG, gr);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != NULL) {
				Vr += pIndVel->GetInducedVelocity(Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR2*Mat3x3(ER_Rot::MatR, gr));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);
			/* Funzione di calcolo delle forze aerodinamiche */
			doublereal Fa0[6];
			Mat6x6 JFa;

			aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, *pvd);

			/* OUTA */
			pvd++;

			doublereal cc = dXds*dsdCsi*PW.dGetWght();

			// rotate force, couple and Jacobian matrix in absolute frame
			Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

			// force and moment about the node
			Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
			Vec3 d(Xr - Xn[iNode]);
			Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + d.Cross(fTmp));

			Mat3x3 Theta1(RR2*Gamma*GammaInv1.MulMT(RR2*dN1));
			Mat3x3 Theta3(RR2*Gamma*GammaInv3.MulMT(RR2*dN3));
			Mat3x3 Theta2(Eye3 - Theta1 - Theta3);

			Vec3 Vrc(Vr*dCoef);
			Mat3x3 Bv1(Vrc.Cross(Theta1) - Mat3x3(Omega1Crossf1*(dN1*dCoef)));
			Mat3x3 Bv2(Vrc.Cross(Theta2) - Mat3x3(Omega2Crossf2*(dN2*dCoef)));
			Mat3x3 Bv3(Vrc.Cross(Theta3) - Mat3x3(Omega3Crossf3*(dN3*dCoef)));

			Vec3 Wrc(Wr*dCoef);
			Mat3x3 Bw1(Wrc.Cross(Theta1) - Mat3x3(Wn1*(dN1*dCoef)));
			Mat3x3 Bw2(Wrc.Cross(Theta2) - Mat3x3(Wn2*(dN2*dCoef)));
			Mat3x3 Bw3(Wrc.Cross(Theta3) - Mat3x3(Wn3*(dN3*dCoef)));

			Mat3x3 WM_F2[6];

			// f <-> x
			WM_F2[DELTAx1] = JFaR.GetMat11()*dN1;

			WM_F2[DELTAx2] = JFaR.GetMat11()*dN2;

			WM_F2[DELTAx3] = JFaR.GetMat11()*dN3;

			doublereal delta;

			// c <-> x
			delta = (iNode == NODE1) ? 1. : 0.;
			WM_M[DELTAx1] += JFaR.GetMat21()*dN1 - Mat3x3(fTmp*(dN1 - delta));

			delta = (iNode == NODE2) ? 1. : 0.;
			WM_M[DELTAx2] += JFaR.GetMat21()*dN2 - Mat3x3(fTmp*(dN2 - delta));

			delta = (iNode == NODE3) ? 1. : 0.;
			WM_M[DELTAx3] += JFaR.GetMat21()*dN3 - Mat3x3(fTmp*(dN3 - delta));

			// f <-> g
			WM_F2[DELTAg1] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f1Tmp))*dN1;
			WM_F2[DELTAg1] += JFaR.GetMat11()*Bv1 + JFaR.GetMat12()*Bw1;
			WM_F2[DELTAg1] -= fTmp.Cross(Theta1);
		
			WM_F2[DELTAg2] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f2Tmp))*dN2;
			WM_F2[DELTAg2] += JFaR.GetMat11()*Bv2 + JFaR.GetMat12()*Bw2;
			WM_F2[DELTAg2] -= fTmp.Cross(Theta2);
		
			WM_F2[DELTAg3] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f3Tmp))*dN3;
			WM_F2[DELTAg3] += JFaR.GetMat11()*Bv3 + JFaR.GetMat12()*Bw3;
			WM_F2[DELTAg3] -= fTmp.Cross(Theta3);

			// c <-> g
			WM_M[DELTAg1] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f1Tmp))*dN1;
			WM_M[DELTAg1] += JFaR.GetMat21()*Bv1 + JFaR.GetMat22()*Bw1;
			WM_M[DELTAg1] -= cTmp.Cross(Theta1);
			WM_M[DELTAg1] += Mat3x3(fTmp, f1Tmp*dN1);

			WM_M[DELTAg2] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f2Tmp))*dN2;
			WM_M[DELTAg2] += JFaR.GetMat21()*Bv2 + JFaR.GetMat22()*Bw2;
			WM_M[DELTAg2] -= cTmp.Cross(Theta2);
			WM_M[DELTAg2] += Mat3x3(fTmp, f2Tmp*dN2);

			WM_M[DELTAg3] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f3Tmp))*dN3;
			WM_M[DELTAg3] += JFaR.GetMat21()*Bv3 + JFaR.GetMat22()*Bw3;
			WM_M[DELTAg3] -= cTmp.Cross(Theta3);
			WM_M[DELTAg3] += Mat3x3(fTmp, f3Tmp*dN3);

			for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
				WM_F[iCnt] += WM_F2[iCnt];
				WM_M[iCnt] += d.Cross(WM_F2[iCnt]);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
			WM.Sub(6*iNode + 1, 3*iCnt + 1, WM_F[iCnt]);
			WM.Sub(6*iNode + 4, 3*iCnt + 1, WM_M[iCnt]);
		}
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
AerodynamicBeam::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssRes");
	WorkVec.ResizeReset(18);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio iniziale residuo */
SubVectorHandler&
AerodynamicBeam::InitialAssRes(SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("AerodynamicBeam::InitialAssRes");
	WorkVec.ResizeReset(18);

	integer iNode1FirstIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstPositionIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WorkVec.PutRowIndex(12 + iCnt, iNode3FirstIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}


/* assemblaggio residuo */
void
AerodynamicBeam::AssVec(SubVectorHandler& WorkVec)
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssVec");

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[3];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Xn[NODE3] = pNode3->GetXCurr();
	const Mat3x3& Rn3(pNode3->GetRCurr());
	const Vec3& Vn3(pNode3->GetVCurr());
	const Vec3& Wn3(pNode3->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);
	Vec3 f3Tmp(Rn3*f3);

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);
	Vec3 X3Tmp(Xn[NODE3] + f3Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));
	Vec3 V3Tmp(Vn3 + Wn3.Cross(f3Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);
	Mat3x3 RR3(Rn3*Ra3);

	/*
	 * Parametri di rotazione dai nodi 1 e 3 al nodo 2 (nell'ipotesi
	 * che tale trasformazione non dia luogo ad una singolarita')
	 */
	Vec3 g1(ER_Rot::Param, RR2.MulTM(RR1));
	Vec3 g3(ER_Rot::Param, RR2.MulTM(RR3));

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	doublereal rho, c, p, T;
	GetAirProps(Xn[NODE2], rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	/* OUTA */
	doublereal** pvd = pvdOuta;
	int iPnt = 0;

	ResetIterator();

	for (int iNode = 0; iNode < LASTNODE; iNode++) {

		/* Resetta le forze */
		F[iNode].Reset();
		M[iNode].Reset();

		doublereal dsi = pdsi3[iNode];
		doublereal dsf = pdsf3[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi3N(ds,
				Xn[NODE1], Xn[NODE2], Xn[NODE3]);

			doublereal dN1 = ShapeFunc3N(ds, 1);
			doublereal dN2 = ShapeFunc3N(ds, 2);
			doublereal dN3 = ShapeFunc3N(ds, 3);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2 + X3Tmp*dN3);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2 + V3Tmp*dN3);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2 + Wn3*dN3);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != NULL) {
				Vr += pIndVel->GetInducedVelocity(Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR2*Mat3x3(ER_Rot::MatR, g1*dN1 + g3*dN3));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			doublereal dW[6];
			doublereal dTng[6];

			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);

			/* Funzione di calcolo delle forze aerodinamiche */
			aerodata->GetForces(iPnt, dW, dTng, *pvd);

			/* OUTA */
			pvd++;

			// specific for Gauss points force output
			if (fToBeOutput() && IsPGAUSS()) {
				SetData(VTmp, dTng);
			}

			/* Dimensionalizza le forze */
			doublereal dWght = dXds*dsdCsi*PW.dGetWght();
			Vec3 FTmp(RRloc*(Vec3(dTng)*dWght));
			F[iNode] += FTmp;
			M[iNode] += RRloc*(Vec3(&dTng[3])*dWght);
			M[iNode] += (Xr - Xn[iNode]).Cross(FTmp);

			iPnt++;

		} while (GDI.fGetNext(PW));

		/*
		 * Se e' definito il rotore, aggiungere il contributo
		 * alla trazione
		 */
		if (pIndVel != NULL && !fPassiveInducedVelocity) {
			pIndVel->AddForce(GetLabel(),
					F[iNode], M[iNode], Xn[iNode]);
		}

		/* Somma il termine al residuo */
		WorkVec.Add(6*iNode + 1, F[iNode]);
		WorkVec.Add(6*iNode + 4, M[iNode]);
	}
}

void
AerodynamicBeam::AfterConvergence(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	/* Memoria in caso di forze instazionarie */
	switch (aerodata->Unsteady()) {
	case AeroData::STEADY:
		break;

	case AeroData::HARRIS:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case AeroData::BIELAWA:
		for (integer i = 0; i < 3*GDI.iGetNum(); i++) {
			aerodata->Update(i);
		}
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBeam::Output(OutputHandler& OH) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Output");

	if (fToBeOutput()) {
		std::ostream& out = OH.Aerodynamic() << std::setw(8) << GetLabel();

		switch (GetOutput()) {
		case AEROD_OUT_NODE:
			out << " " << std::setw(8) << pBeam->GetLabel()
				<< " ", F[NODE1].Write(out, " ") << " ", M[NODE1].Write(out, " ")
				<< " ", F[NODE2].Write(out, " ") << " ", M[NODE2].Write(out, " ")
				<< " ", F[NODE3].Write(out, " ") << " ", M[NODE3].Write(out, " ");
			break;

		case AEROD_OUT_PGAUSS:
			ASSERT(!pOutput.empty());

			for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
				i != pOutput.end(); i++)
			{
				out << " " << i->alpha
					<< " " << i->f;
			}
			break;

		case AEROD_OUT_STD:
			for (int i = 0; i < 3*GDI.iGetNum(); i++) {
				for (int j = 1; j <= 6; j++) {
					out << " " << pvdOuta[i][j];
				}
				out << " " << pvdOuta[i][OUTA_ALF1]
					<< " " << pvdOuta[i][OUTA_ALF2];
			}
			break;

		default:
			ASSERT(0);
			break;
		}

		out << std::endl;
	}
}

/* AerodynamicBeam - end */


/* Legge un elemento aerodinamico di trave */

Elem *
ReadAerodynamicBeam(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBeam");

	/* Trave */
	unsigned int uBeam = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << std::endl);

	/* verifica di esistenza della trave */
	Elem* p = pDM->pFindElem(Elem::BEAM, uBeam);
	if (p == 0) {
		silent_cerr("Beam3(" << uBeam << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	Beam *pBeam = dynamic_cast<Beam *>(p);
	if (pBeam == 0) {
		silent_cerr("Beam(" << uBeam << ") is not a Beam3 "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	ASSERT(pBeam != 0);

	/* Eventuale rotore */
	InducedVelocity* pIndVel = ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBeam3");

	/* Nodo 1: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode1 = pBeam->pGetNode(1);

	ReferenceFrame RF(pNode1);
	Vec3 f1(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

	Mat3x3 Ra1(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 1 rotation matrix: " << std::endl << Ra1 << std::endl);

	/* Nodo 2: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode2 = pBeam->pGetNode(2);

	RF = ReferenceFrame(pNode2);
	Vec3 f2(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);

	Mat3x3 Ra2(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 2 rotation matrix: " << std::endl << Ra2 << std::endl);

	/* Nodo 3: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode3 = pBeam->pGetNode(3);

	RF = ReferenceFrame(pNode3);
	Vec3 f3(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 3 offset: " << f3 << std::endl);

	Mat3x3 Ra3(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		"Node 3 rotation matrix: " << std::endl << Ra3 << std::endl);

	Shape* pChord = NULL;
	Shape* pForce = NULL;
	Shape* pVelocity = NULL;
	Shape* pTwist = NULL;

	integer iNumber = 0;
	DriveCaller* pDC = NULL;
	AeroData* aerodata = NULL;

	ReadAeroData(pDM, HP,
		&pChord, &pForce, &pVelocity, &pTwist,
		&iNumber, &pDC, &aerodata);

	aerodata->SetNumPoints(3*iNumber);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = true;
	}

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBeam3(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = NULL;

	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBeam,
		AerodynamicBeam(uLabel, pDO, pBeam, pIndVel,
			f1, f2, f3, Ra1, Ra2, Ra3,
			pChord, pForce, pVelocity, pTwist,
			iNumber, aerodata, pDC, bUseJacobian, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ostream& out = pDM->GetLogFile();
	out << "aero3: " << uLabel;

	Vec3 ra1 = Ra1.GetVec(1);
	Vec3 ra3 = Ra1.GetVec(3);
	doublereal dC = pChord->dGet(-1.);
	doublereal dP = pForce->dGet(-1.);
	out
		<< " " << pNode1->GetLabel()
		<< " ", (f1 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f1 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra2.GetVec(1);
	ra3 = Ra2.GetVec(3);
	dC = pChord->dGet(0.);
	dP = pForce->dGet(0.);
	out
		<< " " << pNode2->GetLabel()
		<< " ", (f2 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f2 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra3.GetVec(1);
	ra3 = Ra3.GetVec(3);
	dC = pChord->dGet(1.);
	dP = pForce->dGet(1.);
	out
		<< " " << pNode3->GetLabel()
		<< " ", (f3 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f3 + ra1*(dP + dC/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of DataManager::ReadAerodynamicBeam() */


/* AerodynamicBeam2 - begin */

AerodynamicBeam2::AerodynamicBeam2(
	unsigned int uLabel,
	const DofOwner* pDO,
	const Beam2* pB,
	InducedVelocity* pR,
	const Vec3& fTmp1,
	const Vec3& fTmp2,
	const Mat3x3& Ra1Tmp,
	const Mat3x3& Ra2Tmp,
	const Shape* pC,
	const Shape* pF,
	const Shape* pV,
	const Shape* pT,
	integer iN,
	AeroData* a,
	const DriveCaller* pDC,
	bool bUseJacobian,
	flag fOut
)
: Elem(uLabel, fOut),
AerodynamicElem(uLabel, pDO, fOut),
InitialAssemblyElem(uLabel, fOut),
DriveOwner(pDC),
AerodynamicOutput(fOut, 2*iN),
aerodata(a),
pBeam(pB),
pIndVel(pR),
fPassiveInducedVelocity(0),
f1(fTmp1),
f2(fTmp2),
Ra1(Ra1Tmp),
Ra2(Ra2Tmp),
Ra1_3(Ra1Tmp.GetVec(3)),
Ra2_3(Ra2Tmp.GetVec(3)),
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
GDI(iN),
pdOuta(NULL),
pvdOuta(NULL),
bJacobian(bUseJacobian)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AerodynamicBeam2");

	ASSERT(pBeam != NULL);
	ASSERT(pBeam->GetElemType() == Elem::BEAM);

	pNode1 = pBeam->pGetNode(1);
	pNode2 = pBeam->pGetNode(2);

	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);

#ifdef DEBUG
	if (pIndVel != NULL) {
		ASSERT(pIndVel->GetElemType() == Elem::INDUCEDVELOCITY);
	}
#endif /* DEBUG */

	SAFENEWARR(pdOuta, doublereal, 20*2*iN);
	SAFENEWARR(pvdOuta, doublereal*, 2*iN);
	for (integer i = 20*2*iN; i-- > 0; ) {
		pdOuta[i] = 0.;
	}
	for (integer i = 2*iN; i-- > 0; ) {
		pvdOuta[i] = pdOuta+20*i;
	}
}

AerodynamicBeam2::~AerodynamicBeam2(void)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::~AerodynamicBeam2");

	SAFEDELETEARR(pvdOuta);
	SAFEDELETEARR(pdOuta);

	SAFEDELETE(aerodata);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
void
AerodynamicBeam2::SetOutputFlag(flag f)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::SetOutputFlag");

	ToBeOutput::SetOutputFlag(f);
	AerodynamicOutput::SetOutputFlag(f, 2*GDI.iGetNum());
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AerodynamicBeam2::Restart(std::ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBeam2::Restart");
	out << "  aerodynamic beam2: " << GetLabel()
		<< ", " << pBeam->GetLabel();
	if (pIndVel != NULL) {
		out << ", rotor, " << pIndVel->GetLabel();
	}
	out << ", reference, node, ", f1.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra1.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra1.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f2.Write(out, ", ")
		<< ", reference, node, 1, ", (Ra2.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra2.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", "
		<< GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ";" << std::endl;
}

static const doublereal pdsi2[] = { -1., 0. };
static const doublereal pdsf2[] = { 0., 1. };

/* Jacobian assembly */
VariableSubMatrixHandler&
AerodynamicBeam2::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal  dCoef  ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering AerodynamicBeam2::AssJac()" << std::endl);

	if (!bJacobian)	{
		WorkMat.SetNullMatrix();
		return WorkMat;
	}

	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();

	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}

	doublereal dW[6];

	/* array di vettori per via del ciclo sui nodi ... */
	Vec3 Xn[2];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);

	Vec3 fTmp[2];
	fTmp[NODE1] = f1Tmp;
	fTmp[NODE2] = f2Tmp;

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));

	Vec3 Omega1Crossf1(Wn1.Cross(f1Tmp));
	Vec3 Omega2Crossf2(Wn2.Cross(f2Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);

	/*
	 * half of relative rotation vector from node 1 to node 2
	 */
	Vec3 overline_theta(Vec3(ER_Rot::Param, RR1.MulTM(RR2))/2.);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione del nodo 2 perche'
	 * non dovrebbero cambiare "molto")
	 */
	Vec3 Xmid = (Xn[NODE2] + Xn[NODE1])/2.;
	doublereal rho, c, p, T;
	GetAirProps(Xmid, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	/* OUTA */
	doublereal** pvd = pvdOuta;
	int iPnt = 0;

	ResetIterator();

	for (int iNode = 0; iNode < LASTNODE; iNode++) {
		doublereal dsi = pdsi2[iNode];
		doublereal dsf = pdsf2[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		Mat3x3 WM_F[4], WM_M[4];
		WM_F[DELTAx1].Reset();
		WM_F[DELTAg1].Reset();
		WM_F[DELTAx2].Reset();
		WM_F[DELTAg2].Reset();
		WM_M[DELTAx1].Reset();
		WM_M[DELTAg1].Reset();
		WM_M[DELTAx2].Reset();
		WM_M[DELTAg2].Reset();

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi2N(ds,
				Xn[NODE1], Xn[NODE2]);

			doublereal dN1 = ShapeFunc2N(ds, 1);
			doublereal dN2 = ShapeFunc2N(ds, 2);

			// note: identical to dN1 and dN2; see tecman.pdf
#if 0
			doublereal dNN1 = (1. + dN1 - dN2)/2.;
			doublereal dNN2 = (1. + dN2 - dN1)/2.;
#endif

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2);
			Vec3 thetar(overline_theta*dN2);

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != NULL) {
				Vr += pIndVel->GetInducedVelocity(Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			/* Contributo dell'eventuale sup. mobile */
			dTw += dGet();

			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_;
			 * quindi lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR1*Mat3x3(ER_Rot::MatR, thetar));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw(dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g
				 * e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(&dW[0]);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);
			/* Funzione di calcolo delle forze aerodinamiche */
			doublereal Fa0[6];
			Mat6x6 JFa;

			aerodata->GetForcesJac(iPnt, dW, Fa0, JFa, *pvd);

			/* OUTA */
			pvd++;

			doublereal cc = dXds*dsdCsi*PW.dGetWght();

			// rotate force, couple and Jacobian matrix in absolute frame
			Mat6x6 JFaR = MultRMRt(JFa, RRloc, cc);

			// force and moment about the node
			Vec3 fTmp(RRloc*(Vec3(&Fa0[0])*dCoef));
			Vec3 d(Xr - Xn[iNode]);
			Vec3 cTmp(RRloc*(Vec3(&Fa0[3])*dCoef) + d.Cross(fTmp));

			Mat3x3 Bv1((Vr - Omega1Crossf1)*(dN1*dCoef));
			Mat3x3 Bv2((Vr - Omega2Crossf2)*(dN2*dCoef));

			Mat3x3 Bw1((Wr - Wn1)*(dN1*dCoef));
			Mat3x3 Bw2((Wr - Wn2)*(dN2*dCoef));

			Mat3x3 WM_F2[4];

			// f <-> x
			WM_F2[DELTAx1] = JFaR.GetMat11()*dN1;

			WM_F2[DELTAx2] = JFaR.GetMat11()*dN2;

			doublereal delta;

			// c <-> x
			delta = (iNode == NODE1) ? 1. : 0.;
			WM_M[DELTAx1] += JFaR.GetMat21()*dN1 - Mat3x3(fTmp*(dN1 - delta));

			delta = (iNode == NODE2) ? 1. : 0.;
			WM_M[DELTAx2] += JFaR.GetMat21()*dN2 - Mat3x3(fTmp*(dN2 - delta));

			// f <-> g
			WM_F2[DELTAg1] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f1Tmp))*dN1;
			WM_F2[DELTAg1] += JFaR.GetMat11()*Bv1 + JFaR.GetMat12()*Bw1;
			WM_F2[DELTAg1] -= Mat3x3(fTmp*dN1);
		
			WM_F2[DELTAg2] = (JFaR.GetMat12() - JFaR.GetMat11()*Mat3x3(f2Tmp))*dN2;
			WM_F2[DELTAg2] += JFaR.GetMat11()*Bv2 + JFaR.GetMat12()*Bw2;
			WM_F2[DELTAg2] -= Mat3x3(fTmp*dN2);
		
			// c <-> g
			WM_M[DELTAg1] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f1Tmp))*dN1;
			WM_M[DELTAg1] += JFaR.GetMat21()*Bv1 + JFaR.GetMat22()*Bw1;
			WM_M[DELTAg1] -= Mat3x3(cTmp*dN1);
			WM_M[DELTAg1] += Mat3x3(fTmp, f1Tmp*dN1);

			WM_M[DELTAg2] += (JFaR.GetMat22() - JFaR.GetMat21()*Mat3x3(f2Tmp))*dN2;
			WM_M[DELTAg2] += JFaR.GetMat21()*Bv2 + JFaR.GetMat22()*Bw2;
			WM_M[DELTAg2] -= Mat3x3(cTmp*dN2);
			WM_M[DELTAg2] += Mat3x3(fTmp, f2Tmp*dN2);

			for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
				WM_F[iCnt] += WM_F2[iCnt];
				WM_M[iCnt] += d.Cross(WM_F2[iCnt]);
			}

			iPnt++;

		} while (GDI.fGetNext(PW));

		for (int iCnt = 0; iCnt < 2*LASTNODE; iCnt++) {
			WM.Sub(6*iNode + 1, 3*iCnt + 1, WM_F[iCnt]);
			WM.Sub(6*iNode + 4, 3*iCnt + 1, WM_M[iCnt]);
		}
	}

	return WorkMat;
}

/* assemblaggio residuo */
SubVectorHandler&
AerodynamicBeam2::AssRes(
	SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */
)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AssRes");
	WorkVec.ResizeReset(12);

	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio iniziale residuo */
SubVectorHandler&
AerodynamicBeam2::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& /* XCurr */
)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::InitialAssRes");
	WorkVec.ResizeReset(12);

	integer iNode1FirstIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstIndex + iCnt);
	}

	AssVec(WorkVec);

	return WorkVec;
}

/* assemblaggio residuo */
void
AerodynamicBeam2::AssVec(SubVectorHandler& WorkVec)
{
	DEBUGCOUTFNAME("AerodynamicBeam2::AssVec");

	doublereal dTng[6];
	doublereal dW[6];

	Vec3 Xn[LASTNODE];

	/* Dati dei nodi */
	Xn[NODE1] = pNode1->GetXCurr();
	const Mat3x3& Rn1(pNode1->GetRCurr());
	const Vec3& Vn1(pNode1->GetVCurr());
	const Vec3& Wn1(pNode1->GetWCurr());

	Xn[NODE2] = pNode2->GetXCurr();
	const Mat3x3& Rn2(pNode2->GetRCurr());
	const Vec3& Vn2(pNode2->GetVCurr());
	const Vec3& Wn2(pNode2->GetWCurr());

	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);

	Vec3 X1Tmp(Xn[NODE1] + f1Tmp);
	Vec3 X2Tmp(Xn[NODE2] + f2Tmp);

	Vec3 V1Tmp(Vn1 + Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2 + Wn2.Cross(f2Tmp));

	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	Mat3x3 RR2(Rn2*Ra2);

	/*
	 * half of relative rotation vector from node 1 to node 2
	 */
	Vec3 overline_theta(Vec3(ER_Rot::Param, RR1.MulTM(RR2))/2.);

	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pIndVel != NULL) {
		Rotor *pRotor = dynamic_cast<Rotor *>(pIndVel);
		if (pRotor != 0) {
			dOmega = pRotor->dGetOmega();
		}
	}

	/*
	 * Dati "permanenti" (uso solo la posizione di mezzo perche'
	 * non dovrebbero cambiare "molto")
	 */
	Vec3 Xmid = (Xn[NODE2] + Xn[NODE1])/2.;
	doublereal rho, c, p, T;
	GetAirProps(Xmid, rho, c, p, T);	/* p, T no used yet */
	aerodata->SetAirData(rho, c);

	/* OUTA */
	doublereal** pvd = pvdOuta;
	int iPnt = 0;

	ResetIterator();

	for (int iNode = 0; iNode < LASTNODE; iNode++) {

		/* Resetta i dati */
		F[iNode].Reset();
		M[iNode].Reset();

		doublereal dsi = pdsi2[iNode];
		doublereal dsf = pdsf2[iNode];

		doublereal dsm = (dsf + dsi)/2.;
		doublereal dsdCsi = (dsf - dsi)/2.;

		/* Ciclo sui punti di Gauss */
		PntWght PW = GDI.GetFirst();
		do {
			doublereal dCsi = PW.dGetPnt();
			doublereal ds = dsm + dsdCsi*dCsi;
			doublereal dXds = DxDcsi2N(ds, Xn[NODE1], Xn[NODE2]);

			doublereal dN1 = ShapeFunc2N(ds, 1);
			doublereal dN2 = ShapeFunc2N(ds, 2);

			Vec3 Xr(X1Tmp*dN1 + X2Tmp*dN2);
			Vec3 Vr(V1Tmp*dN1 + V2Tmp*dN2);
			Vec3 Wr(Wn1*dN1 + Wn2*dN2);
			Vec3 thetar(overline_theta*((1. + dN2 - dN1)/2.));

			/* Contributo di velocita' del vento */
			Vec3 VTmp(0.);
			if (fGetAirVelocity(VTmp, Xr)) {
				Vr -= VTmp;
			}

			/*
			 * Se l'elemento e' collegato ad un rotore,
			 * aggiunge alla velocita' la velocita' indotta
			 */
			if (pIndVel != NULL) {
				Vr += pIndVel->GetInducedVelocity(Xr);
			}

			/* Copia i dati nel vettore di lavoro dVAM */
			doublereal dTw = Twist.dGet(ds);
			dTw += dGet(); /* Contributo dell'eventuale sup. mobile */
			aerodata->SetSectionData(dCsi,
				Chord.dGet(ds),
				ForcePoint.dGet(ds),
				VelocityPoint.dGet(ds),
				dTw,
				dOmega);

			/*
			 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
			 * lo uso per correggere la matrice di rotazione
			 * dal sistema aerodinamico a quello globale
			 */
			Mat3x3 RRloc(RR1*Mat3x3(ER_Rot::MatR, thetar));
			if (dTw != 0.) {
				doublereal dCosT = cos(dTw);
				doublereal dSinT = sin(dTw);
				/* Assumo lo svergolamento positivo a cabrare */
				Mat3x3 RTw( dCosT, dSinT, 0.,
					-dSinT, dCosT, 0.,
					0.,    0.,    1.);
				/*
				 * Allo stesso tempo interpola le g e aggiunge lo svergolamento
				 */
				RRloc = RRloc*RTw;
			}

			/*
			 * Ruota velocita' e velocita' angolare nel sistema
			 * aerodinamico e li copia nel vettore di lavoro dW
			 */
			VTmp = RRloc.MulTV(Vr);
			VTmp.PutTo(dW);

			Vec3 WTmp = RRloc.MulTV(Wr);
			WTmp.PutTo(&dW[3]);

			/* Funzione di calcolo delle forze aerodinamiche */
			aerodata->GetForces(iPnt, dW, dTng, *pvd);

			/* OUTA */
			pvd++;

			// specific for Gauss points force output
			if (fToBeOutput() && IsPGAUSS()) {
				SetData(VTmp, dTng);
			}

			/* Dimensionalizza le forze */
			doublereal dWght = dXds*dsdCsi*PW.dGetWght();
			Vec3 FTmp(RRloc*(Vec3(&dTng[0])*dWght));
			F[iNode] += FTmp;
			M[iNode] += RRloc*(Vec3(&dTng[3])*dWght);
			M[iNode] += (Xr - Xn[iNode]).Cross(FTmp);

			iPnt++;

		} while (GDI.fGetNext(PW));

		/* Se e' definito il rotore, aggiungere il contributo alla trazione */
		if (pIndVel != NULL && !fPassiveInducedVelocity) {
			pIndVel->AddForce(GetLabel(),
				F[iNode], M[iNode], Xn[iNode]);
		}

		/* Somma il termine al residuo */
		WorkVec.Add(6*iNode + 1, F[iNode]);
		WorkVec.Add(6*iNode + 4, M[iNode]);
	}
}

void
AerodynamicBeam2::AfterConvergence(const VectorHandler& /* X */ ,
	const VectorHandler& /* XP */ )
{
	/* Memoria in caso di forze instazionarie */
	switch (aerodata->Unsteady()) {
	case AeroData::STEADY:
		break;

	case AeroData::HARRIS:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case AeroData::BIELAWA:
		for (integer i = 0; i < 2*GDI.iGetNum(); i++) {
	 		aerodata->Update(i);
		}
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBeam2::Output(OutputHandler& OH ) const
{
	DEBUGCOUTFNAME("AerodynamicBeam2::Output");

	if (fToBeOutput()) {
		std::ostream& out = OH.Aerodynamic() << std::setw(8) << GetLabel();

		switch (GetOutput()) {

		case AEROD_OUT_NODE:
			out << " " << std::setw(8) << pBeam->GetLabel()
				<< " ", F[NODE1].Write(out, " ") << " ", M[NODE1].Write(out, " ")
				<< " ", F[NODE2].Write(out, " ") << " ", M[NODE2].Write(out, " ");
			break;

		case AEROD_OUT_PGAUSS:
			ASSERT(!pOutput.empty());

			for (std::vector<Aero_output>::const_iterator i = pOutput.begin();
				i != pOutput.end(); i++)
			{
				out << " " << i->alpha
					<< " " << i->f;
			}
			break;

		case AEROD_OUT_STD:
			for (int i = 0; i < 2*GDI.iGetNum(); i++) {
				for (int j = 1; j <= 6; j++) {
					out << " " << pvdOuta[i][j];
				}
				out << " " << pvdOuta[i][OUTA_ALF1]
					<< " " << pvdOuta[i][OUTA_ALF2];
			}
			break;

		default:
			ASSERT(0);
			break;
		}

		out << std::endl;
	}
}

/* AerodynamicBeam2 - end */


/* Legge un elemento aerodinamico di trave a due nodi */

Elem *
ReadAerodynamicBeam2(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner *pDO,
	unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBeam2");

	/* Trave */
	unsigned int uBeam = (unsigned int)HP.GetInt();

	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << std::endl);

	/* verifica di esistenza della trave */
	Elem *p = pDM->pFindElem(Elem::BEAM, uBeam);
	if (p == 0) {
		silent_cerr("Beam2(" << uBeam << ") not defined "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	Beam2* pBeam = dynamic_cast<Beam2 *>(p);
	if (pBeam == 0) {
		silent_cerr("Beam(" << uBeam << ") is not a Beam2 "
				"at line " << HP.GetLineData()
				<< std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/* Eventuale rotore */
	InducedVelocity* pIndVel = ReadInducedVelocity(pDM, HP, uLabel, "AerodynamicBeam2");

	/* Nodo 1: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode1 = pBeam->pGetNode(1);

	ReferenceFrame RF(pNode1);
	Vec3 f1(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << std::endl);

	Mat3x3 Ra1(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 1 rotation matrix: " << std::endl << Ra1 << std::endl);

	/* Nodo 2: */

	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode2 = pBeam->pGetNode(2);

	RF = ReferenceFrame(pNode2);
	Vec3 f2(HP.GetPosRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << std::endl);

	Mat3x3 Ra2(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 2 rotation matrix: " << std::endl << Ra2 << std::endl);

	Shape* pChord = NULL;
	Shape* pForce = NULL;
	Shape* pVelocity = NULL;
	Shape* pTwist = NULL;

	integer iNumber = 0;
	DriveCaller* pDC = NULL;
	AeroData* aerodata = NULL;

	ReadAeroData(pDM, HP,
		     &pChord, &pForce, &pVelocity, &pTwist,
		     &iNumber, &pDC, &aerodata);

	aerodata->SetNumPoints(2*iNumber);

	bool bUseJacobian(false);
	if (HP.IsKeyWord("jacobian")) {
		bUseJacobian = true;
	}

	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	if (HP.IsArg()) {
		if (HP.IsKeyWord("std")) {
			fOut |= AerodynamicOutput::AEROD_OUT_STD;
		} else if (HP.IsKeyWord("gauss")) {
			fOut |= AerodynamicOutput::AEROD_OUT_PGAUSS;
		} else if (HP.IsKeyWord("node")) {
			fOut |= AerodynamicOutput::AEROD_OUT_NODE;
		} else {
			silent_cerr("AerodynamicBeam2(" << uLabel << "): "
				"unknown output mode at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else if (fOut) {
		fOut |= AerodynamicOutput::AEROD_OUT_STD;
	}

	Elem* pEl = NULL;

	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBeam2,
		AerodynamicBeam2(uLabel, pDO, pBeam, pIndVel,
				f1, f2, Ra1, Ra2,
				pChord, pForce, pVelocity, pTwist,
				iNumber, aerodata, pDC, bUseJacobian, fOut));

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::ostream& out = pDM->GetLogFile();
	out << "aero2: " << uLabel;

	Vec3 ra1 = Ra1.GetVec(1);
	Vec3 ra3 = Ra1.GetVec(3);
	doublereal dC = pChord->dGet(-1.);
	doublereal dP = pForce->dGet(-1.);
	out
		<< " " << pNode1->GetLabel()
		<< " ", (f1 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f1 + ra1*(dP + dC/4.)).Write(out, " ");

	ra1 = Ra2.GetVec(1);
	ra3 = Ra2.GetVec(3);
	dC = pChord->dGet(1.);
	dP = pForce->dGet(1.);
	out
		<< " " << pNode2->GetLabel()
		<< " ", (f2 + ra1*(dP - dC*3./4.)).Write(out, " ")
		<< " ", (f2 + ra1*(dP + dC/4.)).Write(out, " ")
		<< std::endl;

	return pEl;
} /* End of DataManager::ReadAerodynamicBeam2() */

