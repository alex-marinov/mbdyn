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

/*
 * Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cfloat>
#include <limits>

#include "dataman.h"
#include "constltp.h"
#include "shapefnc.h"
#include "beam.h"
#include "hbeam.h"
#if 0	/* not implemented yet */
#include <pzhbeam.h>
#endif
#include "hbeam_interp.h"

/*
 * Nota: non e' ancora stato implementato il contributo 
 * della ViscoElasticHBeam all'assemblaggio iniziale
 */

/*
 * Nota: la parte viscoelastica va rivista in accordo con la piu' 
 * recente formulazione delle derivate delle deformazioni nel sistema
 * materiale
 */

/* HBeam - begin */

/* Costruttore normale */
HBeam::HBeam(unsigned int uL, 
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& F1,
		const Vec3& F2,
		const Mat3x3& R1, const Mat3x3& R2,
		const ConstitutiveLaw6D* pd,
		flag fOut)
: Elem(uL, fOut), 
ElemGravityOwner(uL, fOut), 
InitialAssemblyElem(uL, fOut),
bFirstRes(true)
{
	/* Validazione dati */
	ASSERT(pN1 != NULL);
	ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN2 != NULL);
	ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
   
	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	f[NODE1] = F1;
	f[NODE2] = F2;

	Rn[NODE1] = R1;
	Rn[NODE2] = R2;
	
	/*
	RRef = R = interpolazione di R1 e R2 ;
	 */
	
	pD = NULL; 
	SAFENEWWITHCONSTRUCTOR(pD,
			ConstitutiveLaw6DOwner,
			ConstitutiveLaw6DOwner(pd));
	
	Omega = Zero3; 
	Az = Zero6;
	AzRef = Zero6;
	AzLoc = Zero6;
	DefLoc = Zero6;
	DefLocRef = Zero6;
	p = Zero3;
	g = Zero3;
	L0 = Zero3;
	L = Zero3;

	DsDxi();
}


HBeam::~HBeam(void) 
{
	/* Distrugge il legame costitutivo */
	ASSERT(pD != NULL);
	if (pD != NULL) {      
		SAFEDELETE(pD);
	}
}

/* Accesso ai dati privati */
unsigned int
HBeam::iGetNumPrivData(void) const
{
	return 12;
}

unsigned int
HBeam::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	/*
	 * {ex|k{x|y|z}}
	 */

	unsigned int idx = 0;

	switch (s[0]) {
	case 'F':
		idx += 6;
	case 'e':
		switch (s[1]) {
		case 'x':
			idx += 1;
			break;

		case 'y':
		case 'z':
			return 0;

		default:
			return 0;
		}
		break;

	case 'M':
		idx += 6;
	case 'k':
		idx += 3;
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
HBeam::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= 6);

	switch (i) {
	case 1:
	case 4:
	case 5:
	case 6:
		return DefLoc.dGet(i);

	case 7:
	case 10:
	case 11:
	case 12:
		return AzLoc.dGet(i);

	case 2:
	case 3:
		silent_cerr("HBeam(" << GetLabel() << "): "
			"not allowed to return shear strain" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	case 8:
	case 9:
		silent_cerr("HBeam(" << GetLabel() << "): "
			"not allowed to return shear force" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);

	default:
		silent_cerr("HBeam(" << GetLabel() << "): "
			"illegal private data " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
HBeam::DsDxi(void)
{
	/* Calcola il ds/dxi e le deformazioni iniziali */
	Mat3x3 RTmp[NUMNODES];
	Vec3 yTmp[NUMNODES];
	Vec3 fTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		RTmp[i] = pNode[i]->GetRCurr()*Rn[i];
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
		yTmp[i] = pNode[i]->GetXCurr() + fTmp[i];
	}
	
	xi = 0.5;
	dsdxi = 1.0;
	/* Calcolo i wder ... */
	ComputeInterpolation(yTmp, RTmp, fTmp,
			xi, dsdxi,
			p, R,
			L, Rho);
	
	doublereal d = L.Dot();
	if (d > std::numeric_limits<doublereal>::epsilon()) {
		d = std::sqrt(d);
	} else {
		silent_cerr("HBeam(" << GetLabel() << ") "
			"has singular metric; aborting..." << std::endl);
		
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}
	
	dsdxi = dsdxi*d;
	dxids = 1./dsdxi;

	/* Calcolo le caratteristiche iniziali ... */
	ComputeInterpolation(yTmp, RTmp, fTmp,
			xi, dxids,
			p, R,
			L, Rho);

	/* Grandezze iniziali e di riferimento */
	/* FIXME: fare un temporaneo per i trasposti ... */
	RRef = R;
	Mat3x3 RT(R.Transpose());
	Rho0 = RT*Rho;
	LRef = L;
	L0 = RT*L;
}


/* Calcola la velocita' angolare delle sezioni a partire da quelle dei nodi */
void 
HBeam::Omega0(void)
{
#if 0
	/* Modo consistente: */      
	Mat3x3 RNod[NUMNODES];
	Vec3 w[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {     
		RNod[i] = pNode[i]->GetRCurr();
		w[i] = pNode[i]->GetWCurr();
	}
	
	/*
	 * Calcolo i parametri di rotazione della rotazione relativa
	 * tra inizio e fine e li dimezzo nell'ipotesi che siano limitati
	 */
	Vec3 gTmp(MatR2gparam(RNod[NODE2].Transpose()*RNod[NODE1]));
	
	/*
	 * Le derivate dei parametri di rotazione si ricavano da omega
	 */
	Vec3 g1P(Mat3x3(MatGm1, gTmp*(-.5))*w[NODE1]);
	Vec3 g2P(Mat3x3(MatGm1, gTmp*.5)*w[NODE2]);

        Vec3 gPTmp(g1P*dN2[NODE1]+g2P*dN2[NODE2]);
        Omega = Mat3x3(MatG, gTmp)*gPTmp;
	
#if 0
	/* Modo brutale: interpolo le velocita' dei nodi */
	Vec3 w[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		w[i] = pNode[i]->GetWCurr();
	}
	Omega[i] = w[NODE1]*dN2[NODE1]+w[NODE2]*dN2[NODE2];
#endif /* 0 */
#endif
}


/* Contributo al file di restart */
std::ostream&
HBeam::Restart(std::ostream& out) const
{
	return Restart_(out)<< ';' << std::endl;
}

std::ostream&
HBeam::Restart_(std::ostream& out) const
{ 
	out << "  beam2: " << GetLabel();
	for (unsigned int i = 0; i < NUMNODES; i++) {
		out << ", " << pNode[i]->GetLabel() << ", reference, node, ", 
		f[i].Write(out, ", ");
	}
	out << ", reference, global,"
		<< "1, ", (R.GetVec(1)).Write(out, ", ") << ", "
		<< "2, ", (R.GetVec(2)).Write(out, ", ") << ", ",
	pD->pGetConstLaw()->Restart(out);
	
	return out;
}

void
HBeam::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	pD->AfterConvergence(DefLoc);
}

/* Assembla la matrice */
void
HBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& /* WMB */ ,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("HBeam::AssStiffnessMat");
	
	/*
	 * La matrice arriva gia' dimensionata
	 * e con gli indici di righe e colonne a posto
	 */
   
	/* Recupera i dati dei nodi */
	Vec3 yTmp[NUMNODES];
	Vec3 fTmp[NUMNODES];
	Mat3x3 RTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		RTmp[i] = pNode[i]->GetRCurr()*Rn[i];
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
		yTmp[i] = pNode[i]->GetXCurr() + fTmp[i];
	}   
	ComputeFullInterpolation(yTmp, RTmp, fTmp,
			xi, dxids,
			p, R,
			RdP, pdP, pdp,
			L, Rho,
			RhodP, LdP, Ldp);
	/* Legame costitutivo (viene generato sempre) */
	DRef = MultRMRt(pD->GetFDE(), R);
	
	/* Derivate delle deformazioni rispetto alle incognite nodali */
	Mat6x6 AzTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Delta - deformazioni */
		AzTmp[i] = Mat6x6(Ldp[i]*dCoef,
				Zero3x3,
				(LdP[i] + L.Cross(RdP[i]))*dCoef,
				(RhodP[i] + Rho.Cross(RdP[i]))*dCoef);
		
		/* Delta - azioni interne */
		AzTmp[i] = DRef*AzTmp[i];
		
		/* Correggo per la rotazione da locale a globale */
		AzTmp[i].SubMat12((AzRef.GetVec1()*dCoef).Cross(RdP[i]));
		AzTmp[i].SubMat22((AzRef.GetVec2()*dCoef).Cross(RdP[i]));
	}
   
	Vec3 bTmp[2];
	
	bTmp[0] = p - pNode[NODE1]->GetXCurr();
	bTmp[1] = p - pNode[NODE2]->GetXCurr();
   
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Equazione all'indietro: */
		WMA.Sub(1, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Sub(1, 6*i + 4, 
			AzTmp[i].GetMat12()
			- (AzRef.GetVec1()*dCoef).Cross(pdP[i]));
		
		WMA.Sub(4, 6*i + 1,
			AzTmp[i].GetMat21()
			- (AzRef.GetVec1()*dCoef).Cross(pdp[i])
			+ bTmp[0].Cross(AzTmp[i].GetMat11()));
		WMA.Sub(4, 6*i + 4, 
			AzTmp[i].GetMat22()
			+ bTmp[0].Cross(AzTmp[i].GetMat12()));
		
		/* Equazione in avanti: */
		WMA.Add(7, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Add(7, 6*i + 4,
			AzTmp[i].GetMat12()
			- (AzRef.GetVec1()*dCoef).Cross(pdP[i]));
		
		WMA.Add(10, 6*i + 1,
			AzTmp[i].GetMat21()
			- (AzRef.GetVec1()*dCoef).Cross(pdp[i])
			+ bTmp[1].Cross(AzTmp[i].GetMat11()));
		WMA.Add(10, 6*i + 4,
			AzTmp[i].GetMat22()
			+ bTmp[1].Cross(AzTmp[i].GetMat12()));
	}
	
	/* correzione alle equazioni */
	Mat3x3 FTmp(MatCross, AzRef.GetVec1()*dCoef);
	WMA.Sub(4, 1, FTmp);
	WMA.Add(10, 7, FTmp);
};


/* Assembla il residuo */
void
HBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("HBeam::AssStiffnessVec");
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei tre nodi
	 */
	
	/* Recupera i dati dei nodi */
	Vec3 yTmp[NUMNODES];
	Vec3 fTmp[NUMNODES];
	Mat3x3 RTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		RTmp[i] = pNode[i]->GetRCurr()*Rn[i];
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
		yTmp[i] = pNode[i]->GetXCurr() + fTmp[i];
	}   
	
	/* Interpolazione generica */
	ComputeInterpolation(yTmp, RTmp, fTmp,
			xi, dxids,
			p, R,
			L, Rho);
	
	Mat3x3 RT(R.Transpose());
	DefLoc = Vec6(RT*L - L0, RT*Rho - Rho0);

	/* Calcola le azioni interne */
	pD->Update(DefLoc);
	AzLoc = pD->GetF();
	
	/* corregge le azioni interne locali (piezo, ecc) */
	AddInternalForces(AzLoc);

	/* Porta le azioni interne nel sistema globale */
	Az = MultRV(AzLoc, R);
	
	WorkVec.Add(1, Az.GetVec1());
	WorkVec.Add(4, (p - pNode[NODE1]->GetXCurr()).Cross(Az.GetVec1()) + Az.GetVec2());
	WorkVec.Sub(7, Az.GetVec1());
	WorkVec.Sub(10, Az.GetVec2() + (p - pNode[NODE2]->GetXCurr()).Cross(Az.GetVec1()));
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler&
HBeam::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("HBeam::AssJac => AssStiffnessMat");
	
	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
	
	FullSubMatrixHandler& WM = WorkMat.SetFull();   
	
	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	WM.ResizeReset(12, 12);
	
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}      
	
	AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
	
	return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler&
HBeam::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("HBeam::AssRes => AssStiffnessVec");
	
	integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
	integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
	
	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	WorkVec.ResizeReset(12);

	for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstMomIndex + iCnt);
	}      
	
	AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
	
	return WorkVec;
}

    
/* Settings iniziali, prima della prima soluzione */
void
HBeam::SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	(Mat3x3&)RRef = R;
	(Vec3&)LRef = L;
	(Vec6&)DefLocRef = DefLoc;
	(Vec6&)AzRef = Az;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	(Mat6x6&)DRef = MultRMRt(pD->GetFDE(), RRef);      
	
	bFirstRes = true;
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void
HBeam::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{  
	/*
	 * Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE
	 */
#if 0	
	/* Recupera i dati dei nodi */  
	Vec3   gNod[NUMNODES];
	Vec3   xTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {            
		gNod[i] = pNode[i]->GetgRef();
		xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRRef()*f[i];
	}
	
	Mat3x3 RDelta;
	Vec3 gGrad;
	
	/* Aggiorna le grandezze della trave nel punto di valutazione */
	
	/* Posizione */
	p = InterpState(xTmp[NODE1], xTmp[NODE2]);
	
	/* Matrici di rotazione */
	g = InterpState(gNod[NODE1], gNod[NODE2]);
	RDelta = Mat3x3(MatR, g);
	R = RRef = RDelta*R;
	
	/* Derivate della posizione */
	L = LRef = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	
	/* Derivate dei parametri di rotazione */
	gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
	
	/* Per le deformazioni nel sistema del materiale */
	Mat3x3 RTmp(R.Transpose());
	
	/*
	 * Calcola le deformazioni nel sistema locale
	 * nei punti di valutazione
	 */
	DefLoc = DefLocRef = Vec6(RTmp*L-L0,
			RTmp*(Mat3x3(MatG, g)*gGrad)+DefLoc.GetVec2());
	
	/* Calcola le azioni interne */
	pD->Update(DefLoc);
	AzLoc = pD->GetF();
	
	/* corregge le azioni interne locali (piezo, ecc) */
	AddInternalForces(AzLoc);
	
	/* Porta le azioni interne nel sistema globale */
	Az = AzRef = MultRV(AzLoc, R);
	
	/* Aggiorna il legame costitutivo di riferimento */
	DRef = MultRMRt(pD->GetFDE(), RRef);

	bFirstRes = true;
#endif
}


/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
HBeam::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Beams() << std::setw(8) << GetLabel() << " " 
			<< AzLoc.GetVec1() << " " << AzLoc.GetVec2() << std::endl;
	}
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
HBeam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) 
{ 
	DEBUGCOUTFNAME("HBeam::InitialAssJac => AssStiffnessMat");
	
	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(12, 12);
	
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
		WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	}      
	
	AssStiffnessMat(WM, WM, 1., XCurr, XCurr);
	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler&
HBeam::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr) 
{ 
	DEBUGCOUTFNAME("HBeam::InitialAssRes => AssStiffnessVec");

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	WorkVec.ResizeReset(12);

	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	}      
	
	AssStiffnessVec(WorkVec, 1., XCurr, XCurr);
	return WorkVec;
}


const StructNode*
HBeam::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 2);
	switch (i) {
	case 1:
	case 2:
		return pNode[i-1];
	default:
		throw HBeam::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}


void
HBeam::GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& r) const
{
	ASSERT(part == 1);
	part--;
   
	x = p;
	r = R;
}

void
HBeam::GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const
{
	ASSERT(part == 1);
	part--;
   
	v = Zero3;
	w = Zero3;
}


#ifdef USE_ADAMS
std::ostream& 
HBeam::WriteAdamsDummyPartCmd(std::ostream& out,
		unsigned int part,
		unsigned int firstId) const
{
	Vec3 xTmp[NUMNODES];
	
	part--;
	
	for (unsigned int i = 0; i <= 1; i++) {
		xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
	}
	
	Mat3x3 RT(R.Transpose());
	
	out << psAdamsElemCode[GetElemType()] << "_" << GetLabel()
		<< "_" << 1 << std::endl
		<< firstId << " "
		<< p << " " 
		<< MatR2EulerAngles(R)*dRaDegr << " "
		<< RT*(xTmp[NODE1]-p) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[part]->GetRCurr())*dRaDegr */ << " "
		<< RT*(xTmp[NODE2]-p) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[1+part]->GetRCurr())*dRaDegr */ << std::endl;
	
	return out;
}
#endif /* USE_ADAMS */

/* HBeam - end */


#ifdef VISCOELASTIC_BEAM
/* ViscoElasticHBeam - begin */

/* Costruttore normale */
ViscoElasticHBeam::ViscoElasticHBeam(unsigned int uL, 
		const StructNode* pN1, 
		const StructNode* pN2, 
		const Vec3& F1,
		const Vec3& F2,
		const Mat3x3& r,
		const ConstitutiveLaw6D* pd, 
		flag fOut)
: Elem(uL, fOut),
HBeam(uL, pN1, pN2, F1, F2, r, pd, fOut)
{
	LPrimeRef = LPrime = Zero3;  
	gPrime = Zero3;
	
	DefPrimeLoc = DefPrimeLocRef = Zero6;
	
	/* Nota: DsDxi() viene chiamata dal costruttore di Beam */
	HBeam::Omega0();
}


/* Assembla la matrice */
void
ViscoElasticHBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticHBeam::AssStiffnessMat");
	
	/*
	 * La matrice arriva gia' dimensionata
	 * e con gli indici di righe e colonne a posto
	 */
	
	/* offset nel riferimento globale */
	Vec3 fTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		fTmp[i] = pNode[i]->GetRCurr()*f[i];
	}
	
	Mat6x6 AzTmp[NUMNODES];
	Mat6x6 AzPrimeTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Delta - deformazioni */
		AzTmp[i] = AzPrimeTmp[i] = Mat6x6(Mat3x3(dN2P[i]*dsdxi),
				Zero3x3,
				Mat3x3(LRef*(dN2[i])-fTmp[i]*(dN2P[i]*dsdxi)),
				Mat3x3(dN2P[i]*dsdxi));
		
		AzTmp[i] = DRef*AzTmp[i]*dCoef;
		
		AzTmp[i] += ERef*Mat6x6(Mat3x3(OmegaRef*(-dN2P[i]*dsdxi*dCoef)),
				Zero3x3, 
				(Mat3x3(LPrimeRef)-Mat3x3(Omega,LRef))
				*(dN2[i]*dCoef)
				+Mat3x3(Omega, fTmp[i]*(dN2P[i]*dsdxi*dCoef))
				+Mat3x3(fTmp[i].Cross(pNode[i]->GetWRef()
						*(dN2P[i]*dsdxi*dCoef))),
				Mat3x3(OmegaRef*(-dN2P[i]*dsdxi*dCoef)));
		
		AzPrimeTmp[i] = ERef*AzPrimeTmp[i];
		
		/* Correggo per la rotazione da locale a globale */
		AzTmp[i].SubMat12(Mat3x3(AzRef.GetVec1()*(dN2[i]*dCoef)));
		AzTmp[i].SubMat22(Mat3x3(AzRef.GetVec2()*(dN2[i]*dCoef)));
	}
	
	Vec3 bTmp[2];
	
	bTmp[0] = p-pNode[NODE1]->GetXCurr();
	bTmp[1] = p-pNode[NODE2]->GetXCurr();
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Equazione all'indietro: */
		WMA.Sub(1, 6*i+1, AzTmp[i].GetMat11());
		WMA.Sub(1, 6*i+4, AzTmp[i].GetMat12());
		
		WMA.Sub(4, 6*i+1,
				AzTmp[i].GetMat21()
				-Mat3x3(AzRef.GetVec1()*(dCoef*dN2[i]))
				+Mat3x3(bTmp[0])*AzTmp[i].GetMat11());
		WMA.Sub(4, 6*i+4, 
				AzTmp[i].GetMat22()
				-Mat3x3(AzRef.GetVec1()*(-dCoef*dN2[i]),
					fTmp[i])
				+Mat3x3(bTmp[0])*AzTmp[i].GetMat12());
		
		/* Equazione in avanti: */
		WMA.Add(7, 6*i+1, AzTmp[i].GetMat11());
		WMA.Add(7, 6*i+4, AzTmp[i].GetMat12());
		
		WMA.Add(10, 6*i+1,
				AzTmp[i].GetMat21()
				-Mat3x3(AzRef.GetVec1()*(dCoef*dN2[i]))
				+Mat3x3(bTmp[1])*AzTmp[i].GetMat11());
		WMA.Add(10, 6*i+4, 
				AzTmp[i].GetMat22()
				+Mat3x3(AzRef.GetVec1()*(dCoef*dN2[i]),
					fTmp[i])
				+Mat3x3(bTmp[1])*AzTmp[i].GetMat12());
		
		/* Equazione viscosa all'indietro: */
		WMB.Sub(1, 6*i+1, AzPrimeTmp[i].GetMat11());
		WMB.Sub(1, 6*i+4, AzPrimeTmp[i].GetMat12());
		
		WMB.Sub(4, 6*i+1,
				AzPrimeTmp[i].GetMat21()
				+Mat3x3(bTmp[0])*AzPrimeTmp[i].GetMat11());
		WMB.Sub(4, 6*i+4,
				AzPrimeTmp[i].GetMat22()
				+Mat3x3(bTmp[0])*AzPrimeTmp[i].GetMat12());
		
		/* Equazione viscosa in avanti: */
		WMB.Add(7, 6*i+1, AzPrimeTmp[i].GetMat11());
		WMB.Add(7, 6*i+4, AzPrimeTmp[i].GetMat12());
		
		WMB.Add(10, 6*i+1,
				AzPrimeTmp[i].GetMat21()	       
				+Mat3x3(bTmp[1])*AzPrimeTmp[i].GetMat11());
		WMB.Add(10, 6*i+4, 
				AzPrimeTmp[i].GetMat22()	     
				+Mat3x3(bTmp[1])*AzPrimeTmp[i].GetMat12());
	}
	
	/* correzione alle equazioni */
	WMA.Add(4, 1, Mat3x3(AzRef.GetVec1()*(-dCoef)));
	WMA.Add(10, 7, Mat3x3(AzRef.GetVec1()*dCoef));
};


/* Assembla il residuo */
void
ViscoElasticHBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticHBeam::AssStiffnessVec");
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei due nodi
	 */
	
	/*
	 * Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
	 * (ora e' superato)
	 */
	
	/* Recupera i dati dei nodi */
	Vec3 xNod[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		xNod[i] = pNode[i]->GetXCurr();
	}
	
	if (bFirstRes) {
		bFirstRes = false; /* AfterPredict ha gia' calcolato tutto */

	} else {
		Vec3 gNod[NUMNODES];    
		Vec3 xTmp[NUMNODES];
		
		Vec3 gPrimeNod[NUMNODES];    
		Vec3 xPrimeTmp[NUMNODES];
		
		for (unsigned int i = 0; i < NUMNODES; i++) {      
			gNod[i] = pNode[i]->GetgCurr();
			Vec3 fTmp = pNode[i]->GetRCurr()*f[i];
			xTmp[i] = xNod[i]+fTmp;
			gPrimeNod[i] = pNode[i]->GetgPCurr();
			xPrimeTmp[i] = pNode[i]->GetVCurr()
				+pNode[i]->GetWCurr().Cross(fTmp);
		}
		
		Mat3x3 RDelta;
		Vec3 gGrad;
		Vec3 gPrimeGrad;
		
		/*
		 * Aggiorna le grandezze della trave nel punto di valutazione
		 */
		
		/* Posizione */
		p = InterpState(xTmp[NODE1], xTmp[NODE2]);
		
		/* Matrici di rotazione */
		g = InterpState(gNod[NODE1], gNod[NODE2]);
		RDelta = Mat3x3(MatR, g);
		R = RDelta*RRef;
		
		/* Velocita' angolare della sezione */	 
		gPrime = InterpState(gPrimeNod[NODE1], gPrimeNod[NODE2]);
		Omega = Mat3x3(MatG, g)*gPrime+RDelta*OmegaRef;
		
		/* Derivate della posizione */
		L = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
		
		/* Derivate della velocita' */
		LPrime = InterpDeriv(xPrimeTmp[NODE1], xPrimeTmp[NODE2]);
		
		/* Derivate dei parametri di rotazione */
		gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
		
		/*
		 * Derivate delle derivate spaziali dei parametri di rotazione
		 */
		gPrimeGrad = InterpDeriv(gPrimeNod[NODE1], gPrimeNod[NODE2]);
		
		/* Per le deformazioni nel sistema del materiale */
		Mat3x3 RTmp(R.Transpose());
		
		/* 
		 * Calcola le deformazioni nel sistema locale nel punto
		 * di valutazione
		 */
		DefLoc = Vec6(RTmp*L-L0, RTmp*(Mat3x3(MatG, g)*gGrad)
				+DefLocRef.GetVec2());
		
		/*
		 * Calcola le velocita' di deformazione nel sistema locale
		 * nel punto di valutazione
		 */
		DefPrimeLoc = Vec6(RTmp*(LPrime+L.Cross(Omega)),
				RTmp*(Mat3x3(MatG, g)*gPrimeGrad
				+(Mat3x3(MatG, g)*gGrad).Cross(Omega))
				+DefPrimeLocRef.GetVec2());
		
		/* Calcola le azioni interne */
		pD->Update(DefLoc, DefPrimeLoc);
		AzLoc = pD->GetF();
		
		/* corregge le azioni interne locali (piezo, ecc) */
		AddInternalForces(AzLoc);
		
		/* Porta le azioni interne nel sistema globale */
		Az = MultRV(AzLoc, R);
	}
	
	WorkVec.Add(1, Az.GetVec1());
	WorkVec.Add(4, (p-xNod[NODE1]).Cross(Az.GetVec1())+Az.GetVec2());
	WorkVec.Sub(7, Az.GetVec1());
	WorkVec.Sub(10, Az.GetVec2()+(p-xNod[NODE2]).Cross(Az.GetVec1()));
}


/* Settings iniziali, prima della prima soluzione */
void
ViscoElasticHBeam::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	HBeam::SetValue(pDM, X, XP. ph);
	
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	(Vec3&)OmegaRef = Omega;
	(Vec3&)LPrimeRef = LPrime;
	(Vec6&)DefPrimeLocRef = DefPrimeLoc;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	(Mat6x6&)ERef = MultRMRt(pD->GetFDEPrime(), RRef);
	
	ASSERT(bFirstRes == true);
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void
ViscoElasticHBeam::AfterPredict(VectorHandler& /* X */ , 
		VectorHandler& /* XP */ )
{
	/*
	 * Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE
	 */
	
	/* Recupera i dati dei nodi */  
	Vec3   gNod[NUMNODES];
	Vec3   xTmp[NUMNODES];
	
	Vec3   gPrimeNod[NUMNODES];
	Vec3   xPrimeTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {            
		gNod[i] = pNode[i]->GetgRef();
		Vec3 fTmp = pNode[i]->GetRRef()*f[i];
		xTmp[i] = pNode[i]->GetXCurr()+fTmp;
		gPrimeNod[i] = pNode[i]->GetgPRef();
		xPrimeTmp[i] = pNode[i]->GetVCurr()
			+pNode[i]->GetWRef().Cross(fTmp);
	}
	
	Mat3x3 RDelta;
	Vec3 gGrad;
	Vec3 gPrimeGrad;
	
	/* Aggiorna le grandezze della trave nel punto di valutazione */
	/* Posizione */
	p = InterpState(xTmp[NODE1], xTmp[NODE2]);
	
	/* Matrici di rotazione */
	g = InterpState(gNod[NODE1], gNod[NODE2]);
	RDelta = Mat3x3(MatR, g);
	R = RRef = RDelta*R;
	
	/* Velocita' angolare della sezione */	 
	gPrime = InterpState(gPrimeNod[NODE1], gPrimeNod[NODE2]);
	Omega = OmegaRef = Mat3x3(MatG, g)*gPrime;
	
	/* Derivate della posizione */
	L = LRef = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	
	/* Derivate della velocita' */
	LPrime = LPrimeRef = InterpDeriv(xPrimeTmp[NODE1], xPrimeTmp[NODE2]);
	
	/* Derivate dei parametri di rotazione */
	gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
	
	/* Derivate delle derivate spaziali dei parametri di rotazione */
	gPrimeGrad = InterpDeriv(gPrimeNod[NODE1], gPrimeNod[NODE2]);
	
	/* Per le deformazioni nel sistema del materiale */
	Mat3x3 RTmp(R.Transpose());
	
	/*
	 * Calcola le deformazioni nel sistema locale nel punto di valutazione
	 */
	DefLoc = DefLocRef = Vec6(RTmp*L-L0,
			RTmp*(Mat3x3(MatG, g)*gGrad)+DefLoc.GetVec2());
	
	/*
	 * Calcola le velocita' di deformazione nel sistema locale
	 * nel punto di valutazione
	 */
	DefPrimeLoc = DefPrimeLocRef = Vec6(RTmp*(LPrime+L.Cross(Omega)),
			RTmp*(Mat3x3(MatG, g)*gPrimeGrad
			+(Mat3x3(MatG, g)*gGrad).Cross(Omega)));
	
	/* Calcola le azioni interne */
	pD->Update(DefLoc, DefPrimeLoc);
	AzLoc = pD->GetF();
	
	/* corregge le azioni interne locali (piezo, ecc) */
	AddInternalForces(AzLoc);
	
	/* Porta le azioni interne nel sistema globale */
	Az = AzRef = MultRV(AzLoc, R);
	
	/* Aggiorna il legame costitutivo */
	DRef = MultRMRt(pD->GetFDE(), R);
	ERef = MultRMRt(pD->GetFDEPrime(), R);
	
	bFirstRes = true;
}

/* ViscoElasticBeam - end */
#endif /* VISCOELASTIC_BEAM */


/* Legge una trave */
Elem*
ReadHBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadHBeam");
	
	const char* sKeyWords[] = {
		"piezoelectric",
		NULL
	};
	
	/* enum delle parole chiave */
	enum KeyWords {
		UNKNOWN = -1,

		PIEZOELECTRIC = 0,

		LASTKEYWORD 
	};
	
	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);
	
	/* Nodo 1 */
	const StructNode* pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
	Mat3x3 R1;
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		R1 = HP.GetRotRel(ReferenceFrame(pNode1));
	} else {
		R1 = pNode1->GetRCurr();
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): " 
			<< f1 << std::endl << "(global frame): "
			<< pNode1->GetXCurr()+pNode1->GetRCurr()*R1*f1 << std::endl);
	
	/* Nodo 2 */
	const StructNode* pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
	Mat3x3 R2;
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		R2 = HP.GetRotRel(ReferenceFrame(pNode2));
	} else {
		R2 = pNode1->GetRCurr();
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "node 2 offset (node reference frame): " 
			<< f2 << std::endl << "(global frame): "
			<< pNode2->GetXCurr()+pNode2->GetRCurr()*R2*f2 << std::endl);
	
	/* Legame costitutivo */
	ConstLawType::Type CLType = ConstLawType::UNKNOWN;
	ConstitutiveLaw6D* pD = HP.GetConstLaw6D(CLType);
	
	if (pD->iGetNumDof() != 0) {
     		silent_cerr("line " << HP.GetLineData()
			<< ": HBeam(" << uLabel << ") does not support "
			"dynamic constitutive laws yet"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
#ifdef DEBUG   
	Mat6x6 MTmp(pD->GetFDE());
	Mat3x3 D11(MTmp.GetMat11());
	Mat3x3 D12(MTmp.GetMat12());
	Mat3x3 D21(MTmp.GetMat21());
	Mat3x3 D22(MTmp.GetMat22());
	
	DEBUGLCOUT(MYDEBUG_INPUT, 
			"First point matrix D11: " << std::endl << D11 << std::endl
			<< "First point matrix D12: " << std::endl << D12 << std::endl
			<< "First point matrix D21: " << std::endl << D21 << std::endl
			<< "First point matrix D22: " << std::endl << D22 << std::endl);
#endif /* DEBUG */

#ifdef PIEZO_BEAM
	flag fPiezo(0);
	Mat3xN PiezoMat[2];
	integer iNumElec = 0;
	ScalarDifferentialNode** pvElecDofs = NULL;
	if (HP.IsKeyWord("piezoelectricactuator")) {
		fPiezo = flag(1);
		DEBUGLCOUT(MYDEBUG_INPUT, 
				"Piezoelectric actuator beam is expected"
				<< std::endl);
		
		iNumElec = HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT, 
				"piezo actuator " << uLabel
				<< " has " << iNumElec 
				<< " electrodes" << std::endl);
		if (iNumElec <= 0) {
			silent_cerr("HBeam(" << uLabel << "): "
				"illegal number of electric nodes "
				<< iNumElec 
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		SAFENEWARR(pvElecDofs, ScalarDifferentialNode*, iNumElec);
		
		for (integer i = 0; i < iNumElec; i++) {
			unsigned int uL = HP.GetInt();
			DEBUGLCOUT(MYDEBUG_INPUT, "linked to abstract node "
					<< uL << std::endl);
			pvElecDofs[i] = (ScalarDifferentialNode*)(pDM->pFindNode(Node::ABSTRACT, uL));
			if (pvElecDofs[i] == NULL) {
				silent_cerr("HBeam(" << uLabel << "): "
					"can't find AbstractNode(" << uL << ") "
					"at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
		
		PiezoMat[0].Resize(iNumElec);
		PiezoMat[1].Resize(iNumElec);
		
		/* leggere le matrici (6xN sez. 1, 6xN sez. 2) */
		HP.GetMat6xN(PiezoMat[0], PiezoMat[1], iNumElec);
		
#if 0
		DEBUGLCOUT(MYDEBUG_INPUT, "Piezo matrix I:" << std::endl
				<< PiezoMat[0][0] << PiezoMat[1][0]);
#endif /* 0 */
	}
#endif /* PIEZO_BEAM */

	flag fOut = pDM->fReadOutput(HP, Elem::BEAM);       
	
	Elem* pEl = NULL;
	
	if (CLType == ConstLawType::ELASTIC) {
#ifdef PIEZO_BEAM
		if (fPiezo == 0) {	 
#endif /* PIEZO_BEAM */
			SAFENEWWITHCONSTRUCTOR(pEl,
					HBeam,
					HBeam(uLabel,
						pNode1, pNode2,
						f1, f2,
						R1, R2,
						pD,
						fOut));
#ifdef PIEZO_BEAM
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
					PiezoActuatorHBeam,
					PiezoActuatorHBeam(uLabel,
						pNode1, pNode2,
						f1, f2,
						R,
						pD,
						iNumElec,
						pvElecDofs,
						PiezoMat[0], PiezoMat[1],
						fOut));
		}
#endif /* PIEZO_BEAM */
		
	} else /* At least one is VISCOUS or VISCOELASTIC */ {
#ifdef VISCOELASTIC_BEAM      
#ifdef PIEZO_BEAM
		if (fPiezo == 0) {	 
#endif /* PIEZO_BEAM */
			SAFENEWWITHCONSTRUCTOR(pEl, 
					ViscoElasticHBeam,
					ViscoElasticHBeam(uLabel,
						pNode1, pNode2,
						f1, f2,
						R,
						pD,
						fOut));
#ifdef PIEZO_BEAM
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
					PiezoActuatorVEHBeam,
					PiezoActuatorVEHBeam(uLabel,
						pNode1, pNode2,
						f1, f2,
						R,
						pD,
						iNumElec,
						pvElecDofs,
						PiezoMat[0], PiezoMat[1],
						fOut));
		}
#endif /* PIEZO_BEAM */

#else /* VISCOELASTIC_BEAM */
		silent_cerr("Sorry, the ViscoElasticHBeam element"
			" is not available yet" << std::endl);
		throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
#endif /* VISCOELASTIC_BEAM */
	}
	
	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	return pEl;
} /* End of ReadHBeam() */

