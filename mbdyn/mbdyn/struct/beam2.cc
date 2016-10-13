/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
#include "beam2.h"
#include "pzbeam2.h"
#include "Rot.hh"

/*
 * Nota: non e' ancora stato implementato il contributo 
 * della ViscoElasticBeam2 all'assemblaggio iniziale
 */

/*
 * Nota: la parte viscoelastica va rivista in accordo con la piu' 
 * recente formulazione delle derivate delle deformazioni nel sistema
 * materiale
 */

/* Beam2 - begin */

/* Costruttore normale */
Beam2::Beam2(unsigned int uL, 
		const StructNode* pN1,
		const StructNode* pN2,
		const Vec3& F1,
		const Vec3& F2,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& r,
		const ConstitutiveLaw6D* pd,
		OrientationDescription ood,
		flag fOut)
: Elem(uL, fOut), 
ElemGravityOwner(uL, fOut), 
InitialAssemblyElem(uL, fOut),
od(ood),
#ifdef USE_NETCDF
Var_X(0),
Var_Phi(0),
Var_F(0),
Var_M(0),
Var_Nu(0),
Var_K(0),
Var_NuP(0),
Var_KP(0),
#endif /* USE_NETCDF */
bFirstRes(false),
bFirstIDRes(true)
{
	/* Validazione dati */
	ASSERT(pN1 != NULL);
	ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN2 != NULL);
	ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
 
	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	const_cast<Vec3&>(f[NODE1]) = F1;
	const_cast<Vec3&>(f[NODE2]) = F2;
	const_cast<Mat3x3&>(RNode[NODE1]) = R1;
	const_cast<Mat3x3&>(RNode[NODE2]) = R2;
	RPrev = RRef = R = (Mat3x3&)r;
	
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
	DefLocPrev = Zero6;
	p = Zero3;
	g = Zero3;
	L0 = Zero3;
	L = Zero3;
	
	DsDxi();
	
	Vec3 xTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {      
		xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRCurr()*f[i];
	}      
	
	/* Aggiorna le grandezze della trave nel punto di valutazione */
	p = InterpState(xTmp[NODE1], xTmp[NODE2]);
}


Beam2::~Beam2(void) 
{
	ASSERT(pD != NULL);
	if (pD != NULL) {      
		SAFEDELETE(pD);
	}
}

/* Accesso ai dati privati */
unsigned int
Beam2::iGetNumPrivData(void) const
{
	return Beam::iNumPrivData;
}

unsigned int
Beam2::iGetPrivDataIdx(const char *s) const
{
	ConstLawType::Type type = ConstLawType::ELASTIC;
	if (dynamic_cast<const ViscoElasticBeam2 *>(this)) {
		type = ConstLawType::VISCOUS;
	}

	return Beam::iGetPrivDataIdx_int(s, type);
}

doublereal
Beam2::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= iGetNumPrivData());

	switch (i) {
	case 1:
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
		return DefLoc.dGet(i);

	case 7:
	case 8:
	case 9:
	case 10:
	case 11:
	case 12:
		return AzLoc.dGet(i - 7);

	case 13:
	case 14:
	case 15:
		return p.dGet(i - 12);

	case 16:
	case 17:
	case 18:
		return RotManip::VecRot(R).dGet(i - 15);

	case 19:
	case 20:
	case 21:
		return Omega.dGet(i - 18);

	default:
		silent_cerr("Beam2(" << GetLabel() << "): "
			"illegal private data " << i << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Vec3 
Beam2::InterpState(const Vec3& v1, const Vec3& v2)
{
	doublereal* pv1 = (doublereal*)v1.pGetVec();
	doublereal* pv2 = (doublereal*)v2.pGetVec();
	return Vec3(pv1[0]*dN2[0] + pv2[0]*dN2[1],
			pv1[1]*dN2[0] + pv2[1]*dN2[1],
			pv1[2]*dN2[0] + pv2[2]*dN2[1]);
}


Vec3
Beam2::InterpDeriv(const Vec3& v1, const Vec3& v2)
{
	doublereal* pv1 = (doublereal*)v1.pGetVec();
	doublereal* pv2 = (doublereal*)v2.pGetVec();
	return Vec3((pv1[0]*dN2P[0] + pv2[0]*dN2P[1])*dsdxi,
			(pv1[1]*dN2P[0] + pv2[1]*dN2P[1])*dsdxi,
			(pv1[2]*dN2P[0] + pv2[2]*dN2P[1])*dsdxi);
}


void
Beam2::DsDxi(void)
{
	/* Calcola il ds/dxi e le deformazioni iniziali */
	Vec3 xNod[NUMNODES];
	Mat3x3 RNod[NUMNODES];
	Vec3 xTmp[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {
		xNod[i] = pNode[i]->GetXCurr();
		RNod[i] = pNode[i]->GetRCurr();
		xTmp[i] = xNod[i] + RNod[i]*f[i];
	}
	
	dsdxi = 1.;

	Vec3 xGrad = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	doublereal d = xGrad.Dot();
	if (d > std::numeric_limits<doublereal>::epsilon()) {
		dsdxi = 1./std::sqrt(d);
	} else {
		silent_cerr("warning, Beam2(" << GetLabel() << ") "
			"has singular metric; aborting..." << std::endl);
		
		throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
	}

	/* Calcola le deformazioni iniziali */
	L0 = R.MulTV(InterpDeriv(xTmp[NODE1], xTmp[NODE2]));
	pD->Update(Zero6);
	DRef = MultRMRt(pD->GetFDE(), R);
}


/* Calcola la velocita' angolare delle sezioni a partire da quelle dei nodi */
void 
Beam2::Omega0(void)
{   
	/* Modo consistente: */      
	Mat3x3 RNod[NUMNODES];
	Vec3 w[NUMNODES];
	for (unsigned int i = 0; i < NUMNODES; i++) {     
		RNod[i] = pNode[i]->GetRCurr()*RNode[i];
		w[i] = pNode[i]->GetWCurr();
	}
	
	/*
	 * Calcolo i parametri di rotazione della rotazione relativa
	 * tra inizio e fine e li dimezzo nell'ipotesi che siano limitati
	 */
	Vec3 gTmp(CGR_Rot::Param, RNod[NODE2].MulTM(RNod[NODE1]));
	
	/*
	 * Le derivate dei parametri di rotazione si ricavano da omega
	 */
	Vec3 g1P(Mat3x3(CGR_Rot::MatGm1, gTmp*(-.5))*w[NODE1]);
	Vec3 g2P(Mat3x3(CGR_Rot::MatGm1, gTmp*.5)*w[NODE2]);

        Vec3 gPTmp(g1P*dN2[NODE1] + g2P*dN2[NODE2]);
        Omega = Mat3x3(CGR_Rot::MatG, gTmp)*gPTmp;
	
#if 0
	/* Modo brutale: interpolo le velocita' dei nodi */
	Omega = pNode[NODE1]->GetWCurr()*dN2[NODE1]
		+ pNode[NODE2]->GetWCurr()*dN2[NODE2];
#endif /* 0 */
}


/* Contributo al file di restart */
std::ostream&
Beam2::Restart(std::ostream& out) const
{
	return Restart_(out)<< ';' << std::endl;
}

std::ostream&
Beam2::Restart_(std::ostream& out) const
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
Beam2::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	RPrev = R;
	DefLocPrev = DefLoc;
	pD->AfterConvergence(DefLoc);
}

/* Assembla la matrice */
void
Beam2::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& /* WMB */ ,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("Beam2::AssStiffnessMat");
	
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
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Delta - deformazioni */
		AzTmp[i] = Mat6x6(mb_deye<Mat3x3>(dN2P[i]*dsdxi*dCoef),
				Zero3x3,
				Mat3x3(MatCross, L*(dN2[i]*dCoef) -fTmp[i]*(dN2P[i]*dsdxi*dCoef)),
				mb_deye<Mat3x3>(dN2P[i]*dsdxi*dCoef));
		
		/* Delta - azioni interne */
		AzTmp[i] = DRef*AzTmp[i];
		
		/* Correggo per la rotazione da locale a globale */
		AzTmp[i].SubMat12(Mat3x3(MatCross, Az.GetVec1()*(dN2[i]*dCoef)));
		AzTmp[i].SubMat22(Mat3x3(MatCross, Az.GetVec2()*(dN2[i]*dCoef)));
	}
   
	Vec3 bTmp[NUMNODES];
	
	bTmp[NODE1] = p - pNode[NODE1]->GetXCurr();
	bTmp[NODE2] = p - pNode[NODE2]->GetXCurr();
   
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Equazione all'indietro: */
		WMA.Sub(1, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Sub(1, 6*i + 4, AzTmp[i].GetMat12());
		
		WMA.Sub(3 + 1, 6*i + 1,
				AzTmp[i].GetMat21()
				- Mat3x3(MatCross, Az.GetVec1()*(dCoef*dN2[i]))
				+ bTmp[NODE1].Cross(AzTmp[i].GetMat11()));
		WMA.Sub(3 + 1, 6*i + 4, 
				AzTmp[i].GetMat22()
				- Mat3x3(MatCrossCross, Az.GetVec1()*(-dCoef*dN2[i]), fTmp[i])
				+ bTmp[NODE1].Cross(AzTmp[i].GetMat12()));
		
		/* Equazione in avanti: */
		WMA.Add(6 + 1, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Add(6 + 1, 6*i + 4, AzTmp[i].GetMat12());
		
		WMA.Add(9 + 1, 6*i + 1,
				AzTmp[i].GetMat21()
				- Mat3x3(MatCross, Az.GetVec1()*(dCoef*dN2[i]))
				+ bTmp[NODE2].Cross(AzTmp[i].GetMat11()));
		WMA.Add(9 + 1, 6*i + 4, 
				AzTmp[i].GetMat22()
				+ Mat3x3(MatCrossCross, Az.GetVec1()*(dCoef*dN2[i]), fTmp[i])
				+ bTmp[NODE2].Cross(AzTmp[i].GetMat12()));
	}
	
	/* correzione alle equazioni */
	Mat3x3 FTmp(MatCross, Az.GetVec1()*dCoef);
	WMA.Sub(3 + 1, 1, FTmp);
	WMA.Add(9 + 1, 6 + 1, FTmp);
};


/* Assembla il residuo */
void
Beam2::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("Beam2::AssStiffnessVec");
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei tre nodi
	 */
	
	/*
	 * Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
	 * (ora e' superato)
	 */
	
	if (bFirstRes) {
		bFirstRes = false; /* AfterPredict ha gia' calcolato tutto */

	} else {
		Vec3 gNod[NUMNODES];    
		Vec3 xTmp[NUMNODES];
		
		for (unsigned int i = 0; i < NUMNODES; i++) {      
			gNod[i] = pNode[i]->GetgCurr();	 
			xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRCurr()*f[i];
		}      
		
		Mat3x3 RDelta;
		Vec3 gGrad;
		
		/*
		 * Aggiorna le grandezze della trave nel punto di valutazione
		 */
		
		/* Posizione */
		p = InterpState(xTmp[NODE1], xTmp[NODE2]);
		
		/* Matrici di rotazione */
		g = InterpState(gNod[NODE1], gNod[NODE2]);
		RDelta = Mat3x3(CGR_Rot::MatR, g);
		R = RDelta*RRef;
		
		/* Derivate della posizione */
		L = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
		
		/* Derivate dei parametri di rotazione */
		gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
		
		/*
		 * Calcola le deformazioni nel sistema locale
		 * nei punti di valutazione
		 */
		DefLoc = Vec6(R.MulTV(L) - L0,
			R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gGrad) + DefLocRef.GetVec2());
		
		/* Calcola le azioni interne */
		pD->Update(DefLoc);
		AzLoc = pD->GetF();
		
		/* corregge le azioni interne locali (piezo, ecc) */
		AddInternalForces(AzLoc);
		
		/* Porta le azioni interne nel sistema globale */
		Az = MultRV(AzLoc, R);
	}

	WorkVec.Add(1, Az.GetVec1());
	WorkVec.Add(4, (p - pNode[NODE1]->GetXCurr()).Cross(Az.GetVec1()) + Az.GetVec2());
	WorkVec.Sub(7, Az.GetVec1());
	WorkVec.Sub(10, Az.GetVec2() + (p - pNode[NODE2]->GetXCurr()).Cross(Az.GetVec1()));
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler&
Beam2::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("Beam2::AssJac => AssStiffnessMat");
	
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
Beam2::AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("Beam2::AssRes => AssStiffnessVec");
	
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
Beam2::SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph)
{
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	RRef = R;
	LRef = L;
	DefLocRef = DefLoc;
	AzRef = Az;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	DRef = MultRMRt(pD->GetFDE(), RRef);      
	
	bFirstRes = true;
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void
Beam2::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{  
	/*
	 * Calcola le deformazioni, aggiorna il legame costitutivo
	 * e crea la FDE
	 */
	
	/* Recupera i dati dei nodi */  
	Vec3   gNod[NUMNODES];
	Vec3   xTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {            
		gNod[i] = pNode[i]->GetgRef();
		xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRRef()*f[i];
	}
	
	Mat3x3 RDelta;
	Vec3 gGrad;
	
	/* Aggiorna le grandezze della trave nel punto di valutazione */
	
	/* Posizione */
	p = InterpState(xTmp[NODE1], xTmp[NODE2]);
	
	/* Matrici di rotazione */
	g = InterpState(gNod[NODE1], gNod[NODE2]);
	RDelta = Mat3x3(CGR_Rot::MatR, g);
	R = RRef = RDelta*RPrev;
	
	/* Derivate della posizione */
	L = LRef = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	
	/* Derivate dei parametri di rotazione */
	gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
	
	/*
	 * Calcola le deformazioni nel sistema locale
	 * nei punti di valutazione
	 */
	DefLoc = DefLocRef = Vec6(R.MulTV(L) - L0,
		R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gGrad) + DefLocPrev.GetVec2());
	
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
}


void
Beam2::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::BEAMS)) {
			ASSERT(OH.IsOpen(OutputHandler::NETCDF));

			const char *type = 0;
			switch (GetBeamType()) {
			case Beam::ELASTIC:
				type = "elastic";
				break;

			case Beam::VISCOELASTIC:
				type = "viscoelastic";
				break;

			case Beam::PIEZOELECTRICELASTIC:
				type = "piezoelectric elastic";
				break;

			case Beam::PIEZOELECTRICVISCOELASTIC:
				type = "piezoelectric viscoelastic";
				break;

			default:
				type = "unknown";
				break;
			}

			std::ostringstream os;
			os << "elem.beam." << GetLabel();

			(void)OH.CreateVar(os.str(), type);

			os << '.';
			std::string name(os.str());

			unsigned uOutputFlags = (fToBeOutput() & ToBeOutput::OUTPUT_PRIVATE_MASK);

			if (uOutputFlags & Beam::OUTPUT_EP_X) {
				Var_X = OH.CreateVar<Vec3>(name + "X", "m",
					"evaluation point global position vector (X, Y, Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_R) {
				Var_Phi = OH.CreateRotationVar(name, "", od,
					" evaluation point global orientation matrix");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_F) {
				Var_F = OH.CreateVar<Vec3>(name + "F", "N",
					"evaluation point internal force in local frame (F_X, F_Y, F_Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_M) {
				Var_M = OH.CreateVar<Vec3>(name + "M", "Nm",
					"evaluation point internal force in local frame (M_X, M_Y, M_Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_NU) {
				Var_Nu = OH.CreateVar<Vec3>(name + "nu", "-",
					"evaluation point linear strain in local frame (nu_X, nu_Y, nu_Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_K) {
				Var_K = OH.CreateVar<Vec3>(name + "k", "1/m",
					"evaluation point angular strain in local frame (K_X, K_Y, K_Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_NUP) {
				Var_NuP = OH.CreateVar<Vec3>(name + "nuP", "1/s",
					"evaluation point linear strain rate in local frame (nuP_X, nuP_Y, nuP_Z)");
			}

			if (uOutputFlags & Beam::OUTPUT_EP_KP) {
				Var_KP = OH.CreateVar<Vec3>(name + "kP", "1/ms",
					"evaluation point angular strain rate in local frame (KP_X, KP_Y, KP_Z)");
			}
		}
#endif // USE_NETCDF
	}
}

/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void
Beam2::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::BEAMS)) {
			if (Var_X) {
				Var_X->put_rec(p.pGetVec(), OH.GetCurrentStep());
			}

			if (Var_Phi) {
				Vec3 E;
				switch (od) {
				case EULER_123:
					E = MatR2EulerAngles123(R)*dRaDegr;
					break;

				case EULER_313:
					E = MatR2EulerAngles313(R)*dRaDegr;
					break;

				case EULER_321:
					E = MatR2EulerAngles321(R)*dRaDegr;
					break;

				case ORIENTATION_VECTOR:
					E = RotManip::VecRot(R);
					break;

				case ORIENTATION_MATRIX:
					break;

				default:
					/* impossible */
					break;
				}

				switch (od) {
				case EULER_123:
				case EULER_313:
				case EULER_321:
				case ORIENTATION_VECTOR:
					Var_Phi->put_rec(E.pGetVec(), OH.GetCurrentStep());
					break;

				case ORIENTATION_MATRIX:
					Var_Phi->put_rec(R.pGetMat(), OH.GetCurrentStep());
					break;

				default:
					/* impossible */
					break;
				}
			}

			if (Var_F) {
				Var_F->put_rec(AzLoc.GetVec1().pGetVec(), OH.GetCurrentStep());
			}

			if (Var_M) {
				Var_M->put_rec(AzLoc.GetVec2().pGetVec(), OH.GetCurrentStep());
			}

			if (Var_Nu) {
				Var_Nu->put_rec(DefLoc.GetVec1().pGetVec(), OH.GetCurrentStep());
			}

			if (Var_K) {
				Var_K->put_rec(DefLoc.GetVec2().pGetVec(), OH.GetCurrentStep());
			}

			if (Var_NuP) {
				Var_NuP->put_rec(DefPrimeLoc.GetVec1().pGetVec(), OH.GetCurrentStep());
			}

			if (Var_KP) {
				Var_KP->put_rec(DefPrimeLoc.GetVec2().pGetVec(), OH.GetCurrentStep());
			}
		}
#endif /* USE_NETCDF */

		if (OH.UseText(OutputHandler::BEAMS)) {
			OH.Beams() << std::setw(8) << GetLabel()
				<< " " << AzLoc.GetVec1()
				<< " " << AzLoc.GetVec2()
				<< std::endl;
		}
	}
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
Beam2::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) 
{ 
	DEBUGCOUTFNAME("Beam2::InitialAssJac => AssStiffnessMat");
	
	/* Dimensiona la matrice, la azzera e pone gli indici corretti */
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeReset(12, 12);
	
	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
		WM.PutColIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}      
	
	AssStiffnessMat(WM, WM, 1., XCurr, XCurr);
	return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler&
Beam2::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr) 
{ 
	DEBUGCOUTFNAME("Beam2::InitialAssRes => AssStiffnessVec");

	/* Dimensiona il vettore, lo azzera e pone gli indici corretti */
	WorkVec.ResizeReset(12);

	integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
	integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iNode2FirstPosIndex + iCnt);
	}      
	
	AssStiffnessVec(WorkVec, 1., XCurr, XCurr);
	return WorkVec;
}


const StructNode*
Beam2::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 2);
	switch (i) {
	case 1:
	case 2:
		return pNode[i-1];
	default:
		throw Beam2::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}


void
Beam2::GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& r) const
{
	ASSERT(part == 1);
	part--;
   
	x = p;
	r = R;
}

void
Beam2::GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const
{
	ASSERT(part == 1);
	part--;
   
	v = Zero3;
	w = Zero3;
}

/* inverse dynamics capable element */
bool
Beam2::bInverseDynamics(void) const
{
	return true;
}


/* Inverse Dynamics Jacobian matrix assembly */
VariableSubMatrixHandler&
Beam2::AssJac(VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// silent_cout("Beam2::IDAssJac()" << std::endl);

	DEBUGCOUT("Entering Beam2::[InverseDynamics]AssJac()" << std::endl);

#if 0
	// iOrder not available
	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS
		|| (iOrder == InverseDynamics::POSITION && bIsErgonomy()));
#endif

	return AssJac(WorkMat, 1., XCurr, XCurr);
}

/* Inverse Dynamics Residual Assembly */
SubVectorHandler&
Beam2::AssRes(SubVectorHandler& WorkVec,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr, 
	const VectorHandler& /* XPrimePrimeCurr */ ,
	InverseDynamics::Order iOrder)
{
	// silent_cout("Beam2::IDAssRes(" << iOrder << ")" << std::endl);

	DEBUGCOUT("Entering Beam2::[InverseDynamics]AssRes()" << std::endl);

	ASSERT(iOrder == InverseDynamics::INVERSE_DYNAMICS
		|| (iOrder == InverseDynamics::POSITION && bIsErgonomy()));

	// if (iOrder == InverseDynamics::POSITION) {
		// ASSERT(bIsErgonomy());

		if (bFirstIDRes) {
			// prepare for new step
			AfterPredict(const_cast<VectorHandler&>(XCurr),
				const_cast<VectorHandler&>(XPrimeCurr));
	
			bFirstIDRes = false;
			bFirstRes = true;
		}
	// }
	
	return AssRes(WorkVec, 1., XCurr, XPrimeCurr);
}

/* Inverse Dynamics update */
void
Beam2::Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder)
{
	NO_OP;
}

/* Inverse Dynamics after convergence */
void
Beam2::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP, const VectorHandler& XPP)
{
	AfterConvergence(X, XP);
	bFirstIDRes = true;
}

#ifdef USE_ADAMS
std::ostream& 
Beam2::WriteAdamsDummyPartCmd(std::ostream& out,
		unsigned int part,
		unsigned int firstId) const
{
	Vec3 xTmp[NUMNODES];
	
	part--;
	
	for (unsigned int i = 0; i <= 1; i++) {
		xTmp[i] = pNode[i]->GetXCurr() + pNode[i]->GetRCurr()*f[i];
	}
	
	out << psAdamsElemCode[GetElemType()] << "_" << GetLabel()
		<< "_" << 1 << std::endl
		<< firstId << " "
		<< p << " " 
		<< MatR2EulerAngles(R)*dRaDegr << " "
		<< R.MulTV(xTmp[NODE1]-p) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[part]->GetRCurr())*dRaDegr */ << " "
		<< R.MulTV(xTmp[NODE2]-p) << " "
		<< Zero3 /* MatR2EulerAngles(pNode[1+part]->GetRCurr())*dRaDegr */ << std::endl;
	
	return out;
}
#endif /* USE_ADAMS */

/* Beam2 - end */


/* ViscoElasticBeam2 - begin */

/* Costruttore normale */
ViscoElasticBeam2::ViscoElasticBeam2(unsigned int uL, 
		const StructNode* pN1, 
		const StructNode* pN2, 
		const Vec3& F1,
		const Vec3& F2,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& r,
		const ConstitutiveLaw6D* pd, 
		OrientationDescription ood,
		flag fOut)
: Elem(uL, fOut),
Beam2(uL, pN1, pN2, F1, F2, R1, R2, r, pd, ood, fOut)
{
	LPrimeRef = LPrime = Zero3;  
	gPrime = Zero3;
	
	DefPrimeLoc = DefPrimeLocRef = Zero6;
	
	/* Nota: DsDxi() viene chiamata dal costruttore di Beam */
	Beam2::Omega0();

	OmegaRef = Omega;
	ERef = MultRMRt(pD->GetFDEPrime(), RRef);
}

void
ViscoElasticBeam2::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	RPrev = R;
	DefLocPrev = DefLoc;
	pD->AfterConvergence(DefLoc, DefPrimeLoc);
}


/* Assembla la matrice */
void
ViscoElasticBeam2::AssStiffnessMat(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticBeam2::AssStiffnessMat");
	
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
		AzTmp[i] = AzPrimeTmp[i] = Mat6x6(mb_deye<Mat3x3>(dN2P[i]*dsdxi),
				Zero3x3,
				Mat3x3(MatCross, L*dN2[i] - fTmp[i]*(dN2P[i]*dsdxi)),
				mb_deye<Mat3x3>(dN2P[i]*dsdxi));
		
		AzTmp[i] = DRef*AzTmp[i]*dCoef;
		
		AzTmp[i] += ERef*Mat6x6(Mat3x3(MatCross, Omega*(-dN2P[i]*dsdxi*dCoef)),
				Zero3x3, 
				(Mat3x3(MatCross, LPrime) - Mat3x3(MatCrossCross, Omega, L))*(dN2[i]*dCoef)
					+ Mat3x3(MatCrossCross, Omega, fTmp[i]*(dN2P[i]*dsdxi*dCoef))
					+ Mat3x3(MatCross, fTmp[i].Cross(pNode[i]->GetWCurr()*(dN2P[i]*dsdxi*dCoef))),
				Mat3x3(MatCross, Omega*(-dN2P[i]*dsdxi*dCoef)));
		
		AzPrimeTmp[i] = ERef*AzPrimeTmp[i];
		
		/* Correggo per la rotazione da locale a globale */
		AzTmp[i].SubMat12(Mat3x3(MatCross, Az.GetVec1()*(dN2[i]*dCoef)));
		AzTmp[i].SubMat22(Mat3x3(MatCross, Az.GetVec2()*(dN2[i]*dCoef)));
	}
	
	Vec3 bTmp[NUMNODES];
	
	bTmp[NODE1] = p-pNode[NODE1]->GetXCurr();
	bTmp[NODE2] = p-pNode[NODE2]->GetXCurr();
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		/* Equazione all'indietro: */
		WMA.Sub(1, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Sub(1, 6*i + 4, AzTmp[i].GetMat12());
		
		WMA.Sub(4, 6*i + 1,
				AzTmp[i].GetMat21()
				- Mat3x3(MatCross, Az.GetVec1()*(dCoef*dN2[i]))
				+ bTmp[NODE1].Cross(AzTmp[i].GetMat11()));
		WMA.Sub(4, 6*i + 4, 
				AzTmp[i].GetMat22()
				- Mat3x3(MatCrossCross, Az.GetVec1()*(-dCoef*dN2[i]), fTmp[i])
				+ bTmp[NODE1].Cross(AzTmp[i].GetMat12()));
		
		/* Equazione in avanti: */
		WMA.Add(7, 6*i + 1, AzTmp[i].GetMat11());
		WMA.Add(7, 6*i + 4, AzTmp[i].GetMat12());
		
		WMA.Add(10, 6*i + 1,
				AzTmp[i].GetMat21()
				- Mat3x3(MatCross, Az.GetVec1()*(dCoef*dN2[i]))
				+ bTmp[NODE2].Cross(AzTmp[i].GetMat11()));
		WMA.Add(10, 6*i + 4, 
				AzTmp[i].GetMat22()
				+ Mat3x3(MatCrossCross, Az.GetVec1()*(dCoef*dN2[i]), fTmp[i])
				+ bTmp[NODE2].Cross(AzTmp[i].GetMat12()));
		
		/* Equazione viscosa all'indietro: */
		WMB.Sub(1, 6*i + 1, AzPrimeTmp[i].GetMat11());
		WMB.Sub(1, 6*i + 4, AzPrimeTmp[i].GetMat12());
		
		WMB.Sub(4, 6*i + 1,
				AzPrimeTmp[i].GetMat21()
				+ bTmp[NODE1].Cross(AzPrimeTmp[i].GetMat11()));
		WMB.Sub(4, 6*i + 4,
				AzPrimeTmp[i].GetMat22()
				+ bTmp[NODE1].Cross(AzPrimeTmp[i].GetMat12()));
		
		/* Equazione viscosa in avanti: */
		WMB.Add(7, 6*i + 1, AzPrimeTmp[i].GetMat11());
		WMB.Add(7, 6*i + 4, AzPrimeTmp[i].GetMat12());
		
		WMB.Add(10, 6*i + 1,
				AzPrimeTmp[i].GetMat21()	       
				+ bTmp[NODE2].Cross(AzPrimeTmp[i].GetMat11()));
		WMB.Add(10, 6*i + 4, 
				AzPrimeTmp[i].GetMat22()	     
				+ bTmp[NODE2].Cross(AzPrimeTmp[i].GetMat12()));
	}
	
	/* correzione alle equazioni */
	Mat3x3 FTmp(MatCross, Az.GetVec1()*dCoef);
	WMA.Sub(4, 1, FTmp);
	WMA.Add(10, 7, FTmp);
};


/* Assembla il residuo */
void
ViscoElasticBeam2::AssStiffnessVec(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("ViscoElasticBeam2::AssStiffnessVec");
	
	/*
	 * Riceve il vettore gia' dimensionato e con gli indici a posto 
	 * per scrivere il residuo delle equazioni di equilibrio dei due nodi
	 */
	
	/*
	 * Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
	 * (ora e' superato)
	 */
	
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
			xTmp[i] = pNode[i]->GetXCurr() + fTmp;
			gPrimeNod[i] = pNode[i]->GetgPCurr();
			xPrimeTmp[i] = pNode[i]->GetVCurr()
				+ pNode[i]->GetWCurr().Cross(fTmp);
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
		RDelta = Mat3x3(CGR_Rot::MatR, g);
		R = RDelta*RRef;
		
		/* Velocita' angolare della sezione */	 
		gPrime = InterpState(gPrimeNod[NODE1], gPrimeNod[NODE2]);
		Omega = Mat3x3(CGR_Rot::MatG, g)*gPrime + RDelta*OmegaRef;
		
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
		
		/* 
		 * Calcola le deformazioni nel sistema locale nel punto
		 * di valutazione
		 */
		DefLoc = Vec6(R.MulTV(L) - L0,
			R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gGrad) + DefLocRef.GetVec2());
		
		/*
		 * Calcola le velocita' di deformazione nel sistema locale
		 * nel punto di valutazione
		 */
		DefPrimeLoc = Vec6(R.MulTV(LPrime + L.Cross(Omega)),
			R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gPrimeGrad
			+ (Mat3x3(CGR_Rot::MatG, g)*gGrad).Cross(Omega))
			+ DefPrimeLocRef.GetVec2());
		
		/* Calcola le azioni interne */
		pD->Update(DefLoc, DefPrimeLoc);
		AzLoc = pD->GetF();

		/* corregge le azioni interne locali (piezo, ecc) */
		AddInternalForces(AzLoc);

		/* Porta le azioni interne nel sistema globale */
		Az = MultRV(AzLoc, R);
	}

	WorkVec.Add(1, Az.GetVec1());
	WorkVec.Add(4, (p - pNode[NODE1]->GetXCurr()).Cross(Az.GetVec1()) + Az.GetVec2());
	WorkVec.Sub(7, Az.GetVec1());
	WorkVec.Sub(10, Az.GetVec2() + (p - pNode[NODE2]->GetXCurr()).Cross(Az.GetVec1()));
}


/* Settings iniziali, prima della prima soluzione */
void
ViscoElasticBeam2::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	Beam2::SetValue(pDM, X, XP, ph);
	
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	OmegaRef = Omega;
	LPrimeRef = LPrime;
	DefPrimeLocRef = DefPrimeLoc;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	ERef = MultRMRt(pD->GetFDEPrime(), RRef);
	
	ASSERT(bFirstRes == true);
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void
ViscoElasticBeam2::AfterPredict(VectorHandler& /* X */ , 
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
		xTmp[i] = pNode[i]->GetXCurr() + fTmp;
		gPrimeNod[i] = pNode[i]->GetgPRef();
		xPrimeTmp[i] = pNode[i]->GetVCurr()
			+ pNode[i]->GetWRef().Cross(fTmp);
	}
	
	Mat3x3 RDelta;
	Vec3 gGrad;
	Vec3 gPrimeGrad;
	
	/* Aggiorna le grandezze della trave nel punto di valutazione */
	/* Posizione */
	p = InterpState(xTmp[NODE1], xTmp[NODE2]);
	
	/* Matrici di rotazione */
	g = InterpState(gNod[NODE1], gNod[NODE2]);
	RDelta = Mat3x3(CGR_Rot::MatR, g);
	R = RRef = RDelta*RPrev;
	
	/* Velocita' angolare della sezione */	 
	gPrime = InterpState(gPrimeNod[NODE1], gPrimeNod[NODE2]);
	Omega = OmegaRef = Mat3x3(CGR_Rot::MatG, g)*gPrime;
	
	/* Derivate della posizione */
	L = LRef = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	
	/* Derivate della velocita' */
	LPrime = LPrimeRef = InterpDeriv(xPrimeTmp[NODE1], xPrimeTmp[NODE2]);
	
	/* Derivate dei parametri di rotazione */
	gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
	
	/* Derivate delle derivate spaziali dei parametri di rotazione */
	gPrimeGrad = InterpDeriv(gPrimeNod[NODE1], gPrimeNod[NODE2]);
	
	/*
	 * Calcola le deformazioni nel sistema locale nel punto di valutazione
	 */
	DefLoc = DefLocRef = Vec6(R.MulTV(L) - L0,
		R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gGrad) + DefLocPrev.GetVec2());
	
	/*
	 * Calcola le velocita' di deformazione nel sistema locale
	 * nel punto di valutazione
	 */
	DefPrimeLoc = DefPrimeLocRef = Vec6(R.MulTV(LPrime + L.Cross(Omega)),
		R.MulTV(Mat3x3(CGR_Rot::MatG, g)*gPrimeGrad
		+ (Mat3x3(CGR_Rot::MatG, g)*gGrad).Cross(Omega)));
	
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

doublereal
ViscoElasticBeam2::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= iGetNumPrivData());

	switch (i) {
	case 22:
	case 23:
	case 24:
	case 25:
	case 26:
	case 27:
		return DefPrimeLoc.dGet(i - 21);

	default:
		return Beam2::dGetPrivData(i);
	}
}

/* ViscoElasticBeam - end */


/* Legge una trave */
Elem*
ReadBeam2(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadBeam2");
	
	/* Nodo 1 */
	const StructNode* pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	
	Mat3x3 R1(pNode1->GetRCurr());   
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
	Mat3x3 Rn1(Eye3);
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		Rn1 = HP.GetRotRel(ReferenceFrame(pNode1));
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): " 
			<< f1 << std::endl << "(global frame): "
			<< pNode1->GetXCurr() + pNode1->GetRCurr()*f1 << std::endl);
	
	/* Nodo 2 */
	const StructNode* pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
	
	Mat3x3 R2(pNode2->GetRCurr());
	if (HP.IsKeyWord("position")) {
		/* just eat it! */
		NO_OP;
	}
	Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
	Mat3x3 Rn2(Eye3);
	if (HP.IsKeyWord("orientation") || HP.IsKeyWord("rot")) {
		Rn2 = HP.GetRotRel(ReferenceFrame(pNode2));
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "node 2 offset (node reference frame): " 
			<< f2 << std::endl << "(global frame): "
			<< pNode2->GetXCurr() + pNode2->GetRCurr()*f2 << std::endl);
	
	/* Matrice R */
	Mat3x3 R;
	flag f(0);
	if (HP.IsKeyWord("from" "nodes") || HP.IsKeyWord("node")) {
		f = flag(1);
	} else {
		R = HP.GetRotAbs(::AbsRefFrame);
	}
	
	/* Legame costitutivo */
	ConstLawType::Type CLType = ConstLawType::UNKNOWN;
	ConstitutiveLaw6D* pD = HP.GetConstLaw6D(CLType);

	if (pD->iGetNumDof() != 0) {
     		silent_cerr("line " << HP.GetLineData()
			<< ": Beam2(" << uLabel << ") does not support "
			"dynamic constitutive laws yet"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	Beam::Type Type;
	if (CLType == ConstLawType::ELASTIC) {
		Type = Beam::ELASTIC;
	} else {
		Type = Beam::VISCOELASTIC;
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
	
	flag fPiezo(0);
	Mat3xN PiezoMat[2];
	integer iNumElec = 0;
	const ScalarDifferentialNode** pvElecDofs = 0;
	if (HP.IsKeyWord("piezoelectric" "actuator")) {
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
			silent_cerr("Beam2(" << uLabel << "): "
				"illegal number of electric nodes "
				<< iNumElec
				<< " at line " << HP.GetLineData()
				<< std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		
		SAFENEWARR(pvElecDofs, const ScalarDifferentialNode *, iNumElec);
		
		for (integer i = 0; i < iNumElec; i++) {
			pvElecDofs[i] = pDM->ReadNode<const ScalarDifferentialNode, Node::ABSTRACT>(HP);
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

	OrientationDescription od = UNKNOWN_ORIENTATION_DESCRIPTION;
	unsigned uFlags = Beam::OUTPUT_NONE;
	ReadOptionalBeamCustomOutput(pDM, HP, uLabel, Type, uFlags, od);

	flag fOut = pDM->fReadOutput(HP, Elem::BEAM);       
	if (fOut) {
		fOut |= uFlags;
	}

	/* Se necessario, interpola i parametri di rotazione delle sezioni */
	if (f) {
		Mat3x3 RR2 = R2*Rn2;
		Vec3 g(Vec3(CGR_Rot::Param, RR2.MulTM(R1*Rn1))/2);
		R = RR2*Mat3x3(CGR_Rot::MatR, g);
	}
	
	std::ostream& out = pDM->GetLogFile();
	out << "beam2: " << uLabel
		<< " " << pNode1->GetLabel()
		<< " ", f1.Write(out, " ")
		<< " " << pNode2->GetLabel()
		<< " ", f2.Write(out, " ")
		<< std::endl;

	Elem* pEl = NULL;
	
	if (CLType == ConstLawType::ELASTIC) {
		if (fPiezo == 0) {	 
			SAFENEWWITHCONSTRUCTOR(pEl,
					Beam2,
					Beam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						od,
						fOut));
		} else {	 
			SAFENEWWITHCONSTRUCTOR(pEl,
					PiezoActuatorBeam2,
					PiezoActuatorBeam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						iNumElec,
						pvElecDofs,
						PiezoMat[0], PiezoMat[1],
						od,
						fOut));
		}
		
	} else /* At least one is VISCOUS or VISCOELASTIC */ {
		if (fPiezo == 0) {	 
			SAFENEWWITHCONSTRUCTOR(pEl, 
					ViscoElasticBeam2,
					ViscoElasticBeam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						od,
						fOut));
		} else {
			SAFENEWWITHCONSTRUCTOR(pEl,
					PiezoActuatorVEBeam2,
					PiezoActuatorVEBeam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						iNumElec,
						pvElecDofs,
						PiezoMat[0], PiezoMat[1],
						od,
						fOut));
		}
	}
	
	// add here inverse dynamics
	bool bIsErgonomy(false);
	bool bIsRightHandSide(true);

	if (HP.IsKeyWord("inverse" "dynamics")) {
		if (HP.IsKeyWord("torque")) {
			silent_cerr("Beam2(" << uLabel << "): \"torque\" meaningless in this context "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (HP.IsKeyWord("prescribed" "motion")) {
			silent_cerr("Beam2(" << uLabel << "): \"prescribed motion\" meaningless in this context "
				"at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (HP.IsKeyWord("right" "hand" "side")) {
			bIsRightHandSide = HP.GetYesNoOrBool(bIsRightHandSide);
		}

		if (HP.IsKeyWord("ergonomy")) {
			bIsErgonomy = HP.GetYesNoOrBool(bIsErgonomy);
			if (bIsErgonomy) {
				ConstLawType::Type type = pD->GetConstLawType();
				if (type != ConstLawType::ELASTIC) {
					silent_cerr("Beam2(" << uLabel << "): invalid constitutive law type (must be ELASTIC)" << std::endl);
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				if (bIsRightHandSide) {
					silent_cerr("warning, Beam2(" << uLabel << ") is both \"ergonomy\" and \"right hand side\"" << std::endl);
				}
			}
		}
	}

	// set flags for inverse dynamics
	if (pDM->bIsInverseDynamics() && pEl->bInverseDynamics()) {
		unsigned flags = 0;
		if (bIsRightHandSide) {
			flags |= InverseDynamics::RIGHT_HAND_SIDE;
		}
		if (bIsErgonomy) {
			flags |= InverseDynamics::ERGONOMY;
		}
		pEl->SetInverseDynamicsFlags(flags);
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		silent_cerr("semicolon expected at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	return pEl;
} /* End of ReadBeam2() */

