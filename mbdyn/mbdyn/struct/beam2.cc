/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/float.h>

#include <dataman.h>
#include <constltp.h>
#include <shapefnc.h>
#include <beam.h>
#include <beam2.h>
#include <pzbeam2.h>

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
		flag fOut)
: Elem(uL, Elem::BEAM, fOut), 
ElemGravityOwner(uL, Elem::BEAM, fOut), 
InitialAssemblyElem(uL, Elem::BEAM, fOut),
BeamT(Beam::ELASTIC),
fFirstRes(1)
{
	/* Validazione dati */
	ASSERT(pN1 != NULL);
	ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pN2 != NULL);
	ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
   
	pNode[NODE1] = pN1;
	pNode[NODE2] = pN2;
	(Vec3&)f[NODE1] = F1;
	(Vec3&)f[NODE2] = F2;
	(Mat3x3&)RNode[NODE1] = R1;
	(Mat3x3&)RNode[NODE2] = R2;
	RRef = R = (Mat3x3&)r;
	
	pD = NULL; 
	SAFENEWWITHCONSTRUCTOR(pD,
			ConstitutiveLaw6DOwner,
			ConstitutiveLaw6DOwner(pd));
	
	Omega = Vec3(0.); 
	Az = Vec6(0.);
	AzRef = Vec6(0.);
	AzLoc = Vec6(0.);
	AzLocRef = Vec6(0.);
	DefLoc = Vec6(0.);
	DefLocRef = Vec6(0.);
	p = Vec3(0.);
	g = Vec3(0.);
	L0 = Vec3(0.);
	L = Vec3(0.);
	
	DsDxi();
	
	Vec3 xTmp[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {      
		xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
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
	return 6;
}

unsigned int
Beam2::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	/*
	 * {ex|k{x|y|z}}
	 */

	unsigned int idx = 0;

	switch (s[0]) {
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
Beam2::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= 6);

	switch (i) {
	case 1:
	case 4:
	case 5:
	case 6:
		return DefLoc.dGet(i);
	case 2:
	case 3:
		std::cerr << "Beam2 " << GetLabel() 
			<< ": not allowed to return shear strain" << std::endl;
		THROW(ErrGeneric());
	default:
		std::cerr << "Beam2 " << GetLabel()
			<< ": illegal private data " << i << std::endl;
		THROW(ErrGeneric());
	}
#ifndef USE_EXCEPTIONS
	return 0.;
#endif /* USE_EXCEPTIONS */
}

Vec3 
Beam2::InterpState(const Vec3& v1, const Vec3& v2)
{
	doublereal* pv1 = (doublereal*)v1.pGetVec();
	doublereal* pv2 = (doublereal*)v2.pGetVec();
	return Vec3(pv1[0]*dN2[0]+pv2[0]*dN2[1],
			pv1[1]*dN2[0]+pv2[1]*dN2[1],
			pv1[2]*dN2[0]+pv2[2]*dN2[1]);
}


Vec3
Beam2::InterpDeriv(const Vec3& v1, const Vec3& v2)
{
	doublereal* pv1 = (doublereal*)v1.pGetVec();
	doublereal* pv2 = (doublereal*)v2.pGetVec();
	return Vec3((pv1[0]*dN2P[0]+pv2[0]*dN2P[1])*dsdxi,
			(pv1[1]*dN2P[0]+pv2[1]*dN2P[1])*dsdxi,
			(pv1[2]*dN2P[0]+pv2[2]*dN2P[1])*dsdxi);
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
		xTmp[i] = xNod[i]+RNod[i]*f[i];
	}
	
	dsdxi = 1.;

	Vec3 xGrad = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	doublereal d = xGrad.Dot();
	if (d > DBL_EPSILON) {
		dsdxi = 1./sqrt(d);
	} else {
		std::cerr << "warning, beam " << GetLabel() 
			<< " has singular metric; aborting ..." << std::endl;
		
		THROW(Beam2::ErrGeneric());
	}

	/* Calcola le deformazioni iniziali */
	L0 = R.Transpose()*InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
	pD->Update(0.);
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
	Vec3 gTmp(gparam(RNod[NODE2].Transpose()*RNod[NODE1]));
	
	/*
	 * Le derivate dei parametri di rotazione si ricavano da omega
	 */
	Vec3 g1P(Mat3x3(MatGm1, gTmp*(-.5))*w[NODE1]);
	Vec3 g2P(Mat3x3(MatGm1, gTmp*.5)*w[NODE2]);

        Vec3 gPTmp(g1P*dN2[NODE1]+g2P*dN2[NODE2]);
        Omega = Mat3x3(MatG, gTmp)*gPTmp;
	
#if 0
	/* Modo brutale: interpolo le velocita' dei nodi */
	Omega = pNode[NODE1]->GetWCurr()*dN2[NODE1]
		+pNode[NODE2]->GetWCurr()*dN2[NODE2];
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
		AzTmp[i] = Mat6x6(Mat3x3(dN2P[i]*dsdxi*dCoef),
				Zero3x3,
				Mat3x3(LRef*(dN2[i]*dCoef)
					-fTmp[i]*(dN2P[i]*dsdxi*dCoef)),
				Mat3x3(dN2P[i]*dsdxi*dCoef));
		
		/* Delta - azioni interne */
		AzTmp[i] = DRef*AzTmp[i];
		
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
	}
	
	/* correzione alle equazioni */
	WMA.Add(4, 1, Mat3x3(AzRef.GetVec1()*(-dCoef)));
	WMA.Add(10, 7, Mat3x3(AzRef.GetVec1()*dCoef));
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
	
	/* Recupera i dati dei nodi */
	Vec3 xNod[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		xNod[i] = pNode[i]->GetXCurr();
	}   
	
	if (fFirstRes) {
		fFirstRes = flag(0); /* AfterPredict ha gia' calcolato tutto */
	} else {
		Vec3 gNod[NUMNODES];    
		Vec3 xTmp[NUMNODES];
		
		for (unsigned int i = 0; i < NUMNODES; i++) {      
			gNod[i] = pNode[i]->GetgCurr();	 
			xTmp[i] = xNod[i]+pNode[i]->GetRCurr()*f[i];
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
		RDelta = Mat3x3(MatR, g);
		R = RDelta*RRef;
		
		/* Derivate della posizione */
		L = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
		
		/* Derivate dei parametri di rotazione */
		gGrad = InterpDeriv(gNod[NODE1], gNod[NODE2]);
		
		/* Per le deformazioni nel sistema del materiale */
		Mat3x3 RTmp(R.Transpose());
		
		/*
		 * Calcola le deformazioni nel sistema locale
		 * nei punti di valutazione
		 */
		DefLoc = Vec6(RTmp*L-L0, RTmp*(Mat3x3(MatG, g)*gGrad)
				+DefLocRef.GetVec2());
		
		/* Calcola le azioni interne */
		pD->Update(DefLoc);
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
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
		WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
		WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
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
		WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
		WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
	}      
	
	AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);
	
	return WorkVec;
}

    
/* Settings iniziali, prima della prima soluzione */
void
Beam2::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	(Mat3x3&)RRef = R;
	(Vec3&)LRef = L;
	(Vec6&)DefLocRef = DefLoc;
	(Vec6&)AzLocRef = AzLoc;
	(Vec6&)AzRef = Az;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	(Mat6x6&)DRef = MultRMRt(pD->GetFDE(), RRef);      
	
	(flag&)fFirstRes = flag(1);
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
	
	AzLocRef = AzLoc;
	
	/* Porta le azioni interne nel sistema globale */
	Az = AzRef = MultRV(AzLoc, R);
	
	/* Aggiorna il legame costitutivo di riferimento */
	DRef = MultRMRt(pD->GetFDE(), RRef);
	
	fFirstRes = flag(1);
}


/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
Beam2::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.Beams() << std::setw(8) << GetLabel() << " " 
			<< AzLoc.GetVec1() << " " << AzLoc.GetVec2() << std::endl;
	}
}

/* 
 * Output di un modello NASTRAN equivalente nella configurazione corrente 
 */
void
Beam2::Output_pch(std::ostream& out) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		unsigned int label = GetLabel();
		const char *name = GetName();
		out << "$ Beam " << GetLabel();
		if (name) {
			out << " (" << name << ")";
		}
		
#define __NASTRAN_FORMAT__ __HACK_NASTRAN_MODES__
		Vec3 F1(pNode[NODE1]->GetRCurr()*f[NODE1]);
		Vec3 F2(pNode[NODE2]->GetRCurr()*f[NODE2]);

#if __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED__
		out << std::endl
			/* PBEAM */
			<< "PBEAM   "
			<< std::setw(8) << label             /* label */
			<< std::setw(8) << 1                 /* material */
			<< std::setw(8) << 1.                /* area */
			<< std::setw(8) << 1.                /* J1 */
			<< std::setw(8) << 1.                /* J2 */
			<< std::setw(8) << " "               /* J12 */
			<< std::setw(8) << 1.                /* Jp */
			<< std::endl
			
			/* CBEAM */
			<< "CBEAM   "
			<< std::setw(8) << label             /* label */
			<< std::setw(8) << label             /* prop */
			<< std::setw(8) << pNode[NODE1]->GetLabel()  /* node 1 */
			<< std::setw(8) << pNode[NODE2]->GetLabel()  /* node 2 */
			<< std::setw(32) << " "
			<< "+" << std::setw(7) << 1
			<< std::endl
			<< "+" << std::setw(7) << 1
			<< std::setw(16) << " "
			<< std::setw(8) << F1.dGet(1)
			<< std::setw(8) << F1.dGet(2)
			<< std::setw(8) << F1.dGet(3)
			<< std::setw(8) << F2.dGet(1)
			<< std::setw(8) << F2.dGet(2)
			<< std::setw(8) << F2.dGet(3)
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED16__
		out << std::endl
			/* PBEAM */
			<< "PBEAM*  "
			<< std::setw(16) << label            /* label */
			<< std::setw(16) << 1                /* material */
			<< std::setw(16) << 1.               /* area */
			<< std::setw(16) << 1.               /* J1 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(16) << 1.               /* J2 */
			<< std::setw(16) << " "              /* J12 */
			<< std::setw(16) << 1.               /* Jp */
			<< std::endl
			
			/* CBEAM */
			<< "CBEAM*  "
			<< std::setw(16) << label            /* label */
			<< std::setw(16) << label            /* prop */
			<< std::setw(16) << pNode[NODE1]->GetLabel() /* node 1 */
			<< std::setw(16) << pNode[NODE2]->GetLabel() /* node 2 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(64) << " "
			<< "*" << std::setw(7) << 2
			<< std::endl
			<< "*" << std::setw(7) << 2
			<< std::setw(32) << " "
			<< std::setw(16) << F1.dGet(1)
			<< std::setw(16) << F1.dGet(2)
			<< "*" << std::setw(7) << 3
			<< std::endl
			<< "*" << std::setw(7) << 3
			<< std::setw(16) << F1.dGet(3)
			<< std::setw(16) << F2.dGet(1)
			<< std::setw(16) << F2.dGet(2)
			<< std::setw(16) << F2.dGet(3)
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FREE__
		out << std::endl
			/* PBEAM */
			<< "PBEAM,"
			<< label << ","
			<< 1 << ","
			<< 1. << ","
			<< 1. << ","
			<< 1. << ","
			<< ","
			<< 1. << std::endl
			
			/* CBEAM */
			<< "CBEAM,"
			<< label << ","
			<< label << ","
			<< pNode[NODE1]->GetLabel() << ","
			<< pNode[NODE2]->GetLabel() << ",,,,"
#if 0
			<< "," 
#endif
			<< std::endl
#if 1
			<< ","
#endif
			<< " ,,", F1.Write(out, ",") << ",", F2.Write(out, ",")
			<< std::endl;
#else
#error "unknown NASTRAN format"
#endif
	}
#endif /* __HACK_NASTRAN_MODES__ */
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
Beam2::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr) 
{ 
	DEBUGCOUTFNAME("Beam2::InitialAssRes => AssStiffnessVec");

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
Beam2::pGetNode(unsigned int i) const
{
	ASSERT(i >= 1 && i <= 2);
	switch (i) {
	case 1:
	case 2:
		return pNode[i-1];
	default:
		THROW(Beam2::ErrGeneric());
	}
#ifndef USE_EXCEPTIONS
	return NULL;
#endif /* USE_EXCEPTIONS */
}


#ifdef USE_ADAMS
void
Beam2::GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& r) const
{
	ASSERT(part == 1);
	part--;
   
	x = p;
	r = R;
}


std::ostream& 
Beam2::WriteAdamsDummyPartCmd(std::ostream& out,
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
		flag fOut)
: Elem(uL, Elem::BEAM, fOut),
Beam2(uL, pN1, pN2, F1, F2, R1, R2, r, pd, fOut)
{
	SetBeamType(Beam::VISCOELASTIC);
	
	LPrimeRef = LPrime = Vec3(0.);  
	gPrime = Vec3(0.);
	
	DefPrimeLoc = DefPrimeLocRef = Vec6(0.);
	
	/* Nota: DsDxi() viene chiamata dal costruttore di Beam */
	Beam2::Omega0();
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
	
	/* Recupera i dati dei nodi */
	Vec3 xNod[NUMNODES];
	
	for (unsigned int i = 0; i < NUMNODES; i++) {
		xNod[i] = pNode[i]->GetXCurr();
	}
	
	if (fFirstRes) {
		fFirstRes = flag(0); /* AfterPredict ha gia' calcolato tutto */
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
ViscoElasticBeam2::SetValue(VectorHandler& X, VectorHandler& XP) const
{
	Beam2::SetValue(X, XP);
	
	/* Aggiorna le grandezze della trave nei punti di valutazione */
	(Vec3&)OmegaRef = Omega;
	(Vec3&)LPrimeRef = LPrime;
	(Vec6&)DefPrimeLocRef = DefPrimeLoc;
	
	/*
	 * Aggiorna il legame costitutivo di riferimento
	 * (la deformazione e' gia' stata aggiornata dall'ultimo residuo)
	 */
	(Mat6x6&)ERef = MultRMRt(pD->GetFDEPrime(), RRef);
	
	ASSERT(fFirstRes == flag(1));
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
	
	AzLocRef = AzLoc;
	
	/* Porta le azioni interne nel sistema globale */
	Az = AzRef = MultRV(AzLoc, R);
	
	/* Aggiorna il legame costitutivo */
	DRef = MultRMRt(pD->GetFDE(), R);
	ERef = MultRMRt(pD->GetFDEPrime(), R);
	
	fFirstRes = flag(1);
}

/* ViscoElasticBeam - end */


/* Legge una trave */
Elem*
ReadBeam2(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadBeam2");
	
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
	StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
	
	Mat3x3 R1(pNode1->GetRCurr());   
	Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
	Mat3x3 Rn1(Eye3);
	if (HP.IsKeyWord("rot")) {
		Rn1 = HP.GetRotRel(ReferenceFrame(pNode1));
	}

	DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): " 
			<< f1 << std::endl << "(global frame): "
			<< pNode1->GetXCurr()+pNode1->GetRCurr()*f1 << std::endl);
	
	/* Nodo 2 */
	StructNode* pNode2 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
	
	Mat3x3 R2(pNode2->GetRCurr());
	Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
	Mat3x3 Rn2(Eye3);
	if (HP.IsKeyWord("rot")) {
		Rn2 = HP.GetRotRel(ReferenceFrame(pNode2));
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "node 2 offset (node reference frame): " 
			<< f2 << std::endl << "(global frame): "
			<< pNode2->GetXCurr()+pNode2->GetRCurr()*f2 << std::endl);
	
	/* Matrice R */
	Mat3x3 R;
	flag f(0);
	if (HP.IsKeyWord("node")) {
		f = flag(1);
	} else {
		R = HP.GetRotAbs(AbsRefFrame);
	}
	
	/* Legame costitutivo */
	ConstLawType::Type CLType = ConstLawType::UNKNOWN;
	ConstitutiveLaw6D* pD = HP.GetConstLaw6D(CLType);

	if (pD->iGetNumDof() != 0) {
     		std::cerr << "line " << HP.GetLineData()
			<< ": beam2 does not support "
			"dynamic constitutive laws yet"
			<< std::endl;
		THROW(ErrGeneric());
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
	
#if defined(USE_ELECTRIC_NODES)
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
			std::cerr << "illegal number of electric nodes "
				<< iNumElec << std::endl;
			THROW(ErrGeneric());
		}
		
		SAFENEWARR(pvElecDofs, ScalarDifferentialNode*, iNumElec);
		
		for (integer i = 0; i < iNumElec; i++) {
			unsigned int uL = HP.GetInt();
			DEBUGLCOUT(MYDEBUG_INPUT, "linked to abstract node "
					<< uL << std::endl);
			pvElecDofs[i] = (ScalarDifferentialNode*)(pDM->pFindNode(Node::ABSTRACT, uL));
			if (pvElecDofs[i] == NULL) {
				std::cerr << "can't find abstract node "
					<< uL << std::endl;
				THROW(ErrGeneric());
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
#endif /* defined(USE_ELECTRIC_NODES) */

	flag fOut = pDM->fReadOutput(HP, Elem::BEAM);       
	
	
	/* Se necessario, interpola i parametri di rotazione delle sezioni */
	if (f) {
		Mat3x3 RR2 = R2*Rn2;
		Vec3 g(gparam(RR2.Transpose()*(R1*Rn1))*.5);
		R = RR2*Mat3x3(MatR, g);
	}
	
	Elem* pEl = NULL;
	
	if (CLType == ConstLawType::ELASTIC) {
#if defined(USE_ELECTRIC_NODES)
		if (fPiezo == 0) {	 
#endif /* defined(USE_ELECTRIC_NODES) */
			SAFENEWWITHCONSTRUCTOR(pEl,
					Beam2,
					Beam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						fOut));
#if defined(USE_ELECTRIC_NODES)
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
						fOut));
		}
#endif /* defined(USE_ELECTRIC_NODES) */
		
	} else /* At least one is VISCOUS or VISCOELASTIC */ {
#if defined(USE_ELECTRIC_NODES)      
		if (fPiezo == 0) {	 
#endif /* defined(USE_ELECTRIC_NODES) */
			SAFENEWWITHCONSTRUCTOR(pEl, 
					ViscoElasticBeam2,
					ViscoElasticBeam2(uLabel,
						pNode1, pNode2,
						f1, f2,
						Rn1, Rn2,
						R,
						pD,
						fOut));
#if defined(USE_ELECTRIC_NODES)
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
						fOut));
		}
#endif /* defined(USE_ELECTRIC_NODES) */
	}
	
	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		std::cerr << "semicolon expected at line "
			<< HP.GetLineData() << std::endl;      
		THROW(DataManager::ErrGeneric());
	}
	
	return pEl;
} /* End of ReadBeam2() */

