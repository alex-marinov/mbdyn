/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#include <constltp.h>
#include <shapefnc.h>
#include <beam.h>
#include <pzbeam.h>
#include <dataman.h>

/*
 * Nota: non e' ancora stato implementato il contributo 
 * della ViscoElasticBeam all'assemblaggio iniziale
 */

/*
 * Nota: la parte viscoelastica va rivista in accordo con la piu' 
 * recente formulazione delle derivate delle deformazioni nel sistema
 * materiale
 */

/* Beam - begin */

/* Costruttore normale */
Beam::Beam(unsigned int uL, 
	   const StructNode* pN1,
	   const StructNode* pN2,
	   const StructNode* pN3,
	   const Vec3& F1,
	   const Vec3& F2,
	   const Vec3& F3,
	   const Mat3x3& R1,
	   const Mat3x3& R2,
	   const Mat3x3& R3,
	   const Mat3x3& r_I,
	   const Mat3x3& rII,
	   const ConstitutiveLaw6D* pD_I,
	   const ConstitutiveLaw6D* pDII,
	   flag fOut)
: Elem(uL, Elem::BEAM, fOut), 
ElemGravityOwner(uL, Elem::BEAM, fOut), 
InitialAssemblyElem(uL, Elem::BEAM, fOut),
BeamT(Beam::ELASTIC),
fConsistentInertia(0), 
dMass_I(0.),
S0_I(0.),
J0_I(0.),
dMassII(0.),
S0II(0.),
J0II(0.),
fFirstRes(1)
{
    ASSERT(pN1 != NULL);
    ASSERT(pN1->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pN2 != NULL);
    ASSERT(pN2->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pN3 != NULL);
    ASSERT(pN3->GetNodeType() == Node::STRUCTURAL);

    pNode[NODE1] = pN1;
    pNode[NODE2] = pN2;
    pNode[NODE3] = pN3;
    
    (Vec3&)f[NODE1] = F1;
    (Vec3&)f[NODE2] = F2;
    (Vec3&)f[NODE3] = F3;
    (Mat3x3&)RNode[NODE1] = R1;
    (Mat3x3&)RNode[NODE2] = R2;
    (Mat3x3&)RNode[NODE3] = R3;
    RRef[S_I] = R[S_I] = (Mat3x3&)r_I;
    RRef[SII] = R[SII] = (Mat3x3&)rII;
  
    pD[S_I] = NULL; 
    SAFENEWWITHCONSTRUCTOR(pD[S_I], 
			   ConstitutiveLaw6DOwner,
			   ConstitutiveLaw6DOwner(pD_I));
    pD[SII] = NULL;   
    SAFENEWWITHCONSTRUCTOR(pD[SII],
			   ConstitutiveLaw6DOwner,
			   ConstitutiveLaw6DOwner(pDII));
   
    Omega[S_I]     = Omega[SII]     = Vec3(0.); 
    Az[S_I]        = Az[SII]        = Vec6(0.);
    AzRef[S_I]     = AzRef[SII]     = Vec6(0.);
    AzLoc[S_I]     = AzLoc[SII]     = Vec6(0.);
    AzLocRef[S_I]  = AzLocRef[SII]  = Vec6(0.);
    DefLoc[S_I]    = DefLoc[SII]    = Vec6(0.);
    DefLocRef[S_I] = DefLocRef[SII] = Vec6(0.);
    p[S_I]         = p[SII]         = Vec3(0.);
    g[S_I]         = g[SII]         = Vec3(0.);
    L0[S_I]        = L0[SII]        = Vec3(0.);
    L[S_I]         = L[SII]         = Vec3(0.);
    
    DsDxi();
   
    Vec3 xTmp[NUMNODES];
      
    for (unsigned int i = 0; i < NUMNODES; i++) {      
        xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
    }      
      
    /* Aggiorna le grandezze della trave nei punti di valutazione */
    for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
        p[iSez] = InterpState(xTmp[NODE1],
	                      xTmp[NODE2],
			      xTmp[NODE3],
			      Beam::Section(iSez));
    }
}


/* Costruttore per la trave con forze d'inerzia consistenti */
Beam::Beam(unsigned int uL,
	   const StructNode* pN1,
	   const StructNode* pN2,
	   const StructNode* pN3,
	   const Vec3& F1,
	   const Vec3& F2,
	   const Vec3& F3,
	   const Mat3x3& R1,
	   const Mat3x3& R2,
	   const Mat3x3& R3,
	   const Mat3x3& r_I,
	   const Mat3x3& rII,
	   const ConstitutiveLaw6D* pD_I,
	   const ConstitutiveLaw6D* pDII,
	   doublereal dM_I,
	   const Mat3x3& s0_I,
	   const Mat3x3& j0_I,
	   doublereal dMII,
	   const Mat3x3& s0II,
	   const Mat3x3& j0II,
	   flag fOut)
: Elem(uL, Elem::BEAM, fOut), 
ElemGravityOwner(uL, Elem::BEAM, fOut), 
InitialAssemblyElem(uL, Elem::BEAM, fOut),
BeamT(Beam::ELASTIC),
fConsistentInertia(1),
dMass_I(dM_I),
S0_I(s0_I),
J0_I(j0_I),
dMassII(dMII),
S0II(s0II),
J0II(j0II),
fFirstRes(1)
{
    pNode[NODE1] = pN1;
    pNode[NODE2] = pN2;
    pNode[NODE3] = pN3;
    (Vec3&)f[NODE1] = F1;
    (Vec3&)f[NODE2] = F2;
    (Vec3&)f[NODE3] = F3;
    (Mat3x3&)RNode[NODE1] = R1;
    (Mat3x3&)RNode[NODE2] = R2;
    (Mat3x3&)RNode[NODE3] = R3;
    RRef[S_I] = R[S_I] = (Mat3x3&)r_I;
    RRef[SII] = R[SII] = (Mat3x3&)rII;

    SAFENEWWITHCONSTRUCTOR(pD[S_I], 
			   ConstitutiveLaw6DOwner,
			   ConstitutiveLaw6DOwner(pD_I));
    SAFENEWWITHCONSTRUCTOR(pD[SII],
			   ConstitutiveLaw6DOwner,
			   ConstitutiveLaw6DOwner(pDII));
   
    Omega[S_I]     = Omega[SII]     = Vec3(0.); 
    Az[S_I]        = Az[SII]        = Vec6(0.);
    AzRef[S_I]     = AzRef[SII]     = Vec6(0.);
    AzLoc[S_I]     = AzLoc[SII]     = Vec6(0.);
    AzLocRef[S_I]  = AzLocRef[SII]  = Vec6(0.);
    DefLoc[S_I]    = DefLoc[SII]    = Vec6(0.);
    DefLocRef[S_I] = DefLocRef[SII] = Vec6(0.);
    p[S_I]         = p[SII]         = Vec3(0.);
    g[S_I]         = g[SII]         = Vec3(0.);
    L0[S_I]        = L0[SII]        = Vec3(0.);
    L[S_I]         = L[SII]         = Vec3(0.);
    
    DsDxi();
}  


Beam::~Beam(void) 
{
    ASSERT(pD[S_I] != NULL);
    if (pD[S_I] != NULL) {      
        SAFEDELETE(pD[S_I]);
    }
   
    ASSERT(pD[SII] != NULL);
    if (pD[SII] != NULL) {      
        SAFEDELETE(pD[SII]);
    }
}


/* Accesso ai dati privati */
unsigned int
Beam::iGetNumPrivData(void) const
{
	return 12;
}

unsigned int
Beam::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	/*
	 * p{I|II}.{ex|k{x|y|z}}
	 */

	unsigned int idx = 0;

	if (strncmp(s, "pII.", sizeof("pII.") - 1) == 0) {
		s += sizeof("pII.") - 1;
		idx += 6;

	} else if (strncmp(s, "pI.", sizeof("pI.") - 1) == 0) {
		s += sizeof("pI.") - 1;

	} else {
		return 0;
	}

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
Beam::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 0 && i <= 12);
	switch (i) {
	case 1:
	case 4:
	case 5:
	case 6:
	case 7:
	case 10:
	case 11:
	case 12:
		return DefLoc[(i-1)/6].dGet((i-1)%6+1);
	case 2:
	case 3:
	case 8:
	case 9:
		std::cerr << "Beam " << GetLabel() 
			<< ": not allowed to return shear strain" << std::endl;
		THROW(ErrGeneric());
	default:
		std::cerr << "Beam " << GetLabel() << ": illegal private data " 
			<< i << std::endl;
		THROW(ErrGeneric());
	}
#ifndef USE_EXCEPTIONS
	return 0.;
#endif /* USE_EXCEPTIONS */
}

Vec3 
Beam::InterpState(const Vec3& v1,
                  const Vec3& v2,
		  const Vec3& v3,
		  enum Section Sec)
{
    doublereal* pv1 = (doublereal*)v1.pGetVec();
    doublereal* pv2 = (doublereal*)v2.pGetVec();
    doublereal* pv3 = (doublereal*)v3.pGetVec();
    return Vec3(pv1[0]*dN3[Sec][0]+pv2[0]*dN3[Sec][1]+pv3[0]*dN3[Sec][2],
	        pv1[1]*dN3[Sec][0]+pv2[1]*dN3[Sec][1]+pv3[1]*dN3[Sec][2],
	        pv1[2]*dN3[Sec][0]+pv2[2]*dN3[Sec][1]+pv3[2]*dN3[Sec][2]);
}


Vec3
Beam::InterpDeriv(const Vec3& v1,
                  const Vec3& v2,
		  const Vec3& v3,
		  enum Section Sec)
{
    doublereal* pv1 = (doublereal*)v1.pGetVec();
    doublereal* pv2 = (doublereal*)v2.pGetVec();
    doublereal* pv3 = (doublereal*)v3.pGetVec();
    return Vec3((pv1[0]*dN3P[Sec][0]+pv2[0]*dN3P[Sec][1]
		 +pv3[0]*dN3P[Sec][2])*dsdxi[Sec],
	        (pv1[1]*dN3P[Sec][0]+pv2[1]*dN3P[Sec][1]
		 +pv3[1]*dN3P[Sec][2])*dsdxi[Sec],
	        (pv1[2]*dN3P[Sec][0]+pv2[2]*dN3P[Sec][1]
		 +pv3[2]*dN3P[Sec][2])*dsdxi[Sec]);
}


void
Beam::DsDxi(void)
{
    /* Validazione dati */
    ASSERT(pNode[NODE1] != NULL);
    ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE2] != NULL);
    ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE3] != NULL);
    ASSERT(pNode[NODE3]->GetNodeType() == Node::STRUCTURAL);     
   
    /* Calcola il ds/dxi e le deformazioni iniziali */
    Vec3 xNod[NUMNODES];
    Mat3x3 RNod[NUMNODES];
    Vec3 xTmp[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {
        xNod[i] = pNode[i]->GetXCurr();
        RNod[i] = pNode[i]->GetRCurr();
        xTmp[i] = xNod[i]+RNod[i]*f[i];
    }

    dsdxi[S_I] = 1.;
    dsdxi[SII] = 1.;

    Vec3 xGrad[NUMSEZ];
    for (unsigned int i = 0; i < NUMSEZ; i++) {
        xGrad[i] = InterpDeriv(xTmp[NODE1],
	                       xTmp[NODE2],
			       xTmp[NODE3],
			       Beam::Section(i));
        doublereal d = xGrad[i].Dot();
        ASSERT(d > DBL_EPSILON);
        if (d > DBL_EPSILON) {
	    dsdxi[i] = 1./sqrt(d);
        } else {
	    std::cerr << "warning, beam " << GetLabel() 
	        << " has singular metric; aborting ..." << std::endl;
	 
	    THROW(Beam::ErrGeneric());
        }
    }

    /* Calcola le deformazioni iniziali */
    for (unsigned int i = 0; i < NUMSEZ; i++) {
        L0[i] = R[i].Transpose()*InterpDeriv(xTmp[NODE1],
	                                     xTmp[NODE2],
					     xTmp[NODE3],
					     Beam::Section(i));
        pD[i]->Update(0.);
        DRef[i] = MultRMRt(pD[i]->GetFDE(), R[i]);
    }
}


/* Calcola la velocita' angolare delle sezioni a partire da quelle dei nodi */
void 
Beam::Omega0(void)
{   
    /* Validazione dati */
    ASSERT(pNode[NODE1] != NULL);
    ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE2] != NULL);
    ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE3] != NULL);
    ASSERT(pNode[NODE3]->GetNodeType() == Node::STRUCTURAL);     
 
    /* Modo consistente: */      
    Mat3x3 RNod[NUMNODES];
    Vec3 w[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {     
        RNod[i] = pNode[i]->GetRCurr()*RNode[i];
        w[i] = pNode[i]->GetWCurr();
    }
   
    /* Calcolo i parametri di rotazione dei nodi di estremo rispetto a quello 
     * centrale, nell'ipotesi (realistica) che siano limitati */
    Mat3x3 RTmp(RNod[NODE2].Transpose());
    Vec3 g1(gparam(RTmp*RNod[NODE1]));
    Vec3 g3(gparam(RTmp*RNod[NODE3]));
   
    /* Calcolo le derivate dei parametri di rotazione ai nodi di estremo; in
     * quello centrale si assume che sia uguale alla velocita' di rotazione */
    Vec3 g1P(Mat3x3(MatGm1, g1)*w[NODE1]);
    Vec3 g3P(Mat3x3(MatGm1, g3)*w[NODE3]);

    for (unsigned int i = 0; i < NUMSEZ; i++) {
        Vec3 gTmp(g1*dN3[i][NODE1]+g3*dN3[i][NODE3]);
        Vec3 gPTmp(g1P*dN3[i][NODE1]+w[NODE2]*dN3[i][NODE2]+g3P*dN3[i][NODE3]);
        Omega[i] = Mat3x3(MatG, gTmp)*gPTmp;
    }

#if 0
    /* Modo brutale: interpolo le velocita' dei nodi */
    Vec3 w[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {
        w[i] = pNode[i]->GetWCurr();
    }
    for (unsigned int i = 0; i < NUMSEZ; i++) {      
        Omega[i] = (w{NODE1]*dN3[i][NODE1]
	            +w[NODE2]*dN3[i][NODE2]
		    +w[NODE3]*dN3[i][NODE3]);
    }
#endif /* 0 */
}


/* Contributo al file di restart */
std::ostream& Beam::Restart(std::ostream& out) const
{
   return Restart_(out)<< ';' << std::endl;
}

std::ostream& Beam::Restart_(std::ostream& out) const
{ 
   out << "  beam: " << GetLabel();
   for (unsigned int i = 0; i < NUMNODES; i++) {
      out << ", " << pNode[i]->GetLabel() << ", reference, node, ", 
	f[i].Write(out, ", ");
   }
   for (unsigned int i = 0; i < NUMSEZ; i++) {
      out << ", reference, global, 1, ",
	(R[i].GetVec(1)).Write(out, ", ")
	  << ", 2, ", (R[i].GetVec(2)).Write(out, ", ") << ", ",
	pD[i]->pGetConstLaw()->Restart(out);
   }
   return out;
}


/* Assembla la matrice */
void Beam::AssStiffnessMat(FullSubMatrixHandler& WMA,
			   FullSubMatrixHandler& /* WMB */ ,
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("Beam::AssStiffnessMat");
   
   /* La matrice arriva gia' dimensionata e con gli indici di righe e colonne
    * a posto */
   
   /* offset nel riferimento globale */
   Vec3 fTmp[NUMNODES];
   for (unsigned int i = 0; i < NUMNODES; i++) {
      fTmp[i] = pNode[i]->GetRCurr()*f[i];
   }
   
   Mat6x6 AzTmp[NUMSEZ][NUMNODES];
   
   /* per ogni punto di valutazione: */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {            
      for (unsigned int i = 0; i < NUMNODES; i++) {
	 /* Delta - deformazioni */
	 AzTmp[iSez][i] = Mat6x6(Mat3x3(dN3P[iSez][i]*dsdxi[iSez]*dCoef),
				 Zero3x3,
				 Mat3x3(LRef[iSez]*(dN3[iSez][i]*dCoef)
					-fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez]*dCoef)),
				 Mat3x3(dN3P[iSez][i]*dsdxi[iSez]*dCoef));
	 
	 /* Delta - azioni interne */
	 AzTmp[iSez][i] = DRef[iSez]*AzTmp[iSez][i];
	 
	 /* Correggo per la rotazione da locale a globale */
	 AzTmp[iSez][i].SubMat12(Mat3x3(AzRef[iSez].GetVec1()*(dN3[iSez][i]*dCoef)));
	 AzTmp[iSez][i].SubMat22(Mat3x3(AzRef[iSez].GetVec2()*(dN3[iSez][i]*dCoef)));
      }
   } /* end ciclo sui punti di valutazione */
   
   
   Vec3 bTmp[2];
   
   /* ciclo sulle equazioni */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      bTmp[0] = p[iSez]-pNode[iSez]->GetXCurr();
      bTmp[1] = p[iSez]-pNode[iSez+1]->GetXCurr();
   
      unsigned int iRow1 = iSez*6+1;
      unsigned int iRow2 = iRow1+6;
      
      for (unsigned int i = 0; i < NUMNODES; i++) {
	 /* Equazione all'indietro: */
	 WMA.Sub(iRow1, 6*i+1, AzTmp[iSez][i].GetMat11());
	 WMA.Sub(iRow1, 6*i+4, AzTmp[iSez][i].GetMat12());
	 
	 WMA.Sub(iRow1+3, 6*i+1,
		 AzTmp[iSez][i].GetMat21()
		 -Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
		 +Mat3x3(bTmp[0])*AzTmp[iSez][i].GetMat11());
	 WMA.Sub(iRow1+3, 6*i+4, 
		 AzTmp[iSez][i].GetMat22()
		 -Mat3x3(AzRef[iSez].GetVec1()*(-dCoef*dN3[iSez][i]), fTmp[i])
		 +Mat3x3(bTmp[0])*AzTmp[iSez][i].GetMat12());
	 
	 /* Equazione in avanti: */
	 WMA.Add(iRow2, 6*i+1, AzTmp[iSez][i].GetMat11());
	 WMA.Add(iRow2, 6*i+4, AzTmp[iSez][i].GetMat12());

	 WMA.Add(iRow2+3, 6*i+1,
		 AzTmp[iSez][i].GetMat21()
		 -Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
		 +Mat3x3(bTmp[1])*AzTmp[iSez][i].GetMat11());
	 WMA.Add(iRow2+3, 6*i+4, 
		 AzTmp[iSez][i].GetMat22()
		 +Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]), fTmp[i])
		 +Mat3x3(bTmp[1])*AzTmp[iSez][i].GetMat12());
      }
      
      /* correzione alle equazioni */
      WMA.Add(iRow1+3, 6*iSez+1, Mat3x3(AzRef[iSez].GetVec1()*(-dCoef)));
      WMA.Add(iRow2+3, 6*iSez+7, Mat3x3(AzRef[iSez].GetVec1()*dCoef));
      
   } /* end ciclo sui punti di valutazione */
};


/* Assembla il residuo */
void Beam::AssStiffnessVec(SubVectorHandler& WorkVec,
			   doublereal /* dCoef */ ,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("Beam::AssStiffnessVec");
   
   /* Riceve il vettore gia' dimensionato e con gli indici a posto 
    * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */
   
   /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
    * (ora e' superato) */
   
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
      
      Mat3x3 RDelta[NUMSEZ];
      Vec3 gGrad[NUMSEZ];
      
      /* Aggiorna le grandezze della trave nei punti di valutazione */
      for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
	 
	 /* Posizione */
	 p[iSez] = InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));
	 
	 /* Matrici di rotazione */
	 g[iSez] = InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
	 RDelta[iSez] = Mat3x3(MatR, g[iSez]);
	 R[iSez] = RDelta[iSez]*RRef[iSez];
	 
	 /* Derivate della posizione */
	 L[iSez] = InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));
	 
	 /* Derivate dei parametri di rotazione */
	 gGrad[iSez] = InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
	 
	 /* Per le deformazioni nel sistema del materiale */
	 Mat3x3 RTmp(R[iSez].Transpose());
	 
	 /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
	 DefLoc[iSez] = Vec6(RTmp*L[iSez]-L0[iSez],
			     RTmp*(Mat3x3(MatG, g[iSez])*gGrad[iSez])+DefLocRef[iSez].GetVec2());
	 
	 /* Calcola le azioni interne */
	 pD[iSez]->Update(DefLoc[iSez]);
	 AzLoc[iSez] = pD[iSez]->GetF();
	 
	 /* corregge le azioni interne locali (piezo, ecc) */
	 AddInternalForces(AzLoc[iSez], iSez);

	 /* Porta le azioni interne nel sistema globale */
	 Az[iSez] = MultRV(AzLoc[iSez], R[iSez]);
      }      
   }
   
   WorkVec.Add(1, Az[S_I].GetVec1());
   WorkVec.Add(4, (p[S_I]-xNod[NODE1]).Cross(Az[S_I].GetVec1())
	          +Az[S_I].GetVec2());
   WorkVec.Add(7, Az[SII].GetVec1()-Az[S_I].GetVec1());
   WorkVec.Add(10, Az[SII].GetVec2()-Az[S_I].GetVec2()
	           +(p[SII]-xNod[NODE2]).Cross(Az[SII].GetVec1())
	           -(p[S_I]-xNod[NODE2]).Cross(Az[S_I].GetVec1()));
   WorkVec.Sub(13, Az[SII].GetVec1());
   WorkVec.Add(16, Az[SII].GetVec1().Cross(p[SII]-xNod[NODE3])
	           -Az[SII].GetVec2());   
}

   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& Beam::AssJac(VariableSubMatrixHandler& WorkMat,
				       doublereal dCoef, 
				       const VectorHandler& XCurr,
				       const VectorHandler& XPrimeCurr)
{	
   DEBUGCOUTFNAME("Beam::AssJac => AssStiffnessMat");
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

   FullSubMatrixHandler& WM = WorkMat.SetFull();   
   
   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   if(fConsistentInertia) {	
      WM.ResizeInit(36, 18, 0.);
   } else {
      WM.ResizeInit(18, 18, 0.);	
   }   
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }      
   
   AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   
   if (fConsistentInertia) {	
      for (int iCnt = 1; iCnt <= 6; iCnt++) {
	 WM.fPutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WM.fPutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WM.fPutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }      
      
      AssInertiaMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   }
   
   return WorkMat;
}


/* assemblaggio matrici per autovalori */
void Beam::AssMats(VariableSubMatrixHandler& WorkMatA,
		  VariableSubMatrixHandler& WorkMatB,		  
		  const VectorHandler& XCurr,
		  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("Beam::AssMats => AssStiffnessMat");
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();

   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   
   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   if (fConsistentInertia) {    
      WMA.ResizeInit(36, 18, 0.);
      WMB.ResizeInit(36, 18, 0.);
   } else  {
      WMA.ResizeInit(18, 18, 0.);	
      WorkMatB.SetNullMatrix();
   }   
      
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WMA.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WMA.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WMA.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WMA.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WMA.fPutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
      WMA.fPutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }

   AssStiffnessMat(WMA, WMA, 1., XCurr, XPrimeCurr);
   
   if (fConsistentInertia) {	
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
	 WMB.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
	 WMB.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
	 WMB.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
	 WMB.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
	 WMB.fPutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
	 WMB.fPutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
      }     
      
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
	 WMA.fPutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WMA.fPutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WMA.fPutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);

	 WMB.fPutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WMB.fPutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WMB.fPutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }      
      
      AssInertiaMat(WMA, WMB, 1., XCurr, XPrimeCurr);
   }   
}


/* assemblaggio residuo */
SubVectorHandler& Beam::AssRes(SubVectorHandler& WorkVec,
			       doublereal dCoef,
			       const VectorHandler& XCurr, 
			       const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("Beam::AssRes => AssStiffnessVec");
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();
   integer iNode3FirstMomIndex = pNode[NODE3]->iGetFirstMomentumIndex();

   /* Dimensiona il vettore, lo azzera e pone gli indici corretti */
   if (fConsistentInertia) {	
      WorkVec.Resize(36);
   } else {	
      WorkVec.Resize(18);
   }
   
   WorkVec.Reset(0.);
         
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode3FirstMomIndex+iCnt);
   }      
   
   AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);

   if (fConsistentInertia) {	
      integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
      integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
      integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
      
      for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
	 WorkVec.fPutRowIndex(18+iCnt, iNode1FirstPosIndex+iCnt);
	 WorkVec.fPutRowIndex(24+iCnt, iNode2FirstPosIndex+iCnt);
	 WorkVec.fPutRowIndex(30+iCnt, iNode3FirstPosIndex+iCnt);
      }      
      
      AssInertiaVec(WorkVec, dCoef, XCurr, XPrimeCurr);
   }
   
   return WorkVec;
}

    
/* Settings iniziali, prima della prima soluzione */
void Beam::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{
   /* Aggiorna le grandezze della trave nei punti di valutazione */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      (Mat3x3&)RRef[iSez] = R[iSez];
      (Vec3&)LRef[iSez] = L[iSez];
      (Vec6&)DefLocRef[iSez] = DefLoc[iSez];
      (Vec6&)AzLocRef[iSez] = AzLoc[iSez];
      (Vec6&)AzRef[iSez] = Az[iSez];
      
      /* Aggiorna il legame costitutivo di riferimento
       * (la deformazione e' gia' stata aggiornata dall'ultimo residuo) */
      (Mat6x6&)DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), RRef[iSez]);      
   }
   (flag&)fFirstRes = flag(1);
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void Beam::AfterPredict(VectorHandler& /* X */ , 
			VectorHandler& /* XP */ )
{  
   /* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */
   
   /* Recupera i dati dei nodi */  
   Vec3   gNod[NUMNODES];
   Vec3   xTmp[NUMNODES];
   
   for (unsigned int i = 0; i < NUMNODES; i++) {            
      gNod[i] = pNode[i]->GetgRef();
      xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRRef()*f[i];
   }      
      
   Mat3x3 RDelta[NUMSEZ];
   Vec3 gGrad[NUMSEZ];
   
   /* Aggiorna le grandezze della trave nei punti di valutazione */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      
      /* Posizione */
      p[iSez] = InterpState(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));
      
      /* Matrici di rotazione */
      g[iSez] = InterpState(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
      RDelta[iSez] = Mat3x3(MatR, g[iSez]);
      R[iSez] = RRef[iSez] = RDelta[iSez]*R[iSez];
      
      /* Derivate della posizione */
      L[iSez] = LRef[iSez] 
	= InterpDeriv(xTmp[NODE1], xTmp[NODE2], xTmp[NODE3], Beam::Section(iSez));
      
      /* Derivate dei parametri di rotazione */
      gGrad[iSez] 
	= InterpDeriv(gNod[NODE1], gNod[NODE2], gNod[NODE3], Beam::Section(iSez));
      
      /* Per le deformazioni nel sistema del materiale */
      Mat3x3 RTmp(R[iSez].Transpose());
      
      /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
      DefLoc[iSez] = DefLocRef[iSez]
	= Vec6(RTmp*L[iSez]-L0[iSez],
	       RTmp*(Mat3x3(MatG, g[iSez])*gGrad[iSez])+DefLoc[iSez].GetVec2());
      
      /* Calcola le azioni interne */
      pD[iSez]->Update(DefLoc[iSez]);
      AzLoc[iSez] = pD[iSez]->GetF();

      /* corregge le azioni interne locali (piezo, ecc) */
      AddInternalForces(AzLoc[iSez], iSez);

      AzLocRef[iSez] = AzLoc[iSez];
      
      /* Porta le azioni interne nel sistema globale */
      Az[iSez] = AzRef[iSez] = MultRV(AzLoc[iSez], R[iSez]);

      /* Aggiorna il legame costitutivo di riferimento */
      DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), RRef[iSez]);
   }

   fFirstRes = flag(1);
}


/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void Beam::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
      OH.Beams() << std::setw(8) << GetLabel() << " " 
	<< AzLoc[S_I].GetVec1() << " " << AzLoc[S_I].GetVec2() << " " 
	<< AzLoc[SII].GetVec1() << " " << AzLoc[SII].GetVec2() << std::endl;
   }   
}

/* Output di un modello NASTRAN equivalente nella configurazione corrente */
void
Beam::Output_pch(std::ostream& out) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		unsigned int label = GetLabel();
		if (label > 9999999) {
			std::cerr << "label of Beam(" << label <<") is too large" << std::endl;
			THROW(ErrGeneric());
		}
		
		const char *name = GetName();
		out << "$ Beam " << GetLabel();
		if (name) {
			out << " (" << name << ")";
		}
		
#define __NASTRAN_FORMAT__ __HACK_NASTRAN_MODES__
		Vec3 F1(pNode[NODE1]->GetRCurr()*f[NODE1]);
		Vec3 F2(pNode[NODE2]->GetRCurr()*f[NODE2]);
		Vec3 F3(pNode[NODE3]->GetRCurr()*f[NODE3]);
		
#if __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED__
		out << std::endl
			/* PBEAM */
			<< "PBEAM   "
			<< std::setw(8) << label    		/* label */
			<< std::setw(8) << 1                 /* material */
			<< std::setw(8) << DRef[S_I].dGet(1, 1);	/* area */
			<< std::setw(8) << DRef[S_I].dGet(5, 5);	/* J1 */
			<< std::setw(8) << DRef[S_I].dGet(6, 6);     /* J2 */
			<< std::setw(8) << DRef[S_I].dGet(4, 5);	/* J12 */
			<< std::setw(8) << DRef[S_I].dGet(4, 4);	/* Jp */
			<< std::endl
			
			/* CBEAM */
			<< "CBEAM   "
			<< std::setw(8) << label   		/* label */
			<< std::setw(8) << label    		/* prop */
			<< std::setw(8) << pNode[NODE1]->GetLabel()	/* node 1 */
			<< std::setw(8) << pNode[NODE2]->GetLabel()	/* node 2 */
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
			<< std::endl

			/* CBEAM */
			<< "CBEAM   "
			<< std::setw(8) << 10000000+label    /* label */
			<< std::setw(8) << label    		/* prop */
			<< std::setw(8) << pNode[NODE2]->GetLabel()	/* node 2 */
			<< std::setw(8) << pNode[NODE3]->GetLabel()	/* node 3 */
			<< std::setw(32) << " "
			<< "+" << std::setw(7) << 1
			<< std::endl
			<< "+" << std::setw(7) << 1
			<< std::setw(16) << " "
			<< std::setw(8) << F2.dGet(1)
			<< std::setw(8) << F2.dGet(2)
			<< std::setw(8) << F2.dGet(3)
			<< std::setw(8) << F3.dGet(1)
			<< std::setw(8) << F3.dGet(2)
			<< std::setw(8) << F3.dGet(3)
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED16__
		out << std::endl
			/* PBEAM */
			<< "PBEAM*  "
			<< std::setw(16) << label   		/* label */
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
			<< std::setw(16) << label   		/* label */
			<< std::setw(16) << label   		/* prop */
			<< std::setw(16) << pNode[NODE1]->GetLabel()	/* node 1 */
			<< std::setw(16) << pNode[NODE2]->GetLabel()	/* node 2 */
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
			<< std::endl
			
			/* CBEAM */
			<< "CBEAM*  "
			<< std::setw(16) << 10000000+label   /* label */
			<< std::setw(16) << label   		/* prop */
			<< std::setw(16) << pNode[NODE2]->GetLabel()	/* node 2 */
			<< std::setw(16) << pNode[NODE3]->GetLabel()	/* node 3 */
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(64) << " "
			<< "*" << std::setw(7) << 2
			<< std::endl
			<< "*" << std::setw(7) << 2
			<< std::setw(32) << " "
			<< std::setw(16) << F2.dGet(1)
			<< std::setw(16) << F2.dGet(2)
			<< "*" << std::setw(7) << 3
			<< std::endl
			<< "*" << std::setw(7) << 3
			<< std::setw(16) << F2.dGet(3)
			<< std::setw(16) << F3.dGet(1)
			<< std::setw(16) << F3.dGet(2)
			<< std::setw(16) << F3.dGet(3)
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
			<< std::endl
			
			/* CBEAM */
			<< "CBEAM,"
			<< 10000000+label << ","
			<< label << ","
			<< pNode[NODE2]->GetLabel() << ","
			<< pNode[NODE3]->GetLabel() << ",,,,"
#if 0
			<< "," 
#endif
			<< std::endl
#if 1
			<< ","
#endif
			<< " ,,", F2.Write(out, ",") << ",", F3.Write(out, ",")
			<< std::endl;
#else
#error "unknown NASTRAN format"
#endif
	}
#endif /* __HACK_NASTRAN_MODES__ */
}

   
/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
Beam::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		    const VectorHandler& XCurr) 
{ 
   DEBUGCOUTFNAME("Beam::InitialAssJac => AssStiffnessMat");

   /* Dimensiona la matrice, la azzera e pone gli indici corretti */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(18, 18, 0.);
   
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutRowIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
      WM.fPutColIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }      
   
   AssStiffnessMat(WM, WM, 1., XCurr, XCurr);
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& Beam::InitialAssRes(SubVectorHandler& WorkVec,
				      const VectorHandler& XCurr) 
{ 
   DEBUGCOUTFNAME("Beam::InitialAssRes => AssStiffnessVec");

   /* Dimensiona il vettore, lo azzera e pone gli indici corretti */
   WorkVec.Resize(18);
   WorkVec.Reset(0.);
   
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   integer iNode3FirstPosIndex = pNode[NODE3]->iGetFirstPositionIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(12+iCnt, iNode3FirstPosIndex+iCnt);
   }      
   
   AssStiffnessVec(WorkVec, 1., XCurr, XCurr);
   return WorkVec;
}


const StructNode* Beam::pGetNode(unsigned int i) const
{
   ASSERT(i >= 1 && i <= 3);
   switch (i) {
    case 1:
    case 2:
    case 3:
      return pNode[i-1];
    default:
      THROW(Beam::ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return NULL;
#endif /* USE_EXCEPTIONS */
}


void
Beam::GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& r) const
{
   ASSERT(part == 1 || part == 2);
   part--;
   
   x = p[part];
   r = R[part];
}


std::ostream& 
Beam::WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part, unsigned int firstId) const
{
   Vec3 xTmp[NUMNODES];

   part--;
   
   for (unsigned int i = part; i <= part+1; i++) {
      xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
   }

   Mat3x3 RT(R[part].Transpose());
   
   out << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << 1+part << std::endl
     << firstId << " "
     << p[part] << " " 
     << MatR2EulerAngles(R[part]) << " "
     << RT*(xTmp[part]-p[part]) << " "
     << Zero3 /* MatR2EulerAngles(pNode[part]->GetRCurr()) */ << " "
     << RT*(xTmp[1+part]-p[part]) << " "
     << Zero3 /* MatR2EulerAngles(pNode[1+part]->GetRCurr()) */ << std::endl;

   return out;
}

/* Beam - end */


#ifdef VISCOELASTIC_BEAM
/* ViscoElasticBeam - begin */

/* Costruttore normale */
ViscoElasticBeam::ViscoElasticBeam(
		unsigned int uL, 
		const StructNode* pN1, 
		const StructNode* pN2, 
		const StructNode* pN3,
		const Vec3& F1,
		const Vec3& F2,
		const Vec3& F3,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& R3,
		const Mat3x3& r_I,
		const Mat3x3& rII,
		const ConstitutiveLaw6D* pD_I, 
		const ConstitutiveLaw6D* pDII,
		flag fOut
) : Elem(uL, Elem::BEAM, fOut),
Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII, fOut)
{   
   SetBeamType(Beam::VISCOELASTIC);

   LPrimeRef[S_I] = LPrime[S_I] = Vec3(0.);  
   gPrime[S_I] = Vec3(0.);
   LPrimeRef[SII] = LPrime[SII] = Vec3(0.); 
   gPrime[SII] = Vec3(0.);
   
   DefPrimeLoc[S_I] = DefPrimeLocRef[S_I] = Vec6(0.);
   DefPrimeLoc[SII] = DefPrimeLocRef[SII] = Vec6(0.);
   
   /* Nota: DsDxi() viene chiamata dal costruttore di Beam */
   Beam::Omega0();
}


/* Costruttore per la trave con forze d'inerzia consistenti */
ViscoElasticBeam::ViscoElasticBeam(
		unsigned int uL, 
		const StructNode* pN1, 
		const StructNode* pN2, 
		const StructNode* pN3,
		const Vec3& F1, 
		const Vec3& F2, 
		const Vec3& F3,
		const Mat3x3& R1,
		const Mat3x3& R2,
		const Mat3x3& R3,
		const Mat3x3& r_I, const Mat3x3& rII,
		const ConstitutiveLaw6D* pD_I, 
		const ConstitutiveLaw6D* pDII,
		doublereal dM_I,
		const Mat3x3& s0_I, const Mat3x3& j0_I,
		doublereal dMII,
		const Mat3x3& s0II, const Mat3x3& j0II,
		flag fOut
) : Elem(uL, Elem::BEAM, fOut),
Beam(uL, pN1, pN2, pN3, F1, F2, F3, R1, R2, R3, r_I, rII, pD_I, pDII,
     dM_I, s0_I, j0_I, dMII, s0II, j0II, fOut)
{
   SetBeamType(Beam::VISCOELASTIC);

   LPrimeRef[S_I] = LPrime[S_I] = Vec3(0.);  
   gPrime[S_I] = Vec3(0.);
   LPrimeRef[SII] = LPrime[SII] = Vec3(0.); 
   gPrime[SII] = Vec3(0.);
   
   DefPrimeLoc[S_I] = DefPrimeLocRef[S_I] = Vec6(0.);
   DefPrimeLoc[SII] = DefPrimeLocRef[SII] = Vec6(0.);

   /* Nota: DsDxi() viene chiamata dal costruttore di Beam */
   Beam::Omega0();
}  


/* Assembla la matrice */
void ViscoElasticBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
				       FullSubMatrixHandler& WMB,
				       doublereal dCoef,
				       const VectorHandler& /* XCurr */ ,
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("ViscoElasticBeam::AssStiffnessMat");
   
   /* La matrice arriva gia' dimensionata e con gli indici di righe e colonne
    * a posto */
   
   /* offset nel riferimento globale */
   Vec3 fTmp[NUMNODES];
   for (unsigned int i = 0; i < NUMNODES; i++) {
      fTmp[i] = pNode[i]->GetRCurr()*f[i];
   }
   
   Mat6x6 AzTmp[NUMSEZ][NUMNODES];
   Mat6x6 AzPrimeTmp[NUMSEZ][NUMNODES];
   
   /* per ogni punto di valutazione: */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {            
      for (unsigned int i = 0; i < NUMNODES; i++) {
	 /* Delta - deformazioni */
	 AzTmp[iSez][i] = AzPrimeTmp[iSez][i]
	   = Mat6x6(Mat3x3(dN3P[iSez][i]*dsdxi[iSez]),
		    Zero3x3,
		    Mat3x3(LRef[iSez]*(dN3[iSez][i])
			   -fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez])),
		    Mat3x3(dN3P[iSez][i]*dsdxi[iSez]));
	 
	 AzTmp[iSez][i] = DRef[iSez]*AzTmp[iSez][i]*dCoef;
	 
	 AzTmp[iSez][i] +=
	   ERef[iSez]*Mat6x6(Mat3x3(OmegaRef[iSez]*(-dN3P[iSez][i]*dsdxi[iSez]*dCoef)),
			     Zero3x3, 
			     (Mat3x3(LPrimeRef[iSez])
			      -Mat3x3(Omega[iSez],LRef[iSez]))*(dN3[iSez][i]*dCoef)
			     +Mat3x3(Omega[iSez], fTmp[i]*(dN3P[iSez][i]*dsdxi[iSez]*dCoef))
			     +Mat3x3(fTmp[i].Cross(pNode[i]->GetWRef()*(dN3P[iSez][i]*dsdxi[iSez]*dCoef))),
			     Mat3x3(OmegaRef[iSez]*(-dN3P[iSez][i]*dsdxi[iSez]*dCoef)));
	 
	 AzPrimeTmp[iSez][i] = ERef[iSez]*AzPrimeTmp[iSez][i];
	 
	 /* Correggo per la rotazione da locale a globale */
	 AzTmp[iSez][i].SubMat12(Mat3x3(AzRef[iSez].GetVec1()*(dN3[iSez][i]*dCoef)));
	 AzTmp[iSez][i].SubMat22(Mat3x3(AzRef[iSez].GetVec2()*(dN3[iSez][i]*dCoef)));
      }
   } /* end ciclo sui punti di valutazione */
   
   
   Vec3 bTmp[2];
   
   /* ciclo sulle equazioni */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      bTmp[0] = p[iSez]-pNode[iSez]->GetXCurr();
      bTmp[1] = p[iSez]-pNode[iSez+1]->GetXCurr();
   
      unsigned int iRow1 = iSez*6+1;
      unsigned int iRow2 = iRow1+6;
      
      for (unsigned int i = 0; i < NUMNODES; i++) {
	 /* Equazione all'indietro: */
	 WMA.Sub(iRow1, 6*i+1, AzTmp[iSez][i].GetMat11());
	 WMA.Sub(iRow1, 6*i+4, AzTmp[iSez][i].GetMat12());
	 
	 WMA.Sub(iRow1+3, 6*i+1,
		 AzTmp[iSez][i].GetMat21()
		 -Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
		 +Mat3x3(bTmp[0])*AzTmp[iSez][i].GetMat11());
	 WMA.Sub(iRow1+3, 6*i+4, 
		 AzTmp[iSez][i].GetMat22()
		 -Mat3x3(AzRef[iSez].GetVec1()*(-dCoef*dN3[iSez][i]), fTmp[i])
		 +Mat3x3(bTmp[0])*AzTmp[iSez][i].GetMat12());
	 
	 /* Equazione in avanti: */
	 WMA.Add(iRow2, 6*i+1, AzTmp[iSez][i].GetMat11());
	 WMA.Add(iRow2, 6*i+4, AzTmp[iSez][i].GetMat12());

	 WMA.Add(iRow2+3, 6*i+1,
		 AzTmp[iSez][i].GetMat21()
		 -Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]))
		 +Mat3x3(bTmp[1])*AzTmp[iSez][i].GetMat11());
	 WMA.Add(iRow2+3, 6*i+4, 
		 AzTmp[iSez][i].GetMat22()
		 +Mat3x3(AzRef[iSez].GetVec1()*(dCoef*dN3[iSez][i]), fTmp[i])
		 +Mat3x3(bTmp[1])*AzTmp[iSez][i].GetMat12());

      
	 /* Equazione viscosa all'indietro: */
	 WMB.Sub(iRow1, 6*i+1, AzPrimeTmp[iSez][i].GetMat11());
	 WMB.Sub(iRow1, 6*i+4, AzPrimeTmp[iSez][i].GetMat12());
	 
	 WMB.Sub(iRow1+3, 6*i+1,
		 AzPrimeTmp[iSez][i].GetMat21()
		 +Mat3x3(bTmp[0])*AzPrimeTmp[iSez][i].GetMat11());
	 WMB.Sub(iRow1+3, 6*i+4,
		 AzPrimeTmp[iSez][i].GetMat22()
		 +Mat3x3(bTmp[0])*AzPrimeTmp[iSez][i].GetMat12());
	 
	 /* Equazione viscosa in avanti: */
	 WMB.Add(iRow2, 6*i+1, AzPrimeTmp[iSez][i].GetMat11());
	 WMB.Add(iRow2, 6*i+4, AzPrimeTmp[iSez][i].GetMat12());

	 WMB.Add(iRow2+3, 6*i+1,
		 AzPrimeTmp[iSez][i].GetMat21()	       
		 +Mat3x3(bTmp[1])*AzPrimeTmp[iSez][i].GetMat11());
	 WMB.Add(iRow2+3, 6*i+4, 
		 AzPrimeTmp[iSez][i].GetMat22()	     
		 +Mat3x3(bTmp[1])*AzPrimeTmp[iSez][i].GetMat12());
      }
      
      /* correzione alle equazioni */
      WMA.Add(iRow1+3, 6*iSez+1, Mat3x3(AzRef[iSez].GetVec1()*(-dCoef)));
      WMA.Add(iRow2+3, 6*iSez+7, Mat3x3(AzRef[iSez].GetVec1()*dCoef));
      
   } /* end ciclo sui punti di valutazione */
};


/* Assembla il residuo */
void ViscoElasticBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
				       doublereal /* dCoef */ ,
				       const VectorHandler& /* XCurr */ ,
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("ViscoElasticBeam::AssStiffnessVec");
   
   /* Riceve il vettore gia' dimensionato e con gli indici a posto 
    * per scrivere il residuo delle equazioni di equilibrio dei tre nodi */
   
   /* Per la trattazione teorica, il riferimento e' il file ul-travi.tex 
    * (ora e' superato) */
   
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
	 xPrimeTmp[i] = pNode[i]->GetVCurr()+pNode[i]->GetWCurr().Cross(fTmp);
      }      
      
      Mat3x3 RDelta[NUMSEZ];
      Vec3 gGrad[NUMSEZ];
      Vec3 gPrimeGrad[NUMSEZ];
      
      /* Aggiorna le grandezze della trave nei punti di valutazione */
      for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
	 
	 /* Posizione */
	 p[iSez] = InterpState(xTmp[NODE1],
			       xTmp[NODE2],
			       xTmp[NODE3], Beam::Section(iSez));
	 
	 /* Matrici di rotazione */
	 g[iSez] = InterpState(gNod[NODE1],
			       gNod[NODE2], 
			       gNod[NODE3], Beam::Section(iSez));
	 RDelta[iSez] = Mat3x3(MatR, g[iSez]);
	 R[iSez] = RDelta[iSez]*RRef[iSez];

	 /* Velocita' angolare della sezione */	 
	 gPrime[iSez] = InterpState(gPrimeNod[NODE1],
				    gPrimeNod[NODE2], 
				    gPrimeNod[NODE3], Beam::Section(iSez));
	 Omega[iSez] = Mat3x3(MatG, g[iSez])*gPrime[iSez]
	   +RDelta[iSez]*OmegaRef[iSez];
	 	 
	 /* Derivate della posizione */
	 L[iSez] = InterpDeriv(xTmp[NODE1],
			       xTmp[NODE2], 
			       xTmp[NODE3], Beam::Section(iSez));
	 
	 /* Derivate della velocita' */
	 LPrime[iSez] = InterpDeriv(xPrimeTmp[NODE1],
				    xPrimeTmp[NODE2], 
				    xPrimeTmp[NODE3], Beam::Section(iSez));
	 
	 /* Derivate dei parametri di rotazione */
	 gGrad[iSez] = InterpDeriv(gNod[NODE1],
				   gNod[NODE2], 
				   gNod[NODE3], Beam::Section(iSez));
	 
	 /* Derivate delle derivate spaziali dei parametri di rotazione */
	 gPrimeGrad[iSez] = InterpDeriv(gPrimeNod[NODE1], 
					gPrimeNod[NODE2], 
					gPrimeNod[NODE3], Beam::Section(iSez));
	 
	 /* Per le deformazioni nel sistema del materiale */
	 Mat3x3 RTmp(R[iSez].Transpose());
	 
	 /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
	 DefLoc[iSez] = Vec6(RTmp*L[iSez]-L0[iSez],
			     RTmp*(Mat3x3(MatG, g[iSez])*gGrad[iSez])+DefLocRef[iSez].GetVec2());

	 /* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
	 DefPrimeLoc[iSez] = Vec6(RTmp*(LPrime[iSez]+L[iSez].Cross(Omega[iSez])),
				  RTmp*(Mat3x3(MatG, g[iSez])*gPrimeGrad[iSez]
					+(Mat3x3(MatG, g[iSez])*gGrad[iSez]).Cross(Omega[iSez]))
				  +DefPrimeLocRef[iSez].GetVec2());
	 
	 /* Calcola le azioni interne */
	 pD[iSez]->Update(DefLoc[iSez], DefPrimeLoc[iSez]);
	 AzLoc[iSez] = pD[iSez]->GetF();
	 
	 /* corregge le azioni interne locali (piezo, ecc) */
	 AddInternalForces(AzLoc[iSez], iSez);

	 /* Porta le azioni interne nel sistema globale */
	 Az[iSez] = MultRV(AzLoc[iSez], R[iSez]);
      }      
   }
   
   WorkVec.Add(1, Az[S_I].GetVec1());
   WorkVec.Add(4, (p[S_I]-xNod[NODE1]).Cross(Az[S_I].GetVec1())
	          +Az[S_I].GetVec2());
   WorkVec.Add(7, Az[SII].GetVec1()-Az[S_I].GetVec1());
   WorkVec.Add(10, Az[SII].GetVec2()-Az[S_I].GetVec2()
	           +(p[SII]-xNod[NODE2]).Cross(Az[SII].GetVec1())
	           -(p[S_I]-xNod[NODE2]).Cross(Az[S_I].GetVec1()));
   WorkVec.Sub(13, Az[SII].GetVec1());
   WorkVec.Add(16, Az[SII].GetVec1().Cross(p[SII]-xNod[NODE3])
	           -Az[SII].GetVec2());   
}


/* Settings iniziali, prima della prima soluzione */
void ViscoElasticBeam::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   Beam::SetValue(X, XP);
   
   /* Aggiorna le grandezze della trave nei punti di valutazione */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      (Vec3&)OmegaRef[iSez] = Omega[iSez];
      (Vec3&)LPrimeRef[iSez] = LPrime[iSez];
      (Vec6&)DefPrimeLocRef[iSez] = DefPrimeLoc[iSez];
      
      /* Aggiorna il legame costitutivo di riferimento
       * (la deformazione e' gia' stata aggiornata dall'ultimo residuo) */
      (Mat6x6&)ERef[iSez] = MultRMRt(pD[iSez]->GetFDEPrime(), RRef[iSez]);
   }
   ASSERT(fFirstRes == flag(1));
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void ViscoElasticBeam::AfterPredict(VectorHandler& /* X */ , 
				    VectorHandler& /* XP */ )
{
   /* Calcola le deformazioni, aggiorna il legame costitutivo e crea la FDE */
   
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
      xPrimeTmp[i] = pNode[i]->GetVCurr()+pNode[i]->GetWRef().Cross(fTmp);
   }      
      
   Mat3x3 RDelta[NUMSEZ];
   Vec3 gGrad[NUMSEZ];
   Vec3 gPrimeGrad[NUMSEZ];
   
   /* Aggiorna le grandezze della trave nei punti di valutazione */
   for (unsigned int iSez = 0; iSez < NUMSEZ; iSez++) {
      /* Posizione */
      p[iSez] = InterpState(xTmp[NODE1],
			    xTmp[NODE2],
			    xTmp[NODE3], Beam::Section(iSez));
      
      /* Matrici di rotazione */
      g[iSez] = InterpState(gNod[NODE1],
			    gNod[NODE2], 
			    gNod[NODE3], Beam::Section(iSez));
      RDelta[iSez] = Mat3x3(MatR, g[iSez]);
      R[iSez] = RRef[iSez] = RDelta[iSez]*R[iSez];
      
      /* Velocita' angolare della sezione */	 
      gPrime[iSez] = InterpState(gPrimeNod[NODE1],
				 gPrimeNod[NODE2], 
				 gPrimeNod[NODE3], Beam::Section(iSez));
      Omega[iSez] = OmegaRef[iSez] = Mat3x3(MatG, g[iSez])*gPrime[iSez];
      
      /* Derivate della posizione */
      L[iSez] = LRef[iSez] = InterpDeriv(xTmp[NODE1],
					 xTmp[NODE2], 
					 xTmp[NODE3], Beam::Section(iSez));
      
      /* Derivate della velocita' */
      LPrime[iSez] = LPrimeRef[iSez] 
	= InterpDeriv(xPrimeTmp[NODE1],
		      xPrimeTmp[NODE2], 
		      xPrimeTmp[NODE3], Beam::Section(iSez));
      
      /* Derivate dei parametri di rotazione */
      gGrad[iSez] = InterpDeriv(gNod[NODE1],
				gNod[NODE2], 
				gNod[NODE3], Beam::Section(iSez));
      
      /* Derivate delle derivate spaziali dei parametri di rotazione */
      gPrimeGrad[iSez] = InterpDeriv(gPrimeNod[NODE1], 
				     gPrimeNod[NODE2], 
				     gPrimeNod[NODE3], Beam::Section(iSez));
      
      /* Per le deformazioni nel sistema del materiale */
      Mat3x3 RTmp(R[iSez].Transpose());
      
      /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
      DefLoc[iSez] = DefLocRef[iSez] 
	= Vec6(RTmp*L[iSez]-L0[iSez],
	       RTmp*(Mat3x3(MatG, g[iSez])*gGrad[iSez])+DefLoc[iSez].GetVec2());
      
      /* Calcola le velocita' di deformazione nel sistema locale nei punti di valutazione */
      DefPrimeLoc[iSez] = DefPrimeLocRef[iSez] 
	= Vec6(RTmp*(LPrime[iSez]+L[iSez].Cross(Omega[iSez])),
	       RTmp*(Mat3x3(MatG, g[iSez])*gPrimeGrad[iSez]
		     +(Mat3x3(MatG, g[iSez])*gGrad[iSez]).Cross(Omega[iSez])));
      
      /* Calcola le azioni interne */
      pD[iSez]->Update(DefLoc[iSez], DefPrimeLoc[iSez]);
      AzLoc[iSez] = pD[iSez]->GetF();
      
      /* corregge le azioni interne locali (piezo, ecc) */
      AddInternalForces(AzLoc[iSez], iSez);
      
      AzLocRef[iSez] = AzLoc[iSez];
      
      /* Porta le azioni interne nel sistema globale */
      Az[iSez] = AzRef[iSez] = MultRV(AzLoc[iSez], R[iSez]);
      
      /* Aggiorna il legame costitutivo */
      DRef[iSez] = MultRMRt(pD[iSez]->GetFDE(), R[iSez]);
      ERef[iSez] = MultRMRt(pD[iSez]->GetFDEPrime(), R[iSez]);
   }
   
   fFirstRes = flag(1);
}

/* ViscoElasticBeam - end */
#endif /* VISCOELASTIC_BEAM */


/* Legge una trave */

Elem* ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadBeam");
   
   const char* sKeyWords[] = {
      "elastic",
	"viscoelastic",
	"piezoelectric"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	ELASTIC = 0,
	VISCOELASTIC,
	PIEZOELECTRIC,
	LASTKEYWORD 
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
         
   /* Per ogni nodo: */
   
   /* Nodo 1 */
   StructNode* pNode1 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
   
   Mat3x3 R1(pNode1->GetRCurr());   
   Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
   Mat3x3 Rn1(Eye3);
   if (HP.IsKeyWord("rot")) {
	   Rn1 = HP.GetRotRel(ReferenceFrame(pNode1));
   }
        
   DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): " 
	      << f1 << std::endl
	      << "(global frame): "
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
	      << f2 << std::endl
	      << "(global frame): "
	      << pNode2->GetXCurr()+pNode2->GetRCurr()*f2 << std::endl);
   
   /* Nodo 3 */
   StructNode* pNode3 = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
   
   Mat3x3 R3(pNode3->GetRCurr());   
   Vec3 f3(HP.GetPosRel(ReferenceFrame(pNode3)));
   Mat3x3 Rn3(Eye3);
   if (HP.IsKeyWord("rot")) {
	   Rn3 = HP.GetRotRel(ReferenceFrame(pNode3));
   }
   
   DEBUGLCOUT(MYDEBUG_INPUT, "node 3 offset (node reference frame): " 
	      << f3 << std::endl
	      << "(global frame): " 
	      << pNode3->GetXCurr()+pNode3->GetRCurr()*f3 << std::endl);
   
   /*   Per ogni punto: */
   
   /* Punto I */
   
   /*     Matrice R */
   Mat3x3 R_I;
   flag f_I(0);
   if (HP.IsKeyWord("node")) {	
      f_I = flag(1);
   } else {     
      R_I = HP.GetRotAbs(AbsRefFrame);
   }   
   
   /*     Legame costitutivo */
   DefHingeType::Type ConstLawType_I = DefHingeType::UNKNOWN;
   ConstitutiveLaw6D* pD_I = pDM->ReadConstLaw6D(HP, ConstLawType_I);
   
   if (pD_I->iGetNumDof() != 0) {
   	   std::cerr << "line " << HP.GetLineData()
		   << ": beam does not support "
		   "dynamic constitutive laws yet"
		   << std::endl;
	   THROW(ErrGeneric());
   }
	
#ifdef DEBUG   
   Mat6x6 MTmp(pD_I->GetFDE());
   Mat3x3 D11(MTmp.GetMat11());
   Mat3x3 D12(MTmp.GetMat12());
   Mat3x3 D21(MTmp.GetMat21());
   Mat3x3 D22(MTmp.GetMat22());
   
   DEBUGLCOUT(MYDEBUG_INPUT, 
	      "First point matrix D11: " << std::endl << D11 << std::endl
	      << "First point matrix D12: " << std::endl << D12 << std::endl
	      << "First point matrix D21: " << std::endl << D21 << std::endl
	      << "First point matrix D22: " << std::endl << D22 << std::endl);
#endif   

      
   /* Punto II */

   /*     Matrice R */
   Mat3x3 RII;
   flag fII(0);
   if (HP.IsKeyWord("same")) {
      if (f_I) {	     
	 fII = flag(1);
      } else {	     
	 RII = R_I;
      }	
   } else {
      if (HP.IsKeyWord("node")) {	     
	 fII = flag(1);
      } else {	            
	 RII = HP.GetRotAbs(AbsRefFrame);
      }	
   }   

   /*     Legame costitutivo */
   DefHingeType::Type ConstLawTypeII = DefHingeType::UNKNOWN;
   ConstitutiveLaw6D* pDII = NULL;
   
   /* Not allowed any more, since there is no simple way implement
    * the duplication of a constitutive law without a "virtual constructor"
    * for the constitutive law, the drivers and so on */
   
   if (HP.IsKeyWord("same")) {
      /*
      std::cerr << "Sorry, 'same' is not supported any more, you must enter the constitutive law of the second section" << std::endl;
      THROW(DataManager::ErrGeneric());
       */
      pDII = pD_I->pCopy();
   } else {
      pDII = pDM->ReadConstLaw6D(HP, ConstLawTypeII);
      
      if (pDII->iGetNumDof() != 0) {
   	      std::cerr << "line " << HP.GetLineData()
		      << ": beam does not support "
		      "dynamic constitutive laws yet"
		      << std::endl;
	      THROW(ErrGeneric());
      }
   }

#ifdef DEBUG   
   MTmp = pDII->GetFDE();
   D11 = MTmp.GetMat11();
   D12 = MTmp.GetMat12();
   D21 = MTmp.GetMat21();
   D22 = MTmp.GetMat22();
   
   DEBUGLCOUT(MYDEBUG_INPUT, 
	      "Second point matrix D11: " << std::endl << D11 << std::endl
	      << "Second point matrix D12: " << std::endl << D12 << std::endl
	      << "Second point matrix D21: " << std::endl << D21 << std::endl
	      << "Second point matrix D22: " << std::endl << D22 << std::endl);
#endif               

#if defined(USE_ELECTRIC_NODES)
   flag fPiezo(0);
   Mat3xN PiezoMat[2][2];
   integer iNumElec = 0;
   ScalarDifferentialNode** pvElecDofs = NULL;
   if (HP.IsKeyWord("piezoelectricactuator")) {
      fPiezo = flag(1);
      DEBUGLCOUT(MYDEBUG_INPUT, 
		 "Piezoelectric actuator beam is expected" << std::endl);
      
      iNumElec = HP.GetInt();
      DEBUGLCOUT(MYDEBUG_INPUT, 
		 "piezo actuator " << uLabel << " has " << iNumElec 
		 << " electrodes" << std::endl);
      if (iNumElec <= 0) {
	 std::cerr << "illegal number of electric nodes " << iNumElec << std::endl;
	 THROW(ErrGeneric());
      }
      
      SAFENEWARR(pvElecDofs, ScalarDifferentialNode*, iNumElec);
      
      for (integer i = 0; i < iNumElec; i++) {
	 unsigned int uL = HP.GetInt();
	 DEBUGLCOUT(MYDEBUG_INPUT, "linked to abstract node " << uL << std::endl);
	 pvElecDofs[i] = (ScalarDifferentialNode*)(pDM->pFindNode(Node::ABSTRACT, uL));
	 if (pvElecDofs[i] == NULL) {
	    std::cerr << "can't find abstract node " << uL << std::endl;
	    THROW(ErrGeneric());
	 }
      }

      PiezoMat[0][0].Resize(iNumElec);
      PiezoMat[1][0].Resize(iNumElec);
      PiezoMat[0][1].Resize(iNumElec);
      PiezoMat[1][1].Resize(iNumElec);
      
      /* leggere le matrici (6xN sez. 1, 6xN sez. 2) */
      HP.GetMat6xN(PiezoMat[0][0], PiezoMat[1][0], iNumElec);
      if (HP.IsKeyWord("same")) {
	 PiezoMat[0][1].Copy(PiezoMat[0][0]);
	 PiezoMat[1][1].Copy(PiezoMat[1][0]);
      } else {	 
	 HP.GetMat6xN(PiezoMat[0][1], PiezoMat[1][1], iNumElec);
      }      
      
      /*
      DEBUGLCOUT(MYDEBUG_INPUT, "Piezo matrix I:" << std::endl << PiezoMat[0][0] << PiezoMat[1][0]);
      DEBUGLCOUT(MYDEBUG_INPUT, "Piezo matrix II:" << std::endl << PiezoMat[0][1] << PiezoMat[1][1]);
       */
   }
#endif /* defined(USE_ELECTRIC_NODES) */
   
   
   flag fOut = pDM->fReadOutput(HP, Elem::BEAM);       
   
   
   /* Se necessario, interpola i parametri di rotazione delle sezioni */
   if (f_I || fII) {
      Mat3x3 RT((R2*Rn2).Transpose());
      Vec3 g1(gparam(RT*(R1*Rn1)));
      Vec3 g3(gparam(RT*(R3*Rn3)));
      if (f_I) {
	 R_I = R2*Mat3x3(Beam::InterpState(g1, 0., g3, Beam::S_I));
      }
      if (fII) {
	 R_I = R2*Mat3x3(Beam::InterpState(g1, 0., g3, Beam::SII));
      }	
   }


   Elem* pEl = NULL;
   
   if ((ConstLawType_I == DefHingeType::ELASTIC) 
       && (ConstLawType_I == DefHingeType::ELASTIC)) {

#if defined(USE_ELECTRIC_NODES)      
      if (fPiezo == 0) {	 
#endif /* defined(USE_ELECTRIC_NODES) */
	 SAFENEWWITHCONSTRUCTOR(pEl,
				Beam,
				Beam(uLabel,
				     pNode1, pNode2, pNode3,
				     f1, f2, f3,
				     Rn1, Rn2, Rn3,
				     R_I, RII,
				     pD_I, pDII,
				     fOut));
#if defined(USE_ELECTRIC_NODES)
      } else {	 
	 SAFENEWWITHCONSTRUCTOR(pEl,
				PiezoActuatorBeam,
				PiezoActuatorBeam(uLabel,
						  pNode1, pNode2, pNode3,
						  f1, f2, f3,
						  Rn1, Rn2, Rn3,
						  R_I, RII,
						  pD_I, pDII,
						  iNumElec,
						  pvElecDofs,
						  PiezoMat[0][0], PiezoMat[1][0],
						  PiezoMat[0][1], PiezoMat[1][1],
						  fOut));
      }
#endif /* defined(USE_ELECTRIC_NODES) */
      
     
   } else /* At least one is VISCOUS or VISCOELASTIC */ {
#ifdef VISCOELASTIC_BEAM      
#if defined(USE_ELECTRIC_NODES)      
      if (fPiezo == 0) {	 
#endif /* defined(USE_ELECTRIC_NODES) */
	 SAFENEWWITHCONSTRUCTOR(pEl, 
				ViscoElasticBeam,
				ViscoElasticBeam(uLabel,
						 pNode1, pNode2, pNode3,
						 f1, f2, f3,
						 Rn1, Rn2, Rn3,
						 R_I, RII,
						 pD_I, pDII,
						 fOut));
#if defined(USE_ELECTRIC_NODES)
      } else {	 
	 SAFENEWWITHCONSTRUCTOR(pEl,
				PiezoActuatorVEBeam,
				PiezoActuatorVEBeam(uLabel,
						    pNode1, pNode2, pNode3,
						    f1, f2, f3,
						    Rn1, Rn2, Rn3,
						    R_I, RII,
						    pD_I, pDII,
						    iNumElec,
						    pvElecDofs,
						    PiezoMat[0][0], PiezoMat[1][0],
						    PiezoMat[0][1], PiezoMat[1][1],
						    fOut));
      }
#endif /* defined(USE_ELECTRIC_NODES) */

#else /* VISCOELASTIC_BEAM */
      std::cerr << "Sorry, the ViscoElasticBeam element is not available yet" << std::endl;
      THROW(ErrNotImplementedYet());
#endif /* VISCOELASTIC_BEAM */
   }
   
   /* Costruttore normale
    * Beam(unsigned int uL,
    *      const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
    *	   const Vec3& X1, const Vec3& X2, const Vec3& X3,
    *	   const Vec3& F1, const Vec3& F2, const Vec3& F3,
    *	   const Mat3x3& r_I, const Mat3x3& rII,
    *	   const Mat3x3& d11_I, const Mat3x3& d12_I,
    *	   const Mat3x3& d21_I, const Mat3x3& d22_I,
    *	   const Mat3x3& d11II, const Mat3x3& d12II,
    *	   const Mat3x3& d21II, const Mat3x3& d22II,
    *	   const Vec3& eps0_I, const Vec3& k0_I,
    *	   const Vec3& eps0II, const Vec3& k0II);
    */
   

   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      std::cerr << std::endl
	<< "semicolon expected at line " << HP.GetLineData() << std::endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of ReadBeam() */

