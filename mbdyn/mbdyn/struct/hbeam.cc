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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <constltp.h>
#include <beam.h>
#include <hbeam.h>
#include <dataman.h>

static doublereal dN2[2] = {
	.5, 
	.5
};

static doublereal dN2P[2] = {
	-.5,
	.5
};

/* Funzioni di interpolazione */
doublereal 
ShapeFunc2N(doublereal d, integer iNode, integer iOrd)
{
    ASSERT(iOrd == 0 || iOrd == 1);
    ASSERT(iNode == 1 || iNode == 2);
   
    switch (iOrd) {
    case 0:
        switch (iNode) {
	case 1:	
	    return .0;
				
	case 2:	
	    return 0.;
		  		
	default:
	    THROW(ErrGeneric());
	}
	    
    case 1:
        switch (iNode) {
	case 1:		
	    return 0.;
	  
	case 2:
	    return 0.;
	  
	default:
	    THROW(ErrGeneric());
	}
    }
   
    /* Per evitare warnings */
    return 0.;
}

doublereal
DxDcsi2N(doublereal d, const Vec3& X1, const Vec3& X2)
{
    doublereal dNp1 = ShapeFunc3N(d, 1, 1);
    doublereal dNp2 = ShapeFunc3N(d, 2, 1);
    Vec3 DXDcsi(X1*dNp1+X2*dNp2);
    doublereal dd = DXDcsi.Dot();
    if (dd > DBL_EPSILON) {
        return sqrt(dd);
    }
    return 0.;
}


/* HBeam - begin */

/* Costruttore normale */
HBeam::HBeam(unsigned int uL, 
	   const StructNode* pN1,
	   const StructNode* pN2,
	   const Vec3& F1,
	   const Vec3& F2,
	   const Mat3x3& r,
	   const ConstitutiveLaw6D* pd,
	   flag fOut)
: Elem(uL, Elem::BEAM, fOut), 
InitialAssemblyElem(uL, Elem::BEAM, fOut),
fFirstRes(1)
{
    pNode[NODE1] = pN1;
    pNode[NODE2] = pN2;
    (Vec3&)f[NODE1] = F1;
    (Vec3&)f[NODE2] = F2;
    RRef = R = (Mat3x3&)r;
  
    pD = NULL; 
    SAFENEWWITHCONSTRUCTOR(pD, 
			   ConstitutiveLaw6DOwner,
			   ConstitutiveLaw6DOwner(pd),
			   DMmm);
   
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
      
    /* Aggiorna le grandezze della trave nei punti di valutazione */
    p = InterpState(xTmp[NODE1], xTmp[NODE2]);
}


HBeam::~HBeam(void) 
{
    ASSERT(pD != NULL);
    if (pD != NULL) {      
        SAFEDELETE(pD, DMmm);
    }
}


Vec3 
HBeam::InterpState(const Vec3& v1, const Vec3& v2)
{
    doublereal* pv1 = (doublereal*)v1.pGetVec();
    doublereal* pv2 = (doublereal*)v2.pGetVec();
    return Vec3(pv1[0]*dN2[0]+pv2[0]*dN2[1],
	        pv1[1]*dN2[0]+pv2[1]*dN2[1],
	        pv1[2]*dN2[0]+pv2[2]*dN2[1]);
}


Vec3
HBeam::InterpDeriv(const Vec3& v1, const Vec3& v2)
{
    doublereal* pv1 = (doublereal*)v1.pGetVec();
    doublereal* pv2 = (doublereal*)v2.pGetVec();
    return Vec3((pv1[0]*dN2P[0]+pv2[0]*dN2P[1])*dsdxi,
	        (pv1[1]*dN2P[0]+pv2[1]*dN2P[1])*dsdxi,
	        (pv1[2]*dN2P[0]+pv2[2]*dN2P[1])*dsdxi);
}


void
HBeam::DsDxi(void)
{
    /* Validazione dati */
    ASSERT(pNode[NODE1] != NULL);
    ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE2] != NULL);
    ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
   
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

    Vec3 xGrad;
    xGrad = InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
    doublereal d = xGrad.Dot();
    ASSERT(d > DBL_EPSILON);
    if (d > DBL_EPSILON) {
       dsdxi = 1./sqrt(d);
    } else {
       cerr << "warning, beam " << GetLabel() 
	 << " has singular metric; aborting ..." << endl;

       THROW(HBeam::ErrGeneric());
    }

    /* Calcola le deformazioni iniziali */
    L0 = R.Transpose()*InterpDeriv(xTmp[NODE1], xTmp[NODE2]);
    pD->Update(0.);
    DRef = MultRMRt(pD->GetFDE(), R);
}


/* Calcola la velocita' angolare delle sezioni a partire da quelle dei nodi */
void 
HBeam::Omega0(void)
{   
    /* Validazione dati */
    ASSERT(pNode[NODE1] != NULL);
    ASSERT(pNode[NODE1]->GetNodeType() == Node::STRUCTURAL);
    ASSERT(pNode[NODE2] != NULL);
    ASSERT(pNode[NODE2]->GetNodeType() == Node::STRUCTURAL);
 
    /* Modo consistente: */      
    Mat3x3 RNod[NUMNODES];
    Vec3 w[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {     
        RNod[i] = pNode[i]->GetRCurr();
        w[i] = pNode[i]->GetWCurr();
    }
   
#if 0
    /* Modo brutale: interpolo le velocita' dei nodia */
    Vec3 w[NUMNODES];
    for (unsigned int i = 0; i < NUMNODES; i++) {
        w[i] = pNode[i]->GetWCurr();
    }
    for (unsigned int i = 0; i < NUMSEZ; i++) {      
        Omega[i] = (w{NODE1]*dN2[i][NODE1]
	            +w[NODE2]*dN2[i][NODE2]
		    +w[NODE3]*dN2[i][NODE3]);
    }
#endif /* 0 */
}


/* Contributo al file di restart */
ostream& HBeam::Restart(ostream& out) const
{
   return Restart_(out)<< ';' << endl;
}

ostream& HBeam::Restart_(ostream& out) const
{ 
   out << "  beam: " << GetLabel();
   for (unsigned int i = 0; i < NUMNODES; i++) {
      out << ", " << pNode[i]->GetLabel() << ", reference, node, ", 
	f[i].Write(out, ", ");
   }
   out << ", reference, global"
     << ", 1, ", (R.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R.GetVec(2)).Write(out, ", ")
     << ", ", pD->pGetConstLaw()->Restart(out);

   return out;
}


/* Assembla la matrice */
void HBeam::AssStiffnessMat(FullSubMatrixHandler& WMA,
			   FullSubMatrixHandler& /* WMB */ ,
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("HBeam::AssStiffnessMat");
   
   /* La matrice arriva gia' dimensionata e con gli indici di righe e colonne
    * a posto */
   
   /* offset nel riferimento globale */
   Vec3 fTmp[NUMNODES];
   for (unsigned int i = 0; i < NUMNODES; i++) {
      fTmp[i] = pNode[i]->GetRCurr()*f[i];
   }
   
   Mat6x6 AzTmp[NUMNODES];
};


/* Assembla il residuo */
void HBeam::AssStiffnessVec(SubVectorHandler& WorkVec,
			   doublereal /* dCoef */ ,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUTFNAME("HBeam::AssStiffnessVec");
   
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
      
      Mat3x3 RDelta;
      Vec3 gGrad;
      
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
      
      /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
      DefLoc = Vec6(RTmp*L-L0,
        RTmp*(Mat3x3(MatG, g)*gGrad)+DefLocRef.GetVec2());
      
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
VariableSubMatrixHandler& HBeam::AssJac(VariableSubMatrixHandler& WorkMat,
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
   WM.ResizeInit(12, 12, 0.);	
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }      
   
   AssStiffnessMat(WM, WM, dCoef, XCurr, XPrimeCurr);
   
   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& HBeam::AssRes(SubVectorHandler& WorkVec,
			       doublereal dCoef,
			       const VectorHandler& XCurr, 
			       const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("HBeam::AssRes => AssStiffnessVec");
   
   integer iNode1FirstMomIndex = pNode[NODE1]->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode[NODE2]->iGetFirstMomentumIndex();

   /* Dimensiona il vettore, lo azzera e pone gli indici corretti */
   WorkVec.Resize(12);
   WorkVec.Reset(0.);
         
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }      
   
   AssStiffnessVec(WorkVec, dCoef, XCurr, XPrimeCurr);

   return WorkVec;
}

    
/* Settings iniziali, prima della prima soluzione */
void HBeam::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const
{
   /* Aggiorna le grandezze della trave nei punti di valutazione */
   (Mat3x3&)RRef = R;
   (Vec3&)LRef = L;
   (Vec6&)DefLocRef = DefLoc;
   (Vec6&)AzLocRef = AzLoc;
   (Vec6&)AzRef = Az;
      
   /* Aggiorna il legame costitutivo di riferimento
    * (la deformazione e' gia' stata aggiornata dall'ultimo residuo) */
   (Mat6x6&)DRef = MultRMRt(pD->GetFDE(), RRef);      
   (flag&)fFirstRes = flag(1);
}
              

/* Prepara i parametri di riferimento dopo la predizione */
void HBeam::AfterPredict(VectorHandler& /* X */ , 
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
      
   Mat3x3 RDelta;
   Vec3 gGrad;
   
   /* Aggiorna le grandezze della trave nei punti di valutazione */
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
   
   /* Calcola le deformazioni nel sistema locale nei punti di valutazione */
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


/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void HBeam::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      OH.Beams() << setw(8) << GetLabel() << " " 
	<< AzLoc.GetVec1() << " " << AzLoc.GetVec2() << endl;
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
   WM.ResizeInit(12, 12, 0.);
   
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.fPutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }      
   
   AssStiffnessMat(WM, WM, 1., XCurr, XCurr);
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& HBeam::InitialAssRes(SubVectorHandler& WorkVec,
				      const VectorHandler& XCurr) 
{ 
   DEBUGCOUTFNAME("HBeam::InitialAssRes => AssStiffnessVec");

   /* Dimensiona il vettore, lo azzera e pone gli indici corretti */
   WorkVec.Resize(12);
   WorkVec.Reset(0.);
   
   integer iNode1FirstPosIndex = pNode[NODE1]->iGetFirstPositionIndex();
   integer iNode2FirstPosIndex = pNode[NODE2]->iGetFirstPositionIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }      
   
   AssStiffnessVec(WorkVec, 1., XCurr, XCurr);
   return WorkVec;
}


const StructNode* HBeam::pGetNode(unsigned int i) const
{
   ASSERT(i >= 1 && i <= 2);
   switch (i) {
    case 1:
    case 2:
      return pNode[i-1];
    default:
      THROW(HBeam::ErrGeneric());
   }
#ifndef USE_EXCEPTIONS
   return NULL;
#endif /* USE_EXCEPTIONS */
}


void
HBeam::GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& r) const
{
   ASSERT(part == 1);
   part--;
   
   x = p;
   r = R;
}


ostream& 
HBeam::WriteAdamsDummyPartCmd(ostream& out, unsigned int part, unsigned int firstId) const
{
   Vec3 xTmp[NUMNODES];

   for (unsigned int i = 0; i <= 1; i++) {
      xTmp[i] = pNode[i]->GetXCurr()+pNode[i]->GetRCurr()*f[i];
   }

   Mat3x3 RT(R.Transpose());
   
   out << psAdamsElemCode[GetElemType()] << "_" << GetLabel() << "_" << 1 << endl
     << firstId << " "
     << p << " " 
     << EulerAngles(R) << " "
     << RT*(xTmp[0]-p) << " "
     << Zero3 /* EulerAngles(pNode[part]->GetRCurr()) */ << " "
     << RT*(xTmp[1]-p) << " "
     << Zero3 /* EulerAngles(pNode[1+part]->GetRCurr()) */ << endl;

   return out;
}

/* Beam - end */


/* Legge una trave */

Elem* ReadHBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadHBeam");
   
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
   unsigned int uNode = (unsigned int)HP.GetInt();
   
   DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << endl);
   
   /* verifica di esistenza del nodo */
   StructNode* pNode1;
   if ((pNode1 = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << endl
	<< " at line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;      
      THROW(DataManager::ErrGeneric());
   }		     
   
   Mat3x3 R1(pNode1->GetRCurr());   
   Vec3 f1(HP.GetPosRel(ReferenceFrame(pNode1)));
        
   DEBUGLCOUT(MYDEBUG_INPUT, "node 1 offset (node reference frame): " 
	      << f1 << endl
	      << "(global frame): "
	      << pNode1->GetXCurr()+pNode1->GetRCurr()*f1 << endl);
   
   
   /* Nodo 2 */
   uNode = (unsigned int)HP.GetInt();
   
   DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << endl);
   
   /* verifica di esistenza del nodo */
   StructNode* pNode2;
   if ((pNode2 = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << endl
	<< "at line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;      
      THROW(DataManager::ErrGeneric());
   }		     
   
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 f2(HP.GetPosRel(ReferenceFrame(pNode2)));
         
   DEBUGLCOUT(MYDEBUG_INPUT, "node 2 offset (node reference frame): " 
	      << f2 << endl
	      << "(global frame): "
	      << pNode2->GetXCurr()+pNode2->GetRCurr()*f2 << endl);
   
   /* Punto */
   
   /*     Matrice R */
   Mat3x3 R;
   flag f(0);
   if (HP.IsKeyWord("node")) {	
      f = flag(1);
   } else {     
      R = HP.GetRotAbs(AbsRefFrame);
   }   
   
   /* Legame costitutivo */
   DefHingeType::Type ConstLawType = DefHingeType::UNKNOWN;
   ConstitutiveLaw6D* pD = pDM->ReadConstLaw6D(HP, ConstLawType);
   
#ifdef DEBUG   
   Mat6x6 MTmp(pD->GetFDE());
   Mat3x3 D11(MTmp.GetMat11());
   Mat3x3 D12(MTmp.GetMat12());
   Mat3x3 D21(MTmp.GetMat21());
   Mat3x3 D22(MTmp.GetMat22());
   
   DEBUGLCOUT(MYDEBUG_INPUT, 
	      "First point matrix D11: " << endl << D11 << endl
	      << "First point matrix D12: " << endl << D12 << endl
	      << "First point matrix D21: " << endl << D21 << endl
	      << "First point matrix D22: " << endl << D22 << endl);
#endif   

   flag fOut = pDM->fReadOutput(HP, Elem::BEAM);       
   
   /* Se necessario, interpola i parametri di rotazione delle sezioni */
   if (f) {
      R = Eye3;
   }

   Elem* pEl = NULL;
   
   if (ConstLawType == DefHingeType::ELASTIC) {

      SAFENEWWITHCONSTRUCTOR(pEl,
				HBeam,
				HBeam(uLabel,
				     pNode1, pNode2,
				     f1, f2,
				     R,
				     pD,
				     fOut),
				DMmm);
   } else /* At least one is VISCOUS or VISCOELASTIC */ {
      cerr << "Sorry, the ViscoElasticHBeam element is not available yet" << endl;
      THROW(ErrNotImplementedYet());
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << endl
	<< "semicolon expected at line " << HP.GetLineData() << endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of ReadHBeam() */

