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

/* elementi di massa */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <body.h>
#include <dataman.h>

Body::Body(unsigned int uL, 
	   const StructNode* pNodeTmp, 
           doublereal dMassTmp, 
	   const Vec3& XgcTmp, 
	   const Mat3x3& JTmp, 
	   flag fOut)
: Elem(uL, Elem::BODY, fOut), 
ElemGravityOwner(uL, Elem::BODY, fOut), 
InitialAssemblyElem(uL, Elem::BODY, fOut),
pNode(pNodeTmp), 
dMass(dMassTmp), 
Xgc(XgcTmp), 
S0(XgcTmp*dMassTmp), 
J0(JTmp)
{
    ASSERT(pNode != NULL);
    ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
    ASSERT(dMassTmp > 0.);
}


/* momento statico */
Vec3 
Body::_GetS(void) const
{
	return pNode->GetXCurr()*dMass+pNode->GetRCurr()*S0;
}


/* momento d'inerzia */
Mat3x3
Body::_GetJ(void) const
{
	Vec3 s = pNode->GetRCurr()*S0;
	const Vec3& x = pNode->GetXCurr();

	return pNode->GetRCurr()*J0*pNode->GetRCurr().Transpose()
		- Mat3x3(x, x*dMass) - Mat3x3(s, x) - Mat3x3(x, s);
}
 

/* Scrive il contributo dell'elemento al file di restart */
std::ostream& 
Body::Restart(std::ostream& out) const
{
    out << "  body: " << GetLabel() << ", " 
        << pNode->GetLabel() << ", " << dMass << ", "
	<< "reference, node, ", Xgc.Write(out, ", ") << ", "
        << "reference, node, ", J0.Write(out, ", ") << ';' << std::endl;
   
    return out;
}


VariableSubMatrixHandler& 
Body::AssJac(VariableSubMatrixHandler& WorkMat,
	     doublereal dCoef,
	     const VectorHandler& XCurr,
	     const VectorHandler& XPrimeCurr)
{
    DEBUGCOUTFNAME("Body::AssJac");

    /* Casting di WorkMat */
    FullSubMatrixHandler& WM = WorkMat.SetFull();
   
    /* Dimensiona e resetta la matrice di lavoro */
    WM.ResizeInit(6, 6, 0.);
      
    /* Setta gli indici della matrice - le incognite sono ordinate come:
     *   - posizione (3)
     *   - parametri di rotazione (3)
     *   - quantita' di moto (3)
     *   - momento della quantita' di moto 
     * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
     * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
     * e' dato da iGetFirstPositionIndex()+i */
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= 6; iCnt++) {
        WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);   
    }
      
    AssMat_(WM, WM, dCoef, XCurr, XPrimeCurr);
    return WorkMat;
}


void
Body::AssMats(VariableSubMatrixHandler& WorkMatA,
	     VariableSubMatrixHandler& WorkMatB,	       
	     const VectorHandler& XCurr,
	     const VectorHandler& XPrimeCurr)
{
    DEBUGCOUTFNAME("Body::AssMats");

    /* Casting di WorkMat */
    FullSubMatrixHandler& WMA = WorkMatA.SetFull();
    FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   
    /* Dimensiona e resetta la matrice di lavoro */
    WMA.ResizeInit(6, 6, 0.);
    WMB.ResizeInit(6, 6, 0.);
      
    /* Setta gli indici della matrice - le incognite sono ordinate come:
     *   - posizione (3)
     *   - parametri di rotazione (3)
     *   - quantita' di moto (3)
     *   - momento della quantita' di moto 
     * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
     * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
     * e' dato da iGetFirstPositionIndex()+i */
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= 6; iCnt++) {
        WMA.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WMA.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      
        WMB.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WMB.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
    }   
   
    AssMat_(WMA, WMB, 1., XCurr, XPrimeCurr);
}


void
Body::AssMat_(FullSubMatrixHandler& WMA,
	      FullSubMatrixHandler& WMB,
	      doublereal dCoef,
	      const VectorHandler& /* XCurr */ ,
	      const VectorHandler& /* XPrimeCurr */ )
{
    DEBUGCOUTFNAME("Body::AssMat_");
   
    Vec3 V(pNode->GetVCurr());
    Vec3 W(pNode->GetWRef());
   
#if 0
    /* ci pensa AfterPredict() */
    S = pNode->GetRRef()*S0;
    J = pNode->GetRRef()*(J0*(pNode->GetRRef()).Transpose());
#endif

    Mat3x3 SWedge(S);			/* S /\ */
    Vec3 Sc(S*dCoef);

    /* 
     * momentum: 
     *
     * m * I DeltaV - S /\ DeltagP + ( S /\ W ) /\ Deltag 
     */
    WMB.IncCoef(1, 1, dMass);
    WMB.IncCoef(2, 2, dMass);
    WMB.IncCoef(3, 3, dMass);
      
    WMB.Sub(1, 4, SWedge);
    WMA.Add(1, 4, Mat3x3(Sc.Cross(W)));
   
    /* 
     * momenta moment: 
     *
     * S /\ DeltaV + J DeltagP + ( V /\ S /\ - ( J * W ) /\ ) Deltag
     */
    WMB.Add(4, 1, SWedge); 

    WMB.Add(4, 4, J);
    WMA.Add(4, 4, Mat3x3(V, Sc)-Mat3x3(J*(W*dCoef)));
}


SubVectorHandler& 
Body::AssRes(SubVectorHandler& WorkVec,
             doublereal /* dCoef */ ,
	     const VectorHandler& /* XCurr */ ,
	     const VectorHandler& /* XPrimeCurr */ )
{
    DEBUGCOUTFNAME("Body::AssRes");

    /* Se e' definita l'accelerazione di gravita',
     * la aggiunge (solo al residuo) */
    Vec3 GravityAcceleration(0.);
    flag g = GravityOwner::fGetAcceleration(pNode->GetXCurr(), 
		    GravityAcceleration);

    integer iNumRows = 6;
    if (g) {
	    iNumRows = 12;
    }

    WorkVec.Resize(iNumRows);
    WorkVec.Reset(0.);

    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= iNumRows; iCnt++) {
        WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
    }
      
    Vec3 V(pNode->GetVCurr());
    Vec3 W(pNode->GetWCurr());
   
    /* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
    Mat3x3 R(pNode->GetRCurr());
    Vec3 STmp = R*S0;
    Mat3x3 JTmp = R*(J0*R.Transpose());

    /* Quantita' di moto: R[1] = Q - M * V - W /\ S */
    WorkVec.Sub(1, V*dMass + W.Cross(STmp));
   
    /* Momento della quantita' di moto: R[2] = G - S /\ V - J * W */
    WorkVec.Sub(4, JTmp*W + STmp.Cross(V));
  
    if (g) {
        WorkVec.Add(7, GravityAcceleration*dMass);
        WorkVec.Add(10, STmp.Cross(GravityAcceleration));
    }
 
    return WorkVec;
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
Body::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		    const VectorHandler& /* XCurr */ )
{   
    DEBUGCOUTFNAME("Body::InitialAssJac");

    /* Casting di WorkMat */
    FullSubMatrixHandler& WM = WorkMat.SetFull();
   
    /* Dimensiona e resetta la matrice di lavoro */
    integer iNumRows = 0;
    integer iNumCols = 0;
    InitialWorkSpaceDim(&iNumRows, &iNumCols);
    WM.ResizeInit(iNumRows, iNumCols, 0.);
   
    /* Setta gli indici della matrice - le incognite sono ordinate come:
     *   - posizione (3)
     *   - parametri di rotazione (3)
     *   - quantita' di moto (3)
     *   - momento della quantita' di moto 
     * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
     * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
     * e' dato da iGetFirstPositionIndex()+i */
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    integer iFirstVelocityIndex = iFirstPositionIndex+6;
    for (integer iCnt = 1; iCnt <= 6; iCnt++) {
        WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
        WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
    }   
   
    /* Prepara matrici e vettori */
   
    /* Velocita' angolare corrente */
    Vec3 W(pNode->GetWRef());
   
    Vec3 SWedgeW(S.Cross(W));   
    Mat3x3 WWedgeSWedge(-W, S);
    Mat3x3 WWedge(W);
    Mat3x3 WWedgeWWedgeSWedge(WWedge*WWedgeSWedge);
    Mat3x3 FDeltaW(Mat3x3(SWedgeW)+WWedgeSWedge);
   
    Vec3 JW(J*W);
    Mat3x3 JWWedge(JW);
    Mat3x3 MDeltag(WWedge*(J*WWedge-JWWedge));
    Mat3x3 MDeltaW(WWedge*J-JWWedge);
   
    /* Forza */
    WM.Add(1, 1, WWedgeWWedgeSWedge);
    WM.Add(1, 4, FDeltaW);
   
    /* Momento */
    WM.Add(4, 1, MDeltag);
    WM.Add(4, 4, MDeltaW);
     
    /* Derivata forza */
    WM.Add(7, 1, Mat3x3(W.Cross(SWedgeW))+WWedge*FDeltaW);
    WM.Add(7, 4, WWedge*WWedgeWWedgeSWedge);
   
    /* Derivata Momento */
    WM.Add(4, 1, WWedge*MDeltag);
    WM.Add(4, 4, WWedge*MDeltaW-Mat3x3(W.Cross(JW)));

    return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
Body::InitialAssRes(SubVectorHandler& WorkVec,
		    const VectorHandler& /* XCurr */ )
{   
    DEBUGCOUTFNAME("Body::InitialAssRes");

    integer iNumRows;
    integer iNumCols;
    InitialWorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.Resize(iNumRows);
    WorkVec.Reset(0.);
   
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= 12; iCnt++) {
        WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
    }
    
    Vec3 X(pNode->GetXCurr());
    Vec3 W(pNode->GetWCurr());
   
    /* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
    Mat3x3 R(pNode->GetRCurr());
    Vec3 STmp = R*S0;
    Mat3x3 JTmp = R*J0*R.Transpose();

    Vec3 FC(-W.Cross(W.Cross(STmp)));
    Vec3 MC(-W.Cross(JTmp*W));
   
    /* Forza */
    WorkVec.Add(1, FC);
   
    /* Momento */
    WorkVec.Add(4, MC);
   
    /* Derivata forza */
    WorkVec.Add(7, W.Cross(FC));
   
    /* Derivata momento */
    WorkVec.Add(10, W.Cross(MC));
  
    /* Se e' definita l'accelerazione di gravita',
     * la aggiunge (solo al residuo) */
    Vec3 GravityAcceleration(0.);
    if (GravityOwner::fGetAcceleration(X, GravityAcceleration)) {
        WorkVec.Add(1, GravityAcceleration*dMass);
        WorkVec.Add(4, STmp.Cross(GravityAcceleration));
        WorkVec.Add(10, (W.Cross(STmp)).Cross(GravityAcceleration));
    }   

    return WorkVec;
}


/* Usata per inizializzare la quantita' di moto */
void 
Body::SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const
{
    integer iFirstIndex = pNode->iGetFirstMomentumIndex();
   
    Vec3 V(pNode->GetVCurr());
    Vec3 W(pNode->GetWCurr());
    Mat3x3 R(pNode->GetRCurr());

    S = R*S0;  
    J = R*(J0*R.Transpose());
   
    X.Add(iFirstIndex+1, V*dMass+W.Cross(S));
    X.Add(iFirstIndex+4, S.Cross(V)+J*W);
}


void 
Body::AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
    Mat3x3 R = pNode->GetRRef();

    S = R*S0;
    J = R*(J0*R.Transpose());
}


/* Legge un corpo rigido */
Elem* ReadBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel)
{
    DEBUGCOUTFNAME("ReadBody");
   
    const char* sKeyWords[] = {
	    NULL
    };
   
    /* enum delle parole chiave */
    enum KeyWords {
        UNKNOWN = -1,
	LASTKEYWORD = 0
    };
   
    /* tabella delle parole chiave */
    KeyTable K(HP, sKeyWords);
   
    /* nodo collegato */
    StructNode* pNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
   
    if (pNode->GetStructNodeType() != StructNode::DYNAMIC) {
        std::cerr << "Illegal structural node type for body " << uLabel << std::endl;
        THROW(DataManager::ErrGeneric());
    }
      
    integer iNumMasses = 1;
    if (HP.IsKeyWord("condense")) {
        iNumMasses = HP.GetInt();
        if (iNumMasses < 1) {
	    std::cerr << "At least one mass is required" << std::endl;
	    THROW(DataManager::ErrGeneric());
        }
        DEBUGLCOUT(MYDEBUG_INPUT, 
	           iNumMasses << " masses will be condensed" << std::endl);
      
        /* The inertia is calculated as follows:               
         *                                                     
         * dm = Sum(dm_i)                                      
         *                                                     
         * Xgc = Sum(Xgc_i*dm_i)/Sum(dm_i)                     
         *                                                     
         * J = Sum(J_i)-Sum(dm_i*(Xgc_i-Xgc)/\(Xgc_i-Xgc)/\)   
         *                                                     
         * and it can be accomplished by accumulating:         
         *                                                     
         * dm = Sum(dm_i)                                      
         *                                                     
         * ~S = Sum(Xgc_i*dm_i)                                
         *                                                     
         * ~J = Sum(J_i)-Sum(dm_i*Xgc_i/\*Xgc_i/\)             
         *                                                     
         * then calculating                                    
         *                                                     
         * Xgc = S/dm                                          
         *                                                     
         * and finally:                                        
         *                                                     
         * J = ~J-Xgc/\(dm*Xgc/\-2*~S)                         
         *                                                     
         */			                          
    }
   
    doublereal dm = 0.;
    ReferenceFrame RF(pNode);
    Vec3 Xgc(0.);
    Vec3 STmp(0.);
    Mat3x3 J(0.);
   
    for (int iCnt = 1; iCnt <= iNumMasses; iCnt++) {
        /* massa */
        doublereal dmTmp = HP.GetReal();
   
        DEBUGLCOUT(MYDEBUG_INPUT, "Mass(" << iCnt << ") = " << dmTmp << std::endl);
        dm += dmTmp;

        /* posiz. c.g. */
        Vec3 XgcTmp(HP.GetPosRel(RF));
        STmp += XgcTmp*dmTmp;
   
        DEBUGLCOUT(MYDEBUG_INPUT, "position of mass(" << iCnt 
		   << ") center of gravity = " << XgcTmp << std::endl);
      
        /*
	 * matrice del mom. d'inerzia
	 *
         * Usa la funzione che legge una matrice qualsiasi con parole chiave
         * per forme abbreviate:
         *   - null: matrice vuota
         *   - eye:  matrice identita'
         *   - diag: matrice diagonale, seguita da 3 reali
         *   - sym:  matrice simmetrica, seguita da 6 reali, 
         *           letta come triangolare superiore ordinata per righe: 
         *           m11, m12, m13,    m22, m23,    m33
         *   - matrice generica, costituita da 9 reali, letta per righe:
         *           m11, m12, m13,    m21, m22, m23,   m31, m32, m33 
         *
         * Si assume inoltre che la matrice dei momenti di inerzia
	 * sia espressa nel centro di massa del corpo, quindi viene
	 * corretta per l'eventuale offset rispetto al nodo
         */
        Mat3x3 JTmp(HP.GetMatRel(RF));
        DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt 
		   << ") =" << std::endl << JTmp << std::endl);
          
        if (HP.IsKeyWord("inertial")) {
	    DEBUGLCOUT(MYDEBUG_INPUT, 
	               "supplied in inertial reference frame" << std::endl);	 
	    if (HP.IsKeyWord("node")) {
	        NO_OP;
	    } else {	    	   
	        Mat3x3 RTmp(HP.GetRotRel(RF));
	        JTmp = RTmp*(JTmp*RTmp.Transpose());
	    }
	    DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt 
		       << ") in current frame =" << std::endl << JTmp << std::endl);
        }      

        J += (JTmp-Mat3x3(XgcTmp, XgcTmp*dmTmp));      
    }
    
    Xgc = STmp/dm;

    DEBUGLCOUT(MYDEBUG_INPUT, "Total mass: " << dm << std::endl
	       << "Center of mass: " << Xgc << std::endl
	       << "Inertia matrix:" << std::endl << J << std::endl);
	           
    flag fOut = pDM->fReadOutput(HP, Elem::BODY);
      
    /* Allocazione e costruzione */
    Elem* pEl = NULL;
    SAFENEWWITHCONSTRUCTOR(pEl, Body, Body(uLabel, pNode, dm, Xgc, J, fOut));

    /* Se non c'e' il punto e virgola finale */
    if (HP.IsArg()) {
        std::cerr << std::endl
	    << "semicolon expected at line " << HP.GetLineData() << std::endl;
        THROW(DataManager::ErrGeneric());
    }   
   
    return pEl;   
} /* End of ReadBody() */

