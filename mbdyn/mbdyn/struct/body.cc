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
: Elem(uL, ElemType::BODY, fOut), 
ElemGravityOwner(uL, ElemType::BODY, fOut), 
InitialAssemblyElem(uL, ElemType::BODY, fOut),
pNode(pNodeTmp), 
dMass(dMassTmp), 
Xgc(XgcTmp), 
S0(XgcTmp*dMassTmp), 
J0(JTmp)
{
    ASSERT(pNode != NULL);
    ASSERT(pNode->GetNodeType() == NodeType::STRUCTURAL);
    ASSERT(dMassTmp > 0.);
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& 
Body::Restart(ostream& out) const
{
    out << "  body: " << GetLabel() << ", " 
        << pNode->GetLabel() << ", " << dMass << ", reference, node, ",
        Xgc.Write(out, ", ") << ", ",
        J0.Write(out, ", ") << ';' << endl;
   
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
    integer iNumRows = 0;
    integer iNumCols = 0;
    this->WorkSpaceDim(&iNumRows, &iNumCols);
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
    integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
    for (integer iCnt = 1; iCnt <= 6; iCnt++) {
        WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WM.fPutRowIndex(6+iCnt, iFirstMomentumIndex+iCnt);
        WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);   
    }   
      
    AssMat_(WM, WM, dCoef, XCurr, XPrimeCurr);
    return WorkMat;
}


void
Body::AssEig(VariableSubMatrixHandler& WorkMatA,
	     VariableSubMatrixHandler& WorkMatB,	       
	     const VectorHandler& XCurr,
	     const VectorHandler& XPrimeCurr)
{
    DEBUGCOUTFNAME("Body::AssEig");

    /* Casting di WorkMat */
    FullSubMatrixHandler& WMA = WorkMatA.SetFull();
    FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   
    /* Dimensiona e resetta la matrice di lavoro */
    integer iNumRows = 0;
    integer iNumCols = 0;
    this->WorkSpaceDim(&iNumRows, &iNumCols);
    WMA.ResizeInit(iNumRows, iNumCols, 0.);
    WMB.ResizeInit(iNumRows, iNumCols, 0.);
      
    /* Setta gli indici della matrice - le incognite sono ordinate come:
     *   - posizione (3)
     *   - parametri di rotazione (3)
     *   - quantita' di moto (3)
     *   - momento della quantita' di moto 
     * e gli indici sono consecutivi. La funzione pGetFirstPositionIndex() 
     * ritorna il valore del primo indice -1, in modo che l'indice i-esimo
     * e' dato da iGetFirstPositionIndex()+i */
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
    for (integer iCnt = 1; iCnt <= 6; iCnt++) {
        WMA.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WMA.fPutRowIndex(6+iCnt, iFirstMomentumIndex+iCnt);
        WMA.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
      
        WMB.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WMB.fPutRowIndex(6+iCnt, iFirstMomentumIndex+iCnt);
        WMB.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
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
   
    S = pNode->GetRRef()*S0;
    J = pNode->GetRRef()*(J0*(pNode->GetRRef()).Transpose());
               
    /* matrice di massa: J[1,1] = M */  
    WMB.fPutCoef(1, 1, dMass);        
    WMB.fPutCoef(2, 2, dMass);        
    WMB.fPutCoef(3, 3, dMass);        
      
    /* matrice momento statico: J[1,2] = (S/\W)/\Deltag-S/\DeltagP,
     *                          J[2,1] = S/\DeltaV */
    Mat3x3 SWedge(S);
    Mat3x3 VWedgeSWedge(V, S);
    Mat3x3 SWedgeWWedge(S.Cross(W));
    WMB.Add(4, 1, SWedge); 
    WMB.Add(1, 4, -SWedge);
    WMA.Add(1, 4, SWedgeWWedge*dCoef);
   
    /* matrice momenti d'inerzia: J[2,2] = J*DeltagP+(V/\S/\-(J*W)/\)Deltag */
    WMA.Add(4, 4, (VWedgeSWedge-J*W)*dCoef);
    WMB.Add(4, 4, J);
      
    /* J[4,1] = (S/\W)/\DeltaV */
    WMB.Add(10, 1, SWedgeWWedge);
   
    /* J[4,2] = V/\(2S/\W/\W-W/\S/\)Deltag-V/\S/\DeltagP */
    WMA.Add(10,4, Mat3x3(V*dCoef)*SWedgeWWedge);
    WMB.Add(10,4, -VWedgeSWedge);			  
}


SubVectorHandler& 
Body::AssRes(SubVectorHandler& WorkVec,
             doublereal /* dCoef */ ,
	     const VectorHandler& /* XCurr */ ,
	     const VectorHandler& /* XPrimeCurr */ )
{
    DEBUGCOUTFNAME("Body::AssRes");

    integer iNumRows;
    integer iNumCols;
    this->WorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.Resize(iNumRows);
    WorkVec.Reset(0.);
   
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= 12; iCnt++) {
        WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
    }
      
    Vec3 V(pNode->GetVCurr());
    Vec3 W(pNode->GetWCurr());
   
    /* Aggiorna i suoi dati (saranno pronti anche per AssJac) */
    Mat3x3 R(pNode->GetRCurr());
    Vec3 STmp = R*S0;
    Mat3x3 JTmp = R*(J0*R.Transpose());

    Vec3 SWedgeW(STmp.Cross(W));
   
    /* Quantita' di moto: R[1] = Q - M V + S/\W */
    WorkVec.Add(1, SWedgeW-V*dMass);
   
    /* Momento della quantita' di moto: R[2] = G - S/\V -J W */
    WorkVec.Add(4, V.Cross(STmp)-(JTmp*W));
   
    /* Derivata del momento della quantita' di moto: R[4] = - GP + V/\(S/\W) */
    WorkVec.Add(10, V.Cross(SWedgeW));
      
    /* Se e' definita l'accelerazione di gravita',
     * la aggiunge (solo al residuo) */
    Vec3 GravityAcceleration(0.);
    if (GravityOwner::fGetAcceleration(pNode->GetXCurr(), 
                                       GravityAcceleration)) {
        WorkVec.Add(7, GravityAcceleration*dMass);
        WorkVec.Add(10, STmp.Cross(GravityAcceleration));
    }   
  
    return WorkVec;
}


/* Output temporaneo, per controllo dei dati */
#ifdef DEBUG
void Body::Output(OutputHandler& OH) const
#else /* !DEBUG */
void Body::Output(OutputHandler& /* OH */ ) const
#endif /* !DEBUG */
{
    if(fToBeOutput()) {      
#ifdef DEBUG
        if (DEBUG_LEVEL_MATCH(DEBUG_LEVEL_OUTPUT)) {
            OH.Output() << "Body Element " << uLabel << ", linked to node " 
	        << pNode->GetLabel() << ':' << endl 
		<< "Mass: " << dMass << endl 
		<< "Mass Center Position: " << endl << Xgc << endl 
		<< "Static Moment Vector: " << endl << S << endl 
		<< "Inertia Moment Matrix: " << endl << J << endl;
	}
#endif
   }
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
    this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
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
        WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
        WM.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
        WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
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
    this->WorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.Resize(iNumRows);
    WorkVec.Reset(0.);
   
    integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
    for (integer iCnt = 1; iCnt <= 12; iCnt++) {
        WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
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
    (Vec3&)S = R*S0;  
    (Mat3x3&)J = R*(J0*R.Transpose());
   
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
   
    const char* sKeyWords[] = { NULL };
   
    /* enum delle parole chiave */
    enum KeyWords {
        UNKNOWN = -1,
	LASTKEYWORD = 0
    };
   
    /* tabella delle parole chiave */
    KeyTable K((int)LASTKEYWORD, sKeyWords);
   
    /* parser del blocco di controllo */
    HP.PutKeyTable(K);
   
    /* nodo collegato */
    unsigned int uNode = (unsigned int)HP.GetInt();
   
    DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Node " << uNode << endl);
   
    /* verifica di esistenza del nodo */
    StructNode* pNode;
    if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
        cerr << endl
	    << "line " << HP.GetLineData() 
	    << ": structural node " << uNode
	    << " not defined" << endl;
        THROW(DataManager::ErrGeneric());
    }
   
    if (pNode->GetStructNodeType() != StructNodeType::DYNAMIC) {
        cerr << "Illegal structural node type for body " << uLabel << endl;
        THROW(DataManager::ErrGeneric());
    }
      
    integer iNumMasses = 1;
    if (HP.IsKeyWord("condense")) {
        iNumMasses = HP.GetInt();
        if (iNumMasses < 1) {
	    cerr << "At least one mass is required" << endl;
	    THROW(DataManager::ErrGeneric());
        }
        DEBUGLCOUT(MYDEBUG_INPUT, 
	           iNumMasses << " masses will be condensed" << endl);
      
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
   
        DEBUGLCOUT(MYDEBUG_INPUT, "Mass(" << iCnt << ") = " << dmTmp << endl);
        dm += dmTmp;

        /* posiz. c.g. */
        Vec3 XgcTmp(HP.GetPosRel(RF));
        STmp += XgcTmp*dmTmp;
   
        DEBUGLCOUT(MYDEBUG_INPUT, "position of mass(" << iCnt 
		   << ") center of gravity = " << XgcTmp << endl);
      
        /* matrice del mom. d'inerzia */
        /* Usa la funzione che legge una matrice qualsiasi con parole chiave
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
        Mat3x3 JTmp(HP.GetMat3x3());
        DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt 
		   << ") =" << endl << JTmp << endl);
          
        if (HP.IsKeyWord("inertial")) {
	    DEBUGLCOUT(MYDEBUG_INPUT, 
	               "supplied in inertial reference frame" << endl);	 
	    if (HP.IsKeyWord("node")) {
	        NO_OP;
	    } else {	    	   
	        Mat3x3 RTmp(HP.GetRotRel(RF));
	        JTmp = RTmp*(JTmp*RTmp.Transpose());
	    }
	    DEBUGLCOUT(MYDEBUG_INPUT, "Inertia matrix of mass(" << iCnt 
		       << ") in current frame =" << endl << JTmp << endl);
        }      

        J += (JTmp-Mat3x3(XgcTmp, XgcTmp*dmTmp));      
    }
    
    Xgc = STmp/dm;

    DEBUGLCOUT(MYDEBUG_INPUT, "Total mass: " << dm << endl
	       << "Center of mass: " << Xgc << endl
	       << "Inertia matrix:" << endl << J << endl);
	           
    flag fOut = pDM->fReadOutput(HP, ElemType::BODY);
      
    /* Allocazione e costruzione */
    Elem* pEl = NULL;
    SAFENEWWITHCONSTRUCTOR(pEl, 
			   Body,
			   Body(uLabel, pNode, dm, Xgc, J, fOut), 
			   DMmm);

    /* Se non c'e' il punto e virgola finale */
    if (HP.fIsArg()) {
        cerr << endl
	    << "semicolon expected at line " << HP.GetLineData() << endl;
        THROW(DataManager::ErrGeneric());
    }   
   
    return pEl;   
} /* End of ReadBody() */

