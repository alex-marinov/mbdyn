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

/* Vincolo prismatico */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "prismj.h"

/* PrismaticJoint - begin */

/* Costruttore non banale */
PrismaticJoint::PrismaticJoint(unsigned int uL, const DofOwner* pDO,
			       const StructNode* pN1, const StructNode* pN2,
			       const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
			       flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
pNode1(pN1), pNode2(pN2),
R1h(R1hTmp), R2h(R2hTmp), M(Zero3)
{
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
}


/* Distruttore banale */
PrismaticJoint::~PrismaticJoint(void)
{
   NO_OP;
};


DofOrder::Order
PrismaticJoint::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(),
		"INDEX ERROR in PrismaticJoint::GetEqType");
	return DofOrder::ALGEBRAIC;
}


/* Contributo al file di restart */
std::ostream& PrismaticJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", prismatic, "
     << pNode1->GetLabel() 
     << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() 
     << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
PrismaticJoint::AssJac(VariableSubMatrixHandler& WorkMat,
		       doublereal dCoef,
		       const VectorHandler& /* XCurr */ ,
		       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PrismaticJoint::AssJac()" << std::endl);
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato */					
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Ridimensiona la sottomatrice in base alle esigenze */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);
   
   /* Recupera gli indici delle varie incognite */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Recupera i dati che servono */
   Mat3x3 R1hTmp(pNode1->GetRRef()*R1h);
   Mat3x3 R2hTmp(pNode2->GetRRef()*R2h);
   /* Suppongo che le reazioni M siano gia' state aggiornate da AssRes;
    * ricordo che la coppia M e' nel sistema locale */
   
      
   /* 
    * Il vincolo prismatico ha tre equazioni che affermano la coincidenza del
    * riferimento solidale con il vincolo visto dai due nodi.
    * 
    *      (R1 * R1h * e2)^T * (R2 * R2h * e3) = 0
    *      (R1 * R1h * e3)^T * (R2 * R2h * e1) = 0
    *      (R1 * R1h * e1)^T * (R2 * R2h * e2) = 0
    * 
    * A queste equazioni corrisponde una reazione di coppia agente attorno
    * agli assi dati dall'errore di allineamento:
    * 
    *       [ (R1 * R1h * e2) /\ (R2 * R2h *e3), ]
    *       [ (R1 * R1h * e3) /\ (R2 * R2h *e1), ]  M
    *       [ (R1 * R1h * e1) /\ (R2 * R2h *e2)  ]
    * 
    */

   /* Setta gli indici delle equazioni */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WM.PutRowIndex(0+iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(0+iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode2FirstPosIndex+iCnt);
	
      WM.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }
   
   Vec3 MTmp = M*dCoef; /* M e' stato aggiornato da AssRes */

   Vec3 e1a(R1hTmp.GetVec(1));
   Vec3 e2a(R1hTmp.GetVec(2));
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   Vec3 e3b(R2hTmp.GetVec(3));
   
   Mat3x3 MWedge(Mat3x3(MatCrossCross, e3b, e2a*MTmp(1))
		 +Mat3x3(MatCrossCross, e1b, e3a*MTmp(2))
		 +Mat3x3(MatCrossCross, e2b, e1a*MTmp(3)));
   Mat3x3 MWedgeT(MWedge.Transpose());
   
   WM.Add(1, 1, MWedge);
   WM.Sub(1, 4, MWedgeT);
 
   WM.Add(4, 1, MWedgeT); 
   WM.Sub(4, 4, MWedge);

   Vec3 v1(e2a.Cross(e3b));
   Vec3 v2(e3a.Cross(e1b));
   Vec3 v3(e1a.Cross(e2b));
   
   MWedge = Mat3x3(v1, v2, v3);
   
   WM.Add(1, 7, MWedge);
   WM.Sub(4, 7, MWedge);
      
   /* Modifica: divido le equazioni di vincolo per dCoef */
      
   /* Equazione di vincolo del momento
    * 
    * Attenzione: bisogna scrivere il vettore trasposto
    *   (Sb[1]^T*(Sa[3]/\))*dCoef
    * Questo pero' e' uguale a:
    *   (-Sa[3]/\*Sb[1])^T*dCoef,
    * che puo' essere ulteriormente semplificato:
    *   (Sa[3].Cross(Sb[1])*(-dCoef))^T;
    */

   MWedge = MWedge.Transpose();
   
   WM.Add(7, 1, MWedge);
   WM.Sub(7, 4, MWedge);
   
   return WorkMat;
}
   
   
/* Assemblaggio residuo */
SubVectorHandler& PrismaticJoint::AssRes(SubVectorHandler& WorkVec,
					 doublereal dCoef,
					 const VectorHandler& XCurr, 
					 const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering PrismaticJoint::AssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
   
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex()+3;
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(0+iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode2FirstMomIndex+iCnt);
      
      WorkVec.PutRowIndex(6+iCnt, iFirstReactionIndex+iCnt);
   }   
   
   /* Aggiorna i dati propri */
   M = Vec3(XCurr, iFirstReactionIndex+1);

   /* Costruisce i dati propri nella configurazione corrente */
   Mat3x3 R1hTmp(pNode1->GetRCurr()*R1h);
   Mat3x3 R2hTmp(pNode2->GetRCurr()*R2h);
   
   Vec3 e1a(R1hTmp.GetVec(1));
   Vec3 e2a(R1hTmp.GetVec(2));
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   Vec3 e3b(R2hTmp.GetVec(3));
   
   Vec3 MTmp(Mat3x3(e2a.Cross(e3b), e3a.Cross(e1b), e1a.Cross(e2b))*M);
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Sub(1, MTmp); 
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(4, MTmp);

   /* Modifica: divido le equazioni di vincolo per dCoef */
   ASSERT(dCoef != 0.);

   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(7, -(e3b.Dot(e2a)/dCoef));
   WorkVec.PutCoef(8, -(e1b.Dot(e3a)/dCoef));
   WorkVec.PutCoef(9, -(e2b.Dot(e1a)/dCoef));

   return WorkVec;
}



/* Output (da mettere a punto) */
void PrismaticJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Prismatic", OH, name);
		}
#endif // USE_NETCDF
	}
}


void PrismaticJoint::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) {
      Mat3x3 R1Tmp(pNode1->GetRCurr()*R1h);

#ifdef USE_NETCDF
		if (OH.UseBinary(OutputHandler::JOINTS)) {
			Joint::NetCDFOutput(OH, Zero3, M, Zero3, R1Tmp*M);
		}
#endif // USE_NETCDF

		if (OH.UseText(OutputHandler::JOINTS)) {
		  Joint::Output(OH.Joints(), "Prismatic", GetLabel(),
				Zero3, M, Zero3, R1Tmp*M) << std::endl;
		}
   }   
}

void
PrismaticJoint::SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
	if (ph) {
		for (unsigned i = 0; i < ph->size(); i++) {
			Joint::JointHint *pjh = dynamic_cast<Joint::JointHint *>((*ph)[i]);

			if (pjh == 0) {
				continue;
			}

			if (dynamic_cast<Joint::HingeHint<1> *>(pjh)) {
				R1h = pNode1->GetRCurr().Transpose()*pNode2->GetRCurr()*R2h;

			} else if (dynamic_cast<Joint::HingeHint<2> *>(pjh)) {
				R2h = pNode2->GetRCurr().Transpose()*pNode1->GetRCurr()*R1h;

			} else if (dynamic_cast<Joint::ReactionsHint *>(pjh)) {
				/* TODO */
			}
		}
	}
}

Hint *
PrismaticJoint::ParseHint(DataManager *pDM, const char *s) const
{
	if (strncasecmp(s, "hinge{" /*}*/, STRLENOF("hinge{" /*}*/)) == 0) {
		s += STRLENOF("hinge{" /*}*/);

		if (strcmp(&s[1], /*{*/ "}") != 0) {
			return 0;
		}

		switch (s[0]) {
		case '1':
			return new Joint::HingeHint<1>;

		case '2':
			return new Joint::HingeHint<2>;
		}
	}

	return 0;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PrismaticJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			       const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering PrismaticJoint::InitialAssJac()" << std::endl);
   
   /* Per ora usa la matrice piena; eventualmente si puo' 
    * passare a quella sparsa quando si ottimizza */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);
        
    
   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Nota: le reazioni vincolari sono: 
    * Momento,     3 incognite, riferimento locale
    */

   /* Setta gli indici dei nodi */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(0+iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutColIndex(0+iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutColIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WM.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
      WM.PutColIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   /* Setta gli indici delle reazioni */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);	
   }   
   
   /* Recupera i dati */
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());   
   /* M e' gia' stato aggiornato da InitialAssRes */
   Vec3 MPrime(XCurr, iReactionPrimeIndex+1);
   
   /* Distanze e matrici di rotazione dai nodi alla cerniera 
    * nel sistema globale */
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e1a(R1hTmp.GetVec(1));
   Vec3 e2a(R1hTmp.GetVec(2));
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   Vec3 e3b(R2hTmp.GetVec(3));
   
   /* */
   Mat3x3 MWedge(Mat3x3(MatCrossCross, e3b, e2a*M(1))
	 + Mat3x3(MatCrossCross, e1b, e3a*M(2))
	 + Mat3x3(MatCrossCross, e2b, e1a*M(3)));
   Mat3x3 MWedgeT(MWedge.Transpose());
   
   /* Equilibrio */
   WM.Add(1, 1, MWedge);
   WM.Sub(1, 7, MWedgeT);
   
   WM.Add(7, 1, MWedgeT);   
   WM.Sub(7, 7, MWedge);

   /* Derivate dell'equilibrio */
   WM.Add(4, 4, MWedge);
   WM.Sub(4, 10, MWedgeT);
   
   WM.Add(10, 4, MWedgeT);   
   WM.Sub(10, 10, MWedge);
   
   
   MWedge = 
     ( (Mat3x3(MatCrossCross, e3b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e3b*M(1))))
      + Mat3x3(MatCross, e3b*MPrime(1)) )*Mat3x3(MatCross, e2a)
     +( (Mat3x3(MatCrossCross, e1b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e1b*M(2))))
       +Mat3x3(MatCross, e1b*MPrime(2)) )*Mat3x3(MatCross, e3a)
     +( (Mat3x3(MatCrossCross, e2b, Omega1) + Mat3x3(MatCross, Omega2.Cross(e2b*M(3))))
       +Mat3x3(MatCross, e2b*MPrime(3)) )*Mat3x3(MatCross, e1a);

   WM.Add(4, 1, MWedge);
   WM.Sub(10, 1, MWedge);
     
   MWedge =
     ( (Mat3x3(MatCrossCross, e2a, Omega2) + Mat3x3(MatCross, Omega1.Cross(e2a*M(1))))
      + Mat3x3(MatCross, e2a*MPrime(1)) )*Mat3x3(MatCross, e3b)
     +( (Mat3x3(MatCrossCross, e3a, Omega2) + Mat3x3(MatCross, Omega1.Cross(e3a*M(2))))
       + Mat3x3(MatCross, e3a*MPrime(2)) )*Mat3x3(MatCross, e1b)
     +( (Mat3x3(MatCrossCross, e1a, Omega2) + Mat3x3(MatCross, Omega1.Cross(e1a*M(3))))
       + Mat3x3(MatCross, e1a*MPrime(3)) )*Mat3x3(MatCross, e2b);
 
   WM.Sub(4, 7, MWedge);
   WM.Add(10, 7, MWedge);

   Vec3 v1(e2a.Cross(e3b)); 
   Vec3 v2(e3a.Cross(e1b)); 
   Vec3 v3(e1a.Cross(e2b));

   /* Error handling: il programma si ferma, pero' segnala dov'e' l'errore */
   if (v1.Dot() < std::numeric_limits<doublereal>::epsilon()
	|| v2.Dot() < std::numeric_limits<doublereal>::epsilon()
	|| v3.Dot() < std::numeric_limits<doublereal>::epsilon())
   {
      silent_cerr("PrismaticJoint(" << GetLabel() << "): "
	      "first and second node hinge axes are (nearly) orthogonal"
	      << std::endl);
      throw ErrNullNorm(MBDYN_EXCEPT_ARGS);
   }
   
   MWedge = Mat3x3(v1, v2, v3);
   
   /* Equilibrio */
   WM.Add(1, 13, MWedge);
   WM.Sub(7, 13, MWedge);

   /* Derivate dell'equilibrio */
   WM.Add(4, 16, MWedge);
   WM.Sub(10, 16, MWedge);
   

   MWedge = MWedge.Transpose();
   
   /* Equaz. di vincolo */
   WM.Add(13, 1, MWedge);
   WM.Sub(13, 7, MWedge);
      
   /* Devivate delle equaz. di vincolo */
   WM.Add(16, 4, MWedge);
   WM.Sub(16, 10, MWedge);
      
   v1 = e3b.Cross(e2a.Cross(Omega1)) + e2a.Cross(Omega2.Cross(e3b));
   v2 = e1b.Cross(e3a.Cross(Omega1)) + e3a.Cross(Omega2.Cross(e1b));
   v3 = e2b.Cross(e1a.Cross(Omega1)) + e1a.Cross(Omega2.Cross(e2b));
   
   MWedge = Mat3x3(v1, v2, v3);
   
   /* Derivate dell'equilibrio */
   WM.Add(4, 13, MWedge);
   WM.Sub(10, 13, MWedge);
   
   /* Devivate delle equaz. di vincolo */
   Omega1 = Omega1 - Omega2;
   
   v1 = e2a.Cross(e3b.Cross(Omega1));
   Vec3 v1p(e3b.Cross(Omega1.Cross(e2a)));
   v2 = e3a.Cross(e1b.Cross(Omega1));
   Vec3 v2p(e1b.Cross(Omega1.Cross(e3a)));
   v3 = e1a.Cross(e2b.Cross(Omega1));
   Vec3 v3p(e2b.Cross(Omega1.Cross(e1a)));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = v1(iCnt);
      WM.PutCoef(16, 0+iCnt, d);
      d = v1p(iCnt);
      WM.PutCoef(16, 6+iCnt, d);

      d = v2(iCnt);
      WM.PutCoef(17, 0+iCnt, d);
      d = v2p(iCnt);
      WM.PutCoef(17, 6+iCnt, d);
      
      d = v3(iCnt);
      WM.PutCoef(18, 0+iCnt, d);
      d = v3p(iCnt);
      WM.PutCoef(18, 6+iCnt, d);
   }    
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
PrismaticJoint::InitialAssRes(SubVectorHandler& WorkVec,
			      const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering PrismaticJoint::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);

   /* Indici */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex()+3;
   integer iNode1FirstVelIndex = iNode1FirstPosIndex+6;
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex()+3;
   integer iNode2FirstVelIndex = iNode2FirstPosIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   integer iReactionPrimeIndex = iFirstReactionIndex+3;
   
   /* Setta gli indici */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(0+iCnt, iNode1FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iNode1FirstVelIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
      WorkVec.PutRowIndex(9+iCnt, iNode2FirstVelIndex+iCnt);
   }
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {      
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }   

   /* Recupera i dati */
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   /* Aggiorna M, che resta anche per InitialAssJac */
   M = Vec3(XCurr, iFirstReactionIndex+1);
   Vec3 MPrime(XCurr, iReactionPrimeIndex+1);   
   
   /* Vincolo nel sistema globale */
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);

   Vec3 e1a(R1hTmp.GetVec(1));
   Vec3 e2a(R1hTmp.GetVec(2));
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));  
   Vec3 e3b(R2hTmp.GetVec(3));  
   
   Vec3 MTmp(e2a.Cross(e3b*M(1))
	     +e3a.Cross(e1b*M(2))
	     +e1a.Cross(e2b*M(3)));
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, -MTmp);
   
   /* Equazioni di equilibrio, nodo 2 */
   WorkVec.Add(7, MTmp);
   
   MTmp = 
     (e2a.Cross(Omega2.Cross(e3b)) - e3b.Cross(Omega1.Cross(e2a)))*M(1)
     +(e3a.Cross(Omega2.Cross(e1b)) - e1b.Cross(Omega1.Cross(e3a)))*M(2)
     +(e1a.Cross(Omega2.Cross(e2b)) - e2b.Cross(Omega1.Cross(e1a)))*M(3)
     +e2a.Cross(e3b*MPrime(1))
     +e3a.Cross(e1b*MPrime(2))
     +e1a.Cross(e2b*MPrime(3));   

   /* Derivate delle equazioni di equilibrio, nodo 1 */
   WorkVec.Add(4, -MTmp);
   
   /* Derivate delle equazioni di equilibrio, nodo 2 */
   WorkVec.Add(10, MTmp);

   /* Equazioni di vincolo di rotazione */
   WorkVec.PutCoef(13, -(e3b.Dot(e2a)));
   WorkVec.PutCoef(14, -(e1b.Dot(e3a)));
   WorkVec.PutCoef(15, -(e2b.Dot(e1a)));
   
   /* Derivate delle equazioni di vincolo di rotazione */
   Omega2 = Omega2-Omega1;
   WorkVec.PutCoef(16, (e2a.Cross(e3b)).Dot(Omega2));
   WorkVec.PutCoef(17, (e3a.Cross(e1b)).Dot(Omega2));
   WorkVec.PutCoef(18, (e1a.Cross(e2b)).Dot(Omega2));
   
   return WorkVec;
}

const MBUnits::Dimensions
PrismaticJoint::GetEquationDimension(integer index) const {
	// DOF == 3
   MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = MBUnits::Dimensions::rad;
			break;
		case 2:
			dimension = MBUnits::Dimensions::rad;
			break;
      case 3:
			dimension = MBUnits::Dimensions::rad;
			break;
	}

	return dimension;
}

std::ostream&
PrismaticJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"relative orientation constraints" << std::endl;

	return out;
}
/* PrismaticJoint - end */
