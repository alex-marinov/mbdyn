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

/* Continuano i vincoli di rotazione piani */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <brake.h>

/* Brake - begin */

const unsigned int Brake::NumSelfDof(0);
const unsigned int Brake::NumDof(12);

/* Costruttore non banale */
Brake::Brake(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2,
		const Vec3& dTmp1, const Vec3& dTmp2,
		const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
		flag fOut, 
		const doublereal rr,
		const doublereal pref,
		BasicShapeCoefficient *const sh,
		BasicFriction *const f,
		DriveCaller *pdc)
: Elem(uL, Elem::JOINT, fOut), 
Joint(uL, Joint::PLANEHINGE, pDO, fOut), 
pNode1(pN1), pNode2(pN2),
d1(dTmp1), R1h(R1hTmp), d2(dTmp2), R2h(R2hTmp), F(0.), M(0.), dTheta(0.),
Sh_c(sh), fc(f), preF(pref), r(rr), brakeForce(pdc)
{
	NO_OP;
}


/* Distruttore banale */
Brake::~Brake(void)
{
   NO_OP;
};

void
Brake::SetValue(VectorHandler& X, VectorHandler& XP) const
{
	Mat3x3 RTmp((pNode1->GetRCurr()*R1h).Transpose()*(pNode2->GetRCurr()*R2h));
	Vec3 v(MatR2EulerAngles(RTmp));

	dTheta = v.dGet(3);

	fc->SetValue(X,XP,iGetFirstIndex()+NumSelfDof);
}

void
Brake::AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP)
{

	Mat3x3 RTmp(((pNode1->GetRCurr()*R1h).Transpose()
			*pNode1->GetRPrev()*R1h).Transpose()
			*((pNode2->GetRCurr()*R2h).Transpose()
			*pNode2->GetRPrev()*R2h));
	Vec3 V(MatR2EulerAngles(RTmp.Transpose()));

	dTheta += V.dGet(3);

	Mat3x3 R1(pNode1->GetRCurr());
	Mat3x3 R1hTmp(R1*R1h);
	Vec3 e3a(R1hTmp.GetVec(3));
	Vec3 Omega1(pNode1->GetWCurr());
	Vec3 Omega2(pNode2->GetWCurr());
	//relative velocity
	doublereal v = (Omega1-Omega2).Dot(e3a)*r;
	//reaction norm
	doublereal modF = std::max(brakeForce.dGet(), preF);
	fc->AfterConvergence(modF, v, X, XP, iGetFirstIndex() + NumSelfDof);
}


/* Contributo al file di restart */
std::ostream& Brake::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", plane hinge, "
     << pNode1->GetLabel() << ", reference, node, ",
     d1.Write(out, ", ")
     << ", hinge, reference, node, 1, ", (R1h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R1h.GetVec(2)).Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, node, ",
     d2.Write(out, ", ")
     << ", hinge, reference, node, 1, ", (R2h.GetVec(1)).Write(out, ", ")
     << ", 2, ", (R2h.GetVec(2)).Write(out, ", ") << ';' << std::endl;
   
   return out;
}


/* Assemblaggio jacobiano */
VariableSubMatrixHandler& 
Brake::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Brake::AssJac()" << std::endl);
   
   /* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
    * e' complicato */					
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Ridimensiona la sottomatrice in base alle esigenze */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeReset(iNumRows, iNumCols);
   
   /* Recupera gli indici delle varie incognite */
   integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstPosIndex = pNode2->iGetFirstPositionIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();

   /* Setta gli indici delle equazioni */
   for (unsigned int iCnt = 1; iCnt <= 6; iCnt++) {	
      WM.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WM.PutColIndex(iCnt, iNode1FirstPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
      WM.PutColIndex(6+iCnt, iNode2FirstPosIndex+iCnt);
   }
   
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {	
      WM.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
      WM.PutColIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Recupera i dati che servono */
   Mat3x3 R1(pNode1->GetRRef());
   Mat3x3 R2(pNode2->GetRRef());   
   Vec3 d1Tmp(R1*d1);
   Vec3 d2Tmp(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   /* Suppongo che le reazioni F, M siano gia' state aggiornate da AssRes;
    * ricordo che la forza F e' nel sistema globale, mentre la coppia M
    * e' nel sistema locale ed il terzo termine, M(3), e' nullo in quanto
    * diretto come l'asse attorno al quale la rotazione e' consentita */
   
      
   /* 
    * La cerniera piana ha le prime 3 equazioni uguali alla cerniera sferica;
    * inoltre ha due equazioni che affermano la coincidenza dell'asse 3 del
    * riferimento solidale con la cerniera visto dai due nodi.
    * 
    *      (R1 * R1h * e1)^T * (R2 * R2h * e3) = 0
    *      (R1 * R1h * e2)^T * (R2 * R2h * e3) = 0
    * 
    * A queste equazioni corrisponde una reazione di coppia agente attorno 
    * agli assi 1 e 2 del riferimento della cerniera. La coppia attorno 
    * all'asse 3 e' nulla per definizione. Quindi la coppia nel sistema 
    * globale e':
    *      -R1 * R1h * M       per il nodo 1,
    *       R2 * R2h * M       per il nodo 2
    * 
    * 
    *       xa   ga                   xb   gb                     F     M 
    * Qa |  0    0                     0    0                     I     0  | | xa |   | -F           |
    * Ga |  0    c*(F/\da/\-(Sa*M)/\)  0    0                     da/\  Sa | | ga |   | -da/\F-Sa*M |
    * Qb |  0    0                     0    0                    -I     0  | | xb | = |  F           |
    * Gb |  0    0                     0   -c*(F/\db/\-(Sb*M)/\) -db/\ -Sb | | gb |   |  db/\F+Sb*M |
    * F  | -c*I  c*da/\                c*I -c*db/\                0     0  | | F  |   |  xa+da-xb-db |
    * M1 |  0    c*Tb1^T*Ta3/\         0    c*Ta3^T*Tb1/\         0     0  | | M  |   |  Sb^T*Ta3    |
    * M2 |  0    c*Tb2^T*Ta3/\         0    c*Ta3^T*Tb2/\         0     0  | 
    * 
    * con da = R1*d01, db = R2*d02, c = dCoef,
    * Sa = R1*R1h*[e1,e2], Sb = R2*R2h*[e1, e2],
    * Ta3 = R1*R1h*e3, Tb1 = R2*R2h*e1, Tb2 = R2*R2h*e2
    */

   
   /* Moltiplica la forza ed il momento per il coefficiente
    * del metodo */
   Vec3 FTmp = F*dCoef;
   Vec3 MTmp = M*dCoef;

   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   MTmp = e2b*MTmp.dGet(1)-e1b*MTmp.dGet(2);
   
   Mat3x3 MWedgee3aWedge(MTmp, e3a);
   Mat3x3 e3aWedgeMWedge(e3a, MTmp);
   
   
   /* Contributo del momento alle equazioni di equilibrio dei nodi */
   Vec3 Tmp1(e2b.Cross(e3a));
   Vec3 Tmp2(e3a.Cross(e1b));
   
   
   /* Modifica: divido le equazioni di vincolo per dCoef */
   
   /* Equazioni di vincolo degli spostamenti */
   
   /* Equazione di vincolo del momento
    * 
    * Attenzione: bisogna scrivere il vettore trasposto
    *   (Sb[1]^T*(Sa[3]/\))*dCoef
    * Questo pero' e' uguale a:
    *   (-Sa[3]/\*Sb[1])^T*dCoef,
    * che puo' essere ulteriormente semplificato:
    *   (Sa[3].Cross(Sb[1])*(-dCoef))^T;
    */
   

      //retrive
          //friction coef
      doublereal f = fc->fc();
          //shape function
      doublereal shc = Sh_c->Sh_c();
          //omega and omega rif
      Vec3 Omega1(pNode1->GetWCurr());
      Vec3 Omega2(pNode2->GetWCurr());
      Vec3 Omega1r(pNode1->GetWRef());
      Vec3 Omega2r(pNode2->GetWRef());   
      //compute 
          //relative velocity
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
          //reaction norm
      doublereal modF = std::max(brakeForce.dGet(), preF);
          //reaction moment
      //doublereal M3 = shc*modF*f;
      
      ExpandableRowVector dfc;
      ExpandableRowVector dF;
      ExpandableRowVector dv;
          //variation of reaction force
      dF.ReDim(0);
/*
      if ((modF == 0.) or (F.Norm() > preF)) {
          dF.Set(0.,1,12+1);
          dF.Set(0.,2,12+2);
          dF.Set(0.,3,12+3);
      } else {
          dF.Set(F.dGet(1)/modF,1,12+1);
          dF.Set(F.dGet(2)/modF,2,12+2);
          dF.Set(F.dGet(3)/modF,3,12+3);
      }
*/
          //variation of relative velocity
      dv.ReDim(6);
      
/* old (wrong?) relative velocity linearization */

//       dv.Set((e3a.dGet(1)*1.-( e3a.dGet(2)*Omega1r.dGet(3)-e3a.dGet(3)*Omega1r.dGet(2))*dCoef)*r,1,0+4);
//       dv.Set((e3a.dGet(2)*1.-(-e3a.dGet(1)*Omega1r.dGet(3)+e3a.dGet(3)*Omega1r.dGet(1))*dCoef)*r,2,0+5);
//       dv.Set((e3a.dGet(3)*1.-( e3a.dGet(1)*Omega1r.dGet(2)-e3a.dGet(2)*Omega1r.dGet(1))*dCoef)*r,3,0+6);
//       
//       dv.Set(-(e3a.dGet(1)*1.-( e3a.dGet(2)*Omega2r.dGet(3)-e3a.dGet(3)*Omega2r.dGet(2))*dCoef)*r,4,6+4);
//       dv.Set(-(e3a.dGet(2)*1.-(-e3a.dGet(1)*Omega2r.dGet(3)+e3a.dGet(3)*Omega2r.dGet(1))*dCoef)*r,5,6+5);
//       dv.Set(-(e3a.dGet(3)*1.-( e3a.dGet(1)*Omega2r.dGet(2)-e3a.dGet(2)*Omega2r.dGet(1))*dCoef)*r,6,6+6);

/* new (exact?) relative velocity linearization */
// 
//       ExpandableRowVector domega11, domega12, domega13;
//       ExpandableRowVector domega21, domega22, domega23;
//       domega11.ReDim(3); domega12.ReDim(3); domega13.ReDim(3);
//       domega21.ReDim(3); domega22.ReDim(3); domega23.ReDim(3);
//       
//       domega11.Set(1., 1, 0+4);
//           domega11.Set( Omega1r.dGet(3)*dCoef, 2, 0+5);
//           domega11.Set(-Omega1r.dGet(2)*dCoef, 3, 0+6);
//       domega21.Set(1., 1, 6+4);
//           domega21.Set( Omega2r.dGet(3)*dCoef, 2, 6+5);
//           domega21.Set(-Omega2r.dGet(2)*dCoef, 3, 6+6);
//       domega12.Set(1., 1, 0+5);
//           domega12.Set(-Omega1r.dGet(3)*dCoef, 2, 0+4);
//           domega12.Set( Omega1r.dGet(1)*dCoef, 3, 0+6);
//       domega22.Set(1., 1, 6+5);
//           domega22.Set(-Omega2r.dGet(3)*dCoef, 2, 6+4);
//           domega22.Set( Omega2r.dGet(1)*dCoef, 3, 6+6);
//       domega13.Set(1., 1, 0+6);
//           domega13.Set( Omega1r.dGet(2)*dCoef, 2, 0+4);
//           domega13.Set(-Omega1r.dGet(1)*dCoef, 3, 0+5);
//       domega23.Set(1., 1, 6+6);
//           domega23.Set( Omega2r.dGet(2)*dCoef, 2, 6+4);
//           domega23.Set(-Omega2r.dGet(1)*dCoef, 3, 6+5);
// 
//       Vec3 domega = Omega1-Omega2;
//       dv.Set((e3a.dGet(1)*1.-( 
//       		e3a.dGet(2)*(Omega1.dGet(3)-Omega2.dGet(3))-
// 		e3a.dGet(3)*(domega.dGet(2)))*dCoef)*r,1);
// 		dv.Link(1, &domega11);
//       dv.Set((e3a.dGet(2)*1.-(
//       		-e3a.dGet(1)*(Omega1.dGet(3)-Omega2.dGet(3))+
// 		e3a.dGet(3)*(domega.dGet(1)))*dCoef)*r,2);
// 		dv.Link(2, &domega12);
//       dv.Set((e3a.dGet(3)*1.-( 
//       		e3a.dGet(1)*(Omega1.dGet(2)-Omega2.dGet(2))-
// 		e3a.dGet(2)*(domega.dGet(1)))*dCoef)*r,3);
// 		dv.Link(3, &domega13);
// 
//       dv.Set(-(e3a.dGet(1)*1.)*r,4,6+4); dv.Link(4, &domega21);
//       dv.Set(-(e3a.dGet(2)*1.)*r,5,6+5); dv.Link(5, &domega22);
//       dv.Set(-(e3a.dGet(3)*1.)*r,6,6+6); dv.Link(6, &domega23);


/* new (approximate: assume constant triads orientations) 
 * relative velocity linearization 
*/

      dv.Set((e3a.dGet(1)*1.)*r,1, 0+4);
      dv.Set((e3a.dGet(2)*1.)*r,2, 0+5);
      dv.Set((e3a.dGet(3)*1.)*r,3, 0+6);
      
      dv.Set(-(e3a.dGet(1)*1.)*r,4, 6+4);
      dv.Set(-(e3a.dGet(2)*1.)*r,5, 6+5);
      dv.Set(-(e3a.dGet(3)*1.)*r,6, 6+6);


      //assemble friction states
      fc->AssJac(WM,dfc,12+NumSelfDof,iFirstReactionIndex+NumSelfDof,dCoef,modF,v,
      		XCurr,XPrimeCurr,dF,dv);
      ExpandableRowVector dM3;
      ExpandableRowVector dShc;
      //compute 
          //variation of shape function
      Sh_c->dSh_c(dShc,f,modF,v,dfc,dF,dv);
          //variation of moment component
      dM3.ReDim(3);
      dM3.Set(shc*f,1); dM3.Link(1,&dF);
      dM3.Set(modF*f,2); dM3.Link(2,&dShc);
      dM3.Set(shc*modF,3); dM3.Link(3,&dfc);
      //assemble first node
          //variation of moment component
      dM3.Add(WM,0+4,e3a.dGet(1));
      dM3.Add(WM,0+5,e3a.dGet(2));
      dM3.Add(WM,0+6,e3a.dGet(3));
      //assemble second node
          //variation of moment component
      dM3.Sub(WM,6+4,e3a.dGet(1));
      dM3.Sub(WM,6+5,e3a.dGet(2));
      dM3.Sub(WM,6+6,e3a.dGet(3));
   
   return WorkMat;
}


/* Assemblaggio residuo */
SubVectorHandler& Brake::AssRes(SubVectorHandler& WorkVec,
					  doublereal dCoef,
					  const VectorHandler& XCurr, 
					  const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering Brake::AssRes()" << std::endl);
      
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
 
   /* Indici */
   integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomIndex = pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici dei nodi */
   for (int iCnt = 1; iCnt <= 6; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iNode1FirstMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iNode2FirstMomIndex+iCnt);
   }
   
   /* Indici del vincolo */
   for (unsigned int iCnt = 1; iCnt <= iGetNumDof(); iCnt++) {
      WorkVec.PutRowIndex(12+iCnt, iFirstReactionIndex+iCnt);
   }
   
   /* Aggiorna i dati propri */
   //FIXME
   F = Vec3(0.);
   M = Vec3(0.);

   /*
    * FIXME: provare a mettere "modificatori" di forza/momento sui gdl
    * residui: attrito, rigidezze e smorzamenti, ecc.
    */
   
   /* Recupera i dati */
   Vec3 x1(pNode1->GetXCurr());
   Vec3 x2(pNode2->GetXCurr());
   Mat3x3 R1(pNode1->GetRCurr());
   Mat3x3 R2(pNode2->GetRCurr());
   
   /* Costruisce i dati propri nella configurazione corrente */
   Vec3 dTmp1(R1*d1);
   Vec3 dTmp2(R2*d2);
   Mat3x3 R1hTmp(R1*R1h);
   Mat3x3 R2hTmp(R2*R2h);
   
   Vec3 e3a(R1hTmp.GetVec(3));
   Vec3 e1b(R2hTmp.GetVec(1));
   Vec3 e2b(R2hTmp.GetVec(2));
   
   Vec3 MTmp(e2b.Cross(e3a)*M.dGet(1)+e3a.Cross(e1b)*M.dGet(2));
   
   /* Equazioni di equilibrio, nodo 1 */
   
   /* Equazioni di equilibrio, nodo 2 */

   /* Modifica: divido le equazioni di vincolo per dCoef */

      bool ChangeJac(false);
      Vec3 Omega1(pNode1->GetWCurr());
      Vec3 Omega2(pNode2->GetWCurr());
      doublereal v = (Omega1-Omega2).Dot(e3a)*r;
      doublereal modF = std::max(brakeForce.dGet(), preF);
      try {
          fc->AssRes(WorkVec, 12 + NumSelfDof,
	  		iFirstReactionIndex + NumSelfDof,
			modF, v, XCurr, XPrimeCurr);
      }
      catch (Elem::ChangedEquationStructure) {
          ChangeJac = true;
      }
      doublereal f = fc->fc();
      doublereal shc = Sh_c->Sh_c(f, modF, v);
      M3 = shc*modF*f;
      WorkVec.Sub(4, e3a*M3);
      WorkVec.Add(10, e3a*M3);
//!!!!!!!!!!!!!!
//      M += e3a*M3;
      if (ChangeJac) {
          throw Elem::ChangedEquationStructure();
      }
   
   return WorkVec;
}

unsigned int Brake::iGetNumDof(void) const {
   unsigned int i = NumSelfDof;
   if (fc) {
       i+=fc->iGetNumDof();
   } 
   return i;
};


DofOrder::Order
Brake::GetDofType(unsigned int i) const {
   ASSERT(i >= 0 && i < iGetNumDof());
   if (i<NumSelfDof) {
       return DofOrder::ALGEBRAIC; 
   } else {
       return fc->GetDofType(i-NumSelfDof);
   }
};

DofOrder::Order
Brake::GetEqType(unsigned int i) const
{
	ASSERTMSGBREAK(i < iGetNumDof(), 
		"INDEX ERROR in Brake::GetEqType");
   if (i<NumSelfDof) {
       return DofOrder::ALGEBRAIC; 
   } else {
       return fc->GetEqType(i-NumSelfDof);
   }
}

/* Output (da mettere a punto) */
void Brake::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {
      Mat3x3 R2Tmp(pNode2->GetRCurr()*R2h);
      Mat3x3 RTmp((pNode1->GetRCurr()*R1h).Transpose()*R2Tmp);
      Mat3x3 R2TmpT(R2Tmp.Transpose());
      
      std::ostream &of = Joint::Output(OH.Joints(), "PlaneHinge", GetLabel(),
		    R2TmpT*F, M, F, R2Tmp*M)
	<< " " << MatR2EulerAngles(RTmp)*dRaDegr
	  << " " << R2TmpT*(pNode2->GetWCurr()-pNode1->GetWCurr());
      if (fc) {
          of << " " << M3 << " " << fc->fc() << " " << brakeForce.dGet();
      }
      of << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
Brake::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			       const VectorHandler& XCurr)
{
   /* Per ora usa la matrice piena; eventualmente si puo' 
    * passare a quella sparsa quando si ottimizza */
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WM.ResizeReset(iNumRows, iNumCols);
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
Brake::InitialAssRes(SubVectorHandler& WorkVec,
			       const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering Brake::InitialAssRes()" << std::endl);
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkVec.ResizeReset(iNumRows);

   return WorkVec;
}


unsigned int
Brake::iGetNumPrivData(void) const
{
	return 2;
}

unsigned int
Brake::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "rz") == 0) {
		return 1;
	}

	if (strcmp(s, "wz") == 0) {
		return 2;
	}


	return 0;
}

doublereal Brake::dGetPrivData(unsigned int i) const
{
   ASSERT(i >= 1 && i <= iGetNumPrivData());
   
   switch (i) {
    case 1: {
	return dTheta;
    }
      
    case 2: {
       Mat3x3 R2TmpT((pNode2->GetRCurr()*R2h).Transpose());
       Vec3 v(R2TmpT*(pNode2->GetWCurr()-pNode1->GetWCurr()));
       
       return v.dGet(3);
    }
      
    default:
      std::cerr << "Illegal private data" << std::endl;
      THROW(ErrGeneric());
   }
}

/* Brake - end */


