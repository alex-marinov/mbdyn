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

/* Vincolo di contatto con un piano */

#ifndef CONTACTJ_H
#define CONTACTJ_H

#include <limits>

#include "joint.h"

/* ContactJoint - begin */

class ContactJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 n;
   doublereal dD;
   doublereal dF;
   
 public:
   /* Costruttore non banale */
   ContactJoint(unsigned int uL, const DofOwner* pDO,
		const StructNode* pN1, const StructNode* pN2, 
		const Vec3& n, flag fOut) 
     : Elem(uL, fOut), 
     Joint(uL, pDO, fOut),
     pNode1(pN1), pNode2(pN2),
     n(n), dD(0.), dF(0.) {
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(n.Dot() > std::numeric_limits<doublereal>::epsilon());       

	Vec3 D(pNode2->GetXCurr()-pNode1->GetXCurr());
	Vec3 N(pNode1->GetRCurr()*n);
	dD = N.Dot(D);
     };
   
   ~ContactJoint(void) {
      NO_OP;
   };
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::INPLANECONTACT;
   };
   
   /* Contributo al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out << "ContactJoint: not implemented yet" << endl;
   };

   virtual unsigned int iGetNumDof(void) const { 
      return 1;
   };
      
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i == 0);
      return DofOrder::ALGEBRAIC;
   };
   
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 7;
      *piNumCols = 10;
   };
   
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr) {
      SparseSubMatrixHandler& WM = WorkMat.SetSparse();
      
      integer iIndex = iGetFirstIndex()+1;
      
      Vec3 D(pNode2->GetXCurr()-pNode1->GetXCurr());
      Vec3 V(pNode2->GetVCurr()-pNode1->GetVCurr());
      Vec3 N(pNode1->GetRCurr()*n);
      dD = N.Dot(D);
      doublereal dV = N.Dot(V);
      
      dF = XCurr(iIndex);
      
      if (dD > 0. 
	  || (dD == 0. && dF >= 0.) 
	  || (dD == 0. && dF == 0. && dV > 0.)) {
	 /* non attivo */
	 WM.Resize(1, 0);
	 WM.PutItem(1, iIndex, iIndex, 1.);
      } else {
	 /* attivo */
	 WM.Resize(27, 0);
	 
	 integer iNode1RowIndex = pNode1->iGetFirstMomentumIndex();
	 integer iNode1ColIndex = pNode1->iGetFirstPositionIndex();
	 integer iNode2RowIndex = pNode2->iGetFirstMomentumIndex();
	 integer iNode2ColIndex = pNode2->iGetFirstPositionIndex();
	       	 
	 Vec3 Tmp(N.Cross(D));	 
	 for (integer i = 1; i <= 3; i++) {
	    WM.PutItem(i, iIndex, iNode1ColIndex+3+i, Tmp.dGet(i));
	    doublereal d = N.dGet(i);
	    WM.PutItem(9+i, iNode1RowIndex+i, iIndex, -d);
	    WM.PutItem(12+i, iNode2RowIndex+i, iIndex, d);	   
	    WM.PutItem(3+i, iIndex, iNode1ColIndex+i, -d);	    
	    WM.PutItem(6+i, iIndex, iNode2ColIndex+i, d);
	 }
	 
	 Tmp = N*(dF*dCoef);
	 WM.PutCross(16, iNode1RowIndex+4, iNode1ColIndex+4, Tmp);
	 WM.PutCross(22, iNode2RowIndex+4, iNode1ColIndex+4, -Tmp);
      }
      return WorkMat;
   };
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr) {
      integer iIndex = iGetFirstIndex()+1;
      
      Vec3 D(pNode2->GetXCurr()-pNode1->GetXCurr());
      Vec3 V(pNode2->GetVCurr()-pNode1->GetVCurr());
      Vec3 N(pNode1->GetRCurr()*n);
      dD = N.Dot(D);
      doublereal dV = N.Dot(V);
      
      dF = XCurr(iIndex);
      
      if (dD > 0. 
	  || (dD == 0. && dF >= 0.) 
	  || (dD == 0. && dF == 0. && dV > 0.)) {
	 /* non attivo */
	 WorkVec.Resize(1);
	 WorkVec.PutItem(1, iIndex, -dF);
      } else {
	 /* attivo */
	 WorkVec.Resize(7);
	 
	 integer iNode1RowIndex = pNode1->iGetFirstMomentumIndex();       
	 integer iNode2RowIndex = pNode2->iGetFirstMomentumIndex();
	 
	 Vec3 Tmp(N*dF);
	 for (integer i = 1; i <= 3; i++) {
	    doublereal d = Tmp.dGet(i);
	    WorkVec.PutItem(i, iNode1RowIndex+i, d);
	    WorkVec.PutItem(3+i, iNode2RowIndex+i, -d);
	 }
	 WorkVec.PutItem(7, iIndex, -dD/dCoef);
      }
      
      return WorkVec;
   };
  
   void OutputPrepare(OutputHandler& OH) {
	   if (bToBeOutput()) {
#ifdef USE_NETCDF
		   if (OH.UseBinary(OutputHandler::JOINTS)) {
			   std::string name;
			   OutputPrepare_int("Contact", OH, name);
		   }
#endif // USE_NETCDF
	   }
   }

   virtual void Output(OutputHandler& OH) const {
	   if (bToBeOutput()) {
		   if (OH.UseText(OutputHandler::JOINTS)) {
			   Vec3 F(pNode1->GetRCurr()*(n*dF));
			   Joint::Output(OH.Joints(), "Contact", GetLabel(),
					   Vec3(dF, 0., 0.), Zero3, F, Zero3);
		   << " " << dD << endl;
		   }
#ifdef USE_NETCDF
		   if (OH.UseBinary(OutputHandler::JOINTS)) {
			   Joint::NetCDFOutput(OH, Vec3(dF, 0., 0., ), Zero3, F, Zero3);
		   }
#endif // USE_NETCDF
	   }
   };

   
   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 2;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 14;
      *piNumCols = 20; 
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr) {
      WorkMat.SetNullMatrix();
      return WorkMat;
   };
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr) {
      WorkVec.Resize(0);
      return WorkVec;
   };
   
   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X) {
      NO_OP;
   };
   
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0)
   {
      X.PutCoef(iGetFirstIndex() + 1, 0.);
   };

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const {
      // DOF == 1
      MBUnits::Dimensions dimension = MBUnits::Dimensions::UnknownDimension;

	   switch (index)
	   {
		   case 1:
			   dimension = MBUnits::Dimensions::Length;
			   break;
	   }

	return dimension;
   };

   /* describes the dimension of components of equation */
   virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const {

           int iIndex = iGetFirstIndex();

            out
               << prefix << iIndex + 1 << ": " 
               << "contact joint compenetration" << std::endl;
            
            return out;
   }
};

/* ContactJoint - end */

#endif /* CONTACTJ_H */

