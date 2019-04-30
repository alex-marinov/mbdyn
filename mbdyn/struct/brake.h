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


#ifndef BRAKEJ_H
#define BRAKEJ_H

#include "joint.h"
#include "drive.h"
#include "friction.h"

/* Brake - begin */

class Brake : virtual public Elem, public Joint {
 private:
	/*
	 * Brake
	 *
	 * rotation axis about local axis 3
	 */

   /* Freno - asse di rotazione dato dall'asse 3 del sistema di 
    * riferimento della cerniera. Tale sistema e' noto relativamente ai due
    * nodi. In particolare rispetto al nodo 1 la trasformazione dal sistema
    * di riferimento della cerniera al sistema globale e': R1*R1h, mentre per
    * il nodo 2 la medesima trasformazion e': R2*R2h.
    * I vettori d1 e d2 esprimono, nel sistema di riferimento dei rispettivi 
    * nodi, la distanza della cerniera dai nodi stessi. 
    * I vettori F, M esprimono le reazioni vincolari di forza e coppia. */
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 d1;
   Mat3x3 R1h;
   Vec3 d2;
   Mat3x3 R2h;
   //Vec3 F;
   Vec3 M;

   /* if the brake generates a force, Dir is the direction
    * with respect to node 1 (supposed to be the fixed one) */
#if 0
   bool isForce;
   Vec3 Dir;
#endif

   mutable doublereal dTheta;

   /* friction related data */
   BasicShapeCoefficient *const Sh_c;
   BasicFriction *const fc;
   const doublereal preF;
   const doublereal r;
   DriveOwner brakeForce;
   static const unsigned int NumSelfDof;
   static const unsigned int NumDof;
   /* end of friction related data */

#ifdef USE_NETCDF
   MBDynNcVar Var_E;
   MBDynNcVar Var_Omega;
   MBDynNcVar Var_fc;
   MBDynNcVar Var_Fb;
#endif // USE_NETCDF

 public:
   /* Costruttore non banale */
   Brake(unsigned int uL, const DofOwner* pDO,
		   const StructNode* pN1, const StructNode* pN2,
		   const Vec3& dTmp1, const Vec3& dTmp2,
		   const Mat3x3& R1hTmp, const Mat3x3& R2hTmp, flag fOut,
		   const doublereal rr,
		   const doublereal pref,
		   BasicShapeCoefficient *const sh,
		   BasicFriction *const f,
#if 0
		   bool isforce,
		   const Vec3& dir,
#endif
		   DriveCaller *pdc);
   
   /* Distruttore */
   ~Brake(void);

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::BRAKE;
   };
   
   virtual unsigned int iGetNumDof(void) const;
   
   DofOrder::Order GetDofType(unsigned int i) const;

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

   virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = NumDof;
      *piNumCols = NumDof;
      if (fc) {
          *piNumRows += fc->iGetNumDof();
          *piNumCols += fc->iGetNumDof();
      } 
   };
   
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
			    
   DofOrder::Order GetEqType(unsigned int i) const;
  
   void OutputPrepare(OutputHandler &OH);
   void Output(OutputHandler& OH) const;
 

   /* funzioni usate nell'assemblaggio iniziale */
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 10;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const {
      *piNumRows = 34; 
      *piNumCols = 34;
   };
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);
   
   /* Dati privati */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual unsigned int iGetPrivDataIdx(const char *s) const;
   virtual doublereal dGetPrivData(unsigned int i) const;
   
   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
   };
   /* ************************************************ */
};

/* Brake - end */



#endif /* PLANEJ_H */

