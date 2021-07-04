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

/* Vincoli generali */


#ifndef GENJ_H
#define GENJ_H

#include "joint.h"
#include "drive.h"
#include "output.h"

#ifndef MBDYN_X_DISTANCE_JOINT
/* DistanceJoint - begin */

class DistanceJoint : virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructDispNode* pNode1;
   const StructDispNode* pNode2;
   mutable Vec3 v;
   doublereal dAlpha;
#ifdef USE_NETCDF
   MBDynNcVar Var_V;
   MBDynNcVar Var_d;
#endif // USE_NETCDF

 public:
   /* Costruttore non banale */
   DistanceJoint(unsigned int uL, const DofOwner* pDO,
		 const StructDispNode* pN1, const StructDispNode* pN2,
		 const DriveCaller* pDC, flag fOut);

   ~DistanceJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const {
      return Joint::DISTANCE;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const {
      return 4;
   };
   virtual DofOrder::Order GetDofType(unsigned int i) const
   {
      ASSERT(i >= 0 && i < 4);
      return DofOrder::ALGEBRAIC;
   };

   virtual DofOrder::Order GetEqType(unsigned int i) const
   {
      return DofOrder::DIFFERENTIAL;
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 10; *piNumCols = 10; };

   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr);

   void OutputPrepare(OutputHandler& OH);
   virtual void Output(OutputHandler& OH) const;


   /* funzioni usate nell'assemblaggio iniziale */

   virtual unsigned int iGetInitialNumDof(void) const { return 8; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const
     { *piNumRows = 20; *piNumCols = 20; };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);

   /* Contributo al residuo durante l'assemblaggio iniziale */
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X);
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

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

/* DistanceJoint - end */


/* DistanceJointWithOffset - begin */

class DistanceJointWithOffset :
virtual public Elem, public Joint, public DriveOwner {
 private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 f1;
   Vec3 f2;
   mutable Vec3 v;
   doublereal dAlpha;
#ifdef USE_NETCDF
   MBDynNcVar Var_V;
   MBDynNcVar Var_d;
#endif // USE_NETCDF

 public:
   /* Costruttore non banale */
   DistanceJointWithOffset(unsigned int uL, const DofOwner* pDO,
			   const StructNode* pN1, const StructNode* pN2,
			   const Vec3& f1Tmp, const Vec3& f2Tmp,
			   const DriveCaller* pDC, flag fOut);

   ~DistanceJointWithOffset(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const {
      return Joint::DISTANCEWITHOFFSET;
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const {
      return 4;
   };

   virtual DofOrder::Order GetDofType(unsigned int i) const
   {
      ASSERT(i >= 0 && i < 4);
      return DofOrder::ALGEBRAIC;
   };

   virtual DofOrder::Order GetEqType(unsigned int i) const
   {
      return DofOrder::DIFFERENTIAL;
   };

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 16;
      *piNumCols = 16;
   };

   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr);

   void OutputPrepare(OutputHandler& OH);
   virtual void Output(OutputHandler& OH) const;


   /* funzioni usate nell'assemblaggio iniziale */

   virtual unsigned int iGetInitialNumDof(void) const { return 8; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const
     { *piNumRows = 32; *piNumCols = 32; };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);

   /* Contributo al residuo durante l'assemblaggio iniziale */
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   /* Setta il valore iniziale delle proprie variabili */
   virtual void SetInitialValue(VectorHandler& X);
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

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

/* DistanceJointWithOffset - end */
#endif


/* ClampJoint - begin */

class ClampJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode;    /* nodo incastrato */
   Vec3 XClamp;                /* posizione imposta */
   Mat3x3 RClamp;              /* assetto imposto */
   Vec3 F;                     /* forza di reazione */
   Vec3 M;                     /* momento di reazione */

 public:
   /* Costruttore definitivo (da mettere a punto) */
   ClampJoint(unsigned int uL, const DofOwner*pD, const StructNode* pN,
	      const Vec3& X0, const Mat3x3& R0, flag fOut);

   /* Distruttore + o - banale */
   virtual ~ClampJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const {
      return Joint::CLAMP;
   };

   /*Funzione che legge lo stato iniziale dal file di input*/
   void ReadInitialState(MBDynParser& HP);

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Funzioni obbligatorie, per la gestione dei dof */
   virtual unsigned int iGetNumDof(void) const {
      return 6;
   };
   virtual std::ostream& DescribeDof(std::ostream& out,
		   const char *prefix = "",
		   bool bInitial = false) const;
   virtual void DescribeDof(std::vector<std::string>& desc,
		   bool bInitial = false, int i = -1) const;
   virtual std::ostream& DescribeEq(std::ostream& out,
		   const char *prefix = "",
		   bool bInitial = false) const;
   virtual void DescribeEq(std::vector<std::string>& desc,
		   bool bInitial = false, int i = -1) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const
   {
      ASSERT(i >= 0 && i < 6);
      return DofOrder::ALGEBRAIC;
   };

   virtual DofOrder::Order GetEqType(unsigned int i) const {
      ASSERT(i >= 0 && i < iGetNumDof());
      return DofOrder::ALGEBRAIC;
   }

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
     { *piNumRows = 12; *piNumCols = 12; };

   /* Assemblaggio matrice jacobiana */
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);


   /* Inverse Dynamics: AssJac() */
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    const VectorHandler& XCurr);

   /* assemblaggio matrici per autovalori */
   void AssMats(VariableSubMatrixHandler& WorkMatA,
	       VariableSubMatrixHandler& WorkMatB,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);

   /* Assemblaggio residuo */
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr);

   /* inverse dynamics capable element */
   virtual bool bInverseDynamics(void) const;

   /* Inverse Dynamics: AssRes */
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    const VectorHandler& XCurr,
			    const VectorHandler& XPrimeCurr,
			    const VectorHandler& XPrimePrimeCurr,
			    InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

   /* Inverse Dynamics update */
   void Update(const VectorHandler& XCurr, InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

   void OutputPrepare(OutputHandler& OH);
   virtual void Output(OutputHandler& OH) const;


   /* funzioni usate nell'assemblaggio iniziale */

   virtual unsigned int iGetInitialNumDof(void) const { return 12; };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const
     { *piNumRows = 24; *piNumCols = 24; };

   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
					   const VectorHandler& XCurr);

   /* Contributo al residuo durante l'assemblaggio iniziale */
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr);

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);
   /* Metodi per l'estrazione di dati "privati".
    * Si suppone che l'estrattore li sappia interpretare.
    * Come default non ci sono dati privati estraibili */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual unsigned int iGetPrivDataIdx(const char *s) const;
   virtual doublereal dGetPrivData(unsigned int i) const;

   /* *******PER IL SOLUTORE PARALLELO******** */
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode;
   };
   /* ************************************************ */

   /* return s the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* ClampJoint - end */

#endif /* GENJ_H */
