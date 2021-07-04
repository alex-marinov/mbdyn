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

/* Giunti piani */

#ifndef PLANEJ_H
#define PLANEJ_H

#include "joint.h"
#include "drive.h"
#include "friction.h"
#include "output.h"

/* PlaneHingeJoint - begin */

class PlaneHingeJoint : virtual public Elem, public Joint {
 private:
   /* Cerniera piana - asse di rotazione dato dall'asse 3 del sistema di 
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
   Vec3 F;
   Vec3 M;
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
	MBDynNcVar Var_Omega;
	MBDynNcVar Var_MFR;
	MBDynNcVar Var_fc;
#endif // USE_NETCDF

   bool calcInitdTheta;
   mutable int NTheta;
   mutable doublereal dTheta, dThetaWrapped;

   /* friction related data */
   BasicShapeCoefficient *const Sh_c;
   BasicFriction *const fc;
   const doublereal preF;
   const doublereal r;
   doublereal M3;
   static const unsigned int NumSelfDof;
   static const unsigned int NumDof;
   /* end of friction related data */

 protected:
	OrientationDescription od;

 public:
   /* Costruttore non banale */
   PlaneHingeJoint(unsigned int uL, const DofOwner* pDO,
		   const StructNode* pN1, const StructNode* pN2,
		   const Vec3& dTmp1, const Vec3& dTmp2,
		   const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
	       const OrientationDescription& od,
	       flag fOut,
		   const bool _calcInitdTheta = true,
		   const doublereal initDTheta = 0.,
		   const doublereal rr = 0.,
		   const doublereal pref = 0.,
		   BasicShapeCoefficient *const sh = 0,
		   BasicFriction *const f = 0);
   
   /* Distruttore */
   ~PlaneHingeJoint(void);

   virtual void ReadInitialState(MBDynParser& HP);
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::PLANEHINGE;
   };
   
   virtual unsigned int iGetNumDof(void) const;

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
   
   DofOrder::Order GetDofType(unsigned int i) const;

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
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

   /* return s the dimension of the component */
   const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* PlaneHingeJoint - end */


/* PlaneRotationJoint - begin */

class PlaneRotationJoint : virtual public Elem, public Joint {
 private:
   /* Cerniera piana - asse di rotazione dato dall'asse 3 del sistema di 
    * riferimento della cerniera. Tale sistema e' noto relativamente ai due
    * nodi. In particolare rispetto al nodo 1 la trasformazione dal sistema
    * di riferimento della cerniera al sistema globale e': R1*R1h, mentre per
    * il nodo 2 la medesima trasformazion e': R2*R2h.
    * I vettori d1 e d2 esprimono, nel sistema di riferimento dei rispettivi 
    * nodi, la distanza della cerniera dai nodi stessi. 
    * I vettori F, M esprimono le reazioni vincolari di forza e coppia. */
   const StructNode* pNode1;
   const StructNode* pNode2;
   Mat3x3 R1h;
   Mat3x3 R2h;
   Vec3 M;
   mutable int NTheta;
   mutable doublereal dTheta, dThetaWrapped;
#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
	MBDynNcVar Var_Omega;
#endif // USE_NETCDF

 protected:
	OrientationDescription od;

 public:
   /* Costruttore non banale */
   PlaneRotationJoint(unsigned int uL, const DofOwner* pDO,
		   const StructNode* pN1, const StructNode* pN2,
		   const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
	       const OrientationDescription& od,
	       flag fOut);
   
   /* Distruttore */
   ~PlaneRotationJoint(void);

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::PLANEROTATION;
   };
   
   virtual unsigned int iGetNumDof(void) const { 
      return 2;
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
   
   DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 2);
      return DofOrder::ALGEBRAIC; 
   };

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

   void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 3+3+2;
      *piNumCols = 3+3+2; 
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
      return 2+2;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const {
      *piNumRows = 3+3+3+3+2+2; 
      *piNumCols = 3+3+3+3+2+2;
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

   /* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* PlaneRotationJoint - end */


/* AxialRotationJoint - begin */

class AxialRotationJoint : virtual public Elem, 
public Joint, public DriveOwner {
 private:
   /* Rotazione assiale attorno ad una cerniera piana - 
    * asse di rotazione dato dall'asse 3 del sistema di 
    * riferimento della cerniera. Tale sistema e' noto relativamente ai due
    * nodi. In particolare rispetto al nodo 1 la trasformazione dal sistema
    * di riferimento della cerniera al sistema globale e': R1*R1h, mentre per
    * il nodo 2 la medesima trasformazione e': R2*R2h.
    * I vettori d1 e d2 esprimono, nel sistema di riferimento dei rispettivi 
    * nodi, la distanza della cerniera dai nodi stessi. 
    * I vettori F, M esprimono le reazioni vincolari di forza e coppia. 
    * La velocita' di rotazione e' imposta attraverso un driver. */
   const StructNode* pNode1;
   const StructNode* pNode2;
   Vec3 d1;
   Mat3x3 R1h;
   Vec3 d2;
   Mat3x3 R2h;
   Vec3 F;
   Vec3 M;
   mutable int NTheta;
   mutable doublereal dTheta, dThetaWrapped;

#ifdef USE_NETCDF
	MBDynNcVar Var_Phi;
	MBDynNcVar Var_Omega;
	MBDynNcVar Var_MFR;
	MBDynNcVar Var_fc;
#endif // USE_NETCDF

   /* friction related data */
   BasicShapeCoefficient *const Sh_c;
   BasicFriction *const fc;
   const doublereal preF;
   const doublereal r;
   doublereal M3;
   static const unsigned int NumSelfDof;
   static const unsigned int NumDof;
   /* end of friction related data */

 protected:
	OrientationDescription od;

 public:
   /* Costruttore non banale */
   AxialRotationJoint(unsigned int uL, const DofOwner* pDO, 
		      const StructNode* pN1, const StructNode* pN2,
		      const Vec3& dTmp1, const Vec3& dTmp2,
		      const Mat3x3& R1hTmp, const Mat3x3& R2hTmp,
		      const DriveCaller* pDC,
		      const OrientationDescription& od,
		      flag fOut,
		      const doublereal rr = 0.,
		      const doublereal pref = 0.,
		      BasicShapeCoefficient *const sh = 0,
		      BasicFriction *const f = 0);
   
   /* Distruttore */
   ~AxialRotationJoint(void);
   
   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::AXIALROTATION; 
   };

   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
       unsigned int i = NumSelfDof;
       if (fc) {
           i+=fc->iGetNumDof();
       } 
       return i;
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
   
   DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < iGetNumDof());
      if (i<NumSelfDof) {
          return DofOrder::ALGEBRAIC; 
      } else {
          return fc->GetDofType(i-NumSelfDof);
      }
   };

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
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
      return 11;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows, 
				    integer* piNumCols) const {
      *piNumRows = 35; 
      *piNumCols = 35; 
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

   /* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* AxialRotationJoint - end */


/* PlanePinJoint - begin */

/* Incastro con liberta' di rotazione su un asse */

class PlanePinJoint : virtual public Elem, public Joint {
 private:
   const StructNode* pNode;
   Vec3 X0;
   Mat3x3 R0;
   Vec3 d;
   Mat3x3 Rh;
   Vec3 F;
   Vec3 M;
   bool calcInitdTheta;
   mutable int NTheta;
   mutable doublereal dTheta, dThetaWrapped;
   
 public:
   /* Costruttore non banale */
   PlanePinJoint(unsigned int uL, const DofOwner* pDO,
		 const StructNode* pN,
		 const Vec3& X0Tmp, const Mat3x3& R0Tmp, 
		 const Vec3& dTmp, const Mat3x3& RhTmp, 
		 flag fOut, const bool _calcInitdTheta,
		 const doublereal initDTheta);
   
   ~PlanePinJoint(void);

   /* Tipo di Joint */
   virtual Joint::Type GetJointType(void) const { 
      return Joint::PIN; 
   };
   
   /* legge lo sato iniziale*/
   virtual void ReadInitialState(MBDynParser& HP);
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;

   virtual unsigned int iGetNumDof(void) const { 
      return 5;
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
   
   virtual DofOrder::Order GetDofType(unsigned int i) const {
      ASSERT(i >= 0 && i < 5);
      return DofOrder::ALGEBRAIC;
   };

	DofOrder::Order GetEqType(unsigned int i) const;

   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
      *piNumRows = 11; 
      *piNumCols = 11;
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
   
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 10;
   };
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const { 
      *piNumRows = 22; 
      *piNumCols = 22; 
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
     connectedNodes.resize(1);
     connectedNodes[0] = pNode;
   };
   /* ************************************************ */

   /* return s the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
	virtual OutputHandler::Dimensions GetEquationDimension(integer index);
	virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* PlanePinJoint - end */

#endif /* PLANEJ_H */

