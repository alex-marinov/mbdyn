/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifndef ACCELEROMETER_H
#define ACCELEROMETER_H

#if defined(USE_ELECTRIC_NODES) && defined(USE_STRUCT_NODES)

#include "elec.h"
      
/* Accelerometer - begin */

class Accelerometer : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   doublereal dOmega;
   doublereal dTau;
   doublereal dCsi;
   doublereal dKappa;
   
 public:
   Accelerometer(unsigned int uL, const DofOwner* pD, 
		 const StructNode* pS, const AbstractNode* pA,
		 const Vec3& TmpDir,
		 doublereal dO, doublereal dT, doublereal dC, doublereal dK,
		 flag fOut);
   ~Accelerometer(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectricType(void) const {
      return Electric::ACCELEROMETER;
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0);

 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = pStrNode;
     connectedNodes[1] = pAbsNode;
   };
   /* ************************************************ */
};


inline void* Accelerometer::pGet(void) const
{
   return (void*)this;
}

/* Accelerometer - end */


/* TraslAccel - begin */

class TraslAccel : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   Vec3 f;
   
 public:
   TraslAccel(unsigned int uL, const DofOwner* pD, 
	      const StructNode* pS, const AbstractNode* pA,
	      const Vec3& TmpDir, const Vec3& Tmpf,
	      flag fOut);
   ~TraslAccel(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectricType(void) const {
      return Electric::ACCELEROMETER;
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = pStrNode;
     connectedNodes[1] = pAbsNode;
   };
   /* ************************************************ */
};


inline void* TraslAccel::pGet(void) const 
{
   return (void*)this;
}
   
/* TraslAccel - end */


/* RotAccel - begin */

class RotAccel : virtual public Elem, public Electric {
 private:
   const StructNode* pStrNode;
   const AbstractNode* pAbsNode;
   Vec3 Dir;
   
 public:
   RotAccel(unsigned int uL, const DofOwner* pD, 
	      const StructNode* pS, const AbstractNode* pA,
	      const Vec3& TmpDir,
	      flag fOut);
   ~RotAccel(void);
   virtual inline void* pGet(void) const;

   virtual Electric::Type GetElectricType(void) const {
      return Electric::ACCELEROMETER;
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int i) const;

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
   virtual void SetInitialValue(VectorHandler& /* X */ ) const;
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
     connectedNodes.resize(2);
     connectedNodes[0] = pStrNode;
     connectedNodes[1] = pAbsNode;
   };
   /* ************************************************ */
};


inline void* RotAccel::pGet(void) const 
{
   return (void*)this;
}
   
/* RotAccel - end */

#endif /* USE_ELECTRIC_NODES && USE_STRUCT_NODES */

#endif /* ACCELEROMETER_H */

