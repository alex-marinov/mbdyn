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

#ifndef PZBEAM2_H
#define PZBEAM2_H

#include "beam.h"
#include "beam2.h"
#include "matvec3n.h"

/* PiezoActuatorBeam2 - begin */

class PiezoActuatorBeam2 : public Beam2 {
   
 private:
   PiezoActuatorBeam2(const PiezoActuatorBeam2&);
   const PiezoActuatorBeam2& operator = (const PiezoActuatorBeam2&);
   
 protected:
   int iNumElec;
   const ScalarDifferentialNode **pvElecDofs;
   VecN V;
   Mat3xN PiezoMat[NUMDEFORM];
   
   /* Funzioni di calcolo delle matrici */
   virtual void AssStiffnessMat(FullSubMatrixHandler& WMA,
				FullSubMatrixHandler& WMB,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr);
   
   virtual void AssStiffnessVec(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr);
   
   virtual void AddInternalForces(Vec6& AzLoc);
   
 public:
   /* Costruttore normale */
   PiezoActuatorBeam2(unsigned int uL,
		     const StructNode* pN1, const StructNode* pN2,
		     const Vec3& F1, const Vec3& F2,
		     const Mat3x3& R1, const Mat3x3& R2,
		     const Mat3x3& r,
		     const ConstitutiveLaw6D* pd,
		     int iEl,
		     const ScalarDifferentialNode **pEDof,
		     const Mat3xN& T_Ie, const Mat3xN& T_Ik,
		     OrientationDescription ood,
		     flag fOut);
   
   /* Distruttore banale */
   virtual ~PiezoActuatorBeam2(void);
    
   /* Tipo di trave */
   virtual Beam::Type GetBeamType(void) const {
      return Beam::PIEZOELECTRICELASTIC; 
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
    
   /* Dimensioni del workspace; sono 36 righe perche' se genera anche le
    * forze d'inerzia consistenti deve avere accesso alle righe di definizione
    * della quantita' di moto */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    
   /* Settings iniziali, prima della prima soluzione */
   void SetValue(DataManager *pDM,
		   VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0);
   
      /* Prepara i parametri di riferimento dopo la predizione */
   virtual void AfterPredict(VectorHandler& /* X */ ,
			     VectorHandler& /* XP */ );
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
      
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */
   virtual SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);

 /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     Beam2::GetConnectedNodes(connectedNodes);
     int NumNodes = connectedNodes.size();
     connectedNodes.resize(NumNodes + iNumElec);
     for (int i = 0; i < iNumElec; i++) {
       connectedNodes[NumNodes + i] = pvElecDofs[i];
     }
   };
   /* ************************************************ */
};

/* PiezoActuatorBeam2 - end */


/* PiezoActuatorVEBeam2 - begin */

class PiezoActuatorVEBeam2 : public ViscoElasticBeam2 {
   
 private:
   PiezoActuatorVEBeam2(const PiezoActuatorVEBeam2&);
   const PiezoActuatorVEBeam2& operator = (const PiezoActuatorVEBeam2&);
   
 protected:
   int iNumElec;
   const ScalarDifferentialNode **pvElecDofs;
   VecN V;
   Mat3xN PiezoMat[NUMDEFORM];
   
   /* Funzioni di calcolo delle matrici */
   virtual void AssStiffnessMat(FullSubMatrixHandler& WMA,
				FullSubMatrixHandler& WMB,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr);
   
   virtual void AssStiffnessVec(SubVectorHandler& WorkVec,
				doublereal dCoef,
				const VectorHandler& XCurr,
				const VectorHandler& XPrimeCurr);
   
   virtual void AddInternalForces(Vec6& AzLoc);
   
 public:
   /* Costruttore normale */
   PiezoActuatorVEBeam2(unsigned int uL,
		       const StructNode* pN1, const StructNode* pN2,
		       const Vec3& F1, const Vec3& F2,
		       const Mat3x3& R1, const Mat3x3& R2,
		       const Mat3x3& r,
		       const ConstitutiveLaw6D* pd,
		       int iEl,
		       const ScalarDifferentialNode **pEDof,
		       const Mat3xN& T_Ie, const Mat3xN& T_Ik,
		       OrientationDescription ood,
		       flag fOut);
   
   /* Distruttore banale */
   virtual ~PiezoActuatorVEBeam2(void);
    
   /* Tipo di trave */
   virtual Beam::Type GetBeamType(void) const {
      return Beam::PIEZOELECTRICVISCOELASTIC; 
   };
   
   /* Contributo al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
    
   /* Dimensioni del workspace; sono 36 righe perche' se genera anche le
    * forze d'inerzia consistenti deve avere accesso alle righe di definizione
    * della quantita' di moto */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    
   /* Settings iniziali, prima della prima soluzione */
   void SetValue(DataManager *pDM,
		   VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0);
   
      /* Prepara i parametri di riferimento dopo la predizione */
   virtual void AfterPredict(VectorHandler& /* X */ ,
			     VectorHandler& /* XP */ );
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef,
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
      
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler&
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr);
   
   /* Contributo al residuo durante l'assemblaggio iniziale */
   virtual SubVectorHandler&
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     ViscoElasticBeam2::GetConnectedNodes(connectedNodes);
     int NumNodes = connectedNodes.size();
     connectedNodes.resize(NumNodes + iNumElec);
     for (int i = 0; i < iNumElec; i++) {
       connectedNodes[NumNodes + i] = pvElecDofs[i];
     }
   };
   /* ************************************************ */
};

/* PiezoActuatorVEBeam2 - end */

#endif /* PZBEAM2_H */

