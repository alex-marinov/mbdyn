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

/* Point to surface contact */


#ifndef POINT_CONTACT_H
#define POINT_CONTACT_H

#include "joint.h"


/* PointSurfaceContact - begin */

class PointSurfaceContact : 
virtual public Elem, public Joint{		
private:      
protected:
	/* Punto di contatto */
	const StructNode* pNode1;	

	/* Superificie */
	const StructNode* pSup;	

	/* Posizione e orientazione della superficie  */
	Vec3 SupDirection;

	/* Distanza tra punto e superficie */
	doublereal dDeltaL;
	
	/* Elastic "stiffness" */
	doublereal ElasticStiffness;

	/* Output --> forza applicata sul punto */
	Vec3 FNode1;
	
	Vec3 Farm;
	/* Output --> reazioni sulla superficie --> forza + momento di trasporto */
	Vec3 FSup;
	Vec3 MSup;
	
	Vec3 n;
	
	void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);	// costruzione Jacobiano locale 
	void AssVec(SubVectorHandler& WorkVec, doublereal dCoef);     // costruzione residuo locale

public:
	/* Costruttore non banale */
	PointSurfaceContact(unsigned int uL,	       
			const DofOwner* pDO,
			const StructNode* pN1, 
			const StructNode* pNs,
			const Vec3& SDir, const doublereal Ek,
			flag fOut);

	/* Distruttore */
	virtual ~PointSurfaceContact(void);
	
	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const {
		return POINT_SURFACE_CONTACT;
	};

	/* Contributo al file di Restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	virtual void Output(OutputHandler& OH) const;

	virtual unsigned int iGetNumDof(void) const {  
		return 0;
	};

/*	virtual DofOrder::Order GetDofType(unsigned int i) const {
		std::cout << "GetDofType 1" << std::endl;
		ASSERT(i >= 0 && i < iGetNumDof());
		std::cout << "GetDofType 2 " << std::endl;
		return DofOrder::ALGEBRAIC;
	};

	virtual DofOrder::Order GetEqType(unsigned int i) const {
		std::cout << "GetEqType 1" << std::endl;
		ASSERT(i >= 0 && i < iGetNumDof());
		std::cout << "GetEqType 2" << std::endl;
		return DofOrder::ALGEBRAIC;
	};*/

	virtual void WorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 9;	
		*piNumCols = 9; 	
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef, 
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr); 
   
	/* assemblaggio residuo */
	virtual SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};
	
	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 9; 
		*piNumCols = 9;
	};
	


	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
			const VectorHandler& XCurr);

   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr);

   
	/* Dati privati (aggiungere magari le reazioni vincolari) */
	/*virtual unsigned int iGetNumPrivData(void) const{};
	virtual unsigned int iGetPrivDataIdx(const char *s) const{};
	virtual doublereal dGetPrivData(unsigned int i = 0) const{};*/


	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
/*	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
*/	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* PointSurfaceContact - end */

#endif /* POINT_CONTACT_H */
