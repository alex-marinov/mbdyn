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

/* Deformable hinges */


#ifndef GIMBAL_H
#define GIMBAL_H

#include "joint.h"


/* GimbalRotationJoint - begin */

class GimbalRotationJoint : 
virtual public Elem, public Joint {
private:      
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	const Mat3x3 R1h;
	const Mat3x3 R2h;

	Vec3 M;
	doublereal dTheta, dPhi;

	OrientationDescription od;

	void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec, doublereal dCoef);
#ifdef USE_NETCDF
	MBDynNcVar Var_Theta;
	MBDynNcVar Var_Phi;
#endif // USE_NETCDF

public:
		/* Costruttore non banale */
	GimbalRotationJoint(unsigned int uL,	       
		const DofOwner* pDO,
		const StructNode* pN1, 
		const StructNode* pN2,
		const Mat3x3& R1,
		const Mat3x3& R2, 
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	virtual ~GimbalRotationJoint(void);

	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const {
		return GIMBAL;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Output(OutputHandler& OH) const;

	virtual unsigned int iGetNumDof(void) const { 
		return 5;
	};

	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < iGetNumDof());
		return DofOrder::ALGEBRAIC;
	};

	virtual DofOrder::Order GetEqType(unsigned int i) const {
		ASSERT(i >= 0 && i < iGetNumDof());
		return DofOrder::ALGEBRAIC;
	};

	virtual void WorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 11;
		*piNumCols = 11; 
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

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);           
	/* funzioni usate nell'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 5;
	};

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 11; 
		*piNumCols = 11;
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
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i = 0) const;

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

/* GimbalRotationJoint - end */

#endif /* GIMBAL_H */

