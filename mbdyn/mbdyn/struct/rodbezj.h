/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/*
 * Author: Andrea Zanoni <andrea.zanoni@polimi.it>
 */

/* Bezier Rod */

#ifndef RODBEZJ_H
#define RODBEZJ_H

#include "joint.h"
#include "constltp.h"
#include "gauss.h"

extern const char* psRodBezNames[];


/* Bezier Rod - begin */

class RodBezier :
virtual public Elem, public Joint, public ConstitutiveLaw1DOwner {
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	const Vec3 fO;
	const Vec3 fA;
	const Vec3 fB;
	const Vec3 fI;
	doublereal dL0;

	Vec3 l1;
	Vec3 l2;
	doublereal dElle;
	doublereal dEpsilon;
	doublereal dEpsilonPrime;

	unsigned int iIntOrd;
	unsigned int iIntSeg;

	GaussDataIterator *gdi;

public:
	/* Costructor */
	RodBezier(unsigned int uL, const DofOwner* pDO,
			const ConstitutiveLaw1D* pCL,
			const StructNode* pN1, const StructNode* pN2,
			const Vec3& fOTmp, const Vec3& fATmp, 
			const Vec3& fBTmp, const Vec3& fITmp,
			doublereal dLength, bool bFromNodes, 
			const unsigned int iIntOrder, const unsigned int iIntSegments,
			flag fOut);

	/* Destructor */
	virtual ~RodBezier(void);

	/* Joint type */
	virtual Joint::Type GetJointType(void) const {
		return Joint::RODBEZIER;
	};

	/* Restart file contribute */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	virtual unsigned int iGetNumDof(void) const { 
		return 0;
	};

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 12;
	};
	
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	void AssVec(SubVectorHandler& WorkVec);

	void Output(OutputHandler& OH) const;

	/* initial assembly functions */

	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	}

	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 24; 
	};

	/* jacobian contribute during initial assembly */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
     			const VectorHandler& XCurr);
   
	/* residual contribue during initial assembly */   
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);   

	/* Inverse Dynamics */
	/* Inverse Dynamics capable element */
	virtual bool bInverseDynamics(void) const;

	/* Inverse Dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr);
	
	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics update */
	void Update(const VectorHandler& XCurr, 
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics after convergence */
	void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP, const VectorHandler& XPP);

	/* ******* PARALLEL SOLVER ******** */        
	 /* Outputs the type and the label of the nodes that are connected to
	  * the element. Useful for dof connection matrix assembly 	   */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	doublereal dGetPrivData(unsigned int i) const;

private:
	inline doublereal Ap1(doublereal u) { return (-3*u*u + 6*u -3); };
	inline doublereal Ap2(doublereal u) { return (9*u*u - 12*u + 3); };
	inline doublereal Ap3(doublereal u) { return (-9*u*u + 6*u); };
	inline doublereal Ap4(doublereal u) { return (3*u*u); };

};

/* RodBezier - end */

#endif /* RODBEZJ_H */ 















