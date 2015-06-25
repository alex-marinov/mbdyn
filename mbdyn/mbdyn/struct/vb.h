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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Deformable hinges */

#ifndef VB_H
#define VB_H

#include "joint.h"
#include "constltp.h"


/* ViscousBody - begin */

class ViscousBody :
virtual public Elem, public Joint, public ConstitutiveLaw6DOwner {
protected:
	const StructNode* pNode;
	mutable Vec3 tilde_f;
	mutable Mat3x3 tilde_Rh;

	OrientationDescription od;

	Vec6 tilde_kPrime;

	bool bFirstRes;

	Vec3 f;
	Mat3x3 Rh;

	Vec6 F;
	Mat6x6 FDEPrime;

	void AssMats(FullSubMatrixHandler& WMA,
		FullSubMatrixHandler& WMB,
		doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec);

public:
	/* Costruttore non banale */
	ViscousBody(unsigned int uL,
		const DofOwner* pDO,
		const ConstitutiveLaw6D* pCL,
		const StructNode* pN,
		const Vec3& tilde_f,
		const Mat3x3& tilde_Rh,
		const OrientationDescription& od,
		flag fOut);

	/* Distruttore */
	virtual ~ViscousBody(void);

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const {
		return Joint::VISCOUSBODY;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void Output(OutputHandler& OH) const;

	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual unsigned int iGetNumDof(void) const {
		return 0;
	};

	virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const {
		return DofOrder::UNKNOWN;
	};

	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const {
		return 0;
	};

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

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;
 
	/* Inverse Dynamics residual assembly */
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler&  XPrimeCurr,
		const VectorHandler&  XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	/* Tipo di DeformableDispHinge */
	virtual ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOUS;
	};

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* Aggiorna le deformazioni ecc. */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

};

/* ViscousBody - end */

#endif // VB_H

