/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Joint regularizations, Elem::Type JOINT_REGULARIZATION */

#ifndef JOINTREG_H
#define JOINTREG_H

#include <cfloat>

#include "joint.h"

/* JointRegularization - begin */

class JointRegularization
: virtual public Elem, public InitialAssemblyElem
{
protected:
	const Joint *pJ;

public:
	/* JointRegularization types */
	enum Type {
		UNKNOWN = -1,

		TIKHONOV_REGULARIZATION = 0,
		DYNAMIC_REGULARIZATION,
		JACOBIAN_REGULARIZATION,

		LASTJOINTREGTYPE
	};

public:
	JointRegularization(unsigned int uL, const Joint *j, flag fOut);
	virtual ~JointRegularization(void);

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::JOINT_REGULARIZATION;
	};

	/* Tipo di joint regularization */
	virtual JointRegularization::Type
	GetJointRegularizationType(void) const = 0;

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "  joint regularization: " << GetLabel();
	};

	virtual unsigned int iGetInitialNumDof(void) const;
};

/* JointRegularization - end */


/* TikhonovRegularization - begin */

class TikhonovRegularization
: virtual public Elem, public JointRegularization
{
protected:
	std::vector<doublereal> dC;

public:
	TikhonovRegularization(unsigned int uL,
		const Joint *j,
		const std::vector<doublereal>& c,
		flag fOut);
	virtual ~TikhonovRegularization(void);

	virtual JointRegularization::Type
	GetJointRegularizationType(void) const;

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;

	/* Inverse Dynamics Jacobian matrix assembly */
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const;
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* TikhonovRegularization - end */


/* DynamicRegularization - begin */

class DynamicRegularization
: virtual public Elem, public JointRegularization
{
protected:
	std::vector<doublereal> dC;
	std::vector<doublereal> dLambda;

public:
	DynamicRegularization(unsigned int uL,
		const Joint *j,
		const std::vector<doublereal>& c,
		flag fOut);
	virtual ~DynamicRegularization(void);

	virtual JointRegularization::Type
	GetJointRegularizationType(void) const;

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	virtual void
	AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const;
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* DynamicRegularization - end */

/* JacobianRegularization - begin */

class JacobianRegularization
: virtual public Elem, public JointRegularization
{
protected:
	std::vector<doublereal> dC;

public:
	JacobianRegularization(unsigned int uL,
		const Joint *j,
		const std::vector<doublereal>& c,
		flag fOut);
	virtual ~JacobianRegularization(void);

	virtual JointRegularization::Type
	GetJointRegularizationType(void) const;

	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr);

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void
	InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const;
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);
};

/* JacobianRegularization - end */

/* Lettura joint regularizations */
class DataManager;
class MBDynParser;

extern Elem *
ReadJointRegularization(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel);

#endif // JOINTREG_H
