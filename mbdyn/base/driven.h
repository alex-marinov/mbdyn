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

/* Driven elements:
 * elements that are used depending on the (boolean) value
 * of a driver. Example: a driven joint is assembled only
 * if the driver is true, otherwise there is no joint and
 * the reaction unknowns are set to zero
 */

#ifndef DRIVEN_H
#define DRIVEN_H

#include "nestedelem.h"
#include "drive.h"

#include "except.h"

class DrivenElem : virtual public Elem,
	public NestedElem, protected DriveOwner {
protected:
	DataManager *pDM;
	SimulationEntity::Hints *pHints;
	bool bActive;

public:
	DrivenElem(DataManager *pDM, const DriveCaller* pDC,
			const Elem* pE, SimulationEntity::Hints *ph = 0);
	~DrivenElem(void);

	virtual bool bIsActive(void) const;

	virtual void Output(OutputHandler& OH) const;

	virtual void SetValue(DataManager *pdm,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* funzioni proprie */

	/*
	 * Elaborazione vettori e dati prima e dopo la predizione
	 * per MultiStepIntegrator */
	virtual void BeforePredict(VectorHandler& X,
			VectorHandler& XP,
			VectorHandler& XPrev,
			VectorHandler& XPPrev) const;

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	virtual void AfterConvergence(const VectorHandler& X,
     			const VectorHandler& XP);

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
	    		const VectorHandler& XCurr,
	    		const VectorHandler& XPrimeCurr);

     	virtual void AssMats(VariableSubMatrixHandler& WorkMatA,
 			VariableSubMatrixHandler& WorkMatB,
 			const VectorHandler& XCurr,
 			const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
     	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* Inverse Dynamics: */
	virtual void Update(const VectorHandler& XCurr,
			InverseDynamics::Order iOrder);

	/* inverse dynamics Jacobian matrix assembly */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* inverse dynamics residual assembly */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

	/* Inverse Dynamics: */
	virtual void AfterConvergence(const VectorHandler& X,
     		const VectorHandler& XP, const VectorHandler& XPP);

	/* InitialAssemblyElem */
public:
	virtual unsigned int iGetInitialNumDof(void) const;

	/* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
	 * conto del numero di dof che l'elemento definisce in questa fase e dei
	 * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
	 * la posizione e la velocita' dei nodi */
	virtual void InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr);

	/* ElemGravityOwner */
protected:
	virtual Vec3 GetS_int(void) const;
	virtual Mat3x3 GetJ_int(void) const;

	virtual Vec3 GetB_int(void) const;

	// NOTE: gravity owners must provide the momenta moment
	// with respect to the origin of the global reference frame!
	virtual Vec3 GetG_int(void) const;

public:
	virtual doublereal dGetM(void) const;
	Vec3 GetS(void) const;
	Mat3x3 GetJ(void) const;

	/* ElemDofOwner */
public:
	virtual void SetInitialValue(VectorHandler& X);
};

#endif /* DRIVEN_H */

