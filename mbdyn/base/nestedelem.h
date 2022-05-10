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

/* Nested elements
 */

#ifndef NESTEDELEM_H
#define NESTEDELEM_H

#include "elem.h"
#include "aerodyn.h"
#include "gravity.h"

#include "except.h"

class NestedElem : virtual public Elem,
public InitialAssemblyElem,
public AerodynamicElem,
public ElemGravityOwner
{
protected:
	Elem* pElem;

public:
	NestedElem(const Elem* pE);
	~NestedElem(void);

	virtual Elem *pGetElem(void) const;

	virtual void OutputPrepare(OutputHandler& OH);
	virtual void Output(OutputHandler& OH) const;

	virtual void SetOutputFlag(flag f);

	virtual void SetValue(DataManager *pdm,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const;

	/* funzioni di servizio */

	/*
	 * Il metodo iGetNumDof() serve a ritornare il numero di gradi
	 * di liberta' propri che l'elemento definisce. Non e' virtuale
	 * in quanto serve a ritornare 0 per gli elementi che non possiedono
	 * gradi di liberta'.
	 * Viene usato nella costruzione dei DofOwner e quindi deve essere
	 * indipendente da essi. In genere non comporta overhead in quanto
	 * il numero di dof aggiunti da un tipo e' una costante e non
	 * richede dati propri.
	 * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner
	 * dell'oggetto.
	 * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
	 * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
	 * dof propri.
	 * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
	 * E' usato per completare i singoli Dof relativi all'elemento.
	 */

	/*
	 * ritorna il numero di Dofs per gli elementi che sono
	 * anche DofOwners
	 */
	virtual unsigned int iGetNumDof(void) const;

	/* inherited from SimulationEntity */
	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "", bool bInitial = false) const;
	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false, int i = -1) const;
	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "", bool bInitial = false) const;
	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false, int i = -1) const;

	/* esegue operazioni sui dof di proprieta' dell'elemento */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	virtual DofOrder::Order GetEqType(unsigned int i) const;

	/* funzioni proprie */

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

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
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const;

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const;

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* *******PER IL SOLUTORE PARALLELO********
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual int GetNumConnectedNodes(void) const;
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;

	// inverse dynamics
	/* inverse dynamics capable element */
	virtual bool bInverseDynamics(void) const;

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

	// end of inverse dynamics

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

	/* AerodynamicElem */
public:
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const;
	virtual bool NeedsAirProperties(void) const;
	virtual const InducedVelocity *pGetInducedVelocity(void) const;
	virtual void PutAirProperties(const AirProperties* pAP);

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
	virtual const DofOwner* pGetDofOwner(void) const;
	virtual integer iGetFirstIndex(void) const;
	virtual void SetInitialValue(VectorHandler& X);

	/* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;
};

#endif // NESTEDELEM_H

