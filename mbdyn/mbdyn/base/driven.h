/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#include "elem.h"
#include "drive.h"

#include "except.h"

class DrivenElem : /* virtual */ public Elem, protected DriveOwner {
protected: 
	Elem* pElem;
 
public:
	DrivenElem(const DriveCaller* pDC, const Elem* pE);
	~DrivenElem(void);

	inline void* pGet(void) const {
		ASSERT(pElem != NULL);
		return pElem->pGet();
	};
 
	virtual Elem* pGetElem(void) const {
		ASSERT(pElem != NULL);
		return pElem->pGetElem();
	};
 
	virtual ElemWithDofs* pGetElemWithDofs(void) const { 
		ASSERT(pElem != NULL);
		return pElem->pGetElemWithDofs();
	};

	virtual ElemGravityOwner* pGetElemGravityOwner(void) const { 
		ASSERT(pElem != NULL);
		return pElem->pGetElemGravityOwner();
	};

	virtual AerodynamicElem* pGetAerodynamicElem(void) const { 
		ASSERT(pElem != NULL);
		return pElem->pGetAerodynamicElem();
	};

	virtual InitialAssemblyElem* pGetInitialAssemblyElem(void) const { 
		ASSERT(pElem != NULL);
		return pElem->pGetInitialAssemblyElem();
	};
   
	virtual void Output(OutputHandler& OH) const;

	virtual void SetOutputFlag(flag f);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		ASSERT(pElem != NULL);
		return pElem->GetElemType(); 
	};

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
 
	/* esegue operazioni sui dof di proprieta' dell'elemento */
	virtual DofOrder::Order GetDofType(unsigned int i) const;
   
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

	/* *******PER IL SOLUTORE PARALLELO********
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps,
			unsigned int* NdLabels);
};

#endif /* DRIVEN_H */

