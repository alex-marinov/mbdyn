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

#ifndef SIMENTITY_H
#define SIMENTITY_H

#include <vector>

#include "myassert.h"

/* include del programma */
#include "output.h"
#include "withlab.h"
#include "dofown.h"
#include "drive.h"
#include "hint.h"
#include "invdyn.h"
#include "solverbase.h"
/* SimulationEntity - begin */

/*
 * Functional class that introduces methods to handle the simulation.
 * 
 * Ancestor of nodes and elements.
 *
 * Usage pattern:
 * 
 * pre-operation:
 *	iGetNumDof()
 *	bIsValidIndex()
 *	GetDofType()
 *	GetEqType()
 *	SetValue()
 *	iGetNumPrivData()
 *	iGetPrivDataIdx()
 *
 * runtime:
 * 	BeforePredict()		: prepare for prediction
 *	AfterPredict()		: account for predicted state
 *	Update()		: use converged solution
 *	AfterConvergence()	: account for conveged state
 *	dGetPrivData()		: get an internal state
 */

class MBDynParser;
class DataManager;

class SimulationEntity {
protected:
#if 0 
	/* punta a un vettore di due elementi che sono il valore
	 * iniziale dello stato e la sua derivata prima */
	const VectorHandler const** ppX0_Xp0;
#endif

public:
	SimulationEntity(void);
	virtual ~SimulationEntity(void);

	/* used to pass hints to SetValue */
	typedef std::vector<Hint *> Hints;

	/* 
	 * Ritorna il numero di DoFs.
	 * Non usa il DofOwner in quanto viene usata per generale il DofOwner
	 * stesso (per compatibilita' con gli elementi che generano gradi di 
	 * liberta' ed in previsione di nodi con un numero variabile di DoF) 
	 */
	virtual unsigned int iGetNumDof(void) const = 0;

	/*
	 * Describe the degrees of freedom
	 */
	virtual std::ostream& DescribeDof(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const = 0;
	virtual void DescribeDof(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const = 0;

	/*
	 * Describe the degrees of freedom
	 */
	virtual std::ostream& DescribeEq(std::ostream& out,
			const char *prefix = "",
			bool bInitial = false) const = 0;
	virtual void DescribeEq(std::vector<std::string>& desc,
			bool bInitial = false, int i = -1) const = 0;
	 
	/*
	 * Test di validita' di un indice. 
	 * Nota: gli indici vanno da 1 a iGetNumDofs()
	 */
	virtual bool bIsValidIndex(unsigned int i) const;
	
	/*
	 * Esegue operazioni sui DoF di proprieta' dell'elemento.
	 * In particolare ritorna il tipo di DoF in base all'indice i.
	 * Di default i DoF dei nodi sono assunti differenziali.
	 * Il tipo e' preso dall'enum DofOrder.
	 * Nota: gli indici sono in base 0, ovvero deve essere
	 * 0 < i < iGetNumDof()
	 * @see DofOrder
	 */   
	virtual DofOrder::Order GetDofType(unsigned int i) const = 0;

	/*
	 * Complementare di GetDofType(); dice che tipo di equazione
	 * corrisponde al dof i (ALGEBRAIC o DIFFERENTIAL).
	 */
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual DofOrder::Equality GetEqualityType(unsigned int i) const;

	/* Metodi legati all'integrazione */
	
	/*
	 * Setta i valori iniziali dei DoF.
	 * Puo' essere usata per altre inizializzazioni prima di 
	 * iniziare l'integrazione 
	 */
	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints* h = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	/*
	 * Elaborazione vettori e dati prima della predizione.
	 * Per MultiStepIntegrator
	 */
	virtual void BeforePredict(VectorHandler& /* X */ ,
	   		      VectorHandler& /* XP */ ,
		std::deque<VectorHandler*>& /* qXPr */ ,
		std::deque<VectorHandler*>& /* qXPPr */ ) const;
	
	/*
	 * Elaborazione vettori e dati dopo la predizione.
	 * Per MultiStepIntegrator
	 */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	/*
	 * Aggiorna dati in base alla soluzione. 
	 * Usata per operazioni aggiuntive al semplice aggiornamento additivo,
	 * effettuato gia' dall'integratore.
	 */
	virtual void Update(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr);
	
	virtual void DerivativesUpdate(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr);
	
	/* Inverse Dynamics:*/
	virtual void Update(const VectorHandler& XCurr, 
		       InverseDynamics::Order iOrder);
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
	
	/* Inverse Dynamics:*/
	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP,
			const VectorHandler& XPP);

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

	/*
	 * Contributes to output on a stream and/or NetCDF
	 */
	virtual std::ostream& OutputAppend(std::ostream& out) const;
	virtual void NetCDFOutputAppend(OutputHandler& OH) const;
	virtual void OutputAppendPrepare(OutputHandler& OH, const std::string& name);

	virtual void ReadInitialState(MBDynParser& HP);

};

/* SimulationEntity - end */

#endif /* SIMENTITY_H */

