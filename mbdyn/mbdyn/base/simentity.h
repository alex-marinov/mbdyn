/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2002
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


#include <myassert.h>

/* include del programma */
#include <output.h>
#include <withlab.h>
#include <dofown.h>

/* SimulationEntity - begin */

/*
 * Functional class that introduces methods to handle the simulation.
 * 
 * Ancestor of nodes and elements.
 */

class SimulationEntity {
public:
	SimulationEntity(void);
	virtual ~SimulationEntity(void);

	/* 
	 * Ritorna il numero di DoFs.
	 * Non usa il DofOwner in quanto viene usata per generale il DofOwner
	 * stesso (per compatibilita' con gli elementi che generano gradi di 
	 * liberta' ed in previsione di nodi con un numero variabile di DoF) 
	 */
	virtual unsigned int iGetNumDof(void) const = 0;
	 
	/*
	 * Test di validita' di un indice. 
	 * Nota: gli indici vanno da 1 a iGetNumDofs()
	 */
	virtual flag fIsValidIndex(unsigned int i) const;
	
	/*
	 * Esegue operazioni sui DoF di proprieta' dell'elemento.
	 * In particolare ritorna il tipo di DoF in base all'indice i.
	 * Di default i DoF dei nodi sono assunti differenziali.
	 * Il tipo e' preso dall'enum DofOrder.
	 * Nota: gli indici sono in base 0, ovvero deve essere
	 * 0 < i < iGetNumDof()
	 * @see DofOrder
	 */   
	virtual DofOrder::Order SetDof(unsigned int i) const = 0;

	/* Metodi legati all'integrazione */
	
	/*
	 * Setta i valori iniziali dei DoF.
	 * Puo' essere usata per altre inizializzazioni prima di 
	 * iniziare l'integrazione 
	 */
	virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
	         
	/*
	 * Elaborazione vettori e dati prima della predizione.
	 * Per MultiStepIntegrator
	 */
	virtual void BeforePredict(VectorHandler& /* X */ ,
	   		      VectorHandler& /* XP */ ,
			      VectorHandler& /* XPrev */ ,
			      VectorHandler& /* XPPrev */ ) const;
	
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

	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(VectorHandler& X, VectorHandler& XP);
};

/* SimulationEntity - end */

#endif /* SIMENTITY_H */

